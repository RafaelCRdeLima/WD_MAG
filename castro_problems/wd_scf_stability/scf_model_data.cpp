#include <scf_model_data.H>

#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

namespace scf_model
{
    AMREX_GPU_MANAGED int n_varpi = 0;
    AMREX_GPU_MANAGED int n_z = 0;

    AMREX_GPU_MANAGED amrex::Real varpi_lo = 0.0;
    AMREX_GPU_MANAGED amrex::Real varpi_hi = 0.0;
    AMREX_GPU_MANAGED amrex::Real z_lo = 0.0;
    AMREX_GPU_MANAGED amrex::Real z_hi = 0.0;
    AMREX_GPU_MANAGED amrex::Real dvarpi = 0.0;
    AMREX_GPU_MANAGED amrex::Real dz = 0.0;

    AMREX_GPU_MANAGED amrex::Real* dens = nullptr;
    AMREX_GPU_MANAGED amrex::Real* a_phi = nullptr;
    AMREX_GPU_MANAGED amrex::Real* a_z = nullptr;

    AMREX_GPU_MANAGED amrex::Real dens_max = 0.0;

    namespace
    {
        // The writer emits a uniform grid; this checks it instead of trusting
        // it, because the interpolation indexes arithmetically and a
        // non-uniform axis would be read as uniform and silently give the
        // wrong field rather than failing.
        void require_uniform(const std::vector<amrex::Real>& v,
                             const char* name,
                             amrex::Real& lo, amrex::Real& hi, amrex::Real& d)
        {
            if (v.size() < 2) {
                amrex::Abort(std::string("scf_model: axis ") + name
                             + " has fewer than two points");
            }
            lo = v.front();
            hi = v.back();
            d = (hi - lo) / static_cast<amrex::Real>(v.size() - 1);
            amrex::Real worst = 0.0;
            for (std::size_t i = 0; i < v.size(); ++i) {
                amrex::Real expect = lo + d * static_cast<amrex::Real>(i);
                worst = std::max(worst, std::abs(v[i] - expect));
            }
            if (worst > 1.0e-8 * std::abs(hi - lo)) {
                amrex::Abort(std::string("scf_model: axis ") + name
                             + " is not uniform (worst deviation "
                             + std::to_string(worst) + " cm)");
            }
        }
    }

    void read(const std::string& filename)
    {
        std::ifstream f(filename);
        if (!f) {
            amrex::Abort("scf_model: cannot open model file " + filename);
        }

        // Header comments carry the provenance written by the exporter. They
        // are echoed, not parsed: the point is that the log of a run records
        // which configuration it was, including the verification results the
        // writer stored there.
        std::string line;
        std::streampos last = f.tellg();
        while (std::getline(f, line)) {
            if (line.empty()) { last = f.tellg(); continue; }
            if (line[0] != '#') { break; }
            amrex::Print() << "  [scf_model] " << line << "\n";
            last = f.tellg();
        }

        std::istringstream dims(line);
        if (!(dims >> n_varpi >> n_z) || n_varpi < 2 || n_z < 2) {
            amrex::Abort("scf_model: bad dimension line in " + filename);
        }

        std::vector<amrex::Real> vp(n_varpi), zz(n_z);
        for (int i = 0; i < n_varpi; ++i) {
            if (!(f >> vp[i])) { amrex::Abort("scf_model: truncated varpi axis"); }
        }
        for (int j = 0; j < n_z; ++j) {
            if (!(f >> zz[j])) { amrex::Abort("scf_model: truncated z axis"); }
        }
        require_uniform(vp, "varpi", varpi_lo, varpi_hi, dvarpi);
        require_uniform(zz, "z", z_lo, z_hi, dz);

        const std::size_t n = static_cast<std::size_t>(n_varpi)
                            * static_cast<std::size_t>(n_z);
        auto* arena = amrex::The_Managed_Arena();
        dens = static_cast<amrex::Real*>(arena->alloc(n * sizeof(amrex::Real)));
        a_phi = static_cast<amrex::Real*>(arena->alloc(n * sizeof(amrex::Real)));
        a_z = static_cast<amrex::Real*>(arena->alloc(n * sizeof(amrex::Real)));

        dens_max = 0.0;
        for (std::size_t idx = 0; idx < n; ++idx) {
            if (!(f >> dens[idx] >> a_phi[idx] >> a_z[idx])) {
                amrex::Abort("scf_model: truncated data block in " + filename);
            }
            dens_max = std::max(dens_max, dens[idx]);
        }

        amrex::Print() << "  [scf_model] " << n_varpi << " x " << n_z
                       << " grid, varpi to " << varpi_hi << " cm, z in ["
                       << z_lo << ", " << z_hi << "] cm\n"
                       << "  [scf_model] max density " << dens_max
                       << " g/cm^3\n";
    }

    void deallocate()
    {
        auto* arena = amrex::The_Managed_Arena();
        if (dens != nullptr) { arena->free(dens); dens = nullptr; }
        if (a_phi != nullptr) { arena->free(a_phi); a_phi = nullptr; }
        if (a_z != nullptr) { arena->free(a_z); a_z = nullptr; }
    }
}
