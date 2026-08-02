// Azimuthal mode content of the field: is the instability m = 0 or m = 1?
//
// This is the discriminator between the two candidates for what destroys the
// toroidal-dominated configuration.
//
//   MRI    feeds on the rotational shear, its fastest mode is AXISYMMETRIC
//          (m = 0) with vertical structure, and its growth rate is set by
//          the shear alone: sigma = (1/2)|dOmega/dln varpi|, which for the
//          j-constant law of these models peaks at 2.05/s at varpi = R_eq.
//   Tayler feeds on the toroidal field's own energy, needs no rotation, and
//          its fastest mode is the m = 1 kink.
//
// At 192^3 the m = 0 and m = 1 components of the poloidal field grew at the
// same rate and the test was inconclusive. The rate also RISES with
// resolution -- 1.369/s at 192^3 against 2.007/s at 256^3, climbing toward
// the analytic MRI value -- which is the signature of an under-resolved MRI,
// lambda_MRI/dx going 13.4 to 17.9 against the ~25 it needs.
//
// The projection is done as a direct sum over cells rather than by
// interpolating onto rings, which needs no interpolation and no slice file:
//
//     a_m(bin) = |sum_cells f exp(-i m phi)| / N_cells(bin)
//
// binned in cylindrical radius, then averaged over the bins that hold enough
// cells. Reported for B_phi and for B_z, the latter being the poloidal
// component that an m = 0 MRI channel flow shows up in most cleanly.
//
// Fields are read in Castro's Heaviside-Lorentz state and converted to gauss,
// as in fbtbp -- the ratios do not care, but the absolute levels do.
//
// Usage:
//     fmodes [--rho-cut RHO] [--nbins N] plotfile [plotfile ...]

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>

#include <cmath>
#include <array>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace amrex;

namespace {

constexpr Real castro_to_gauss = 3.5449077018110318_rt;   // sqrt(4 pi)
constexpr int MMAX = 3;

bool has_var (const Vector<std::string>& names, const std::string& want)
{
    for (auto const& n : names) { if (n == want) { return true; } }
    return false;
}

void process (const std::string& pltfile, Real rho_cut, int nbins, Real rmax)
{
    PlotFileData pf(pltfile);
    if (pf.spaceDim() != 3 || pf.coordSys() != 0) {
        amrex::Print() << "# " << pltfile << ": skipped, needs 3-d Cartesian\n";
        return;
    }
    for (auto const& want : {"density", "B_x", "B_y", "B_z"}) {
        if (! has_var(pf.varNames(), want)) {
            amrex::Print() << "# " << pltfile << ": skipped, no " << want << '\n';
            return;
        }
    }

    const Array<Real,AMREX_SPACEDIM> problo = pf.probLo();
    const int lev = pf.finestLevel();
    const Array<Real,AMREX_SPACEDIM> dx = pf.cellSize(lev);

    // accumulators: [bin][m], for B_phi and B_z
    std::vector<std::complex<Real>> acc_t(std::size_t(nbins) * (MMAX + 1), {0.0, 0.0});
    std::vector<std::complex<Real>> acc_z(std::size_t(nbins) * (MMAX + 1), {0.0, 0.0});
    std::vector<long> count(nbins, 0);

    const MultiFab& rho = pf.get(lev, "density");
    const MultiFab& bxf = pf.get(lev, "B_x");
    const MultiFab& byf = pf.get(lev, "B_y");
    const MultiFab& bzf = pf.get(lev, "B_z");

    const Real dr = rmax / Real(nbins);

    for (MFIter mfi(rho); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& ra = rho.const_array(mfi);
        auto const& xa = bxf.const_array(mfi);
        auto const& ya = byf.const_array(mfi);
        auto const& za = bzf.const_array(mfi);

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if (ra(i,j,k) <= rho_cut) { continue; }

            const Real x = problo[0] + (Real(i) + 0.5_rt) * dx[0];
            const Real y = problo[1] + (Real(j) + 0.5_rt) * dx[1];
            const Real varpi = std::sqrt(x*x + y*y);
            if (varpi <= 0.0 || varpi >= rmax) { continue; }

            const int b = static_cast<int>(varpi / dr);
            if (b < 0 || b >= nbins) { continue; }

            const Real phi = std::atan2(y, x);
            const Real btor = (-y * xa(i,j,k) + x * ya(i,j,k)) / varpi * castro_to_gauss;
            const Real bz = za(i,j,k) * castro_to_gauss;

            for (int m = 0; m <= MMAX; ++m) {
                const std::complex<Real> e(std::cos(m * phi), -std::sin(m * phi));
                acc_t[std::size_t(b) * (MMAX + 1) + m] += btor * e;
                acc_z[std::size_t(b) * (MMAX + 1) + m] += bz * e;
            }
            ++count[b];
        }}}
    }

    // Average |a_m| over the bins that hold enough cells for the projection to
    // mean anything. A bin with a handful of cells cannot resolve m = 3.
    std::array<Real,MMAX+1> at{}, az{};
    int used = 0;
    for (int b = 0; b < nbins; ++b) {
        if (count[b] < 200) { continue; }
        ++used;
        for (int m = 0; m <= MMAX; ++m) {
            at[m] += std::abs(acc_t[std::size_t(b)*(MMAX+1)+m]) / Real(count[b]);
            az[m] += std::abs(acc_z[std::size_t(b)*(MMAX+1)+m]) / Real(count[b]);
        }
    }
    if (used == 0) { amrex::Print() << "# " << pltfile << ": no bin with enough cells\n"; return; }
    for (int m = 0; m <= MMAX; ++m) { at[m] /= Real(used); az[m] /= Real(used); }

    std::cout << std::setw(11) << std::fixed << std::setprecision(5) << pf.time()
              << std::scientific << std::setprecision(4);
    for (int m = 0; m <= MMAX; ++m) { std::cout << std::setw(12) << at[m]; }
    for (int m = 0; m <= MMAX; ++m) { std::cout << std::setw(12) << az[m]; }
    std::cout << std::setw(12) << (at[0] > 0 ? at[1] / at[0] : 0.0)
              << std::setw(12) << (az[0] > 0 ? az[1] / az[0] : 0.0)
              << std::setw(7) << used
              << '\n' << std::flush;
}

} // namespace

void main_main ()
{
    const int narg = amrex::command_argument_count();

    Real rho_cut = 1.0e5;
    Real rmax = 5.0e8;      // just outside R_eq = 3.917e8 cm
    int nbins = 24;
    std::vector<std::string> plotfiles;

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "--rho-cut")   { rho_cut = std::stod(amrex::get_command_argument(++farg)); }
        else if (name == "--nbins"){ nbins   = std::stoi(amrex::get_command_argument(++farg)); }
        else if (name == "--rmax") { rmax    = std::stod(amrex::get_command_argument(++farg)); }
        else { break; }
        ++farg;
    }
    for (; farg <= narg; ++farg) { plotfiles.push_back(amrex::get_command_argument(farg)); }

    if (plotfiles.empty()) {
        amrex::Print()
            << "\n Azimuthal mode content of B_phi and B_z: m = 0 (MRI) or m = 1 (Tayler).\n\n"
            << " Usage:\n    fmodes [--rho-cut RHO] [--nbins N] [--rmax R] plotfile [...]\n\n";
        return;
    }

    amrex::Print() << "# azimuthal projection over cells above rho_cut = " << rho_cut
                   << ", binned in varpi out to " << rmax << " cm in " << nbins << " bins\n";
    amrex::Print() << "# amplitudes in GAUSS; m1/m0 is the ratio that discriminates the mode\n";
    amrex::Print() << "#          t"
                   << "  Bphi_m0     Bphi_m1     Bphi_m2     Bphi_m3"
                   << "     Bz_m0       Bz_m1       Bz_m2       Bz_m3"
                   << "   Bphi_m1/m0    Bz_m1/m0  nbin\n";

    for (auto const& p : plotfiles) { process(p, rho_cut, nbins, rmax); }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
