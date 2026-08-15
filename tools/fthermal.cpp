// Thermal diagnostic for campaign HZ: where the heat is, and how much of it.
//
// Nothing existing measures temperature. fbtbp.cpp gives mass, radii, field
// and rotation; the run log gives MAXIMUM DENSITY and mass but no T. Under
// ztwd that was correct -- temperature was thermodynamically inert and the
// column would have been noise. Under helmholtz it is the measurement.
//
// WHAT THIS HAS TO SETTLE
//
// The 256^3 helmholtz window (DIARIO 10.1) found a shell at rho ~ 2.5e8,
// about 0.3 R_eq out, at 3.37e9 K while the core sat at 2.25e7 and the
// ambient at 8.6e6. Two readings, and a single run cannot separate them:
//
//   (a) it is the field settling out of magnetohydrostatic equilibrium, which
//       ztwd discarded and helmholtz keeps;
//   (b) it is numerical, through dual_energy_eta2 redefining UEINT as E - K
//       in cells where that subtraction is ill-conditioned.
//
// dir_hz192ctl is the same star with field_scale = 0, so (a) has nothing to
// draw on there. The comparison is only as good as the numbers it is made
// from, and one number -- T_max over the box -- is not enough: a single hot
// cell at the surface would look identical to a heated shell. So this reports
// WHERE the heat sits and HOW MUCH MASS carries it.
//
// THE ERROR THIS IS BUILT TO PREVENT
//
// DIARIO 10.1 records writing that the heating was "eight times more energy
// than the field had to give". That number came from spreading the shell
// temperature over ~30% of the mass. The heat is not spread -- it is in a
// shell, and at ~5% of the mass the same temperature costs 2e49 erg, which
// the field's 6.06e49 can afford. A mean over the star gets this wrong in the
// direction that discards the result. Hence the mass fractions above
// thresholds and the integrated ion thermal energy, which cannot be fooled
// that way: E_ion is the volume integral itself, not a mean times a guessed
// mass.
//
// THE ENERGY, AND ITS CONVENTION
//
// E_ion = int (3/2) rho k T / (Abar m_u) dV, over cells above the density cut.
//
// Abar = 4 from mu2.net -- one species, A = 4, Z = 2, mu_e = 2, matching the
// SCF's own assumption. This is the IONS only. Degenerate electrons carry
// far more energy but almost none of it is thermal, and it is the thermal
// part that the field can have paid for. Comparing E_ion(t) - E_ion(0)
// against the 6.06e49 erg the field loses is the whole point of the column.
//
// Usage:
//     fthermal [--rho-cut RHO] [--abar A] plotfile [plotfile ...]

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace amrex;

namespace {

constexpr Real K_BOLTZ = 1.380649e-16;    // erg/K
constexpr Real M_U     = 1.66053907e-24;  // g
constexpr Real MSUN    = 1.98892e33;      // g, as in fbtbp.cpp

// Temperature thresholds for the mass fractions. 1e8 is "warmed at all"
// against the 1e7 the run starts at; 1e9 and 2e9 bracket the 3.37e9 shell
// seen at 256^3, so a shell of that kind cannot hide inside a single bin.
constexpr int NTHRESH = 4;
constexpr Real THRESH[NTHRESH] = {1.0e8, 5.0e8, 1.0e9, 2.0e9};

// Returns false if the plotfile lacks the fields we need, so that a sweep
// over many files reports and moves on rather than aborting on the first.
bool has_var (const Vector<std::string>& names, const std::string& want)
{
    for (auto const& n : names) {
        if (n == want) { return true; }
    }
    return false;
}

void process (const std::string& pltfile, Real rho_cut, Real abar)
{
    PlotFileData pf(pltfile);

    if (pf.spaceDim() != 3 || pf.coordSys() != 0) {
        amrex::Print() << "# " << pltfile << ": skipped, needs 3-d Cartesian\n";
        return;
    }

    const Vector<std::string>& names = pf.varNames();
    for (auto const& want : {"density", "Temp"}) {
        if (! has_var(names, want)) {
            // A ztwd plotfile has no Temp, or carries an inert one. Say which
            // file was skipped rather than printing a row of zeros that would
            // read as a measurement.
            amrex::Print() << "# " << pltfile << ": skipped, no " << want << '\n';
            return;
        }
    }

    const Array<Real,AMREX_SPACEDIM> problo = pf.probLo();
    const int fine_level = pf.finestLevel();

    // Extrema over the star. rho_max stands in for the central density: the
    // star is centrally condensed and the peak has stayed at the centre in
    // every snapshot checked, while an actual centre-of-the-box sample would
    // move with the star if it drifts.
    Real rho_max = 0.0;
    Real t_max   = 0.0;
    // The ambient is diagnosed separately rather than dropped. At 256^3 it
    // COOLED slightly while the shell heated, and that asymmetry is evidence
    // the heat is not junk leaking in from the boundary.
    Real t_amb_max = 0.0;

    // Sums over the star.
    Real mass    = 0.0;   // g
    Real mt      = 0.0;   // int rho T dV, for the mass-weighted mean
    Real e_ion   = 0.0;   // erg
    Real m_above[NTHRESH] = {0.0, 0.0, 0.0, 0.0};

    // Where the hottest cell is. Carried through the host loop rather than a
    // reduction because the position has to travel with the value, and a
    // ParReduce on T alone would give the temperature and lose the location --
    // which is exactly the information that separates a heated shell from a
    // single bad surface cell.
    Real hot_rho = 0.0, hot_varpi = 0.0, hot_z = 0.0, hot_r = 0.0;

    // Mass-weighted T inside the core, on the same 0.15 R_eq window fbtbp.cpp
    // uses for Om_core, so the two diagnostics refer to the same region.
    constexpr Real R_EQ = 3.917259e8;    // of the initial model
    Real t_core_num = 0.0, t_core_den = 0.0;

    for (int ilev = 0; ilev <= fine_level; ++ilev) {

        const Array<Real,AMREX_SPACEDIM> dx = pf.cellSize(ilev);
        const Real dv = dx[0] * dx[1] * dx[2];

        const MultiFab& rho = pf.get(ilev, "density");
        const MultiFab& tmp = pf.get(ilev, "Temp");

        const Real xlo = problo[0];
        const Real ylo = problo[1];
        const Real zlo = problo[2];
        const Real dx0 = dx[0];
        const Real dx1 = dx[1];
        const Real dx2 = dx[2];

        // A cell covered by a finer level must not be counted twice. The
        // finest level gets an all-zero mask so the loop below stays a single
        // code path.
        iMultiFab mask;
        if (ilev < fine_level) {
            IntVect ratio{pf.refRatio(ilev)};
            mask = makeFineMask(pf.boxArray(ilev), pf.DistributionMap(ilev),
                                pf.boxArray(ilev+1), ratio);
        } else {
            mask.define(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, 0);
            mask.setVal(0);
        }

        const Real e_ion_coeff = 1.5_rt * K_BOLTZ / (abar * M_U);

        // A plain host loop, for the reason fbtbp.cpp gives for its rotation
        // block: this is past what a tuple reduction is comfortable with --
        // eleven accumulators plus a position that has to follow a maximum --
        // and it is not the hot path. A 192^3 plotfile is seconds either way.
        for (MFIter mfi(rho); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const& ra = rho.const_array(mfi);
            auto const& ta = tmp.const_array(mfi);
            auto const& mk = mask.const_array(mfi);
            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (mk(i,j,k) != 0) { continue; }
                const Real r = ra(i,j,k);
                const Real T = ta(i,j,k);

                if (r <= rho_cut) {
                    // Ambient: the maximum only, and no contribution to any
                    // integral over the star.
                    t_amb_max = amrex::max(t_amb_max, T);
                    continue;
                }

                const Real x = xlo + (Real(i) + 0.5_rt) * dx0;
                const Real y = ylo + (Real(j) + 0.5_rt) * dx1;
                const Real z = zlo + (Real(k) + 0.5_rt) * dx2;
                const Real varpi = std::sqrt(x*x + y*y);
                const Real dm = r * dv;

                rho_max = amrex::max(rho_max, r);
                if (T > t_max) {
                    t_max     = T;
                    hot_rho   = r;
                    hot_varpi = varpi;
                    hot_z     = z;
                    hot_r     = std::sqrt(x*x + y*y + z*z);
                }

                mass  += dm;
                mt    += dm * T;
                e_ion += e_ion_coeff * dm * T;
                for (int n = 0; n < NTHRESH; ++n) {
                    if (T > THRESH[n]) { m_above[n] += dm; }
                }
                if (varpi < 0.15_rt * R_EQ) { t_core_num += dm * T; t_core_den += dm; }
            }}}
        }
    }

    ParallelDescriptor::ReduceRealMax(rho_max);
    ParallelDescriptor::ReduceRealMax(t_amb_max);
    ParallelDescriptor::ReduceRealSum(mass);
    ParallelDescriptor::ReduceRealSum(mt);
    ParallelDescriptor::ReduceRealSum(e_ion);
    ParallelDescriptor::ReduceRealSum(m_above, NTHRESH);
    ParallelDescriptor::ReduceRealSum(t_core_num);
    ParallelDescriptor::ReduceRealSum(t_core_den);

    // The location of the maximum has to come from the rank that HOLDS the
    // maximum, so a plain ReduceRealMax on t_max would leave hot_rho and the
    // rest belonging to whichever rank happened to be the I/O processor.
    // Gathering the (T, position) pairs and picking the winner on rank 0 keeps
    // them together. Five doubles per rank, once per plotfile.
    {
        const int nproc = ParallelDescriptor::NProcs();
        const int ioproc = ParallelDescriptor::IOProcessorNumber();
        Real mine[5] = {t_max, hot_rho, hot_varpi, hot_z, hot_r};
        std::vector<Real> all(static_cast<std::size_t>(5 * nproc), 0.0);
        ParallelDescriptor::Gather(mine, 5, all.data(), 5, ioproc);
        if (ParallelDescriptor::IOProcessor()) {
            for (int p = 0; p < nproc; ++p) {
                if (all[5*p] > t_max) {
                    t_max     = all[5*p + 0];
                    hot_rho   = all[5*p + 1];
                    hot_varpi = all[5*p + 2];
                    hot_z     = all[5*p + 3];
                    hot_r     = all[5*p + 4];
                }
            }
        }
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << std::setw(11) << std::fixed << std::setprecision(5) << pf.time()
                  << std::scientific << std::setprecision(5)
                  << std::setw(14) << rho_max
                  << std::setw(14) << (t_core_den > 0.0 ? t_core_num / t_core_den : 0.0)
                  << std::setw(14) << (mass > 0.0 ? mt / mass : 0.0)
                  << std::setw(14) << t_max
                  << std::setw(14) << t_amb_max
                  // Where the hottest cell is. rho there separates "inside the
                  // star" from "at the surface"; varpi/R_eq and |z| say whether
                  // it sits in a shell or on the axis.
                  << std::setw(14) << hot_rho
                  << std::fixed << std::setprecision(4)
                  << std::setw(11) << hot_varpi / R_EQ
                  << std::setw(11) << hot_z / R_EQ
                  << std::setw(11) << hot_r / R_EQ
                  // Mass fractions, which is what the "eight times too much
                  // energy" error came down to: a mean over the wrong mass.
                  << std::scientific << std::setprecision(4)
                  << std::setw(13) << (mass > 0.0 ? m_above[0] / mass : 0.0)
                  << std::setw(13) << (mass > 0.0 ? m_above[1] / mass : 0.0)
                  << std::setw(13) << (mass > 0.0 ? m_above[2] / mass : 0.0)
                  << std::setw(13) << (mass > 0.0 ? m_above[3] / mass : 0.0)
                  << std::setprecision(5)
                  << std::setw(14) << e_ion
                  << std::setw(14) << mass / MSUN
                  << '\n' << std::flush;
    }
}

} // namespace

void main_main ()
{
    const int narg = amrex::command_argument_count();

    Real rho_cut = 1.0e5;   // ambient is 2.0e4, as in fbtbp.cpp
    Real abar    = 4.0;     // mu2.net: one species, A = 4, Z = 2
    std::vector<std::string> plotfiles;

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "--rho-cut") {
            rho_cut = std::stod(amrex::get_command_argument(++farg));
        } else if (name == "--abar") {
            abar = std::stod(amrex::get_command_argument(++farg));
        } else {
            break;
        }
        ++farg;
    }
    for (; farg <= narg; ++farg) {
        plotfiles.push_back(amrex::get_command_argument(farg));
    }

    if (plotfiles.empty()) {
        amrex::Print()
            << "\n"
            << " Thermal diagnostic for campaign HZ: where the heat is, and how much.\n"
            << "\n"
            << " Usage:\n"
            << "    fthermal [--rho-cut RHO] [--abar A] plotfile [plotfile ...]\n"
            << "\n"
            << " args --rho-cut RHO : skip cells at or below this density (default 1.0e5)\n"
            << "      --abar A      : mean ion mass for E_ion (default 4.0, from mu2.net)\n"
            << "\n";
        return;
    }

    amrex::Print() << "# rho_cut = " << rho_cut << ", Abar = " << abar << '\n';
    amrex::Print() << "# temperatures in K; densities in g/cm^3; E_ion = int (3/2) rho k T /"
                      " (Abar m_u) dV in erg, over cells above rho_cut\n";
    amrex::Print() << "# positions of the hottest cell in units of the initial R_eq ="
                      " 3.917259e8 cm; f_T are MASS fractions above each threshold\n";
    amrex::Print() << "# the field has 6.06e49 erg to give: compare E_ion(t) - E_ion(0)\n";
    amrex::Print() << "#          t       rho_max        T_core       T_mean         T_max"
                   << "     T_amb_max   rho_at_Tmax  varpi_hot     z_hot     r_hot"
                   << "     f_T>1e8     f_T>5e8     f_T>1e9     f_T>2e9"
                   << "         E_ion       M/Msun\n";

    for (auto const& pltfile : plotfiles) {
        process(pltfile, rho_cut, abar);
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
