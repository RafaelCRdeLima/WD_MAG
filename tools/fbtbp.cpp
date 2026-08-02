// Toroidal/poloidal magnetic energy split from a Castro MHD plotfile.
//
// The emag_density and etor_density derives in the plotfile are corrupt.
// They, and Castro's own Div_B, register several face-centred components
// (Mag_Type_x/y/z) into one cell-centred FAB and then read dat(i+1,...),
// which runs off the end of the data box. Measured on dir_rot96/plt00000,
// before a single step: Div_B reached 3.1e146 where it should be ~0, and
// emag_density reached 4.2e306 where B^2/2 is ~4e25. Squaring that overflows
// to inf, the inf lands in the plotfile metadata, and every AMReX plotfile
// tool then aborts in the header parse -- all 70 plotfiles of the 192^3 run
// are unreadable until the metadata is patched (see tools/patch_plotfile_inf.sh).
//
// B_x, B_y and B_z are written correctly at every time checked, so the split
// is rebuilt here from them:
//
//     B_tor   = (-y B_x + x B_y) / varpi
//     B_pol^2 = B_x^2 + B_y^2 + B_z^2 - B_tor^2
//
// with the rotation axis taken to be z, which it is for these models. The
// energy convention is 0.5*B^2, matching ca_deremag, so the numbers are
// directly comparable to what the derives were meant to produce.
//
// Cells at or below the density cut are skipped so that the ambient, which
// fills the box out to the boundary, does not dominate the volume integral.
//
// Usage:
//     fbtbp [--rho-cut RHO] plotfile [plotfile ...]
//
// Prints one row per plotfile: time, E_tor, E_pol, E_tor/E_pol, E_tor/E_mag.

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace amrex;

namespace {

// Returns false if the plotfile lacks the fields we need, so that a sweep
// over many files reports and moves on rather than aborting on the first.
bool has_var (const Vector<std::string>& names, const std::string& want)
{
    for (auto const& n : names) {
        if (n == want) { return true; }
    }
    return false;
}

void process (const std::string& pltfile, Real rho_cut)
{
    PlotFileData pf(pltfile);

    if (pf.spaceDim() != 3 || pf.coordSys() != 0) {
        amrex::Print() << "# " << pltfile << ": skipped, needs 3-d Cartesian\n";
        return;
    }

    const Vector<std::string>& names = pf.varNames();
    for (auto const& want : {"density", "B_x", "B_y", "B_z"}) {
        if (! has_var(names, want)) {
            amrex::Print() << "# " << pltfile << ": skipped, no " << want << '\n';
            return;
        }
    }

    const Array<Real,AMREX_SPACEDIM> problo = pf.probLo();
    const int fine_level = pf.finestLevel();

    Real e_tor = 0.0;
    Real e_pol = 0.0;
    // Peak field strengths as well as energies: the model manifest reports
    // Bt_over_Bp_amplitude, a ratio of maxima, and that is what it has to be
    // compared against. The energy ratio weights the whole volume and is a
    // different number.
    //
    // Reported in GAUSS. Castro's MHD state is in Heaviside-Lorentz units,
    // B' = B/sqrt(4 pi) -- see problem_initialize_mhd_data.H, which applies
    // the factor once on the way in. Reading the raw state as gauss
    // understates every field by 3.5449 and made the reconstructed peak look
    // like 27% of the model's max|B_phi| when it is really 99%. The energies
    // need no conversion: 0.5*B'^2 is already B^2/(8 pi) in erg/cm^3.
    constexpr Real castro_to_gauss = 3.5449077018110318_rt;  // sqrt(4 pi)
    Real b_tor_max = 0.0;
    Real b_pol_max = 0.0;
    Real b_max = 0.0;

    for (int ilev = 0; ilev <= fine_level; ++ilev) {

        const Array<Real,AMREX_SPACEDIM> dx = pf.cellSize(ilev);
        const Real dv = dx[0] * dx[1] * dx[2];

        const MultiFab& rho = pf.get(ilev, "density");
        const MultiFab& bxf = pf.get(ilev, "B_x");
        const MultiFab& byf = pf.get(ilev, "B_y");
        const MultiFab& bzf = pf.get(ilev, "B_z");

        auto const& ra = rho.const_arrays();
        auto const& xa = bxf.const_arrays();
        auto const& ya = byf.const_arrays();
        auto const& za = bzf.const_arrays();

        const Real xlo = problo[0];
        const Real ylo = problo[1];
        const Real dx0 = dx[0];
        const Real dx1 = dx[1];

        // A cell covered by a finer level must not be counted twice. The
        // finest level gets an all-zero mask so the reduction below stays a
        // single code path.
        iMultiFab mask;
        if (ilev < fine_level) {
            IntVect ratio{pf.refRatio(ilev)};
            mask = makeFineMask(pf.boxArray(ilev), pf.DistributionMap(ilev),
                                pf.boxArray(ilev+1), ratio);
        } else {
            mask.define(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, 0);
            mask.setVal(0);
        }
        auto const& ma = mask.const_arrays();

        auto rr = ParReduce(TypeList<ReduceOpSum,ReduceOpSum,ReduceOpMax,ReduceOpMax,ReduceOpMax>{},
                            TypeList<Real,Real,Real,Real,Real>{}, rho,
                    [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                        -> GpuTuple<Real,Real,Real,Real,Real>
                    {
                        if (ma[bno](i,j,k) != 0)       { return {0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt}; }
                        if (ra[bno](i,j,k) <= rho_cut) { return {0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt}; }

                        const Real x = xlo + (Real(i) + 0.5_rt) * dx0;
                        const Real y = ylo + (Real(j) + 0.5_rt) * dx1;
                        const Real varpi = std::sqrt(x*x + y*y);

                        const Real bx = xa[bno](i,j,k);
                        const Real by = ya[bno](i,j,k);
                        const Real bz = za[bno](i,j,k);

                        const Real btor = (varpi > 0.0_rt) ? (-y*bx + x*by) / varpi : 0.0_rt;
                        const Real b2 = bx*bx + by*by + bz*bz;
                        const Real bpol2 = amrex::max(b2 - btor*btor, 0.0_rt);

                        return {0.5_rt * btor * btor * dv, 0.5_rt * bpol2 * dv,
                                std::abs(btor), std::sqrt(bpol2), std::sqrt(b2)};
                    });

        e_tor += amrex::get<0>(rr);
        e_pol += amrex::get<1>(rr);
        b_tor_max = amrex::max(b_tor_max, castro_to_gauss * amrex::get<2>(rr));
        b_pol_max = amrex::max(b_pol_max, castro_to_gauss * amrex::get<3>(rr));
        b_max     = amrex::max(b_max,     castro_to_gauss * amrex::get<4>(rr));
    }

    ParallelDescriptor::ReduceRealSum(e_tor);
    ParallelDescriptor::ReduceRealSum(e_pol);
    ParallelDescriptor::ReduceRealMax(b_tor_max);
    ParallelDescriptor::ReduceRealMax(b_pol_max);
    ParallelDescriptor::ReduceRealMax(b_max);

    const Real e_mag = e_tor + e_pol;

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << std::setw(11) << std::fixed << std::setprecision(5) << pf.time()
                  << std::scientific << std::setprecision(5)
                  << std::setw(14) << e_tor
                  << std::setw(14) << e_pol
                  << std::setw(14) << (e_pol > 0.0 ? e_tor / e_pol
                                                   : std::numeric_limits<Real>::infinity())
                  << std::fixed << std::setprecision(7)
                  << std::setw(12) << (e_mag > 0.0 ? e_tor / e_mag : 0.0)
                  << std::scientific << std::setprecision(5)
                  << std::setw(14) << b_tor_max
                  << std::setw(14) << b_pol_max
                  << std::setw(14) << b_max
                  << std::setw(14) << (b_pol_max > 0.0 ? b_tor_max / b_pol_max
                                                       : std::numeric_limits<Real>::infinity())
                  // The ztwd EOS assumes an unquantised electron gas, which
                  // needs B below the critical field. The model starts at
                  // 0.73 B_c, so this column says when the run leaves the
                  // range its own equation of state is valid in.
                  << std::setw(13) << b_max / 4.4140e13
                  << '\n' << std::flush;
    }
}

} // namespace

void main_main ()
{
    const int narg = amrex::command_argument_count();

    Real rho_cut = 1.0e5;   // ambient is 2.0e4
    std::vector<std::string> plotfiles;

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "--rho-cut") {
            rho_cut = std::stod(amrex::get_command_argument(++farg));
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
            << " Toroidal/poloidal magnetic energy split, rebuilt from B_x, B_y, B_z\n"
            << " because the emag_density and etor_density derives are corrupt.\n"
            << "\n"
            << " Usage:\n"
            << "    fbtbp [--rho-cut RHO] plotfile [plotfile ...]\n"
            << "\n"
            << " args --rho-cut RHO : skip cells at or below this density (default 1.0e5)\n"
            << "\n";
        return;
    }

    amrex::Print() << "# rho_cut = " << rho_cut << '\n';
    amrex::Print() << "# energies in erg; peak fields in GAUSS (state is Heaviside-Lorentz,"
                      " converted by sqrt(4 pi)); B_c = 4.414e13 G\n";
    amrex::Print() << "#          t         E_tor         E_pol         Et/Ep     Et/Emag"
                   << "   Btor_max_G    Bpol_max_G       B_max_G     Bt/Bp_amp        B/B_c\n";

    for (auto const& pltfile : plotfiles) {
        process(pltfile, rho_cut);
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
