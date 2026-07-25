#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParReduce.H>
#include <AMReX_PlotFileUtil.H>
#include <algorithm>
#include <limits>
#include <cstdlib>

// Reports the extrema of a variable, computed ONLY over cells that are
// at least `margin` cells away from any box boundary (whether a real
// domain boundary or an interior box-to-box boundary from the level's
// own grid decomposition) -- i.e., cells whose stencil (up to +-margin)
// never reaches outside the box's own valid region into ghost data.
// Single level only (level 0) -- built for the wd_braithwaite
// investigation (docs/teoria.md Sec 6.6): isolating whether Div_B
// extrema come from the interior field (real physics) or specifically
// from box-edge cells (a suspected ghost-fill gap in Castro's own
// ca_derdivb, registered "the_same_box" while its stencil reads i+1).
//
// usage: finterior -v varname -m margin plotfile

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string var_name;
    int margin = 2;

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "-v" || name == "--variable") {
            var_name = amrex::get_command_argument(++farg);
        } else if (name == "-m" || name == "--margin") {
            margin = std::stoi(amrex::get_command_argument(++farg));
        } else {
            break;
        }
        ++farg;
    }

    if (farg > narg || var_name.empty()) {
        amrex::Print() << "\n"
                       << " Report extrema of a variable restricted to cells >= margin\n"
                       << " cells away from any box boundary (interior-only check).\n"
                       << " usage: finterior -v varname [-m margin] plotfile\n\n";
        return;
    }

    const std::string& filename = amrex::get_command_argument(farg);
    PlotFileData pf(filename);

    int ilev = 0;  // level 0 only

    const MultiFab& mf = pf.get(ilev, var_name);
    const BoxArray& ba = pf.boxArray(ilev);

    Real vmin_interior = std::numeric_limits<Real>::max();
    Real vmax_interior = std::numeric_limits<Real>::lowest();
    Real vmin_edge = std::numeric_limits<Real>::max();
    Real vmax_edge = std::numeric_limits<Real>::lowest();
    long n_interior = 0;
    long n_edge = 0;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Box interior_box = bx;
        interior_box.grow(-margin);
        auto const& arr = mf.const_array(mfi);

        amrex::LoopOnCpu(bx, [&] (int i, int j, int k)
        {
            Real val = arr(i, j, k);
            if (interior_box.contains(IntVect(AMREX_D_DECL(i, j, k)))) {
                vmin_interior = std::min(vmin_interior, val);
                vmax_interior = std::max(vmax_interior, val);
                n_interior++;
            } else {
                vmin_edge = std::min(vmin_edge, val);
                vmax_edge = std::max(vmax_edge, val);
                n_edge++;
            }
        });
    }

    amrex::Print() << " plotfile = " << filename << "\n"
                   << " variable = " << var_name << " margin = " << margin << " cells\n"
                   << " time = " << std::setprecision(17) << pf.time() << "\n"
                   << " INTERIOR (n=" << n_interior << "): min = " << vmin_interior
                   << "  max = " << vmax_interior << "\n"
                   << " EDGE     (n=" << n_edge << "): min = " << vmin_edge
                   << "  max = " << vmax_edge << "\n";
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
