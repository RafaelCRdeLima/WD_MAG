#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <cstdlib>

// Dumps a variable along a line of cells (fixed j,k, i ranging over
// [ilo,ihi]) for a quick sanity check of a gravity/field profile near
// a specific point -- built for the wd_braithwaite r=0 gravity-patch
// verification (docs/teoria.md Sec 6.7/6.8).
//
// usage: fline -v varname -j jidx -k kidx -ilo ilo -ihi ihi plotfile

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string var_name;
    int jidx = -1, kidx = -1, ilo = 0, ihi = -1;

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "-v") {
            var_name = amrex::get_command_argument(++farg);
        } else if (name == "-j") {
            jidx = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-k") {
            kidx = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-ilo") {
            ilo = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-ihi") {
            ihi = std::stoi(amrex::get_command_argument(++farg));
        } else {
            break;
        }
        ++farg;
    }

    if (farg > narg || var_name.empty() || jidx < 0 || kidx < 0 || ihi < ilo) {
        amrex::Print() << "usage: fline -v varname -j jidx -k kidx -ilo ilo -ihi ihi plotfile\n";
        return;
    }

    const std::string& filename = amrex::get_command_argument(farg);
    PlotFileData pf(filename);

    int ilev = 0;
    const MultiFab& mf = pf.get(ilev, var_name);
    auto dx = pf.cellSize(ilev);
    auto problo = pf.probLo();

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& arr = mf.const_array(mfi);
        for (int i = ilo; i <= ihi; ++i) {
            IntVect iv(AMREX_D_DECL(i, jidx, kidx));
            if (bx.contains(iv)) {
                Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5);
                amrex::AllPrint() << "i=" << i << " x=" << x
                                  << " " << var_name << "=" << arr(i, jidx, kidx) << "\n";
            }
        }
    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
