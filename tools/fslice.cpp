// Dump one 2-d plane of a Castro MHD plotfile as text.
//
// The point is to get the evolved field off the cluster. A 192^3 plotfile is
// 1.2 GB and there are 70 of them; a single 192^2 plane with four fields is
// about 600 kB, which fits in the repository and can be plotted on a laptop.
//
// Why this and not the existing figures: field_meridional.pdf and
// field_lines_3d.pdf are drawn from the Grad-Shafranov solver, where the
// configuration is axisymmetric and a poloidal field line IS exactly a
// contour of the flux function. The evolved field has no flux function --
// from t ~ 1.3 s the m = 1 kink grows and axisymmetry is gone, which is the
// result being looked for -- so the lines have to come from B itself.
//
// Two planes matter:
//   --normal y   the meridional cut. At y = 0 and x > 0, B_y IS B_phi and
//                (B_x, B_z) is the in-plane poloidal field, so no rotation
//                of components is needed to read it.
//   --normal z   the equatorial cut, where an m = 1 displacement of the
//                toroidal column is directly visible.
//
// The half-shift geometry of these inputs puts a cell CENTRE exactly on each
// axis, so a cut at 0 lands on cell centres rather than between them.
//
// Fields are written in GAUSS: Castro's MHD state is Heaviside-Lorentz,
// B' = B/sqrt(4 pi) (problem_initialize_mhd_data.H), and reading it raw
// understates every field by 3.5449.
//
// Usage:
//     fslice [--normal x|y|z] [--at COORD] [-o FILE] plotfile
//
// Writes to FILE, defaulting to <plotfile>_<axis><at>.txt, rather than to
// stdout: AMReX prints its own banner there and it would land in the middle
// of the table.
//
// The file carries a commented header and then rows of
//     coord1  coord2  density  B_x  B_y  B_z

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace amrex;

namespace {

constexpr Real castro_to_gauss = 3.5449077018110318_rt;   // sqrt(4 pi)

const char* const AXIS = "xyz";

bool has_var (const Vector<std::string>& names, const std::string& want)
{
    for (auto const& n : names) {
        if (n == want) { return true; }
    }
    return false;
}

} // namespace

void main_main ()
{
    const int narg = amrex::command_argument_count();

    int normal = 1;          // y: the meridional cut
    Real at = 0.0;
    std::string pltfile;
    std::string outfile;

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "--normal") {
            const std::string a = amrex::get_command_argument(++farg);
            if      (a == "x") { normal = 0; }
            else if (a == "y") { normal = 1; }
            else if (a == "z") { normal = 2; }
            else { amrex::Abort("--normal must be x, y or z"); }
        } else if (name == "--at") {
            at = std::stod(amrex::get_command_argument(++farg));
        } else if (name == "-o" || name == "--out") {
            outfile = amrex::get_command_argument(++farg);
        } else {
            break;
        }
        ++farg;
    }
    if (farg <= narg) { pltfile = amrex::get_command_argument(farg); }

    if (pltfile.empty()) {
        amrex::Print()
            << "\n"
            << " Dump one 2-d plane of a Castro MHD plotfile as text.\n"
            << "\n"
            << " Usage:\n"
            << "    fslice [--normal x|y|z] [--at COORD] [-o FILE] plotfile\n"
            << "\n"
            << " args --normal x|y|z : plane normal (default y, the meridional cut)\n"
            << " args --at COORD     : where along that normal, in cm (default 0)\n"
            << " args -o|--out FILE  : output file (default <plotfile>_<axis><at>.txt)\n"
            << "\n";
        return;
    }

    PlotFileData pf(pltfile);

    if (pf.spaceDim() != 3 || pf.coordSys() != 0) {
        amrex::Abort("fslice needs a 3-d Cartesian plotfile");
    }

    const Vector<std::string> vars = {"density", "B_x", "B_y", "B_z"};
    for (auto const& v : vars) {
        if (! has_var(pf.varNames(), v)) { amrex::Abort("plotfile has no " + v); }
    }

    // Finest level only. These runs are single-level -- Castro's constrained
    // transport does not do AMR -- and silently flattening a refined grid
    // would put cells of two sizes in one table.
    const int lev = pf.finestLevel();
    if (lev != 0) {
        amrex::Print() << "# WARNING: taking level " << lev << " only\n";
    }

    const Array<Real,AMREX_SPACEDIM> problo = pf.probLo();
    const Array<Real,AMREX_SPACEDIM> dx = pf.cellSize(lev);
    const Box domain = pf.probDomain(lev);

    const int islice = static_cast<int>(std::floor((at - problo[normal]) / dx[normal]));
    if (islice < domain.smallEnd(normal) || islice > domain.bigEnd(normal)) {
        amrex::Abort("--at is outside the domain");
    }
    const Real at_actual = problo[normal] + (Real(islice) + 0.5_rt) * dx[normal];

    // The two in-plane axes, in ascending order, so the output is always
    // (x,y), (x,z) or (y,z) and never depends on which normal was asked for.
    const int a0 = (normal == 0) ? 1 : 0;
    const int a1 = (normal == 2) ? 1 : 2;

    const int n0 = domain.length(a0);
    const int n1 = domain.length(a1);
    const int nv = static_cast<int>(vars.size());

    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    std::vector<Real> buf(std::size_t(n0) * n1 * nv, nan);

    for (int v = 0; v < nv; ++v) {
        const MultiFab& mf = pf.get(lev, vars[v]);
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            if (islice < bx.smallEnd(normal) || islice > bx.bigEnd(normal)) { continue; }
            auto const& a = mf.const_array(mfi);

            IntVect lo = bx.smallEnd();
            IntVect hi = bx.bigEnd();
            lo[normal] = islice;
            hi[normal] = islice;

            for (int j = lo[a1]; j <= hi[a1]; ++j) {
                for (int i = lo[a0]; i <= hi[a0]; ++i) {
                    IntVect iv;
                    iv[normal] = islice;
                    iv[a0] = i;
                    iv[a1] = j;
                    const std::size_t k =
                        (std::size_t(j - domain.smallEnd(a1)) * n0
                         + (i - domain.smallEnd(a0))) * nv + v;
                    buf[k] = a(iv[0], iv[1], iv[2]);
                }
            }
        }
    }

    if (outfile.empty()) {
        outfile = pltfile + "_" + AXIS[normal] + std::to_string(static_cast<long>(at)) + ".txt";
    }
    std::ofstream out(outfile);
    if (! out) { amrex::Abort("cannot write " + outfile); }

    out << "# plotfile = " << pltfile << "\n"
        << "# time = " << std::setprecision(10) << pf.time() << " s\n"
        << "# plane: " << AXIS[normal] << " = " << at_actual
        << " cm (requested " << at << ", cell " << islice << ")\n"
        << "# grid: " << n0 << " x " << n1
        << ", d" << AXIS[a0] << " = " << dx[a0]
        << ", d" << AXIS[a1] << " = " << dx[a1] << " cm\n"
        << "# fields in GAUSS (Heaviside-Lorentz state x sqrt(4 pi));"
           " density in g/cm^3\n"
        << "# " << AXIS[a0] << " " << AXIS[a1] << " density B_x B_y B_z\n";

    out << std::scientific << std::setprecision(7);
    for (int j = 0; j < n1; ++j) {
        const Real c1 = problo[a1] + (Real(j + domain.smallEnd(a1)) + 0.5_rt) * dx[a1];
        for (int i = 0; i < n0; ++i) {
            const Real c0 = problo[a0] + (Real(i + domain.smallEnd(a0)) + 0.5_rt) * dx[a0];
            const std::size_t k = (std::size_t(j) * n0 + i) * nv;
            out << c0 << ' ' << c1
                << ' ' << buf[k + 0]
                << ' ' << castro_to_gauss * buf[k + 1]
                << ' ' << castro_to_gauss * buf[k + 2]
                << ' ' << castro_to_gauss * buf[k + 3] << '\n';
        }
    }
    amrex::Print() << "wrote " << outfile << " (" << n0 << " x " << n1 << ")\n";
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
