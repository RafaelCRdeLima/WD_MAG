// Evaluate the compiled-in EOS along a density profile and print P.
//
// Built twice, once per EOS_DIR, so the same densities go through the same
// code paths the runs use. Reads "rho T" pairs from the file named by
// eos_profile.infile (default profile.in) and writes
//     rho  T  P  e  cs  (and the EOS name)
// to stdout, one line per point.
//
// The question (DIARIO 11): the SCF balanced the star against P_ztwd(rho).
// The HZ campaign then evolved it against P_helmholtz(rho, 1e7 K). If those
// pressures differ by a percent, every cell starts out of hydrostatic
// balance by a percent, and the surface -- smallest scale height -- goes
// first. That is where the heating appeared, and no run has ever tested it.

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <AMReX_ParmParse.H>
#include <eos.H>
#include <extern_parameters.H>
#include <network.H>
#include <unit_test.H>

int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    {
        init_unit_test();
        eos_init(unit_test_rp::small_temp, unit_test_rp::small_dens);
        network_init();

        std::string infile = "profile.in";
        {
            amrex::ParmParse pp("eos_profile");
            pp.query("infile", infile);
        }

        std::ifstream f(infile);
        if (!f) {
            amrex::Abort("could not open " + infile);
        }

        // One species, mass fraction 1 -- mu2.net, abar = 4, zbar = 2.
        std::cout << "# EOS = " << EOS_NAME << ", NumSpec = "
                  << NumSpec << "\n";
        std::cout << "# rho  T  P  e  cs\n";
        std::cout << std::scientific << std::setprecision(10);

        amrex::Real rho, T;
        while (f >> rho >> T) {
            eos_t state;
            state.rho = rho;
            state.T = T;
            for (int n = 0; n < NumSpec; ++n) {
                state.xn[n] = (n == 0) ? 1.0_rt : 0.0_rt;
            }
            eos(eos_input_rt, state);
            std::cout << rho << "  " << T << "  " << state.p << "  "
                      << state.e << "  " << state.cs << "\n";
        }
    }
    amrex::Finalize();
}
