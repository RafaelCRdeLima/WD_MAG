#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <model_parser.H>
#include <seed_field.H>

AMREX_INLINE
void problem_initialize ()
{
    // 1D radial background density profile: the field-free, non-rotating
    // SCF equilibrium selected as the Braithwaite background star (plan
    // Step 1), written by scf/castro_model_writer.py. Read the same way
    // Castro's own StarGrav/toy_convect problems read an external model.
    read_model_file(problem::model_name);

    amrex::Real rho_c = model::profile(0).state(0, model::idens);
    amrex::Real surface_floor = rho_c * problem::surface_density_frac;

    int i_surface = model::npts - 1;
    for (int i = 0; i < model::npts; ++i) {
        if (model::profile(0).state(i, model::idens) < surface_floor) {
            i_surface = i;
            break;
        }
    }
    problem::R_star_cm = model::profile(0).r(i_surface);

    // |W| = G integral_0^Rstar [M(r)/r] dM, dM = 4 pi r^2 rho(r) dr --
    // exact for a spherically symmetric star (Gauss's law shell
    // reduction), computed directly from the 1D model via the
    // trapezoidal rule. Independent of Castro's own self-gravity module
    // (USE_GRAV=FALSE here -- not needed for a static IC check).
    amrex::Real M_enc = 0.0_rt;
    amrex::Real W = 0.0_rt;
    for (int i = 1; i <= i_surface; ++i) {
        amrex::Real r0 = model::profile(0).r(i - 1);
        amrex::Real r1 = model::profile(0).r(i);
        amrex::Real rho0 = model::profile(0).state(i - 1, model::idens);
        amrex::Real rho1 = model::profile(0).state(i, model::idens);
        amrex::Real dr = r1 - r0;

        amrex::Real M0 = M_enc;
        amrex::Real dM = 0.5_rt * (rho0 * r0 * r0 + rho1 * r1 * r1) * 4.0_rt * M_PI * dr;
        M_enc += dM;
        amrex::Real M1 = M_enc;

        amrex::Real integrand0 = M0 * rho0 * r0;
        amrex::Real integrand1 = M1 * rho1 * r1;
        W -= 4.0_rt * M_PI * C::Gconst * 0.5_rt * (integrand0 + integrand1) * dr;
    }
    problem::W_abs_erg = std::abs(W);

    amrex::Print() << "wd_braithwaite: R_star = " << problem::R_star_cm
                    << " cm, M_enclosed(R_star) = " << M_enc / 1.989e33
                    << " Msun, |W| = " << problem::W_abs_erg << " erg\n";

    seed_field::generate_modes();
    seed_field::calibrate_and_diagnose();
}

#endif
