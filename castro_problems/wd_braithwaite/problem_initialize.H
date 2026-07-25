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

    int i_surface = -1;
    for (int i = 0; i < model::npts; ++i) {
        if (model::profile(0).state(i, model::idens) < surface_floor) {
            i_surface = i;
            break;
        }
    }
    if (i_surface < 0) {
        // Never found a point below threshold -- surface_density_frac
        // is at or below the model file's own density floor, so this
        // would silently fall back to the last model row (the domain
        // edge, not the stellar surface) instead of erroring. Caught
        // exactly this way once already (see inputs.evolve's FLUFF
        // HANDLING / surface_density_frac comments) -- fail loudly
        // instead of a second time.
        amrex::Error("wd_braithwaite: surface_density_frac never crossed -- "
                     "it must be set above the model file's own density "
                     "floor (density_floor_fraction in "
                     "scf/castro_model_writer.py), or R_star silently "
                     "becomes the domain edge instead of the stellar surface.");
    }
    problem::R_star_cm = model::profile(0).r(i_surface);

    // |W| = G integral_0^Rstar [M(r)/r] dM, dM = 4 pi r^2 rho(r) dr --
    // exact for a spherically symmetric star (Gauss's law shell
    // reduction), computed directly from the 1D model via the
    // trapezoidal rule. Independent of Castro's own self-gravity module
    // (USE_GRAV=TRUE as of Step 3 -- Castro's Poisson/Monopole solve
    // handles the actual dynamical gravity; this 1D integral is only
    // ever used for the E_mag/|W| calibration and diagnostics below).
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

    // t_dyn = sqrt(R_star^3 / (G M)) -- same formula as
    // dashboard/units.py's dynamical_time(), so the damping window
    // below is sized from the actual star, not a guessed constant.
    // See docs/teoria.md Sec 6.4: the Step 2 IC is hydrostatic against
    // gas pressure + gravity only (the seed field's magnetic pressure
    // was never in that balance) -- releasing it cold let rho_c climb,
    // accelerating, with no sign of turning over (unlike the
    // temperature transient, which visibly peaked and decayed). Global
    // velocity damping (problem_source.H, ALL cells, not just the
    // sponge's low-density exterior) for the first few t_dyn lets the
    // transient radiate and the star find magnetohydrostatic
    // equilibrium before the actual Braithwaite reorganization runs
    // freely -- exactly the accommodation phase the original plan
    // called for and the first attempt skipped.
    problem::t_dyn_s = std::sqrt(std::pow(problem::R_star_cm, 3) /
                                  (C::Gconst * M_enc));
    problem::damping_rate = 1.0_rt / (problem::damping_timescale_in_tdyn * problem::t_dyn_s);
    problem::damping_end_time_s = problem::damping_end_time_in_tdyn * problem::t_dyn_s;
    problem::damping_ramp_start_s = problem::damping_end_time_s -
        problem::damping_ramp_in_tdyn * problem::t_dyn_s;

    amrex::Print() << "wd_braithwaite: R_star = " << problem::R_star_cm
                    << " cm, M_enclosed(R_star) = " << M_enc / 1.989e33
                    << " Msun, |W| = " << problem::W_abs_erg << " erg\n"
                    << "wd_braithwaite: t_dyn = " << problem::t_dyn_s
                    << " s, global damping active [0, " << problem::damping_end_time_s
                    << " s], ramping off from " << problem::damping_ramp_start_s << " s\n";

    seed_field::generate_modes();
    seed_field::calibrate_and_diagnose();
}

#endif
