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
    //
    // Computes R_star, |W|, t_dyn, and the global-damping window (see
    // seed_field::compute_background_star_and_damping() -- shared with
    // problem_restart.H, since Castro only calls problem_initialize()
    // on a cold start, never on restart).
    seed_field::compute_background_star_and_damping();

    seed_field::generate_modes();
    seed_field::calibrate_and_diagnose();
}

#endif
