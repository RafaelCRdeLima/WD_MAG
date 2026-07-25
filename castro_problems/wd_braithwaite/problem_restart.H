#ifndef problem_restart_H
#define problem_restart_H

#include <prob_parameters.H>
#include <seed_field.H>

// Castro calls problem_initialize() only on a cold start
// (Castro::initData()), never on restart (Castro::restart() reads
// state directly from the checkpoint). Our damping window and R_star/
// |W|/t_dyn are computed once in problem_initialize() and stored in
// problem:: -- without also recomputing them here, they would sit at
// their _prob_params defaults (0.0) after any restart, silently
// disabling problem_source.H's global damping (damping_rate=0.0)
// instead of erroring. Found while adding restart support (see
// docs/teoria.md Sec 6.5).
//
// Deliberately does NOT call seed_field::generate_modes() or
// calibrate_and_diagnose() -- the B field itself is restored from the
// checkpoint (already evolved/relaxed to whatever point the run
// reached), and regenerating or recalibrating the seed here would
// silently discard that and substitute a fresh, uncalibrated field.
AMREX_INLINE
void problem_restart (const std::string& checkpoint_dir)
{
    amrex::ignore_unused(checkpoint_dir);

    seed_field::compute_background_star_and_damping();
}

#endif
