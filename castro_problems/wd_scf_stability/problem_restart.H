#ifndef problem_restart_H
#define problem_restart_H

#include <prob_parameters.H>

#include <string>
#include <scf_model.H>
#include <scf_setup.H>

// Castro calls problem_initialize() only on a cold start, never on restart,
// so anything derived from the model has to be recomputed here or the
// restarted run silently carries zeros for R_star and t_dyn.
AMREX_INLINE
void problem_restart (const std::string& /*checkpoint_dir*/)
{
    scf_model::read(problem::model_name);
    scf_model::compute_background_star();
}

#endif
