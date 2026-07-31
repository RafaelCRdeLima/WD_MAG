#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <scf_model.H>
#include <scf_setup.H>

AMREX_INLINE
void problem_initialize ()
{
    // The magnetized, axisymmetric 2 Msun equilibrium written by
    // investigations/export_phase1_models.py: density and the vector
    // potential (A_phi, A_z) on a meridional grid. Phase 1 of
    // docs/PLAN_ULTRAMASSIVE_SNIA.md evolves this to ask whether a
    // toroidal-dominated field survives an Alfven time.
    scf_model::read(problem::model_name);
    scf_model::init_check_covers_domain();
    scf_model::compute_background_star();
    scf_model::diagnose_field();
}

#endif
