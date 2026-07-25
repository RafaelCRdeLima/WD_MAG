#ifndef SEED_FIELD_DATA_H
#define SEED_FIELD_DATA_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

// Fixed, GPU-managed storage for the random multi-scale mode set used by
// the Braithwaite seed-field generator (Step 2 of the plan). Mirrors the
// idiom in Castro's own Util/model_parser/model_parser_data.H: a plain
// extern AMREX_GPU_MANAGED array, sized at compile time, filled once on
// the host (seed_field::generate_modes(), called from
// problem_initialize()) and read on host or device thereafter
// (problem_initialize_mhd_data.H runs per grid point, possibly on GPU).
//
// The modes MUST be generated once and shared across every cell/rank/
// tile -- unlike amrex::RandomEngine-based per-cell perturbations used
// elsewhere in Castro, a spatially-coherent multi-scale field cannot be
// regenerated independently per grid point.
namespace seed_modes
{
    constexpr int MAX_MODES = 64;  // headroom above the dashboard's n_modes slider max (50)
    constexpr int NCOMP = 3;       // A_x, A_y, A_z

    extern AMREX_GPU_MANAGED int n_modes_actual;

    // indexed [comp][mode]
    extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 0, NCOMP - 1, 0, MAX_MODES - 1> kx;
    extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 0, NCOMP - 1, 0, MAX_MODES - 1> ky;
    extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 0, NCOMP - 1, 0, MAX_MODES - 1> kz;
    extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 0, NCOMP - 1, 0, MAX_MODES - 1> phase;
    extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 0, NCOMP - 1, 0, MAX_MODES - 1> amp;
}

#endif
