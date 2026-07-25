#include <seed_field_data.H>

namespace seed_modes
{
    AMREX_GPU_MANAGED int n_modes_actual;

    AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 0, NCOMP - 1, 0, MAX_MODES - 1> kx;
    AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 0, NCOMP - 1, 0, MAX_MODES - 1> ky;
    AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 0, NCOMP - 1, 0, MAX_MODES - 1> kz;
    AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 0, NCOMP - 1, 0, MAX_MODES - 1> phase;
    AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 0, NCOMP - 1, 0, MAX_MODES - 1> amp;
}
