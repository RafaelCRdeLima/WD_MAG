#ifndef SCF_MODEL_DATA_H
#define SCF_MODEL_DATA_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>
#include <string>

// Storage for the 2D axisymmetric magnetized SCF model written by
// scf/axisym_model_writer.py (format "axisym-vector-potential-v1").
//
// Why 2D and not Castro's own model_parser. model_parser reads a 1D radial
// profile and interpolates on |r|. That is right for a spherical star and
// wrong for this one: the configuration is prolate (R_pol/R_eq = 1.10) and
// its field is a function of (varpi, z), not of r. Reducing it to 1D is the
// operation that destroys what is being tested here.
//
// Why the vector potential and not B. Interpolating B onto faces leaves a
// divergence that constrained transport carries forever. The table stores
// A_phi and A_z; the problem takes the DISCRETE curl on cell edges, so
// div B = 0 holds as an identity of the discretisation. Verified in Python
// before any of this was written: investigations/vector_potential_export.py
// and the gates in investigations/export_phase1_models.py both measure
// |div B| ~ 2e-16 normalised.
//
// The grid is required to be uniform in both varpi and z, which the writer
// guarantees (np.linspace) and read() verifies rather than assumes -- that
// lets the interpolation index arithmetically instead of searching.
namespace scf_model
{
    extern AMREX_GPU_MANAGED int n_varpi;
    extern AMREX_GPU_MANAGED int n_z;

    extern AMREX_GPU_MANAGED amrex::Real varpi_lo;
    extern AMREX_GPU_MANAGED amrex::Real varpi_hi;
    extern AMREX_GPU_MANAGED amrex::Real z_lo;
    extern AMREX_GPU_MANAGED amrex::Real z_hi;
    extern AMREX_GPU_MANAGED amrex::Real dvarpi;
    extern AMREX_GPU_MANAGED amrex::Real dz;

    // flattened [i_varpi * n_z + i_z], managed so device code can read them
    extern AMREX_GPU_MANAGED amrex::Real* dens;
    extern AMREX_GPU_MANAGED amrex::Real* a_phi;
    extern AMREX_GPU_MANAGED amrex::Real* a_z;
    // v2 models only; a null pointer means the model is v1 and
    // the star is not rotating.
    extern AMREX_GPU_MANAGED amrex::Real* v_phi;
    extern AMREX_GPU_MANAGED int has_velocity;

    // Peak values, filled by read(). Reported at initialisation so a run
    // whose model file is not the one intended is caught in the log rather
    // than in the physics.
    extern AMREX_GPU_MANAGED amrex::Real dens_max;

    // Host-side. Aborts on a malformed file, on a non-uniform grid, or on a
    // size that would not fit; never returns a partially filled table.
    void read(const std::string& filename);

    void deallocate();
}

#endif
