#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <scf_model.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{
    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);
    Real y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
    Real z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);

    // Density from the 2D table, NOT from a 1D radial profile: this star is
    // prolate (R_pol/R_eq = 1.10) and interpolating on |r| would round it.
    Real rho = scf_model::density_at(x, y, z);

    // Ambient floor, which has to stay above castro.small_dens (1.0e4 for
    // these runs, docs/teoria.md Sec 6.5) or Castro aborts at init.
    state(i,j,k,URHO) = amrex::max(rho, problem::ambient_density);

    // With EOS_DIR = ztwd, the same zero-temperature degenerate equation of
    // state the SCF solved with, the density profile IS an equilibrium under
    // it -- that is what makes this a test of the FIELD's stability rather
    // than of a structural transient, and temperature is carried but
    // thermodynamically inert, so the default 1.0e5 is a placeholder.
    //
    // With EOS_DIR = helmholtz it is NOT inert. Set problem.initial_temperature
    // from the inputs file to a value at which the ion and radiation pressures
    // are still negligible against the degenerate one, or the star is balanced
    // against one EOS and evolved under another. See _prob_params.
    state(i,j,k,UTEMP) = problem::initial_temperature;
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = (n == 0) ? 1.0_rt : 0.0_rt;
    }

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n);
    }
    eos(eos_input_rt, eos_state);

    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,URHO) * eos_state.e;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
    }

    // Rotation, when the model carries it. Castro's own rotation support is
    // a rotating frame at constant Omega, which a j-constant law is not, so
    // differential rotation has to enter as an initial velocity field and
    // then evolve. That it evolves is the point: a differentially rotating
    // magnetised star is MRI-unstable, and whether the shear that supplies
    // the mass survives its own instability is what this run measures.
    //
    // No perturbation is added on top. The Tayler and MRI modes have to grow
    // from the discretisation's own asymmetry, which is the honest test:
    // seeding them would measure the seed's growth, not the configuration's.
    Real vx = scf_model::velocity_at(0, x, y, z);
    Real vy = scf_model::velocity_at(1, x, y, z);

    state(i,j,k,UMX) = state(i,j,k,URHO) * vx;
    state(i,j,k,UMY) = state(i,j,k,URHO) * vy;
    state(i,j,k,UMZ) = 0.0_rt;

    // UEDEN carries the kinetic energy; UEINT must not.
    state(i,j,k,UEDEN) += 0.5_rt * state(i,j,k,URHO) * (vx * vx + vy * vy);
}

#endif
