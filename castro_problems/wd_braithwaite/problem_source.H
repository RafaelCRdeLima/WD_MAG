#ifndef problem_source_H
#define problem_source_H

// Global velocity damping for the initial accommodation phase (see
// docs/teoria.md Sec 6.4): the Step 2 IC is hydrostatic against gas
// pressure + gravity only -- the seed field's magnetic pressure was
// never in that balance. Releasing the field cold (no damping) lets
// the resulting transient (observed: rho_c rising, accelerating, not
// turning over like the temperature spike did) run unchecked. This
// damps velocity EVERYWHERE (unlike Castro's built-in sponge, which
// Castro_sponge.cpp confirms only targets the low-density exterior via
// density/radius/pressure thresholds) for the first
// damping_end_time_s, smoothly ramped off over the preceding
// damping_ramp_start_s->damping_end_time_s window, so the transient
// radiates and the star finds magnetohydrostatic equilibrium before the
// actual Braithwaite m=1 reorganization is allowed to run freely.

#include <prob_parameters.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_source (int i, int j, int k,
                     GeometryData const& geomdata,
                     Array4<const Real> const& state,
                     Array4<Real> const& src,
                     const Real dt, const Real time)
{
    amrex::ignore_unused(i);
    amrex::ignore_unused(j);
    amrex::ignore_unused(k);
    amrex::ignore_unused(geomdata);
    amrex::ignore_unused(dt);

    if (time >= problem::damping_end_time_s) {
        return;
    }

    Real ramp = 1.0_rt;
    if (time > problem::damping_ramp_start_s) {
        Real x = (time - problem::damping_ramp_start_s) /
                 (problem::damping_end_time_s - problem::damping_ramp_start_s);
        ramp = 0.5_rt * (1.0_rt + std::cos(M_PI * x));  // 1 -> 0, smooth cutoff
    }

    Real fac = -problem::damping_rate * ramp;

    Real rho = state(i,j,k,URHO);
    Real rhoInv = 1.0_rt / rho;

    Real Sr[3];
    for (int n = 0; n < 3; n++) {
        Sr[n] = fac * state(i,j,k,UMX+n);
        src(i,j,k,UMX+n) = Sr[n];
    }

    // Same kinetic-energy-consistent treatment as Castro_sponge.cpp's
    // own momentum source (see apply_sponge): d(KE)/dt = u . d(rho u)/dt
    // when the source has no rho contribution.
    Real SrE = 0.0_rt;
    for (int n = 0; n < 3; n++) {
        SrE += state(i,j,k,UMX+n) * rhoInv * Sr[n];
    }
    src(i,j,k,UEDEN) = SrE;
}

#endif
