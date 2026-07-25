#ifndef problem_initialize_mhd_data_H
#define problem_initialize_mhd_data_H

#include <seed_field.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_mhd_data (int i, int j, int k,
                                  Array4<Real> const& B, const int ibdir,
                                  const GeometryData& geomdata)
{
    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    // Castro/AMReX MHD uses Heaviside-Lorentz units, B' = B/sqrt(4 pi)
    // (dashboard/units.py's gauss_to_castro -- same factor, same reason,
    // the R4 "sqrt(4 pi) source of error" rule). The curl below is
    // evaluated in Gaussian units (seed_field::A_component /
    // Ax_at/Ay_at/Az_at) to match the SCF's E_mag = integral B^2/(8 pi)
    // dV convention used to calibrate the amplitude in
    // seed_field::calibrate_and_diagnose(); this factor is applied ONLY
    // here, at the last step before the value is written into Castro's
    // own state.
    const Real gaussian_to_castro = 1.0_rt / std::sqrt(4.0_rt * M_PI);

    if (ibdir == 0) {
        Real Bx = (seed_field::Az_at(i, j + 1, k, dx, problo) - seed_field::Az_at(i, j, k, dx, problo)) / dx[1]
                - (seed_field::Ay_at(i, j, k + 1, dx, problo) - seed_field::Ay_at(i, j, k, dx, problo)) / dx[2];
        B(i,j,k) = Bx * gaussian_to_castro;
    } else if (ibdir == 1) {
        Real By = (seed_field::Ax_at(i, j, k + 1, dx, problo) - seed_field::Ax_at(i, j, k, dx, problo)) / dx[2]
                - (seed_field::Az_at(i + 1, j, k, dx, problo) - seed_field::Az_at(i, j, k, dx, problo)) / dx[0];
        B(i,j,k) = By * gaussian_to_castro;
    } else {
        Real Bz = (seed_field::Ay_at(i + 1, j, k, dx, problo) - seed_field::Ay_at(i, j, k, dx, problo)) / dx[0]
                - (seed_field::Ax_at(i, j + 1, k, dx, problo) - seed_field::Ax_at(i, j, k, dx, problo)) / dx[1];
        B(i,j,k) = Bz * gaussian_to_castro;
    }
}

#endif
