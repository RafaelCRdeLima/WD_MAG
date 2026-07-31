#ifndef problem_initialize_mhd_data_H
#define problem_initialize_mhd_data_H

#include <scf_model.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_mhd_data (int i, int j, int k,
                                  Array4<Real> const& B, const int ibdir,
                                  const GeometryData& geomdata)
{
    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    // Castro/AMReX MHD uses Heaviside-Lorentz units, B' = B/sqrt(4 pi). The
    // curl is taken in Gaussian units, matching the SCF's
    // E_mag = int B^2/(8 pi) dV convention, and the factor is applied ONLY
    // here, at the last step before the value enters Castro's state.
    const Real gaussian_to_castro = 1.0_rt / std::sqrt(4.0_rt * M_PI);

    Real Bx, By, Bz;
    scf_model::B_gaussian_at(i, j, k, dx, problo, Bx, By, Bz);

    if (ibdir == 0) {
        B(i,j,k) = Bx * gaussian_to_castro;
    } else if (ibdir == 1) {
        B(i,j,k) = By * gaussian_to_castro;
    } else {
        B(i,j,k) = Bz * gaussian_to_castro;
    }
}

#endif
