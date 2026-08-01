#ifndef SCF_MODEL_H
#define SCF_MODEL_H

#include <AMReX_REAL.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

#include <scf_model_data.H>
#include <prob_parameters.H>

#include <cmath>

namespace scf_model
{

// Bilinear interpolation on the uniform (varpi, z) table. Indices are
// computed arithmetically, which is why read() insists the axes are uniform.
// Points outside the table are CLAMPED to the boundary value, not zeroed:
// zeroing puts a step at the table edge, and the discrete curl of a step is
// a spurious field. That failure mode is not hypothetical -- the Python
// exporter hit it when the meridional grid was sized to the face of the
// Cartesian box instead of its corners, and produced a reconstructed field
// four times the star's. init_check_covers_domain() below makes sure the
// clamp is never actually needed.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real table_at(const amrex::Real* tab, amrex::Real varpi, amrex::Real z)
{
    amrex::Real fi = (varpi - varpi_lo) / dvarpi;
    amrex::Real fj = (z - z_lo) / dz;

    fi = amrex::max(0.0_rt, amrex::min(fi, static_cast<amrex::Real>(n_varpi - 1)));
    fj = amrex::max(0.0_rt, amrex::min(fj, static_cast<amrex::Real>(n_z - 1)));

    int i0 = amrex::min(static_cast<int>(fi), n_varpi - 2);
    int j0 = amrex::min(static_cast<int>(fj), n_z - 2);
    amrex::Real ti = fi - static_cast<amrex::Real>(i0);
    amrex::Real tj = fj - static_cast<amrex::Real>(j0);

    const amrex::Real v00 = tab[i0 * n_z + j0];
    const amrex::Real v01 = tab[i0 * n_z + (j0 + 1)];
    const amrex::Real v10 = tab[(i0 + 1) * n_z + j0];
    const amrex::Real v11 = tab[(i0 + 1) * n_z + (j0 + 1)];

    return (1.0_rt - ti) * ((1.0_rt - tj) * v00 + tj * v01)
         + ti * ((1.0_rt - tj) * v10 + tj * v11);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real density_at(amrex::Real x, amrex::Real y, amrex::Real z)
{
    return table_at(dens, std::sqrt(x * x + y * y), z);
}

// Azimuthal velocity of a v2 (rotating) model, as Cartesian components.
// v = v_phi e_phi with e_phi = (-y, x, 0)/varpi, so v_x = -v_phi y/varpi
// and v_y = +v_phi x/varpi. A v1 model has no velocity table and returns
// zero, which is the non-rotating case and not a fallback.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real velocity_at(int comp, amrex::Real x, amrex::Real y, amrex::Real z)
{
    if (has_velocity == 0 || comp == 2) {
        return 0.0_rt;
    }
    amrex::Real varpi = std::sqrt(x * x + y * y);
    // v_phi vanishes on the axis as varpi (Omega finite), so the Cartesian
    // components go to zero there; this is the limit, not a guard.
    if (varpi <= 0.0_rt) {
        return 0.0_rt;
    }
    amrex::Real vphi = table_at(v_phi, varpi, z);
    return (comp == 0) ? -vphi * y / varpi : vphi * x / varpi;
}

// Cartesian components of A from the axisymmetric pair (A_phi, A_z):
//   A = A_phi e_phi + A_z e_z,   e_phi = (-y, x, 0)/varpi
// A_varpi is zero by construction of the potential (see the writer).
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real A_component(int comp, amrex::Real x, amrex::Real y, amrex::Real z)
{
    amrex::Real varpi = std::sqrt(x * x + y * y);

    if (comp == 2) {
        return problem::field_scale * table_at(a_z, varpi, z);
    }

    // On the axis e_phi is undefined, but A_phi vanishes there as varpi
    // (A_phi = u/varpi with u ~ varpi^2), so the Cartesian components go to
    // zero. Returning zero is the limit, not a patch over a singularity.
    if (varpi <= 0.0_rt) {
        return 0.0_rt;
    }

    const amrex::Real ap = problem::field_scale * table_at(a_phi, varpi, z);
    return (comp == 0) ? (-ap * y / varpi) : (ap * x / varpi);
}

// Edge-centred evaluation, matching the convention of Castro's own
// Exec/mhd_tests/LoopAdvection and of the wd_braithwaite problem: each
// component is node-centred in the two directions it is NOT aligned with and
// cell-centred in the one it IS aligned with. The discrete curl of A sampled
// this way lands on faces and is divergence-free identically.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real Ax_at(int i, int j, int k, const amrex::Real* dx, const amrex::Real* problo)
{
    amrex::Real x = problo[0] + dx[0] * (static_cast<amrex::Real>(i) + 0.5_rt);
    amrex::Real y = problo[1] + dx[1] * static_cast<amrex::Real>(j);
    amrex::Real z = problo[2] + dx[2] * static_cast<amrex::Real>(k);
    return A_component(0, x, y, z);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real Ay_at(int i, int j, int k, const amrex::Real* dx, const amrex::Real* problo)
{
    amrex::Real x = problo[0] + dx[0] * static_cast<amrex::Real>(i);
    amrex::Real y = problo[1] + dx[1] * (static_cast<amrex::Real>(j) + 0.5_rt);
    amrex::Real z = problo[2] + dx[2] * static_cast<amrex::Real>(k);
    return A_component(1, x, y, z);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real Az_at(int i, int j, int k, const amrex::Real* dx, const amrex::Real* problo)
{
    amrex::Real x = problo[0] + dx[0] * static_cast<amrex::Real>(i);
    amrex::Real y = problo[1] + dx[1] * static_cast<amrex::Real>(j);
    amrex::Real z = problo[2] + dx[2] * (static_cast<amrex::Real>(k) + 0.5_rt);
    return A_component(2, x, y, z);
}

// B = curl A in GAUSSIAN units. The Heaviside-Lorentz conversion Castro's
// state expects (B' = B / sqrt(4 pi)) is applied once, at the last step
// before a value is written into the state array, in
// problem_initialize_mhd_data.H -- not here, so that the diagnostic
// quadratures below share the SCF's own E_mag = int B^2/(8 pi) dV convention.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void B_gaussian_at(int i, int j, int k, const amrex::Real* dx, const amrex::Real* problo,
                   amrex::Real& Bx, amrex::Real& By, amrex::Real& Bz)
{
    Bx = (Az_at(i, j + 1, k, dx, problo) - Az_at(i, j, k, dx, problo)) / dx[1]
       - (Ay_at(i, j, k + 1, dx, problo) - Ay_at(i, j, k, dx, problo)) / dx[2];

    By = (Ax_at(i, j, k + 1, dx, problo) - Ax_at(i, j, k, dx, problo)) / dx[2]
       - (Az_at(i + 1, j, k, dx, problo) - Az_at(i, j, k, dx, problo)) / dx[0];

    Bz = (Ay_at(i + 1, j, k, dx, problo) - Ay_at(i, j, k, dx, problo)) / dx[0]
       - (Ax_at(i, j + 1, k, dx, problo) - Ax_at(i, j, k, dx, problo)) / dx[1];
}

// Host-side, called once from problem_initialize(). The table must reach the
// CORNERS of the domain, at sqrt(3) times its half-width for a cube, or the
// clamp in table_at() starts standing in for real data and the curl of the
// resulting plateau is wrong. Checked rather than assumed.
AMREX_INLINE
void init_check_covers_domain()
{
    amrex::ParmParse pp("geometry");
    amrex::Vector<amrex::Real> lo, hi;
    pp.getarr("prob_lo", lo);
    pp.getarr("prob_hi", hi);

    amrex::Real corner = 0.0;
    amrex::Real z_min = 0.0, z_max = 0.0;
    for (int d = 0; d < 3; ++d) {
        amrex::Real m = amrex::max(std::abs(lo[d]), std::abs(hi[d]));
        corner += m * m;
        if (d == 2) { z_min = lo[d]; z_max = hi[d]; }
    }
    corner = std::sqrt(corner);

    const amrex::Real varpi_corner =
        std::sqrt(amrex::max(std::abs(lo[0]), std::abs(hi[0]))
                  * amrex::max(std::abs(lo[0]), std::abs(hi[0]))
                + amrex::max(std::abs(lo[1]), std::abs(hi[1]))
                  * amrex::max(std::abs(lo[1]), std::abs(hi[1])));

    amrex::Print() << "  [scf_model] domain corner at r = " << corner
                   << " cm, needs varpi <= " << varpi_hi
                   << " and z in [" << z_lo << ", " << z_hi << "]\n";

    if (varpi_corner > varpi_hi * (1.0 + 1.0e-10)) {
        amrex::Abort("scf_model: the model table does not reach the domain's "
                     "varpi corner -- regenerate it with a larger grid "
                     "(see investigations/export_phase1_models.py)");
    }
    if (z_min < z_lo * (1.0 + 1.0e-10) || z_max > z_hi * (1.0 + 1.0e-10)) {
        amrex::Abort("scf_model: the model table does not span the domain in "
                     "z -- regenerate it with a larger grid");
    }
}

}  // namespace scf_model

#endif
