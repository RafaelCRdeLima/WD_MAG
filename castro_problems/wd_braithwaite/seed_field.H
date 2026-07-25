#ifndef SEED_FIELD_H
#define SEED_FIELD_H

// Braithwaite Step 2: random multi-scale vector potential, confined to
// the star, edge-centered so the discrete curl gives div B = 0 to
// machine precision. This step delivers a STATIC initial condition
// only -- fluff treatment (v_A limiter, sponge) is Step 3's job, not
// this one (see calibrate_and_diagnose()'s printed banner).

#include <cmath>
#include <random>

#include <AMReX_REAL.H>
#include <AMReX_Print.H>

#include <prob_parameters.H>
#include <seed_field_data.H>
#include <model_parser.H>

namespace seed_field {

// Reads the 1D background model and (re)computes R_star, |W|, t_dyn,
// and the global-damping time window into the problem:: namespace.
// Shared by problem_initialize() (cold start) and problem_restart()
// (restart from checkpoint) -- both need these, but only cold start
// should also call generate_modes()/calibrate_and_diagnose(), since on
// restart the B field is already restored from the checkpoint and
// must not be regenerated or re-calibrated. Without this being called
// again on restart, problem::damping_rate/damping_end_time_s stay at
// their _prob_params defaults (0.0) because problem_initialize() is
// cold-start-only in Castro (Castro::initData(), never
// Castro::restart()) -- silently disabling global damping after any
// restart, a real bug found and fixed while adding restart support
// (docs/teoria.md Sec 6.5).
inline void compute_background_star_and_damping()
{
    read_model_file(problem::model_name);

    amrex::Real rho_c = model::profile(0).state(0, model::idens);
    amrex::Real surface_floor = rho_c * problem::surface_density_frac;

    int i_surface = -1;
    for (int i = 0; i < model::npts; ++i) {
        if (model::profile(0).state(i, model::idens) < surface_floor) {
            i_surface = i;
            break;
        }
    }
    if (i_surface < 0) {
        amrex::Error("wd_braithwaite: surface_density_frac never crossed -- "
                     "it must be set above the model file's own density "
                     "floor (density_floor_fraction in "
                     "scf/castro_model_writer.py), or R_star silently "
                     "becomes the domain edge instead of the stellar surface.");
    }
    problem::R_star_cm = model::profile(0).r(i_surface);

    // |W| = G integral_0^Rstar [M(r)/r] dM, dM = 4 pi r^2 rho(r) dr --
    // exact for a spherically symmetric star (Gauss's law shell
    // reduction), computed directly from the 1D model via the
    // trapezoidal rule. Independent of Castro's own self-gravity
    // module (USE_GRAV=TRUE as of Step 3 -- Castro's Poisson/Monopole
    // solve handles the actual dynamical gravity; this 1D integral is
    // only ever used for the E_mag/|W| calibration and diagnostics).
    amrex::Real M_enc = 0.0_rt;
    amrex::Real W = 0.0_rt;
    for (int i = 1; i <= i_surface; ++i) {
        amrex::Real r0 = model::profile(0).r(i - 1);
        amrex::Real r1 = model::profile(0).r(i);
        amrex::Real rho0 = model::profile(0).state(i - 1, model::idens);
        amrex::Real rho1 = model::profile(0).state(i, model::idens);
        amrex::Real dr = r1 - r0;

        amrex::Real M0 = M_enc;
        amrex::Real dM = 0.5_rt * (rho0 * r0 * r0 + rho1 * r1 * r1) * 4.0_rt * M_PI * dr;
        M_enc += dM;
        amrex::Real M1 = M_enc;

        amrex::Real integrand0 = M0 * rho0 * r0;
        amrex::Real integrand1 = M1 * rho1 * r1;
        W -= 4.0_rt * M_PI * C::Gconst * 0.5_rt * (integrand0 + integrand1) * dr;
    }
    problem::W_abs_erg = std::abs(W);

    // t_dyn = sqrt(R_star^3 / (G M)) -- same formula as
    // dashboard/units.py's dynamical_time().
    problem::t_dyn_s = std::sqrt(std::pow(problem::R_star_cm, 3) /
                                  (C::Gconst * M_enc));
    problem::damping_rate = 1.0_rt / (problem::damping_timescale_in_tdyn * problem::t_dyn_s);
    problem::damping_end_time_s = problem::damping_end_time_in_tdyn * problem::t_dyn_s;
    problem::damping_ramp_start_s = problem::damping_end_time_s -
        problem::damping_ramp_in_tdyn * problem::t_dyn_s;

    amrex::Print() << "wd_braithwaite: R_star = " << problem::R_star_cm
                    << " cm, M_enclosed(R_star) = " << M_enc / 1.989e33
                    << " Msun, |W| = " << problem::W_abs_erg << " erg\n"
                    << "wd_braithwaite: t_dyn = " << problem::t_dyn_s
                    << " s, global damping active [0, " << problem::damping_end_time_s
                    << " s], ramping off from " << problem::damping_ramp_start_s << " s\n";
}

// Smoothstep taper: 1 inside (1-transition_frac)*R_star, smoothly to 0
// at R_star, identically 0 outside. Confinement radius (R_star) comes
// from the real SCF background star (Step 1); the taper SHAPE is a
// clean analytic smoothstep rather than the model's own (placeholder-
// EOS) density profile -- documented choice, not hidden: Step 2 cares
// about a real, physically-set confinement radius, not exactly
// reproducing rho(r) inside Castro's different EOS.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real envelope(amrex::Real r, amrex::Real R_star, amrex::Real transition_frac)
{
    using namespace amrex;
    if (r >= R_star) {
        return 0.0_rt;
    }
    Real r_taper_start = R_star * (1.0_rt - transition_frac);
    if (transition_frac <= 0.0_rt || r <= r_taper_start) {
        return 1.0_rt;
    }
    Real x = (r - r_taper_start) / (R_star - r_taper_start);
    return 1.0_rt - (3.0_rt * x * x - 2.0_rt * x * x * x);
}

// Sum of n_modes plane waves with random direction, phase, and
// amplitude, for one Cartesian component of A. The wavenumbers span a
// real band (see generate_modes()) -- spectral richness, not n_modes
// copies of one dominant_scale with different phases.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real mode_sum(int comp, amrex::Real x, amrex::Real y, amrex::Real z)
{
    amrex::Real total = 0.0_rt;
    for (int m = 0; m < seed_modes::n_modes_actual; ++m) {
        amrex::Real kdotx = seed_modes::kx(comp, m) * x
                           + seed_modes::ky(comp, m) * y
                           + seed_modes::kz(comp, m) * z;
        total += seed_modes::amp(comp, m) * std::cos(kdotx + seed_modes::phase(comp, m));
    }
    return total;
}

// The envelope is applied to A HERE, before any curl is taken. B is
// always computed as a discrete difference of THIS function (never of
// mode_sum() alone with the envelope multiplied in afterward) -- see
// B_gaussian_at() and problem_initialize_mhd_data.H. Multiplying a B
// field by a confinement window after curl breaks div B = 0, because
// curl(env * A) != env * curl(A) in general; multiplying A before a
// LINEAR difference operator does not introduce that error, since the
// difference is of the same (already windowed) function everywhere.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real A_component(int comp, amrex::Real x, amrex::Real y, amrex::Real z)
{
    amrex::Real r = std::sqrt(x * x + y * y + z * z);
    return problem::A_amplitude_scale
         * envelope(r, problem::R_star_cm, problem::envelope_transition_frac)
         * mode_sum(comp, x, y, z);
}

// Edge-centered evaluation points, matching Castro's own
// Exec/mhd_tests/LoopAdvection/problem_initialize_mhd_data.H convention
// generalized from a single A_z (2D) to A_x, A_y, A_z (3D): each
// component is node-centered in the two directions it is NOT aligned
// with, and cell-centered in the direction it IS aligned with.
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

// Raw discrete curl in GAUSSIAN units (B = curl A, no sqrt(4 pi)
// factor). Used only for the calibration/diagnostic quadrature below,
// which must match the SCF's E_mag = integral B^2/(8 pi) dV convention.
// The Castro-native Heaviside-Lorentz conversion (B' = B/sqrt(4 pi),
// see dashboard/units.py's gauss_to_castro -- same factor, same reason)
// is applied only in problem_initialize_mhd_data.H, at the last step
// before a value is written into Castro's own state array.
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

// Generates the fixed mode set from the seeded RNG. Host-only, called
// once from problem_initialize() before any grid is filled -- the modes
// must be identical for every cell/rank/tile, so (unlike the per-cell
// amrex::RandomEngine perturbations used elsewhere in Castro, e.g.
// Exec/hydro_tests/simple_convection) they cannot be generated per-cell.
inline void generate_modes()
{
    if (problem::n_modes > seed_modes::MAX_MODES) {
        amrex::Error("wd_braithwaite: n_modes exceeds seed_modes::MAX_MODES");
    }
    seed_modes::n_modes_actual = problem::n_modes;

    std::mt19937_64 rng(static_cast<std::uint64_t>(problem::rng_seed));
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);

    amrex::Real R = problem::R_star_cm;
    amrex::Real scale_center = problem::dominant_scale_over_R * R;
    amrex::Real scale_lo = scale_center / problem::mode_scale_band_factor;
    amrex::Real scale_hi = scale_center * problem::mode_scale_band_factor;

    // Real spectral band, not n_modes copies of one scale: k sampled
    // log-uniformly between the scales the band covers (short and long
    // wavelength modes both appear), direction sampled isotropically
    // (no imposed box symmetry), amplitude varies mode-to-mode.
    amrex::Real log_k_lo = std::log(2.0 * M_PI / scale_hi);
    amrex::Real log_k_hi = std::log(2.0 * M_PI / scale_lo);

    for (int comp = 0; comp < seed_modes::NCOMP; ++comp) {
        for (int m = 0; m < seed_modes::n_modes_actual; ++m) {
            amrex::Real k_mag = static_cast<amrex::Real>(
                std::exp(log_k_lo + uniform01(rng) * (log_k_hi - log_k_lo)));
            amrex::Real cos_th = static_cast<amrex::Real>(2.0 * uniform01(rng) - 1.0);
            amrex::Real sin_th = std::sqrt(amrex::max(0.0_rt, 1.0_rt - cos_th * cos_th));
            amrex::Real phi = static_cast<amrex::Real>(2.0 * M_PI * uniform01(rng));

            seed_modes::kx(comp, m) = k_mag * sin_th * std::cos(phi);
            seed_modes::ky(comp, m) = k_mag * sin_th * std::sin(phi);
            seed_modes::kz(comp, m) = k_mag * cos_th;
            seed_modes::phase(comp, m) = static_cast<amrex::Real>(2.0 * M_PI * uniform01(rng));
            // mode-to-mode amplitude variation (0.5-1.5), independent of
            // the global A_amplitude_scale calibrated below
            seed_modes::amp(comp, m) = static_cast<amrex::Real>(0.5 + uniform01(rng));
        }
    }
}

// One-time host-side quadrature over a sample cube covering the star.
// This is a diagnostic/calibration ESTIMATE at a fixed sample
// resolution independent of the actual Castro production grid
// (n_cell in inputs) -- both the achieved E_mag/|W| and the div B /
// exterior-leakage diagnostics are printed so that this is visible,
// not hidden.
//
// Sets problem::A_amplitude_scale so that a LINEAR rescaling of A
// (curl is linear, so this preserves div B = 0 exactly; rescaling B
// directly after the fact would not, and would also undo the
// envelope-before-curl guarantee above) delivers the requested
// E_mag/|W| (E_mag is quadratic in amplitude).
inline void calibrate_and_diagnose(int n_sample = 128)
{
    using namespace amrex;

    Real R = problem::R_star_cm;
    Real L = 2.4_rt * R;
    Real dx_s = L / static_cast<Real>(n_sample);
    Real dxa[3] = {dx_s, dx_s, dx_s};
    Real problo[3] = {-1.2_rt * R, -1.2_rt * R, -1.2_rt * R};
    Real cell_vol = dx_s * dx_s * dx_s;

    // Pass 1: raw energy at A_amplitude_scale = 1, to get the linear
    // rescaling factor.
    problem::A_amplitude_scale = 1.0_rt;
    Real E_mag_raw = 0.0_rt;

    for (int k = 0; k < n_sample; ++k) {
        for (int j = 0; j < n_sample; ++j) {
            for (int i = 0; i < n_sample; ++i) {
                Real Bx, By, Bz;
                B_gaussian_at(i, j, k, dxa, problo, Bx, By, Bz);
                E_mag_raw += (Bx * Bx + By * By + Bz * Bz) * cell_vol;
            }
        }
    }
    E_mag_raw /= (8.0_rt * M_PI);

    Real E_mag_target = problem::E_mag_over_W_target * problem::W_abs_erg;
    problem::A_amplitude_scale = (E_mag_raw > 0.0_rt)
        ? std::sqrt(E_mag_target / E_mag_raw) : 0.0_rt;

    // Pass 2: at the calibrated scale, measure achieved E_mag/|W|, div
    // B, and exterior leakage. The two ratio diagnostics are scale-
    // invariant (a common amplitude factor cancels), so pass 1 and
    // pass 2 would agree on those; recomputed anyway at the final scale
    // so the printed numbers are the ones that will actually go into
    // Castro's B array (up to the Gaussian->Castro sqrt(4 pi) factor,
    // which also cancels in both ratios).
    Real E_mag_achieved = 0.0_rt;
    Real max_divB_h_over_B = 0.0_rt;
    Real max_B_inside = 0.0_rt;
    Real max_B_outside = 0.0_rt;

    for (int k = 0; k < n_sample; ++k) {
        for (int j = 0; j < n_sample; ++j) {
            for (int i = 0; i < n_sample; ++i) {
                Real Bx, By, Bz;
                B_gaussian_at(i, j, k, dxa, problo, Bx, By, Bz);
                Real Bmag = std::sqrt(Bx * Bx + By * By + Bz * Bz);
                E_mag_achieved += (Bx * Bx + By * By + Bz * Bz) * cell_vol;

                Real Bx_r, By_r, Bz_r;
                B_gaussian_at(i + 1, j, k, dxa, problo, Bx_r, By_r, Bz_r);
                Real Bx_j, By_j, Bz_j;
                B_gaussian_at(i, j + 1, k, dxa, problo, Bx_j, By_j, Bz_j);
                Real Bx_k, By_k, Bz_k;
                B_gaussian_at(i, j, k + 1, dxa, problo, Bx_k, By_k, Bz_k);

                Real divB = (Bx_r - Bx) / dx_s + (By_j - By) / dx_s + (Bz_k - Bz) / dx_s;

                Real x = problo[0] + dx_s * (static_cast<Real>(i) + 0.5_rt);
                Real y = problo[1] + dx_s * (static_cast<Real>(j) + 0.5_rt);
                Real z = problo[2] + dx_s * (static_cast<Real>(k) + 0.5_rt);
                Real r = std::sqrt(x * x + y * y + z * z);

                if (Bmag > 0.0_rt) {
                    max_divB_h_over_B = amrex::max(max_divB_h_over_B, std::abs(divB) * dx_s / Bmag);
                }
                if (r < R) {
                    max_B_inside = amrex::max(max_B_inside, Bmag);
                } else {
                    max_B_outside = amrex::max(max_B_outside, Bmag);
                }
            }
        }
    }
    E_mag_achieved /= (8.0_rt * M_PI);

    problem::diag_E_mag_over_W_achieved = E_mag_achieved / problem::W_abs_erg;
    problem::diag_max_divB_h_over_B = max_divB_h_over_B;
    problem::diag_Bmax_outside_over_inside =
        (max_B_inside > 0.0_rt) ? max_B_outside / max_B_inside : 0.0_rt;

    amrex::Print()
        << "\n=== wd_braithwaite seed field diagnostics (Step 2 acceptance criteria) ===\n"
        << "  sample resolution: " << n_sample << "^3 over [-1.2,1.2] R_star "
        << "(host-side estimate, NOT the production Castro grid)\n"
        << "  R_star = " << R << " cm, |W| = " << problem::W_abs_erg << " erg\n"
        << "  E_mag/|W| requested = " << problem::E_mag_over_W_target << "\n"
        << "  E_mag/|W| achieved  = " << problem::diag_E_mag_over_W_achieved << "\n"
        << "  max|div B|*h/|B|    = " << problem::diag_max_divB_h_over_B
        << "   (acceptance: ~1e-14, machine precision)\n"
        << "  max|B|_outside / max|B|_inside = " << problem::diag_Bmax_outside_over_inside
        << "   (acceptance: << 1, confirms envelope confinement)\n"
        << "  Fluff/sponge treatment (v_A limiter etc.) is NOT part of this step -- "
        << "Step 3's job. This is a static initial condition only.\n"
        << "===========================================================================\n\n";
}

} // namespace seed_field

#endif
