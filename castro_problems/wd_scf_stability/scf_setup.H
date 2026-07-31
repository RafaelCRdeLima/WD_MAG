#ifndef SCF_SETUP_H
#define SCF_SETUP_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

#include <scf_model.H>
#include <prob_parameters.H>
#include <fundamental_constants.H>

#include <cmath>

namespace scf_model
{

// Mass, radius and dynamical time by quadrature over the model table, so the
// run reports times in units of the star it is actually evolving rather than
// of a nominal one. Axisymmetric, so dV = 2 pi varpi dvarpi dz.
AMREX_INLINE
void compute_background_star()
{
    amrex::Real mass = 0.0;
    amrex::Real r_star = 0.0;
    const amrex::Real rho_surface = problem::surface_density_frac * dens_max;

    for (int i = 0; i < n_varpi; ++i) {
        const amrex::Real vp = varpi_lo + dvarpi * static_cast<amrex::Real>(i);
        // trapezoid weights in both directions
        const amrex::Real wi = (i == 0 || i == n_varpi - 1) ? 0.5_rt : 1.0_rt;
        for (int j = 0; j < n_z; ++j) {
            const amrex::Real zc = z_lo + dz * static_cast<amrex::Real>(j);
            const amrex::Real wj = (j == 0 || j == n_z - 1) ? 0.5_rt : 1.0_rt;
            const amrex::Real rho = dens[i * n_z + j];
            mass += wi * wj * rho * vp;
            if (rho > rho_surface) {
                r_star = amrex::max(r_star, std::sqrt(vp * vp + zc * zc));
            }
        }
    }
    mass *= 2.0_rt * M_PI * dvarpi * dz;

    problem::R_star_cm = r_star;
    problem::t_dyn_s = std::sqrt(r_star * r_star * r_star
                                 / (C::Gconst * amrex::max(mass, 1.0_rt)));

    // Relax-then-measure. The initial condition is an SCF equilibrium sampled
    // onto a Cartesian mesh, and the sampling alone costs it: the field's peak
    // amplitude is retained only to 84.5% at 64^3 and 95.6% at 128^3, and the
    // magnetic force goes as the square of that. So the star arrives out of
    // balance and rings. Damping every velocity for a fixed window lets it
    // settle onto the mesh, and the window then CLOSES -- after
    // damping_end_time_s nothing is damped, and the measurement runs on an
    // undamped star.
    //
    // This is the opposite of what damping did in the Braithwaite problem,
    // where it ran throughout to hold a random seed field still. Damping
    // during the measurement would suppress the very instability being looked
    // for, so the window must end well before the Alfven time the instability
    // grows on. Here t_dyn = 0.84 s against t_A = 1.26 s, so a window of a
    // couple of t_dyn leaves most of the run clean.
    problem::damping_rate = (problem::damping_timescale_in_tdyn > 0.0_rt)
        ? 1.0_rt / (problem::damping_timescale_in_tdyn * problem::t_dyn_s)
        : 0.0_rt;
    problem::damping_end_time_s =
        problem::damping_end_time_in_tdyn * problem::t_dyn_s;
    problem::damping_ramp_start_s = problem::damping_end_time_s
        - problem::damping_ramp_in_tdyn * problem::t_dyn_s;

    amrex::Print() << "  [scf_model] M = " << mass / 1.98892e33
                   << " Msun, R_star = " << r_star
                   << " cm, t_dyn = " << problem::t_dyn_s << " s\n";
    if (problem::damping_end_time_s > 0.0_rt) {
        amrex::Print() << "  [scf_model] damping until t = "
                       << problem::damping_end_time_s << " s ("
                       << problem::damping_end_time_in_tdyn << " t_dyn), "
                       << "ramping from " << problem::damping_ramp_start_s
                       << " s; rate " << problem::damping_rate << " /s\n"
                       << "  [scf_model] MEASUREMENT WINDOW starts there and "
                       << "is undamped\n";
    } else {
        amrex::Print() << "  [scf_model] no damping\n";
    }
}

// Sample the discrete curl on the actual domain mesh and report what the
// initial condition really contains: the divergence it carries, the peak
// field, and how much of the field sits outside the star. These are the same
// three numbers the Python exporter gates on, recomputed here on Castro's own
// mesh so a mismatch between the two is visible rather than assumed away.
AMREX_INLINE
void diagnose_field()
{
    amrex::ParmParse pp("geometry");
    amrex::Vector<amrex::Real> lo, hi;
    pp.getarr("prob_lo", lo);
    pp.getarr("prob_hi", hi);

    amrex::ParmParse pa("amr");
    amrex::Vector<int> ncell;
    pa.getarr("n_cell", ncell);

    amrex::Real dx[3], problo[3];
    for (int d = 0; d < 3; ++d) {
        dx[d] = (hi[d] - lo[d]) / static_cast<amrex::Real>(ncell[d]);
        problo[d] = lo[d];
    }

    amrex::Real bx_max = 0.0, by_max = 0.0, bz_max = 0.0;
    amrex::Real b_max = 0.0, b_max_outside = 0.0, div_max = 0.0;
    amrex::Real e_mag = 0.0;
    const amrex::Real cell_vol = dx[0] * dx[1] * dx[2];

    for (int i = 0; i < ncell[0]; ++i) {
        for (int j = 0; j < ncell[1]; ++j) {
            for (int k = 0; k < ncell[2]; ++k) {
                amrex::Real bx, by, bz;
                B_gaussian_at(i, j, k, dx, problo, bx, by, bz);
                const amrex::Real b2 = bx * bx + by * by + bz * bz;
                const amrex::Real b = std::sqrt(b2);
                // NOTE: bx, by, bz live on three DIFFERENT faces, so their
                // vector magnitude is not a field magnitude at a point and
                // runs high. Reported for the outside/inside ratio only; the
                // honest scalar comparisons are the per-component maxima and
                // E_mag, which is what the SCF's E_tor is compared against.
                bx_max = amrex::max(bx_max, std::abs(bx));
                by_max = amrex::max(by_max, std::abs(by));
                bz_max = amrex::max(bz_max, std::abs(bz));
                b_max = amrex::max(b_max, b);
                e_mag += b2 / (8.0_rt * M_PI) * cell_vol;

                const amrex::Real x = problo[0] + dx[0] * (i + 0.5_rt);
                const amrex::Real y = problo[1] + dx[1] * (j + 0.5_rt);
                const amrex::Real z = problo[2] + dx[2] * (k + 0.5_rt);
                if (density_at(x, y, z)
                    < problem::surface_density_frac * dens_max) {
                    b_max_outside = amrex::max(b_max_outside, b);
                }

                // div B from the same faces the curl produced, which must
                // vanish to round-off; anything else means the staggering
                // convention here and in the writer have drifted apart.
                amrex::Real bx1, by1, bz1, dum1, dum2;
                B_gaussian_at(i + 1, j, k, dx, problo, bx1, dum1, dum2);
                B_gaussian_at(i, j + 1, k, dx, problo, dum1, by1, dum2);
                B_gaussian_at(i, j, k + 1, dx, problo, dum1, dum2, bz1);
                const amrex::Real divb = (bx1 - bx) / dx[0]
                                       + (by1 - by) / dx[1]
                                       + (bz1 - bz) / dx[2];
                div_max = amrex::max(div_max, std::abs(divb));
            }
        }
    }

    problem::diag_max_divB_h_over_B =
        (b_max > 0.0) ? div_max * dx[0] / b_max : 0.0;
    problem::diag_Bmax_outside_over_inside =
        (b_max > 0.0) ? b_max_outside / b_max : 0.0;

    // Alfven crossing time at the peak field and central density -- the unit
    // the Tayler question is asked in, so it belongs in the run log.
    const amrex::Real v_a = b_max / std::sqrt(4.0_rt * M_PI * dens_max);
    const amrex::Real t_alfven = problem::R_star_cm / amrex::max(v_a, 1.0_rt);

    amrex::Print() << "  [scf_model] initial field on the run mesh:\n"
                   << "    max|Bx|,|By|,|Bz| = " << bx_max << ", " << by_max
                   << ", " << bz_max << " G\n"
                   << "    E_mag = " << e_mag << " erg\n"
                   << "    max|div B| h / |B| = "
                   << problem::diag_max_divB_h_over_B << "\n"
                   << "    max|B| outside / inside = "
                   << problem::diag_Bmax_outside_over_inside << "\n"
                   << "    t_alfven (peak B, central rho) = " << t_alfven
                   << " s = " << t_alfven / problem::t_dyn_s << " t_dyn\n";

    if (problem::field_scale != 0.0_rt
        && problem::diag_max_divB_h_over_B > 1.0e-10) {
        amrex::Abort("scf_model: the initial field is not divergence free on "
                     "the run mesh -- the staggering convention in "
                     "scf_model.H and in the writer have drifted apart");
    }
}

}  // namespace scf_model

#endif
