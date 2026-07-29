"""Magnetic dipole moment of the relaxed field, measured directly from a
plotfile's B grid via a surface integral outside the star -- a number,
not a visual impression of whether field lines "look like" they close
outside the star. This is what tells Jorge whether the Braithwaite
relaxation produces a usable exterior dipole for magnetic braking, and
how strong it is compared to his ~1e9 G target.

Derivation. Outside all currents, any localized source's field has a
multipole expansion; the leading (dipole) term, for a moment vector m at
position r = r*n (n a unit vector), is (Gaussian CGS, matching this
project's B convention throughout -- docs/teoria.md Sec 1.10):

    B(r) = [3 n (n.m) - m] / r^3        (Jackson, Classical
                                          Electrodynamics, eq. 5.56)

so the radial component is B_r = B.n = 2(n.m)/r^3. Multiplying by n_j
and integrating over the full sphere, and using the standard identity
integral(n_i n_j dOmega) = (4*pi/3)*delta_ij (average of n_i*n_j over a
sphere), gives a clean projection formula -- no derivative of B (hence
no numerical-curl noise) is needed, only its value on a sphere:

    m_j = (3 r^3 / 8*pi) * integral( B_r(theta,phi) * n_j(theta,phi) dOmega )

Validity. This formula assumes the sampling sphere sits in a region
where the dipole term actually dominates (current-free exterior, higher
multipoles small). It is not assumed here -- measure_dipole_moment()
evaluates at several radii and reports whether |m| is stable across
them; if it drifts with radius, the exterior isn't cleanly dipolar at
that point and the number should be reported as such, not averaged into
a false single answer. Verified against a synthetic analytic dipole
field before trusting it on real data -- see tests/test_dipole.py.
"""

from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from core.field_reader import estimate_star_radius, load_vector_grid

DEFAULT_RADII_FRAC = (1.2, 1.5, 2.0)  # multiples of R_star -- outside the star, inside the domain
DEFAULT_N_POINTS = 1024


def fibonacci_sphere_directions(n: int) -> np.ndarray:
    """n unit vectors spread evenly over a sphere (golden-angle spiral)
    -- deterministic, so a re-measurement of the same plotfile gives the
    same number, not a slightly different one from a different random
    sample. Equal-area enough that weighting each point by 4*pi/n is a
    good approximation to true solid-angle integration (this is the
    same construction ui/field_view.py uses for streamline seeding,
    duplicated here rather than imported -- this module has no UI
    dependency and shouldn't gain one just to share six lines).
    """
    if n <= 1:
        return np.array([[0.0, 0.0, 1.0]])
    i = np.arange(n)
    y = 1 - 2 * i / (n - 1)
    r_xy = np.sqrt(np.clip(1 - y**2, 0, None))
    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * i
    return np.stack([np.cos(theta) * r_xy, y, np.sin(theta) * r_xy], axis=1)


def _grid_interpolators(g: dict):
    """Trilinear interpolators for B_x/B_y/B_z from the covering-grid
    arrays load_vector_grid() returns -- cell-centered samples on a
    regular grid, so the interpolation nodes are cell centers, not the
    domain corners.
    """
    nx, ny, nz = g["dims"]
    ox, oy, oz = g["origin_cm"]
    sx, sy, sz = g["spacing_cm"]
    xs = ox + sx * (np.arange(nx) + 0.5)
    ys = oy + sy * (np.arange(ny) + 0.5)
    zs = oz + sz * (np.arange(nz) + 0.5)
    return {
        comp: RegularGridInterpolator((xs, ys, zs), g[comp], bounds_error=False, fill_value=0.0)
        for comp in ("B_x", "B_y", "B_z")
    }


def dipole_moment_from_grid(
    g: dict, center_cm: tuple, radius_cm: float, n_points: int = DEFAULT_N_POINTS,
) -> np.ndarray:
    """m = (mx, my, mz) in G*cm^3, from a single sampling sphere of the
    given radius. Low-level entry point -- takes an already-loaded grid
    (from load_vector_grid) so a caller measuring several radii on the
    same plotfile only builds the interpolators once.
    """
    interp = _grid_interpolators(g)
    directions = fibonacci_sphere_directions(n_points)
    points = np.asarray(center_cm) + directions * radius_cm
    b_x = interp["B_x"](points)
    b_y = interp["B_y"](points)
    b_z = interp["B_z"](points)
    b_r = b_x * directions[:, 0] + b_y * directions[:, 1] + b_z * directions[:, 2]
    # dOmega ~= 4*pi/n_points per (near-equal-area) Fibonacci-sphere point
    return (3 * radius_cm**3 / (2 * n_points)) * (b_r[:, None] * directions).sum(axis=0)


def measure_dipole_moment(
    plotfile: Path, radii_frac: tuple = DEFAULT_RADII_FRAC, n_points: int = DEFAULT_N_POINTS,
) -> dict:
    """The full measurement for one plotfile: R_star (same boundary
    convention used everywhere else in the app), the dipole moment
    vector at each of `radii_frac` (multiples of R_star), and a
    stability verdict across them -- consistent |m| across radii is what
    lets a "there is a dipole" claim be made from a number instead of a
    picture.
    """
    g = load_vector_grid(plotfile)
    r_star = estimate_star_radius(g["density"], g["spacing_cm"])
    nx, ny, nz = g["dims"]
    ox, oy, oz = g["origin_cm"]
    sx, sy, sz = g["spacing_cm"]
    center = (ox + sx * nx / 2, oy + sy * ny / 2, oz + sz * nz / 2)

    measurements = []
    for frac in radii_frac:
        radius_cm = r_star * frac
        m = dipole_moment_from_grid(g, center, radius_cm, n_points)
        measurements.append({
            "radius_over_rstar": frac,
            "radius_cm": radius_cm,
            "m_vector_G_cm3": m,
            "m_mag_G_cm3": float(np.linalg.norm(m)),
            "b_pole_at_rstar_G": 2.0 * float(np.linalg.norm(m)) / r_star**3 if r_star > 0 else float("nan"),
        })

    mags = np.array([meas["m_mag_G_cm3"] for meas in measurements])
    # stable = the largest and smallest |m| across the tested radii agree
    # to within 30% -- a loose but honest bar; a clean dipole in vacuum
    # should be flat to a few percent, a contaminated one (higher
    # multipoles or real exterior current) will fail this by a lot, not
    # marginally
    stable = bool(mags.max() <= 1.3 * mags.min()) if mags.min() > 0 else False

    return {
        "plotfile": str(plotfile),
        "time_s": g["time_s"],
        "r_star_cm": r_star,
        "measurements": measurements,
        "stable_across_radii": stable,
    }
