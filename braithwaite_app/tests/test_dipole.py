"""Proves the dipole-moment surface-integral formula (core/dipole.py)
recovers a KNOWN moment from a hand-built analytic dipole field before
it is ever trusted on real plotfile data -- the same "verify against a
known case first" discipline this project used for the Ampere-law test,
the sound-speed finite difference, etc. A general (non-axis-aligned)
moment is used deliberately, so the test exercises the full 3D
projection formula, not just the easier axisymmetric special case.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.dipole import dipole_moment_from_grid


def _analytic_dipole_grid(m_true, n=64, half_width_cm=5e8):
    """B(r) = [3 n(n.m) - m] / r^3 evaluated on a regular grid -- the
    exact formula core/dipole.py's docstring derives the measurement
    from, built independently here (not by calling anything in
    core/dipole.py) so the test is a real check, not a tautology.
    """
    spacing = 2 * half_width_cm / n
    xs = -half_width_cm + spacing * (np.arange(n) + 0.5)
    x, y, z = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(x**2 + y**2 + z**2)
    r_safe = np.where(r < spacing, spacing, r)  # avoid the r=0 singularity; unused far from center anyway
    nx, ny, nz = x / r_safe, y / r_safe, z / r_safe
    n_dot_m = nx * m_true[0] + ny * m_true[1] + nz * m_true[2]
    b_x = (3 * nx * n_dot_m - m_true[0]) / r_safe**3
    b_y = (3 * ny * n_dot_m - m_true[1]) / r_safe**3
    b_z = (3 * nz * n_dot_m - m_true[2]) / r_safe**3
    return {
        "B_x": b_x, "B_y": b_y, "B_z": b_z,
        "dims": (n, n, n),
        "origin_cm": (-half_width_cm, -half_width_cm, -half_width_cm),
        "spacing_cm": (spacing, spacing, spacing),
    }


def test_dipole_moment_recovers_synthetic_field():
    m_true = np.array([1.3e30, -0.7e30, 2.1e30])  # off-axis on purpose
    g = _analytic_dipole_grid(m_true)

    m_measured = dipole_moment_from_grid(g, center_cm=(0.0, 0.0, 0.0), radius_cm=2.5e8, n_points=2048)

    rel_err = np.linalg.norm(m_measured - m_true) / np.linalg.norm(m_true)
    assert rel_err < 0.02, f"recovered {m_measured}, expected {m_true}, rel_err={rel_err:.4f}"


def test_dipole_moment_is_stable_across_radii_for_a_pure_dipole():
    """A real check that the multi-radius stability criterion
    measure_dipole_moment() relies on actually behaves as designed: a
    pure analytic dipole (no higher multipoles at all) must pass it."""
    m_true = np.array([0.0, 0.0, 4.0e30])
    g = _analytic_dipole_grid(m_true)

    mags = [
        np.linalg.norm(dipole_moment_from_grid(g, (0.0, 0.0, 0.0), radius_cm=frac * 1e8, n_points=1024))
        for frac in (2.0, 3.0, 4.0)
    ]
    assert max(mags) <= 1.05 * min(mags)  # tighter than the 1.3x production bar -- this case has zero contamination
