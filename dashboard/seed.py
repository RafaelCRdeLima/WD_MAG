"""Radius guess to seed the SCF initial profile — not physics (the
(rho_c,k0) SCF lets the real radius emerge from convergence), it just needs
to have the right order of magnitude so the initial guess is not absurd.
Shared across pages to avoid duplicating the fit."""


def r_guess(rho_c):
    """Rough fit from the M-R relation of degenerate white dwarfs (mu_e=2),
    calibrated against tests/test_chandrasekhar_shooting.py."""
    return 1.09e9 * (rho_c / 1.0e6) ** (-0.2436)
