"""Unit conversions — single source of truth (dashboard rule R4).

The SCF (scf.py, gradshafranov.py, diagnostics.py) works entirely in
Gaussian CGS: B in Gauss directly from B=curl(A), magnetic energy as
B^2/(8 pi). The dashboard ALWAYS displays field in gauss.

Castro loads the field as B' = B/sqrt(4 pi) (unit-permeability /
Heaviside-Lorentz convention). This is the number-1 source of the
sqrt(4 pi) ~ 3.5449 factor error mentioned in the plan — the conversion
lives only here.
"""

import numpy as np

M_SUN = 1.989e33       # g
G_CONST = 6.674e-8     # cm^3 g^-1 s^-2
C_LIGHT = 2.998e10     # cm/s
KM = 1.0e5             # cm


def gauss_to_castro(B_gauss):
    """B' = B/sqrt(4 pi) — what Castro expects in problem_initialize."""
    return np.asarray(B_gauss) / np.sqrt(4 * np.pi)


def castro_to_gauss(B_prime):
    """Inverse: B = B' * sqrt(4 pi)."""
    return np.asarray(B_prime) * np.sqrt(4 * np.pi)


def cm_to_km(x_cm):
    return np.asarray(x_cm) / KM


def g_to_msun(m_g):
    return np.asarray(m_g) / M_SUN


def msun_to_g(m_msun):
    return np.asarray(m_msun) * M_SUN


def dynamical_time(M_g, R_cm):
    """t_dyn = sqrt(R^3 / (G M))."""
    return np.sqrt(R_cm**3 / (G_CONST * M_g))


def alfven_speed(B_gauss_mean, rho_mean):
    """v_A = <B> / sqrt(4 pi mean_rho), in cm/s. B in gauss (pure CGS, not B')."""
    return np.asarray(B_gauss_mean) / np.sqrt(4 * np.pi * np.asarray(rho_mean))


def alfven_time(R_cm, v_alfven_cms):
    """t_Alfven = R / v_A."""
    return R_cm / v_alfven_cms


# ---------------------------------------------------------------------------
# Display formatting — project rule: field in gauss always in scientific
# notation, radii in km always with 1-2 decimal places. Single source of
# truth for both the NUMBER and the displayed STRING, so no page invents
# its own formatting.
# ---------------------------------------------------------------------------

def format_gauss(B_gauss, sig=3):
    """Field in gauss, scientific notation (e.g. 7.98e+13 G)."""
    return f"{B_gauss:.{sig}e} G"


def format_km(x_cm, decimals=2):
    """Radius (already in cm) displayed in km, with fixed decimal places —
    converts AND formats."""
    return f"{cm_to_km(x_cm):.{decimals}f} km"


def format_km_value(x_km, decimals=2):
    """Formats a value ALREADY converted to km (does not convert again) —
    for pages that store the scalar in km from the moment it is computed
    (Tabs 1/2/4, to stay consistent with what goes into the index/registry)."""
    return f"{x_km:.{decimals}f} km"
