"""Solver de Poisson por expansao em polinomios de Legendre (Hachisu 1986).

Malha esferica (r, theta), r em [0, r_max], theta em [0, pi] (sem simetria
equatorial, ver D3 em plano_wd_magnetizada.md).
"""

import numpy as np
from numpy.polynomial import legendre as npleg

G_CONST = 6.674e-8  # cm^3 g^-1 s^-2


def legendre_matrix(mu, lmax):
    """P_l(mu) para l = 0..lmax, shape (lmax+1, len(mu))."""
    n = len(mu)
    P = np.zeros((lmax + 1, n))
    for l in range(lmax + 1):
        c = np.zeros(l + 1)
        c[l] = 1.0
        P[l] = npleg.legval(mu, c)
    return P


def solve_poisson(rho, r, theta, lmax=16):
    """Resolve nabla^2 Phi = 4 pi G rho por funcao de Green radial + Legendre.

    rho: array (nr, ntheta)
    r: array (nr,) crescente, r[0] pode ser 0
    theta: array (ntheta,), theta in [0, pi]

    Retorna Phi: array (nr, ntheta), com Phi -> -GM/r no infinito (Phi<0).
    """
    nr = len(r)
    mu = np.cos(theta)
    P = legendre_matrix(mu, lmax)  # (lmax+1, ntheta)

    # coeficientes angulares de rho: rho_l(r) = (2l+1)/2 * int rho(r,mu) P_l(mu) dmu
    rho_l = np.zeros((lmax + 1, nr))
    for l in range(lmax + 1):
        integrand = rho * P[l][None, :]
        rho_l[l] = (2 * l + 1) / 2 * np.trapezoid(integrand, mu, axis=1)
        # trapz com mu decrescente (theta de 0 a pi -> mu de 1 a -1) da sinal trocado
    rho_l *= -1.0 if mu[0] > mu[-1] else 1.0

    Phi_l = np.zeros((lmax + 1, nr))
    for l in range(lmax + 1):
        f = rho_l[l]
        # D_l(r) = int_0^r f(r') r'^{l+2} dr' ;  E_l(r) = int_r^rmax f(r') r'^{1-l} dr'
        inner = f * r ** (l + 2)
        D_l = np.concatenate(([0.0], np.cumsum(
            0.5 * (inner[1:] + inner[:-1]) * np.diff(r))))
        with np.errstate(divide="ignore", invalid="ignore"):
            r_pow_outer = np.where(r > 0, r ** (1 - l), 0.0)
        outer_int = f * r_pow_outer
        rev_cum = np.concatenate(([0.0], np.cumsum(
            0.5 * (outer_int[1:] + outer_int[:-1]) * np.diff(r))))
        E_l = rev_cum[-1] - rev_cum  # int_r^rmax

        with np.errstate(divide="ignore", invalid="ignore"):
            term1 = np.where(r > 0, D_l / r ** (l + 1), 0.0)
        term2 = E_l * r ** l
        Phi_l[l] = -4 * np.pi * G_CONST / (2 * l + 1) * (term1 + term2)

    Phi = np.einsum("lr,lt->rt", Phi_l, P)
    return Phi
