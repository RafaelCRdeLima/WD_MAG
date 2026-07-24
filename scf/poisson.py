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


def _inner_terms(f, r, lmax):
    """r^-(l+1) * int_0^r f(r') r'^(l+2) dr', para l=0..lmax simultaneamente.

    REVISADO (transbordo de ponto flutuante achado investigando VE em K alto
    no ramo toroidal autoconsistente, dominio grande -- ver docs/teoria.md
    Sec 6): a forma original formava r'^(l+2) e r^-(l+1) SEPARADAMENTE (por
    exemplo r'~9e8 cm, l=48 -> r'^50 ~ 1e447, transborda float64 muito antes
    de dividir; mesmo em l=16 os dois fatores absolutos ficam separados por
    ~270 ordens de grandeza, perdendo algarismos silenciosamente na soma
    cumulativa mesmo sem chegar a NaN). Reescrito como recursao avancando em
    r, carregando o termo ja' normalizado por r_i^-(l+1) e avancando por uma
    razao (r_{i-1}/r_i)^(l+1) <= 1 a cada passo -- nenhum fator jamais excede
    1 em modulo, nenhuma potencia absoluta de r e' formada. Matematicamente
    identico a' forma antiga (mesma integral), so' a ordem das operacoes
    muda. Custo O(nr) por l, igual a antes (vetorizado em l, loop em r)."""
    nr = len(r)
    l_arr = np.arange(lmax + 1)
    T1 = np.zeros((lmax + 1, nr))  # T1[:, i] = r[i]^-(l+1) * D_l(r[i])
    for i in range(1, nr):
        ratio = r[i - 1] / r[i] if r[i] > 0 else 0.0
        ratio_pow = ratio ** (l_arr + 1)  # (lmax+1,), sempre <= 1
        dr_i = r[i] - r[i - 1]
        increment = 0.5 * (f[:, i - 1] * r[i - 1] * ratio_pow + f[:, i] * r[i]) * dr_i
        T1[:, i] = T1[:, i - 1] * ratio_pow + increment
    return T1


def _outer_terms(f, r, lmax):
    """r^l * int_r^rmax f(r') r'^(1-l) dr', para l=0..lmax simultaneamente.

    Espelho de _inner_terms (ver docstring la'): recursao retrocedendo em r
    a partir de r_max, razao (r_i/r_{i+1})^l <= 1 a cada passo."""
    nr = len(r)
    l_arr = np.arange(lmax + 1)
    T2 = np.zeros((lmax + 1, nr))  # T2[:, i] = r[i]^l * E_l(r[i])
    for i in range(nr - 2, -1, -1):
        ratio = r[i] / r[i + 1] if r[i + 1] > 0 else 0.0
        ratio_pow = ratio ** l_arr  # (lmax+1,), sempre <= 1 (0**0 = 1, ok)
        dr_i = r[i + 1] - r[i]
        increment = 0.5 * (f[:, i] * r[i] + f[:, i + 1] * r[i + 1] * ratio_pow) * dr_i
        T2[:, i] = T2[:, i + 1] * ratio_pow + increment
    return T2


def solve_poisson(rho, r, theta, lmax=16):
    """Resolve nabla^2 Phi = 4 pi G rho por funcao de Green radial + Legendre.

    rho: array (nr, ntheta)
    r: array (nr,) crescente, r[0] pode ser 0
    theta: array (ntheta,), theta in [0, pi]

    Retorna Phi: array (nr, ntheta), com Phi -> -GM/r no infinito (Phi<0).
    """
    mu = np.cos(theta)
    P = legendre_matrix(mu, lmax)  # (lmax+1, ntheta)

    # coeficientes angulares de rho: rho_l(r) = (2l+1)/2 * int rho(r,mu) P_l(mu) dmu
    rho_l = np.zeros((lmax + 1, len(r)))
    for l in range(lmax + 1):
        integrand = rho * P[l][None, :]
        rho_l[l] = (2 * l + 1) / 2 * np.trapezoid(integrand, mu, axis=1)
        # trapz com mu decrescente (theta de 0 a pi -> mu de 1 a -1) da sinal trocado
    rho_l *= -1.0 if mu[0] > mu[-1] else 1.0

    l_arr = np.arange(lmax + 1)
    term1 = _inner_terms(rho_l, r, lmax)
    term2 = _outer_terms(rho_l, r, lmax)
    Phi_l = -4 * np.pi * G_CONST / (2 * l_arr[:, None] + 1) * (term1 + term2)

    Phi = np.einsum("lr,lt->rt", Phi_l, P)
    return Phi
