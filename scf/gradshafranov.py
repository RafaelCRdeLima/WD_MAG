"""Solver de Grad-Shafranov por funcao de Green radial, base angular P_l^1.

Delta* u = -4 pi omega^2 rho f(u) - beta beta'(u),  Delta* = d^2/domega^2 - (1/omega) d/domega + d^2/dz^2

Malha esferica (r, theta), mesma malha do poisson.py. omega = r sin(theta).

A base angular correta para Delta* e' P_l^1(cos theta) (Legendre associada de
ordem 1, l >= 1), NAO P_l(cos theta) — ver D_pontos e secao 4 do
plano_wd_magnetizada.md ("Este e' um erro facil de cometer").

As solucoes homogeneas de Delta* u = 0 separaveis como
u = R(r) sin(theta) P_l^1(cos theta) exigem R(r) = r^{l+1} ou r^{-l}
(equacao indicial r^2 R'' - l(l+1) R = 0), um grau deslocado em relacao ao
caso escalar de Poisson (r^l, r^{-(l+1)}).

Com a fonte expandida como S(r,theta) = sum_l S_l(r) sin(theta) P_l^1(cos theta)
e a ortogonalidade  int_{-1}^1 P_l^1(mu) P_m^1(mu) dmu = [2 l(l+1)/(2l+1)] delta_lm,
a equacao radial e' R''(r) - l(l+1)/r^2 R(r) = S_l(r) (verificado por
substituicao simbolica, nao so' pela equacao indicial). Por variacao de
parametros com y1=r^{l+1}, y2=r^{-l} e Wronskiano W=-(2l+1):

u_l(r) = -(1/(2l+1)) [ r^{-l} int_0^r S_l(r') r'^{l+1} dr'
                        + r^{l+1} int_r^{rmax} S_l(r') r'^{-l} dr' ]

ATENCAO: uma versao anterior deste modulo usava r'^{l+2} e r'^{1-l} nessas
integrais (copiado por engano do padrao de poisson.py, que tem uma potencia
de r a mais por causa da substituicao chi=r*Phi_l usada na reducao do
Laplaciano escalar — Delta* nao precisa dessa substituicao). O erro so' foi
achado comparando com a lei de Ampere integral (nenhuma derivada nova,
independente da propria funcao de Green), porque a auto-consistencia
"resolve com a formula, confere com a mesma formula" nao pega esse tipo de
erro — nem testes contra formas fechadas derivadas usando a mesma equacao
indicial (r^{l+1}, r^{-l} continuam corretos; so' os EXPOENTES DENTRO das
integrais estavam errados).
"""

import numpy as np
from scipy.special import lpmv


def assoc_legendre_matrix(mu, lmax):
    """P_l^1(mu) para l = 1..lmax, shape (lmax, len(mu))."""
    n = len(mu)
    P1 = np.zeros((lmax, n))
    for l in range(1, lmax + 1):
        P1[l - 1] = lpmv(1, l, mu)
    return P1


def solve_gradshafranov(source, r, theta, lmax=16):
    """Resolve Delta* u = source por funcao de Green radial + P_l^1.

    source: array (nr, ntheta), o lado direito -4 pi omega^2 rho f(u) - beta beta'(u)
    r, theta: mesma malha do poisson.solve_poisson

    Retorna u: array (nr, ntheta).
    """
    nr = len(r)
    mu = np.cos(theta)
    P1 = assoc_legendre_matrix(mu, lmax)  # (lmax, ntheta), l = 1..lmax

    # projecao angular: S_l(r) = (2l+1)/(2 l(l+1)) * int_0^pi S(r,theta) P_l^1(mu) dtheta
    S_l = np.zeros((lmax, nr))
    for idx, l in enumerate(range(1, lmax + 1)):
        norm = (2 * l + 1) / (2 * l * (l + 1))
        integrand = source * P1[idx][None, :]
        S_l[idx] = norm * np.trapezoid(integrand, theta, axis=1)

    u_l = np.zeros((lmax, nr))
    for idx, l in enumerate(range(1, lmax + 1)):
        f = S_l[idx]
        inner = f * r ** (l + 1)
        D_l = np.concatenate(([0.0], np.cumsum(
            0.5 * (inner[1:] + inner[:-1]) * np.diff(r))))
        with np.errstate(divide="ignore", invalid="ignore"):
            r_pow_outer = np.where(r > 0, r ** (-l), 0.0)
        outer_int = f * r_pow_outer
        rev_cum = np.concatenate(([0.0], np.cumsum(
            0.5 * (outer_int[1:] + outer_int[:-1]) * np.diff(r))))
        E_l = rev_cum[-1] - rev_cum

        with np.errstate(divide="ignore", invalid="ignore"):
            term1 = np.where(r > 0, D_l / r ** l, 0.0)
        term2 = E_l * r ** (l + 1)
        u_l[idx] = -(term1 + term2) / (2 * l + 1)

    sin_theta = np.sin(theta)
    u = np.einsum("lr,lt->rt", u_l, P1) * sin_theta[None, :]
    return u
