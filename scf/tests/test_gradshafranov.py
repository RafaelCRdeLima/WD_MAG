"""Valida gradshafranov.solve_gradshafranov contra solucao analitica de modo l=1 puro,
e contra a lei de Ampere integral (independente da forma fechada).

Fonte: S(r,theta) = S0 sin(theta) P_1^1(cos theta) para r <= R, 0 fora.

NOTA DE PROJETO: uma versao anterior deste solver tinha um bug de expoente
(r'^{l+2} em vez de r'^{l+1}, e r'^{1-l} em vez de r'^{-l} nas integrais de
Green — copiado por engano da estrutura de poisson.py, que tem uma potencia
de r a mais por causa de uma substituicao que Delta* nao usa). O bug so' foi
achado comparando com a lei de Ampere integral (oint B.dl = -int int
source/sin(theta) dr dtheta numa regiao do plano meridional): testar contra
uma "forma fechada" derivada com a MESMA equacao indicial nao pega esse tipo
de erro, porque ambas herdam a mesma normalizacao errada. A licao: uma
checagem que resolve e confere com formulas relacionadas nao e' independente
so' porque o codigo esta em arquivos diferentes.

Para l=1 com fonte constante, a equacao radial R''-2R/r^2=S0 e' ressonante
(r^2 ja e' solucao homogenea), entao a solucao particular tem um termo
log(r), nao e' um polinomio simples.
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gradshafranov import solve_gradshafranov
from scipy.special import lpmv


def u1_exact(r, R, S0):
    """Solucao fechada de R''(r) - 2R(r)/r^2 = S0 (r<=R) verificada por
    substituicao direta na EDO (nao so' pela equacao indicial)."""
    r_safe = np.where(r > 0, r, np.nan)
    inside = S0 * r**2 / 9.0 + (S0 / 3.0) * r**2 * np.log(R / r_safe)
    outside = S0 * R**3 / (9.0 * r_safe)
    return -np.where(r <= R, inside, outside)


def test_single_mode_l1():
    R = 1.0e8
    S0 = 1.0e3
    rmax = 4 * R

    nr, ntheta = 1601, 65
    r = np.linspace(0, rmax, nr)
    theta = np.linspace(0, np.pi, ntheta)
    mu = np.cos(theta)

    P11 = lpmv(1, 1, mu)  # = -sin(theta)
    source = np.where(r[:, None] <= R, S0, 0.0) * (np.sin(theta) * P11)[None, :]

    u = solve_gradshafranov(source, r, theta, lmax=4)

    u1 = u1_exact(r, R, S0)
    u_exact = u1[:, None] * (np.sin(theta) * P11)[None, :]

    mask = np.abs(np.sin(theta)) > 0.1
    rel_err = np.abs(u[1:][:, mask] - u_exact[1:][:, mask]) / (
        np.abs(u_exact[1:][:, mask]) + 1e-30
    )
    max_err = np.nanmax(rel_err)
    print(f"erro relativo maximo (Grad-Shafranov, modo l=1 puro): {max_err:.3e}")
    # tolerancia frouxa (~4%, converge devagar com a resolucao): S0 constante
    # e' um caso l=1 ressonante (fonte nao vai a zero em r=0 como fontes
    # fisicas reais, que tem o fator omega^2 embutido) — a integral "externa"
    # tem um nucleo log(r) perto da origem que o trapezio simples resolve mal.
    # test_ampere_law() abaixo e' a validacao robusta (nao depende de resolver
    # esse nucleo); mantido para nao perder cobertura do caso fechado.
    assert max_err < 0.05, f"erro relativo {max_err} acima da tolerancia"


def test_ampere_law():
    """Lei de Ampere integral numa regiao do plano meridional: oint B.dl deve
    bater com -int int [source/sin(theta)] dr dtheta na mesma regiao. So' usa
    B (uma derivada de u, ja validado em outro teste) e integrais — nao
    envolve a formula de Green nem uma segunda derivada, entao e' um teste
    genuinamente independente da normalizacao interna do solver."""
    from diagnostics import poloidal_field

    R = 1.0e8
    S0 = 1.0e3
    rmax = 4 * R
    nr, ntheta = 801, 81
    r = np.linspace(0, rmax, nr)
    theta = np.linspace(0, np.pi, ntheta)
    mu = np.cos(theta)
    P11 = lpmv(1, 1, mu)
    source = np.where(r[:, None] <= R, S0, 0.0) * (np.sin(theta) * P11)[None, :]

    u = solve_gradshafranov(source, r, theta, lmax=4)
    Br, Bth = poloidal_field(u, r, theta)

    i1, i2 = nr // 4, nr // 4 + 80
    j1, j2 = ntheta // 4, 3 * ntheta // 4
    r1, r2 = r[i1], r[i2]

    seg1 = np.trapezoid(Br[i1:i2 + 1, j1], r[i1:i2 + 1])
    seg2 = np.trapezoid(Bth[i2, j1:j2 + 1] * r2, theta[j1:j2 + 1])
    seg3 = np.trapezoid(Br[i1:i2 + 1, j2], r[i1:i2 + 1])
    seg4 = np.trapezoid(Bth[i1, j1:j2 + 1] * r1, theta[j1:j2 + 1])
    oint_B = seg1 + seg2 - seg3 - seg4

    sin_t = np.sin(theta)
    integrand = source[i1:i2 + 1, j1:j2 + 1] / sin_t[None, j1:j2 + 1]
    rhs = -np.trapezoid(np.trapezoid(integrand, theta[j1:j2 + 1], axis=1), r[i1:i2 + 1])

    rel_err = abs(oint_B / rhs - 1.0)
    print(f"lei de Ampere: oint B.dl={oint_B:.4e}  RHS={rhs:.4e}  erro relativo={rel_err:.3e}")
    assert rel_err < 0.05, f"erro relativo {rel_err} acima da tolerancia"


if __name__ == "__main__":
    test_single_mode_l1()
    test_ampere_law()
    print("OK")
