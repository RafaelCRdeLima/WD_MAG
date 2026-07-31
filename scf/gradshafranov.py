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

# numpy renamed trapz to trapezoid in 2.0. This project develops against 2.x
# but runs on clusters with older numpy -- CENAPAD's system python3 is 3.6.8 --
# so bind whichever exists. The two are the same function.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz

from scipy.special import lpmv


def assoc_legendre_matrix(mu, lmax):
    """P_l^1(mu) para l = 1..lmax, shape (lmax, len(mu))."""
    n = len(mu)
    P1 = np.zeros((lmax, n))
    for l in range(1, lmax + 1):
        P1[l - 1] = lpmv(1, l, mu)
    return P1


def _inner_terms(S, r, lmax):
    """r^-l * int_0^r S_l(r') r'^(l+1) dr', para l = 1..lmax simultaneamente.

    REVISADO (mesmo transbordo de ponto flutuante que poisson.py levou na
    Sec 6.2c, achado aqui ao varrer k0 para cima: o ramo poloidal parava de
    convergir em E_pol/|W| ~ 0.018 com "overflow encountered in multiply"
    nesta linha). A forma original montava r'^(l+1) e r^-l SEPARADAMENTE
    antes de dividir -- em raios estelares (~1e8-1e9 cm) r'^17 ja' e' ~1e146
    e a soma cumulativa transborda, e mesmo quando nao transborda os dois
    fatores absolutos ficam separados por centenas de ordens de grandeza e
    perdem algarismos em silencio.

    Reescrito como recursao avancando em r, carregando o termo ja'
    normalizado por r_i^-l e avancando por uma razao (r_{i-1}/r_i)^l <= 1 --
    nenhuma potencia absoluta de r e' formada. Matematicamente identico,
    so' muda a ordem das operacoes.

    NAO e' copia literal de poisson.py: la' o expoente e' l+1 e aqui e' l,
    porque Delta* nao usa a substituicao chi = r*Phi_l do Laplaciano escalar
    (ver o ATENCAO no topo deste modulo -- essa confusao de expoentes ja'
    produziu um bug real neste arquivo)."""
    nr = len(r)
    l_arr = np.arange(1, lmax + 1)
    T1 = np.zeros((lmax, nr))          # T1[:, i] = r[i]^-l * D_l(r[i])
    for i in range(1, nr):
        ratio = r[i - 1] / r[i] if r[i] > 0 else 0.0
        ratio_pow = ratio ** l_arr     # (lmax,), sempre <= 1
        dr_i = r[i] - r[i - 1]
        increment = 0.5 * (S[:, i - 1] * r[i - 1] * ratio_pow + S[:, i] * r[i]) * dr_i
        T1[:, i] = T1[:, i - 1] * ratio_pow + increment
    return T1


def _outer_terms(S, r, lmax):
    """r^(l+1) * int_r^rmax S_l(r') r'^(-l) dr', para l = 1..lmax.

    Espelho de _inner_terms, recuando em r a partir de r_max. Uma diferenca
    em relacao ao caso escalar: aqui o integrando e' S(r') * r * (r/r')^l,
    com o prefator r sendo o PONTO DE AVALIACAO e nao r'. Trocar a base de
    r_{i+1} para r_i multiplica esse prefator por r_i/r_{i+1} alem da razao
    (r_i/r_{i+1})^l do integrando, de modo que o fator de transporte leva
    expoente l+1 enquanto o incremento leva l. Em poisson.py os dois levam
    o mesmo expoente, porque la' o r' fica dentro da integral."""
    nr = len(r)
    l_arr = np.arange(1, lmax + 1)
    T2 = np.zeros((lmax, nr))          # T2[:, i] = r[i]^(l+1) * E_l(r[i])
    for i in range(nr - 2, -1, -1):
        ratio = r[i] / r[i + 1] if r[i + 1] > 0 else 0.0
        ratio_pow = ratio ** l_arr          # <= 1, para o integrando
        carry = ratio_pow * ratio           # (r_i/r_{i+1})^(l+1), <= 1
        dr_i = r[i + 1] - r[i]
        increment = 0.5 * r[i] * (S[:, i] + S[:, i + 1] * ratio_pow) * dr_i
        T2[:, i] = T2[:, i + 1] * carry + increment
    return T2


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
        S_l[idx] = norm * _trapezoid(integrand, theta, axis=1)

    T1 = _inner_terms(S_l, r, lmax)
    T2 = _outer_terms(S_l, r, lmax)
    l_arr = np.arange(1, lmax + 1)[:, None]
    u_l = -(T1 + T2) / (2 * l_arr + 1)

    sin_theta = np.sin(theta)
    u = np.einsum("lr,lt->rt", u_l, P1) * sin_theta[None, :]
    return u
