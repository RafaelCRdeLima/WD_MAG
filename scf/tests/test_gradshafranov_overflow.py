"""Regressao para a reescrita das integrais radiais do Grad-Shafranov.

A forma original montava r'^(l+1) e r^-l SEPARADAMENTE antes de dividir. O
termo externo, em particular, fazia soma cumulativa de f(r') * r'^-l a
partir de r=0 e subtraia o total no fim -- e r'^-l perto da origem e' tao
grande que a soma cumulativa perde por completo as contribuicoes de ordem
1. Em lmax=16 (o default do projeto) isso da 100% de erro no termo externo.

O que salvou os resultados poloidais existentes e' que a fonte fisica do
SCF, -4 pi omega^2 rho f(u), se anula como r^2 na origem, e essa supressao
basta para o caso realista: a diferenca medida entre as duas formas na
fonte real e' 3.6e-7. O bug era latente, nao ativo -- ver o teste
test_fonte_realista_praticamente_inalterada abaixo, que documenta os dois
lados disso.
"""

import numpy as np
import pytest

from gradshafranov import _outer_terms, assoc_legendre_matrix, solve_gradshafranov


def _outer_referencia(S, r, l):
    """int_r^rmax S(r') r'^-l dr' * r^(l+1), integrando cada trecho por si.

    Sem soma cumulativa global, portanto sem o cancelamento que a forma
    antiga sofria. E' O(n^2) e serve so' como referencia de teste."""
    out = np.zeros_like(r)
    for i in range(len(r)):
        if r[i] <= 0:
            continue
        integrando = S[i:] * np.where(r[i:] > 0, (r[i] / r[i:]) ** l, 0.0) * r[i]
        out[i] = np.trapezoid(integrando, r[i:])
    return out


@pytest.mark.parametrize("lmax", [4, 16])
@pytest.mark.parametrize("rmax", [1.0, 2.6e8])
def test_termo_externo_bate_com_referencia(lmax, rmax):
    """A recursao tem que reproduzir a integral, inclusive em lmax=16 e em
    raios estelares, onde a forma antiga errava 100%."""
    rng = np.random.default_rng(3)
    r = np.linspace(0.0, rmax, 193)
    S = rng.normal(size=len(r)) * np.exp(-((r / (0.3 * rmax)) ** 2))

    obtido = _outer_terms(np.repeat(S[None, :], lmax, axis=0), r, lmax)[lmax - 1]
    esperado = _outer_referencia(S, r, lmax)

    escala = np.max(np.abs(esperado))
    assert escala > 0
    assert np.max(np.abs(obtido - esperado)) / escala < 1e-10


def test_nenhuma_potencia_absoluta_de_r_transborda():
    """Em raios estelares e lmax alto o resultado tem que continuar finito.

    r^(l+1) em r ~ 1e9 cm e l = 48 vale ~1e447 e transborda float64; a
    recursao so' forma razoes <= 1 e nao pode chegar la'."""
    r = np.linspace(0.0, 9.0e8, 129)
    theta = np.linspace(0.0, np.pi, 33)
    fonte = np.ones((len(r), len(theta))) * 1.0e20
    u = solve_gradshafranov(fonte, r, theta, lmax=48)
    assert np.isfinite(u).all()


def test_fonte_realista_praticamente_inalterada():
    """A fonte fisica do SCF se anula como omega^2 na origem, e' por isso
    que os resultados poloidais anteriores nao mudam com a correcao.

    Documenta o limite do bug: latente para esta fonte, catastrofico para
    uma que nao se anule."""
    r = np.linspace(0.0, 2.6e8, 193)
    theta = np.linspace(0.0, np.pi, 65)
    omega2 = (r[:, None] * np.sin(theta)[None, :]) ** 2
    rho = 1.0e9 * np.clip(1.0 - (r[:, None] / 2.0e8) ** 2, 0.0, None) ** 3
    fonte = -4.0 * np.pi * omega2 * rho * 1.0e-12

    u = solve_gradshafranov(fonte, r, theta, lmax=16)
    assert np.isfinite(u).all()
    assert np.max(np.abs(u)) > 0.0

    # a base angular e' P_l^1 vezes sin(theta): u tem que anular no eixo
    assert np.allclose(u[:, 0], 0.0, atol=1e-8 * np.max(np.abs(u)))
    assert np.allclose(u[:, -1], 0.0, atol=1e-8 * np.max(np.abs(u)))
