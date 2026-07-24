"""EOS degenerada (gás de elétrons, T=0) para a anã branca. Ver plano_wd_magnetizada.md secao 4."""

import numpy as np

A_CONST = 6.01e22   # dyn cm^-2


def B_of_mu_e(mu_e):
    """B(mu_e) = 9.82e5 * mu_e g/cm^3, com Y_e = 1/mu_e."""
    return 9.82e5 * mu_e


B_CONST = B_of_mu_e(2.0)  # valor padrao historico do modulo (Y_e=0.5)


def pressure(x):
    """P(x), x = p_F / (m_e c). Independe de mu_e."""
    x = np.asarray(x, dtype=float)
    return A_CONST * (x * (2 * x**2 - 3) * np.sqrt(x**2 + 1) + 3 * np.arcsinh(x))


def density(x, mu_e=2.0):
    """rho(x) = B(mu_e) x^3."""
    x = np.asarray(x, dtype=float)
    return B_of_mu_e(mu_e) * x**3


def enthalpy(x, mu_e=2.0):
    """H(x) = (8A/B) [sqrt(1+x^2) - 1], normalizada a H=0 na superficie (x=0)."""
    x = np.asarray(x, dtype=float)
    B = B_of_mu_e(mu_e)
    return (8 * A_CONST / B) * (np.sqrt(1 + x**2) - 1)


def x_of_enthalpy(H, mu_e=2.0):
    """Inversa de enthalpy(x): x(H) = sqrt[(1 + H B/(8A))^2 - 1]. H < 0 -> x = 0 (vacuo)."""
    H = np.asarray(H, dtype=float)
    B = B_of_mu_e(mu_e)
    arg = (1 + H * B / (8 * A_CONST)) ** 2 - 1
    x = np.sqrt(np.clip(arg, 0.0, None))
    return np.where(H < 0, 0.0, x)


def density_of_enthalpy(H, mu_e=2.0):
    """rho(H) = EOS^{-1}(H), com rho=0 para H<0."""
    return density(x_of_enthalpy(H, mu_e), mu_e)


def sound_speed(x, mu_e=2.0):
    """c_s = sqrt(dP/drho). dP/dx = 8A x^4/sqrt(1+x^2) (derivada analitica
    de pressure(x)); drho/dx = 3B x^2. c_s^2 = (dP/dx)/(drho/dx)."""
    x = np.asarray(x, dtype=float)
    B = B_of_mu_e(mu_e)
    cs2 = (8 * A_CONST / (3 * B)) * x**2 / np.sqrt(1 + x**2)
    return np.sqrt(cs2)


# Limiares de decaimento beta inverso (neutronizacao), Tabela 2 de
# Boshkayev, Rueda, Ruffini & Siutsou 2013, ApJ 762, 117 (energias de
# limiar de Audi et al. 2003, densidade critica via EOS RFMT) --
# CONFERIDOS contra o artigo publicado (nao re-derivados aqui). mu_e=A/Z.
#   nuclideo   mu_e        rho_threshold (g/cm^3)
#   4He        2.0         1.39e11
#   12C        2.0         3.97e10
#   16O        2.0         1.94e10
#   56Fe       56/26=2.1538 1.18e9
_NEUTRONIZATION_TABLE = [
    (2.0, 1.94e10, "16O"),
    (56.0 / 26.0, 1.18e9, "56Fe"),
]


def neutronization_threshold_rho_c(mu_e=2.0):
    """Densidade acima da qual o decaimento beta inverso (captura
    eletronica) passa a importar -- a EOS usada em todo este projeto
    (mu_e=const, gas degenerado a T=0) deixa de descrever uma anã branca
    real dali em diante: a composicao comeca a neutronizar, mudando mu_e
    localmente, e o objeto deixa de ser uma anã branca em qualquer sentido
    reconhecivel (ver a nota de contexto em plano_wd_magnetizada.md: o
    modelo de 2.58 Msun de Das & Mukhopadhyay tem raio de 70 km — "nessa
    altura o objeto nao e' uma ana branca em nenhum sentido reconhecivel").

    mu_e=A/Z NAO determina a composicao de forma unica: 4He, 12C e 16O tem
    mu_e=2 exatamente, mas limiares diferentes (1.39e11, 3.97e10, 1.94e10
    g/cm^3 — ver _NEUTRONIZATION_TABLE) porque o limiar de captura
    eletronica depende da energia de ligacao do NUCLIDEO especifico, nao
    so' de mu_e. Este projeto nao tem seletor de composicao alem de mu_e,
    entao esta funcao devolve o valor CONSERVADOR (o mais baixo, tipo
    16O) para mu_e=2 -- melhor avisar cedo demais que tarde demais. Para
    mu_e entre 2.0 e 56/26 (~2.154, 56Fe), interpola linearmente entre as
    duas referencias tabeladas; fora dessa faixa, devolve o extremo mais
    proximo (extrapolacao PLANA, nao de tendencia -- nao ha dado tabelado
    para extrapolar uma tendencia).
    """
    tbl = sorted(_NEUTRONIZATION_TABLE)
    if mu_e <= tbl[0][0]:
        return tbl[0][1]
    if mu_e >= tbl[-1][0]:
        return tbl[-1][1]
    for (mu_lo, rho_lo, _), (mu_hi, rho_hi, _) in zip(tbl, tbl[1:]):
        if mu_lo <= mu_e <= mu_hi:
            frac = (mu_e - mu_lo) / (mu_hi - mu_lo)
            return rho_lo + frac * (rho_hi - rho_lo)
    return tbl[0][1]  # inalcancavel, defensivo
