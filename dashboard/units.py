"""Conversoes de unidade — ponto unico de verdade (regra R4 do dashboard).

O SCF (scf.py, gradshafranov.py, diagnostics.py) trabalha inteiramente em
CGS gaussiano: B em Gauss diretamente de B=rot(A), energia magnetica em
B^2/(8 pi). O dashboard SEMPRE exibe campo em gauss.

O Castro carrega o campo como B' = B/sqrt(4 pi) (convencao de permeabilidade
unitaria / Heaviside-Lorentz). Essa e' a fonte no. 1 de erro de fator
sqrt(4 pi) ~ 3.5449 mencionada no plano — a conversao mora so' aqui.
"""

import numpy as np

M_SUN = 1.989e33       # g
G_CONST = 6.674e-8     # cm^3 g^-1 s^-2
C_LIGHT = 2.998e10     # cm/s
KM = 1.0e5             # cm


def gauss_to_castro(B_gauss):
    """B' = B/sqrt(4 pi) — o que o Castro espera no problem_initialize."""
    return np.asarray(B_gauss) / np.sqrt(4 * np.pi)


def castro_to_gauss(B_prime):
    """Inversa: B = B' * sqrt(4 pi)."""
    return np.asarray(B_prime) * np.sqrt(4 * np.pi)


def cm_to_km(x_cm):
    return np.asarray(x_cm) / KM


def g_to_msun(m_g):
    return np.asarray(m_g) / M_SUN


def msun_to_g(m_msun):
    return np.asarray(m_msun) * M_SUN


def dynamical_time(M_g, R_cm):
    """t_din = sqrt(R^3 / (G M))."""
    return np.sqrt(R_cm**3 / (G_CONST * M_g))


def alfven_speed(B_gauss_mean, rho_mean):
    """v_A = <B> / sqrt(4 pi rho_medio), em cm/s. B em gauss (CGS puro, nao B')."""
    return np.asarray(B_gauss_mean) / np.sqrt(4 * np.pi * np.asarray(rho_mean))


def alfven_time(R_cm, v_alfven_cms):
    """t_Alfven = R / v_A."""
    return R_cm / v_alfven_cms


# ---------------------------------------------------------------------------
# Formatacao de exibicao — regra de projeto: campo em gauss sempre em notacao
# cientifica, raios em km sempre com 1-2 casas decimais. Ponto unico de
# verdade tanto para o NUMERO quanto para a STRING exibida, para nenhuma
# pagina inventar sua propria formatacao.
# ---------------------------------------------------------------------------

def format_gauss(B_gauss, sig=3):
    """Campo em gauss, notacao cientifica (ex: 7.98e+13 G)."""
    return f"{B_gauss:.{sig}e} G"


def format_km(x_cm, decimals=2):
    """Raio (ja em cm) exibido em km, com casas decimais fixas — converte E formata."""
    return f"{cm_to_km(x_cm):.{decimals}f} km"


def format_km_value(x_km, decimals=2):
    """Formata um valor JA convertido para km (nao converte de novo) — para
    paginas que guardam o escalar em km desde o calculo (Abas 1/2/4, para
    ficarem consistentes com o que vai pro indice/registro)."""
    return f"{x_km:.{decimals}f} km"
