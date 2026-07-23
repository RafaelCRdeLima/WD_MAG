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
