"""Diagnosticos fisicos derivados de um equilibrio SCF convergido.

Regra R1 do dashboard: toda fisica mora em scf/, nao no dashboard. Este
modulo e' a interface de leitura (grandezas derivadas) sobre o resultado
de scf.hachisu_scf — nao resolve nada, so' integra/deriva.
"""

import numpy as np
from eos import pressure, x_of_enthalpy


def volume_integral(field, r, theta):
    """integral de field(r,theta) dV, dV = r^2 sin(theta) dr dtheta dphi (phi->2pi)."""
    integrand = field * r[:, None] ** 2 * np.sin(theta)[None, :]
    over_theta = np.trapezoid(integrand, theta, axis=1)
    return 2 * np.pi * np.trapezoid(over_theta, r)


def gravitational_energy(rho, Phi, r, theta):
    """W = (1/2) integral rho * Phi dV (auto-energia gravitacional; Phi<0 => W<0)."""
    return 0.5 * volume_integral(rho * Phi, r, theta)


def pressure_integral(H, r, theta, mu_e=2.0):
    """integral P dV, com P(x) da EOS (eos.py), x obtido de H."""
    x = x_of_enthalpy(H, mu_e)
    P = pressure(x)
    return volume_integral(P, r, theta)


def poloidal_field(u, r, theta):
    """B_r, B_theta a partir da funcao de fluxo u=omega*A_phi (diferencas finitas).

    Perto do eixo (theta=0 ou pi) e da origem, r^2*sin(theta) fica pequeno e
    o quociente amplifica ruido de discretizacao (theta=pi da sin(pi)~1e-16
    em ponto flutuante, nao exatamente 0 — uma comparacao "!= 0" NAO pega
    isso e deixa a divisao quase-por-zero passar). Usa um limiar em vez de
    igualdade exata, e zera explicitamente no eixo: fisicamente correto la'
    (dV ~ sin(theta) dtheta -> 0), e evita o artefato.
    """
    sin_theta = np.sin(theta)
    dudtheta = np.gradient(u, theta, axis=1)
    dudr = np.gradient(u, r, axis=0)

    r_col = r[:, None]
    sin_row = sin_theta[None, :]
    eps = 1e-8
    safe = (r_col > 0) & (np.abs(sin_row) > eps)
    with np.errstate(divide="ignore", invalid="ignore"):
        Br = np.where(safe, dudtheta / (r_col ** 2 * sin_row), 0.0)
        Btheta = np.where(safe, -dudr / (r_col * sin_row), 0.0)
    return Br, Btheta


def magnetic_energies(Br, Btheta, Bphi, r, theta):
    """E_pol, E_tor, E_mag = integral B^2/(8 pi) dV, separado pol/tor."""
    E_pol = volume_integral((Br ** 2 + Btheta ** 2) / (8 * np.pi), r, theta)
    E_tor = volume_integral(Bphi ** 2 / (8 * np.pi), r, theta)
    return E_pol, E_tor, E_pol + E_tor


def virial_error(rho, Phi, H, Br, Btheta, Bphi, r, theta, mu_e=2.0):
    """VE = |W + 3 int P dV + E_mag| / |W|  (T=0, sem rotacao, Chandrasekhar-Fermi)."""
    W = gravitational_energy(rho, Phi, r, theta)
    Pi = pressure_integral(H, r, theta, mu_e)
    _, _, E_mag = magnetic_energies(Br, Btheta, Bphi, r, theta)
    residual = W + 3 * Pi + E_mag
    VE = abs(residual) / abs(W)
    return VE, W, Pi, E_mag


def surface_radius(rho, r, j_index):
    """raio onde rho cai a zero ao longo de theta[j_index], por interpolacao linear."""
    col = rho[:, j_index]
    nz = np.nonzero(col > 0)[0]
    if len(nz) == 0:
        return 0.0
    i_last = nz[-1]
    if i_last >= len(r) - 1 or col[i_last + 1] > 0:
        return r[i_last]
    r0, r1 = r[i_last], r[i_last + 1]
    rho0, rho1 = col[i_last], col[i_last + 1]
    frac = rho0 / (rho0 - rho1) if rho0 != rho1 else 0.0
    return r0 + frac * (r1 - r0)


def equatorial_polar_radii(rho, r, theta):
    j_eq = len(theta) // 2
    j_pole = 0
    return surface_radius(rho, r, j_eq), surface_radius(rho, r, j_pole)


def density_peak_location(rho, r, theta):
    """(r, theta) do maximo global de rho — deve ficar em r=0 se o centro
    continua sendo o ponto mais denso; usado para checar se a ancoragem em
    rho_c (fixada em r=0) ainda faz sentido fisico em campo forte."""
    i, j = np.unravel_index(np.argmax(rho), rho.shape)
    return r[i], theta[j], i, j
