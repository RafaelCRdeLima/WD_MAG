"""Imposicao do campo toroidal pos-SCF (D6). Ver plano_wd_magnetizada.md secao 4.

O toroidal NAO e' extraido da barotropia (D6: isso so' da' Bt/Bp de poucos
por cento). E' imposto por cima do poloidal ja convergido:

    B_phi = zeta * (u - u_c)^(m+1) / omega     para u > u_c, m >= 1

confinado a regiao das ultimas linhas de fluxo fechadas (u > u_c). O
expoente m+1 >= 2 mantem B_phi e sua derivada continuas em u=u_c.
"""

import numpy as np

import diagnostics as diag


def find_uc(u, rho, r, theta):
    """u_c = valor maximo de u ao longo da superficie estelar (H=0).

    Contornos com u > u_c nao tocam a superficie em lugar nenhum — ficam
    inteiramente internos (linhas de fluxo fechadas). Esta e' a fronteira
    do toro (D6)."""
    ntheta = len(theta)
    u_surf = np.zeros(ntheta)
    for j in range(ntheta):
        R_surf_j = diag.surface_radius(rho, r, j)
        u_surf[j] = np.interp(R_surf_j, r, u[:, j])
    return np.max(u_surf)


def impose_toroidal(u, rho, r, theta, zeta, m_tor=1):
    """B_phi = zeta*(u-u_c)^(m_tor+1)/omega para u>u_c, 0 fora. Retorna (B_phi, u_c)."""
    if m_tor < 1:
        raise ValueError("m_tor >= 1 exigido para B_phi e sua derivada serem continuas em u_c")
    u_c = find_uc(u, rho, r, theta)
    omega = r[:, None] * np.sin(theta)[None, :]
    mask = u > u_c
    with np.errstate(divide="ignore", invalid="ignore"):
        omega_safe = np.where(omega > 0, omega, 1.0)
        Bphi = np.where(mask & (omega > 0),
                         zeta * (u - u_c) ** (m_tor + 1) / omega_safe, 0.0)
    return Bphi, u_c


def solve_zeta_for_energy_ratio(u, rho, r, theta, target_ratio, m_tor=1):
    """Ajusta zeta para atingir Bt/Bp = target_ratio (razao de ENERGIAS
    E_tor/E_pol — ver dashboard, que mostra tambem a razao de amplitudes).

    B_phi e' linear em zeta, logo E_tor e' quadratico em zeta: basta uma
    corrida de referencia com zeta=1 e reescalar.
    """
    Bphi_unit, u_c = impose_toroidal(u, rho, r, theta, zeta=1.0, m_tor=m_tor)
    Br, Bth = diag.poloidal_field(u, r, theta)
    E_pol, _, _ = diag.magnetic_energies(Br, Bth, np.zeros_like(Br), r, theta)
    _, E_tor_unit, _ = diag.magnetic_energies(Br, Bth, Bphi_unit, r, theta)

    if E_tor_unit <= 0:
        raise ValueError(
            "regiao de linhas fechadas vazia (u_c >= max(u) no dominio) — "
            "nao ha onde impor o toroidal com este equilibrio poloidal"
        )

    zeta = np.sqrt(target_ratio * E_pol / E_tor_unit)
    Bphi = zeta * Bphi_unit
    return Bphi, zeta, u_c


def bt_bp_ratios(Br, Bth, Bphi, r, theta):
    """As duas definicoes de Bt/Bp que a literatura mistura — mostrar as
    duas sempre (ver prompt do dashboard, Aba 1)."""
    E_pol, E_tor, E_mag = diag.magnetic_energies(Br, Bth, Bphi, r, theta)
    ratio_energy = E_tor / E_pol if E_pol > 0 else float("nan")
    B_pol_max = np.max(np.sqrt(Br**2 + Bth**2))
    B_tor_max = np.max(np.abs(Bphi))
    ratio_amplitude = B_tor_max / B_pol_max if B_pol_max > 0 else float("nan")
    return ratio_energy, ratio_amplitude


def torus_radial_extent(u, rho, r, theta, u_c, j_index=None):
    """Espessura radial do toro (regiao u>u_c) ao longo de um theta fixo
    (equador por padrao) — usado na Aba 3 para checar quantas celulas do
    Castro atravessam o toro."""
    j = j_index if j_index is not None else len(theta) // 2
    in_torus = (u[:, j] > u_c) & (rho[:, j] > 0)
    if not np.any(in_torus):
        return 0.0, 0.0, 0.0
    idx = np.nonzero(in_torus)[0]
    r_inner, r_outer = r[idx[0]], r[idx[-1]]
    return r_inner, r_outer, r_outer - r_inner


def closed_torus_volume_fraction(u, rho, r, theta, u_c):
    """fracao do volume estelar com u > u_c (o toro de linhas fechadas)."""
    inside_star = rho > 0
    in_torus = (u > u_c) & inside_star
    vol_torus = diag.volume_integral(in_torus.astype(float), r, theta)
    vol_star = diag.volume_integral(inside_star.astype(float), r, theta)
    return vol_torus / vol_star if vol_star > 0 else 0.0
