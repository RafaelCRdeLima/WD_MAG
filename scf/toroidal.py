"""Imposicao do campo toroidal pos-SCF (D6). Ver plano_wd_magnetizada.md secao 4.

O toroidal NAO e' extraido da barotropia (D6: isso so' da' Bt/Bp de poucos
por cento). E' imposto por cima do poloidal ja convergido:

    B_phi = zeta * (u - u_c)^(m+1) / omega     para u > u_c, m >= 1

confinado a regiao das ultimas linhas de fluxo fechadas (u > u_c). O
expoente m+1 >= 2 mantem B_phi e sua derivada continuas em u=u_c.
"""

import numpy as np

import diagnostics as diag


def _surface_u_profile(u, H, r, theta):
    """u(r,theta) interpolado no raio da superficie (H=0), para cada
    theta — usado por find_uc() e pela verificacao de tangencia
    (check_uc_tangency). Recebe H, NAO rho -- ver diagnostics.surface_radius
    para o motivo (rho e' clipado exatamente a 0 alem da superficie, o que
    faz a interpolacao degenerar para o ponto de grade; H e' continuo e
    cruza zero de verdade entre dois pontos de grade)."""
    ntheta = len(theta)
    u_surf = np.zeros(ntheta)
    r_surf = np.zeros(ntheta)
    for j in range(ntheta):
        r_surf[j] = diag.surface_radius(H, r, j)
        u_surf[j] = np.interp(r_surf[j], r, u[:, j])
    return u_surf, r_surf


def find_uc(u, H, r, theta):
    """u_c = valor maximo de u ao longo da superficie estelar (H=0).

    Contornos com u > u_c nao tocam a superficie em lugar nenhum — ficam
    inteiramente internos (linhas de fluxo fechadas). Esta e' a fronteira
    do toro (D6)."""
    u_surf, _ = _surface_u_profile(u, H, r, theta)
    return np.max(u_surf)


def check_uc_tangency(u, rho, H, r, theta, u_c):
    """Verificacao NUMERICA de que u=u_c e' tangente a' superficie estelar
    em exatamente um ponto (D6) — em vez de comparar visualmente os
    contornos u=u_c e H=0 num plot. O contour() do matplotlib usa
    interpolacao 2D (marching squares) na malha curva (r,theta)->(varpi,z),
    que pode deslocar visualmente duas curvas matematicamente tangentes em
    ate' ~1 celula de malha (ver plots.plot_flux_contours) — a checagem
    real precisa ser feita aqui, com a mesma interpolacao 1D por linha de
    theta usada em find_uc(), nao no desenho. Recebe rho E H: H para achar
    a superficie (surface_radius), rho para a checagem pontual de vacuo
    abaixo.

    REVISADO (varredura da classe de bug de surface_radius, ver
    docs/teoria.md Sec 1.11): `vacuum = rho <= 0` e' ADEQUADO, nao sofre
    do bug. rho<=0 <=> H<=0 exatamente em CADA ponto de malha (pela
    propria construcao do clip da EOS, eos.density_of_enthalpy — os dois
    campos nunca discordam sobre qual ponto e' vacuo), e o teste e' so'
    booleano por ponto de malha ja' existente, nao uma busca de posicao
    sub-malha (nao ha' "interpolacao" para degenerar aqui).

    Retorna dict:
      theta_tangent, r_tangent : localizacao do ponto de tangencia (rad, cm)
      unique                   : True se so' um theta atinge u_c na superficie
                                  (tangencia genuina, nao um platô degenerado)
      margin                   : (u_c - segundo maior u_surf) / |u_c| —
                                  grande => tangencia limpa; perto de 0 =>
                                  quase-degenerado
      vacuum_leak               : True se algum ponto fora da estrela
                                  (rho<=0) tem u > u_c — indicaria que o
                                  toro extrapola a superficie fisica
                                  (bug real, nao artefato de plot)
    """
    u_surf, r_surf = _surface_u_profile(u, H, r, theta)
    jmax = int(np.argmax(u_surf))
    sorted_surf = np.sort(u_surf)
    second = sorted_surf[-2] if len(sorted_surf) > 1 else -np.inf
    margin = (u_c - second) / abs(u_c) if u_c != 0 else float("nan")
    vacuum = rho <= 0
    vacuum_leak = bool(np.any(u[vacuum] > u_c)) if np.any(vacuum) else False
    return {
        "theta_tangent": theta[jmax],
        "r_tangent": r_surf[jmax],
        "unique": margin > 1e-6,
        "margin": margin,
        "vacuum_leak": vacuum_leak,
    }


def impose_toroidal(u, H, r, theta, zeta, m_tor=1):
    """B_phi = zeta*(u-u_c)^(m_tor+1)/omega para u>u_c, 0 fora. Retorna
    (B_phi, u_c). Recebe H (nao rho) -- so' usado para localizar a
    superficie via find_uc(); ver diagnostics.surface_radius."""
    if m_tor < 1:
        raise ValueError("m_tor >= 1 exigido para B_phi e sua derivada serem continuas em u_c")
    u_c = find_uc(u, H, r, theta)
    omega = r[:, None] * np.sin(theta)[None, :]
    mask = u > u_c
    with np.errstate(divide="ignore", invalid="ignore"):
        omega_safe = np.where(omega > 0, omega, 1.0)
        Bphi = np.where(mask & (omega > 0),
                         zeta * (u - u_c) ** (m_tor + 1) / omega_safe, 0.0)
    return Bphi, u_c


def solve_zeta_for_energy_ratio(u, H, r, theta, target_ratio, m_tor=1):
    """Ajusta zeta para atingir Bt/Bp = target_ratio (razao de ENERGIAS
    E_tor/E_pol — ver dashboard, que mostra tambem a razao de amplitudes).

    B_phi e' linear em zeta, logo E_tor e' quadratico em zeta: basta uma
    corrida de referencia com zeta=1 e reescalar.
    """
    Bphi_unit, u_c = impose_toroidal(u, H, r, theta, zeta=1.0, m_tor=m_tor)
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
    Castro atravessam o toro.

    REVISADO (varredura da classe de bug de surface_radius, ver
    diagnostics.surface_radius e docs/teoria.md Sec 1.11): `rho[:, j] > 0`
    aqui e' ADEQUADO, nao precisa de H. Este e' um teste booleano por
    ponto de malha (rho<=0 <=> H<=0 exatamente, sempre, pela propria
    construcao do clip da EOS — nao ha' divergencia possivel entre os
    dois), e o proposito explicito da funcao e' CONTAR CELULAS de malha
    (r_inner/r_outer sao raios de grade, nao posicoes sub-malha) — a
    quantizacao e' o comportamento correto aqui, nao um artefato."""
    j = j_index if j_index is not None else len(theta) // 2
    in_torus = (u[:, j] > u_c) & (rho[:, j] > 0)
    if not np.any(in_torus):
        return 0.0, 0.0, 0.0
    idx = np.nonzero(in_torus)[0]
    r_inner, r_outer = r[idx[0]], r[idx[-1]]
    return r_inner, r_outer, r_outer - r_inner


def closed_torus_volume_fraction(u, rho, r, theta, u_c):
    """fracao do volume estelar com u > u_c (o toro de linhas fechadas).

    REVISADO (mesma varredura acima): `inside_star = rho > 0` e' ADEQUADO
    aqui. Alimenta uma integral de VOLUME (soma sobre celulas de malha,
    via volume_integral/trapezoid) — precisao sub-celula na fronteira nao
    faz parte do desenho desta quadratura de jeito nenhum, com qualquer
    criterio de fronteira (mesma limitacao inerente de W, Pi, E_mag, todas
    somadas do mesmo jeito). Nao e' a mesma categoria de erro que
    surface_radius (que buscava uma POSICAO precisa, nao uma soma)."""
    inside_star = rho > 0
    in_torus = (u > u_c) & inside_star
    vol_torus = diag.volume_integral(in_torus.astype(float), r, theta)
    vol_star = diag.volume_integral(inside_star.astype(float), r, theta)
    return vol_torus / vol_star if vol_star > 0 else 0.0
