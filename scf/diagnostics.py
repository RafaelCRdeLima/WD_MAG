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
    igualdade exata para blindar a divisao.

    NO EIXO (theta=0 ou pi), r>0: B_theta e' zero de verdade por
    axissimetria (du/dr tambem vai a zero la', entao o zero abaixo esta'
    correto). B_r **nao** e' zero em geral — e' o campo polar (du/dtheta
    e sin(theta) vao a zero juntos; por L'Hopital, du/dtheta/sin(theta) ->
    d2u/dtheta2 = valor finito, tipicamente da mesma ordem de B_pol,max,
    nao zero). Usar o ponto interior mais proximo do eixo como estimativa
    (a razao dudtheta/sin(theta) ja e' suave perto do eixo, entao o vizinho
    e' uma boa aproximacao regularizada — conferido contra uma estimativa
    analitica por L'Hopital de 2 pontos, concordancia a ~0.02% em
    nr=ntheta=129). Zerar Br no eixo (comportamento antigo) e' inofensivo
    para integrais de volume (dV~sin(theta)dtheta->0 ali, entao nao muda
    W/E_pol/VE), mas errado para leituras pontuais do campo — por exemplo
    o escalar "B_central" do dashboard, que ficava reportando 0 G onde o
    campo real e' ~B_pol,max.
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

    on_r = r_col[:, 0] > 0
    if len(theta) > 2:
        if abs(sin_theta[0]) <= eps:
            Br[on_r, 0] = Br[on_r, 1]
        if abs(sin_theta[-1]) <= eps:
            Br[on_r, -1] = Br[on_r, -2]

    return Br, Btheta


def magnetic_energies(Br, Btheta, Bphi, r, theta):
    """E_pol, E_tor, E_mag = integral B^2/(8 pi) dV, separado pol/tor."""
    E_pol = volume_integral((Br ** 2 + Btheta ** 2) / (8 * np.pi), r, theta)
    E_tor = volume_integral(Bphi ** 2 / (8 * np.pi), r, theta)
    return E_pol, E_tor, E_pol + E_tor


def virial_error(rho, Phi, H, Br, Btheta, Bphi, r, theta, mu_e=2.0, T=0.0):
    """VE = |2T + W + 3 int P dV + E_mag| / |W|. T=0 (padrao) reduz ao caso
    sem rotacao (Chandrasekhar-Fermi), igual a antes desta generalizacao —
    o termo 2T ja' estava previsto no plano original (secao 4, "Controle
    de qualidade obrigatorio"), so' nao havia rotacao implementada ainda."""
    W = gravitational_energy(rho, Phi, r, theta)
    Pi = pressure_integral(H, r, theta, mu_e)
    _, _, E_mag = magnetic_energies(Br, Btheta, Bphi, r, theta)
    residual = 2 * T + W + 3 * Pi + E_mag
    VE = abs(residual) / abs(W)
    return VE, W, Pi, E_mag, T


def virial_error_terms(rho, Phi, H, r, theta, mu_e=2.0, rotation=None, poloidal=None, toroidal=None):
    """Wrapper de virial_error() para a arquitetura de termos plugaveis
    (scf/terms/) — reune T/Br/Btheta/Bphi/E_pol/E_tor dos termos ativos
    (None => contribuicao zero) e chama virial_error() (fonte UNICA da
    formula do residuo 2T+W+3Pi+E_mag; nao reimplementada aqui, para nao
    repetir o erro de normalizacao duplicada que ja mordeu este projeto
    uma vez — ver docs/teoria.md secao 1.4)."""
    T = rotation.energy(rho, r, theta)["T"] if rotation is not None else 0.0

    if poloidal is not None:
        pol = poloidal.energy(rho, r, theta)
        Br, Btheta, E_pol = pol["Br"], pol["Bth"], pol["E_pol"]
    else:
        Br = np.zeros_like(rho)
        Btheta = np.zeros_like(rho)
        E_pol = 0.0

    if toroidal is not None:
        tor = toroidal.energy(rho, r, theta)
        Bphi, E_tor = tor["Bphi"], tor["E_tor"]
    else:
        Bphi = np.zeros_like(rho)
        E_tor = 0.0

    VE, W, Pi, E_mag, T = virial_error(rho, Phi, H, Br, Btheta, Bphi, r, theta, mu_e, T=T)
    return {"VE": VE, "W": W, "Pi": Pi, "T": T, "E_pol": E_pol, "E_tor": E_tor,
            "E_mag": E_mag, "Br": Br, "Btheta": Btheta, "Bphi": Bphi}


def equatorial_mass_loss_ratio(Phi, rotation, r, theta, R_eq):
    """Omega^2(R_eq)*R_eq / (dPhi/dr em R_eq, equador) — razao critica de
    perda de massa por rotacao. O numerador e' a forca centrifuga por
    massa na superficie equatorial (para fora); o denominador e' a forca
    gravitacional por massa ali (dPhi/dr > 0, ja' que Phi<0 cresce em
    direcao a 0 para fora — a forca de gravidade em si e' -dPhi/dr, para
    dentro, entao seu MODULO e' dPhi/dr). Razao -> 1 significa que a
    gravidade efetiva se anula no equador (breakup Kepleriano, a estrela
    perde massa ali); razao > 1 e' configuracao que nao pode existir em
    equilibrio. rotation=None ou R_eq<=0 -> 0.0 (sem rotacao, sem risco).

    REVISADO (varredura da classe de bug de surface_radius, ver
    docs/teoria.md Sec 1.11): esta funcao NAO usa rho diretamente -- recebe
    R_eq como argumento (do chamador, que ja' vem de
    equatorial_polar_radii(H,...), corrigido) e interpola dPhi/dr (campo
    continuo, nunca clipado) nesse raio via np.interp(). Ja' estava correta
    antes desta varredura, por construcao."""
    if rotation is None or R_eq <= 0:
        return 0.0
    j_eq = len(theta) // 2
    dPhidr = np.gradient(Phi, r, axis=0)
    g_grav = float(np.interp(R_eq, r, dPhidr[:, j_eq]))
    if g_grav <= 0:
        return float("inf")
    Om_eq = float(rotation.Omega(np.array([R_eq]))[0])
    centrifugal = Om_eq**2 * R_eq
    return centrifugal / g_grav


def surface_radius(H, r, j_index):
    """raio onde H cai a zero (a fronteira fisica H=0 da estrela) ao longo
    de theta[j_index], por interpolacao linear -- em H, NAO em rho.

    BUG HISTORICO (corrigido): esta funcao antes recebia rho e interpolava
    nele. rho = EOS^{-1}(H) e' CLIPADO exatamente a 0.0 para H<=0
    (eos.density_of_enthalpy) -- ou seja, o primeiro ponto de grade alem da
    superficie tem rho EXATAMENTE 0.0, sempre, por construcao. Substituindo
    isso na formula frac = rho0/(rho0-rho1) com rho1=0.0 da' frac=1.0
    SEMPRE, entao a "interpolacao linear" degenerava e devolvia sempre
    r[i_last+1] -- um ponto de grade cru, nunca uma posicao sub-malha real.
    Isso so' foi percebido ao investigar por que uma sequencia de rotacao
    parecia ter um "ponto de virada" ficticio em Omega_c: R_pol/R_eq estava
    caindo em degraus de razoes de inteiros pequenos (37/39, 38/39, ...),
    a assinatura classica de quantizacao de malha, nao de fisica.

    H, em contraste, NUNCA e' clipado -- vem direto da formula continua do
    Bernoulli (H = C - Phi + ...) em todo o dominio, e cruza zero de forma
    suave entre dois pontos de grade adjacentes (nao exatamente EM um deles,
    em geral). Interpolar em H da' a posicao real da superficie -- testado:
    ~46 km de diferenca em relacao ao valor antigo (grid-preso) numa malha
    nr=65, um efeito grande, nao um refinamento cosmetico.
    """
    col = H[:, j_index]
    nz = np.nonzero(col > 0)[0]
    if len(nz) == 0:
        return 0.0
    i_last = nz[-1]
    if i_last >= len(r) - 1 or col[i_last + 1] > 0:
        return r[i_last]
    r0, r1 = r[i_last], r[i_last + 1]
    H0, H1 = col[i_last], col[i_last + 1]
    frac = H0 / (H0 - H1) if H0 != H1 else 0.0
    return r0 + frac * (r1 - r0)


def surface_field(field, H, r, theta):
    """field(r,theta) interpolado no raio da superficie estelar (H=0,
    via surface_radius), para cada theta — usado para reportar campo (ou
    outra grandeza) NA SUPERFICIE, em vez de no maximo/centro interior.
    Fisicamente e' a grandeza comparavel a medidas observacionais (p.ex.
    campo polar/equatorial de uma ana branca magnetizada real), diferente
    de B_pol,max (que ocorre tipicamente no interior, nao na superficie)."""
    ntheta = len(theta)
    out = np.zeros(ntheta)
    for j in range(ntheta):
        R_surf_j = surface_radius(H, r, j)
        out[j] = np.interp(R_surf_j, r, field[:, j])
    return out


def surface_dipolarity(Bpol, H, r, theta):
    """B_pol na superficie estelar (via surface_field/surface_radius — a
    MESMA interpolacao 1D por linha de theta que toroidal.find_uc() usa
    para localizar H=0; nenhum localizador de superficie novo aqui).

    Retorna B_pole (theta=0), B_eq (theta=pi/2), B_surf_max (maximo sobre
    toda a curva da superficie) e a razao B_pole/B_eq ("dipolaridade":
    exatamente 2 para um dipolo puro; desvio mede contaminacao
    multipolar).

    B_pole cai exatamente no eixo de simetria — a mesma singularidade de
    coordenada (B_theta ~ 1/sin(theta)) corrigida em poloidal_field() via
    valor do vizinho interior mais proximo. Como checagem INDEPENDENTE
    (nao so reusar o mesmo valor regularizado), B_pole tambem e' estimado
    por extrapolacao quadratica a partir do campo de superficie nos dois
    pontos de grade mais proximos MAS FORA do eixo (simetria par perto do
    eixo: B_pol,surf(theta) ~ A + B*theta^2) — grande desacordo entre os
    dois indicaria que a correcao de poloidal_field() nao e' suficiente na
    resolucao atual. Conferido numericamente: concordancia a ~0.03% em
    ntheta=129, ~0.01% em ntheta=257 (converge com a malha).
    """
    Bpol_surf = surface_field(Bpol, H, r, theta)
    ntheta = len(theta)
    j_eq = ntheta // 2
    B_pole = float(Bpol_surf[0])
    B_eq = float(Bpol_surf[j_eq])
    B_surf_max = float(np.max(Bpol_surf))
    if ntheta > 2:
        B_pole_extrapolated = float((4 * Bpol_surf[1] - Bpol_surf[2]) / 3.0)
    else:
        B_pole_extrapolated = float("nan")
    dipolarity = B_pole / B_eq if B_eq != 0 else float("nan")
    return {
        "B_pole": B_pole,
        "B_pole_extrapolated": B_pole_extrapolated,
        "B_eq": B_eq,
        "B_surf_max": B_surf_max,
        "dipolarity": dipolarity,
    }


def equatorial_polar_radii(H, r, theta):
    j_eq = len(theta) // 2
    j_pole = 0
    return surface_radius(H, r, j_eq), surface_radius(H, r, j_pole)


def domain_overflow_check(R_eq, R_pol, r_max, tol=0.1):
    """Flags whether the computed surface radius has landed suspiciously
    close to the outer edge of the computational domain (r_max = r[-1])
    in EITHER direction, independently. This is the item-11 gap (docs/
    teoria.md Sec 1.13/Sec 8): a radius near r_max is not trustworthy --
    the true surface may be truncated by too small a box, not genuinely
    located there.

    Checks R_eq and R_pol separately, on purpose: an oblate
    (rotation-dominated) star overflows at the EQUATOR first (R_eq grows,
    R_pol shrinks); a prolate (toroidal-dominated) star overflows at the
    POLE first (R_pol grows, R_eq shrinks) -- the two failure modes are
    each other's mirror image in shape, so a check that only looks at one
    of the two radii misses exactly the family of configurations where
    the OTHER one is the one running away. Neither direction is checked
    by any existing code path before this function."""
    frac_eq = R_eq / r_max
    frac_pol = R_pol / r_max
    overflow_eq = frac_eq > (1.0 - tol)
    overflow_pol = frac_pol > (1.0 - tol)
    return {
        "frac_eq": frac_eq,
        "frac_pol": frac_pol,
        "overflow_eq": overflow_eq,
        "overflow_pol": overflow_pol,
        "overflow": overflow_eq or overflow_pol,
    }


def density_peak_location(rho, r, theta):
    """(r, theta) do maximo global de rho — deve ficar em r=0 se o centro
    continua sendo o ponto mais denso; usado para checar se a ancoragem em
    rho_c (fixada em r=0) ainda faz sentido fisico em campo forte.

    REVISADO (varredura da classe de bug de surface_radius, ver
    docs/teoria.md Sec 1.11): np.argmax(rho) e' ADEQUADO aqui, categoria
    diferente de erro. Isto busca um MAXIMO interior (discreto, exato por
    construcao — argmax nao interpola nada), nao um cruzamento de
    fronteira; nao ha' clip da EOS envolvido (rho e' estritamente positivo
    no interior, o clip so' afeta o exterior/vacuo) nem posicao sub-malha
    sendo procurada.

    NOTA DE ESTADO: esta funcao nao e' chamada por nenhuma pagina do
    dashboard ainda (gap ja' registrado antes desta varredura) — a
    migracao de pico de densidade em campo/rotacao fortes hoje so' e'
    inspecionada via o grafico plots.plot_density_profile() (Aba 1)."""
    i, j = np.unravel_index(np.argmax(rho), rho.shape)
    return r[i], theta[j], i, j
