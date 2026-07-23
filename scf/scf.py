"""Loop de iteracao SCF. Ver plano_wd_magnetizada.md secao 4.

NOTA DE PROJETO (nao estava no plano original, descoberta ao implementar):
a receita original do plano fixa a constante de Bernoulli C (e, no caso com
campo, a amplitude do campo) impondo H=0 em pontos da SUPERFICIE (polo e
equador) — o metodo classico de Hachisu para estrelas em rotacao. Ao testar
isso para o caso esferico sem campo, a iteracao de Picard (substituicao
direta rho -> Phi -> H -> rho_novo) e' LINEARMENTE INSTAVEL para esta EOS
perto do limite de Chandrasekhar: o indice politropico efetivo tende a n=3
no nucleo, onde a massa fica quase independente da densidade central e o
raio despenca: resolver para C a partir do raio inverte um mapa quase
singular. Sub-relaxacao (em rho ou em H) NAO resolve isso — testado com
nr=150/300/600 e rho_c de 1e6 a 1e10 g/cm^3, instavel em toda a faixa.

A parametrizacao usada aqui evita o problema por construcao: (rho_c, k0)
sao as duas entradas INDEPENDENTES (densidade central e amplitude do campo
poloidal f(u)=k0), nenhuma condicao de superficie entra.

    C = H(rho_c) + Phi_c - M(u_c)

onde Phi_c, u_c sao os valores no centro (r=0). Como omega=r*sin(theta) se
anula no centro, u_c=0 sempre e M(u_c)=0 — o termo fica na formula por
generalidade/clareza, nao porque contribua.

k0=0 reduz exatamente ao caso sem campo (Fase 1, V1, validado a 0.78% do
limite de Chandrasekhar em tests/test_scf_v1.py).

NOTA (poloidal, k0!=0): gradshafranov.py tinha um bug de expoente na funcao
de Green (corrigido — ver nota de projeto no topo daquele modulo) que
inflava o campo por um fator ~7000x. Com o fix, para rho_c=1e9 g/cm^3,
R~3e8 cm, nr=161, ntheta=65, lmax=16, a sequencia em k0 e' fisicamente sa'
ate' k0~1.6e-12 (M~1.50 Msun, VE~6.6e-4). Acima disso (testado ate' k0~2.3e-12,
M~2.02 Msun) o SCF ainda converge numericamente mas VE ultrapassa 1e-3 (o
V3 do plano) e NAO melhora com resolucao (testado lmax 16->48 e malha
129^2->385^2, VE fica em ~1.2-1.6e-3, nao cai) — e' terminacao fisica da
sequencia (evacuacao polar: o pico de rho migra do polo para fora do eixo,
Rpol/Req cai a 0.61), nao artefato numerico. Nao compare esses M contra o
M_max~1.9 Msun de Bera & Bhattacharya (2014): aquele numero e' o maximo
sobre TODO o plano (rho_c, k0), nao ao longo de uma fatia de rho_c fixo —
rho_c=1e9 aqui da' M(k0=0)=1.39 Msun, ja abaixo do limite de Chandrasekhar
sem campo (1.44), entao essa fatia especifica nem comeca no lugar certo.
A varredura de verdade (V2) precisa escanear rho_c tambem — Aba 2 do
dashboard.
"""

import numpy as np

from eos import enthalpy, density_of_enthalpy, B_of_mu_e
from poisson import solve_poisson
from gradshafranov import solve_gradshafranov


def initial_guess(r, theta, rho_c, r_target):
    """Perfil inicial tipo politropo n=3: rho_c (1 - (r/r_target)^2)^3, clipado em 0."""
    profile = np.clip(1 - (r / r_target) ** 2, 0.0, None) ** 3
    return rho_c * profile[:, None] * np.ones((1, len(theta)))


def total_mass(rho, r, theta):
    """M = integral rho dV em coordenadas esfericas (r, theta), simetria azimutal."""
    integrand = rho * r[:, None] ** 2 * np.sin(theta)[None, :]
    over_theta = np.trapezoid(integrand, theta, axis=1)
    return 2 * np.pi * np.trapezoid(over_theta, r)


def hachisu_scf(rho, r, theta, rho_c, k0=0.0, mu_e=2.0, lmax=16, tol=1e-8,
                 max_iter=200, verbose=False, track_virial=False):
    """SCF parametrizado por (rho_c, k0) — ver nota de projeto no topo do modulo.

    rho: chute inicial (nr, ntheta) — ver initial_guess()
    rho_c: densidade central alvo (g/cm^3), fixa H(r=0) = H(EOS(rho_c))
    k0: amplitude do campo poloidal f(u)=k0 (D6). k0=0 desliga o campo.
    mu_e: peso molecular medio por eletron (Y_e = 1/mu_e)
    track_virial: se True, calcula o erro virial (diagnostics.virial_error)
        a cada iteracao e retorna em "ve_history" — custa integrais extras
        por iteracao, desligado por padrao. Usado pelo dashboard (Aba 1,
        grafico de convergencia do virial) — mora aqui, nao la, por R1.

    Retorna dict com rho, Phi, u, H, C, iterations, converged, history,
    ve_history (vazio se track_virial=False).
    """
    x_c = (rho_c / B_of_mu_e(mu_e)) ** (1.0 / 3.0)
    H_c = enthalpy(x_c, mu_e)

    u = np.zeros_like(rho)
    history = []
    ve_history = []
    converged = False
    it = 0
    C = 0.0
    H = np.zeros_like(rho)

    for it in range(max_iter):
        Phi = solve_poisson(rho, r, theta, lmax=lmax)

        if k0 != 0.0:
            omega2 = (r[:, None] * np.sin(theta)[None, :]) ** 2
            source = -4 * np.pi * omega2 * rho * k0
            u = solve_gradshafranov(source, r, theta, lmax=lmax)
        else:
            u = np.zeros_like(rho)

        M_u = k0 * u
        C = H_c + Phi[0, 0] - M_u[0, 0]
        H = C - Phi + M_u
        rho_new = density_of_enthalpy(H, mu_e)

        rel_delta = np.max(np.abs(rho_new - rho)) / rho_c
        history.append(rel_delta)

        if track_virial:
            import diagnostics as _diag
            Br, Bth = _diag.poloidal_field(u, r, theta)
            Bphi_zero = np.zeros_like(Br)
            VE, _, _, _ = _diag.virial_error(rho_new, Phi, H, Br, Bth, Bphi_zero, r, theta, mu_e)
            ve_history.append(VE)

        rho = rho_new

        if verbose:
            print(f"iter {it:3d}  rel_delta={rel_delta:.3e}")

        if rel_delta < tol:
            converged = True
            break

    return {"rho": rho, "Phi": Phi, "u": u, "H": H, "C": C,
            "iterations": it + 1, "converged": converged, "history": history,
            "ve_history": ve_history}
