"""Loop de iteracao SCF. Ver plano_wd_magnetizada.md secao 4 e docs/teoria.md.

REFATORACAO (rotacao + toroidal autoconsistente, D0): UM UNICO loop, com
termos plugaveis (scf/terms/) em vez de uma funcao por combinacao de
fisica -- cada ingrediente (rotacao, poloidal, toroidal autoconsistente)
contribui um termo aditivo ao MESMO Bernoulli:

    H + Phi - C_rot(varpi) - M_pol(u) + M_tor(rho*varpi^2) = C

rearranjado para o loop:

    H = C - Phi + C_rot(varpi) + M_pol(u) - M_tor(rho*varpi^2)

rotation/poloidal/toroidal = None desliga o termo correspondente
(contribuicao e energia identicamente zero). Ver scf/terms/__init__.py
para o protocolo completo. poloidal e toroidal sao MUTUAMENTE EXCLUSIVOS
(fora de escopo: campo misto autoconsistente -- D6 continua sendo a unica
forma de ter poloidal+toroidal juntos, via imposicao a posteriori).

Rotacao e poloidal sao termos EXPLICITOS: seu potencial, numa dada
iteracao, nao depende do rho_novo sendo resolvido (rotacao depende so' da
malha; poloidal depende de u, resolvido a partir do rho da iteracao
ANTERIOR via Grad-Shafranov -- o mesmo defasamento que ja existia antes
desta refatoracao). O toroidal autoconsistente e' IMPLICITO: seu potencial
depende do proprio rho_novo, no MESMO ponto de malha -- o passo 9 deixa de
ser inversao direta da EOS e vira uma busca de raiz por ponto
(_solve_rho_implicit).

NOTA DE PROJETO (historico, nao estava no plano original, descoberta ao
implementar): a receita original do plano fixa a constante de Bernoulli C
(e, no caso com campo, a amplitude do campo) impondo H=0 em pontos da
SUPERFICIE (polo e equador) — o metodo classico de Hachisu para estrelas em
rotacao. Ao testar isso para o caso esferico sem campo, a iteracao de
Picard (substituicao direta rho -> Phi -> H -> rho_novo) e' LINEARMENTE
INSTAVEL para esta EOS perto do limite de Chandrasekhar: o indice
politropico efetivo tende a n=3 no nucleo, onde a massa fica quase
independente da densidade central e o raio despenca: resolver para C a
partir do raio inverte um mapa quase singular. Sub-relaxacao (em rho ou em
H) NAO resolve isso — testado com nr=150/300/600 e rho_c de 1e6 a 1e10
g/cm^3, instavel em toda a faixa.

A parametrizacao usada aqui evita o problema por construcao: rho_c e os
TERMOS (rotation, poloidal, toroidal) sao as entradas INDEPENDENTES,
nenhuma condicao de superficie entra:

    C = H(rho_c) + Phi_c - C_rot(0) - M_pol(u_c) + M_tor(0)

onde Phi_c, u_c sao os valores no centro (r=0). C_rot(0)=0 (varpi=0),
M_tor(0)=0 (varpi=0) e M_pol(u_c)=0 (u_c=0 no centro, ver nota antiga) --
TODOS os termos se anulam no centro por construcao, entao C reduz sempre a
H(rho_c)+Phi_c, com ou sem qualquer termo ligado.

**Ressalva (rotacao ou campo fortes):** a configuracao pode migrar o pico
de densidade para fora do centro (evacuacao polar/equatorial) -- a mesma
patologia ja observada em campo poloidal forte, pelo mesmo motivo: a
ancora em rho_c(r=0) perde sentido fisico ali. Ver
diagnostics.density_peak_location() e o aviso correspondente no dashboard.

rotation=poloidal=toroidal=None reduz exatamente ao caso sem campo/rotacao
(Fase 1, V1, validado a 0.78% do limite de Chandrasekhar em
tests/test_scf_v1.py) -- e essa reducao tem que ser BIT A BIT identica ao
comportamento anterior a' refatoracao (V-R0, tests/test_regression_v0.py).

NOTA (poloidal): gradshafranov.py tinha um bug de expoente na funcao de
Green (corrigido — ver nota de projeto no topo daquele modulo) que inflava
o campo por um fator ~7000x. Com o fix, para rho_c=1e9 g/cm^3, R~3e8 cm,
nr=161, ntheta=65, lmax=16, a sequencia em k0 e' fisicamente sa' ate'
k0~1.6e-12 (M~1.50 Msun, VE~6.6e-4). Acima disso (testado ate' k0~2.3e-12,
M~2.02 Msun) o SCF ainda converge numericamente mas VE ultrapassa 1e-3 (o
V3 do plano) e NAO melhora com resolucao — e' terminacao fisica da
sequencia (evacuacao polar), nao artefato numerico. Nao compare esses M
contra o M_max~1.9 Msun de Bera & Bhattacharya (2014): aquele numero e' o
maximo sobre TODO o plano (rho_c, k0), nao ao longo de uma fatia de rho_c
fixo. A varredura de verdade precisa escanear rho_c tambem — Aba 2 do
dashboard.
"""

import numpy as np
from scipy.optimize import brentq

from eos import enthalpy, density_of_enthalpy, B_of_mu_e
from poisson import solve_poisson


def initial_guess(r, theta, rho_c, r_target):
    """Perfil inicial tipo politropo n=3: rho_c (1 - (r/r_target)^2)^3, clipado em 0."""
    profile = np.clip(1 - (r / r_target) ** 2, 0.0, None) ** 3
    return rho_c * profile[:, None] * np.ones((1, len(theta)))


def total_mass(rho, r, theta):
    """M = integral rho dV em coordenadas esfericas (r, theta), simetria azimutal."""
    integrand = rho * r[:, None] ** 2 * np.sin(theta)[None, :]
    over_theta = np.trapezoid(integrand, theta, axis=1)
    return 2 * np.pi * np.trapezoid(over_theta, r)


def _solve_rho_implicit(RHS, toroidal, r, theta, mu_e, rho_c):
    """Resolve, PONTO A PONTO, H(rho) + M_tor(rho*varpi^2) = RHS(r,theta)
    para rho, via busca de raiz protegida (scipy.optimize.brentq).

    O lado esquerdo e' monotonico CRESCENTE em rho: H(rho) e' crescente
    (propriedade padrao da EOS, ja' usada em toda parte deste projeto) e
    M_tor(rho*varpi^2) tambem e' crescente em rho para varpi>0, m>=1 (ver
    scf/terms/toroidal_sc.py) -- portanto a soma e' crescente, e a busca de
    raiz e' segura (no maximo uma raiz, colchetavel).

    RHS <= 0 (o valor do lado esquerdo em rho=0, sempre 0 por
    normalizacao: H(0)=0, M_tor(0)=0) => vacuo (rho=0), sem busca de raiz
    -- ver o prompt de fisica ("se o lado direito ficar abaixo do valor em
    rho=0, a solucao e' rho=0").

    Nao vetorizado (loop Python por ponto de malha) -- aceito
    deliberadamente (custo extra ainda cabe em segundos por iteracao, ver
    prompt de fisica: "vetorize ou aceite o custo").
    """
    varpi2 = (r[:, None] * np.sin(theta)[None, :]) ** 2
    B = B_of_mu_e(mu_e)

    def f(rho_val, rhs_val, varpi2_val):
        if rho_val <= 0.0:
            H_val = 0.0
        else:
            x = (rho_val / B) ** (1.0 / 3.0)
            H_val = enthalpy(x, mu_e)
        Mtor_val = toroidal.M_tor(rho_val * varpi2_val)
        return H_val + Mtor_val - rhs_val

    nr, ntheta = RHS.shape
    rho_out = np.zeros((nr, ntheta))
    rho_hi_seed = max(rho_c * 10.0, 1.0)

    for i in range(nr):
        for j in range(ntheta):
            rhs_ij = RHS[i, j]
            if rhs_ij <= 0.0:
                continue  # vacuo, rho_out ja' e' 0.0
            varpi2_ij = varpi2[i, j]
            hi = rho_hi_seed
            f_hi = f(hi, rhs_ij, varpi2_ij)
            n_expand = 0
            while f_hi < 0.0 and n_expand < 60:
                hi *= 4.0
                f_hi = f(hi, rhs_ij, varpi2_ij)
                n_expand += 1
            if f_hi < 0.0:
                # nao colchetou depois de expandir bastante -- defensivo,
                # nao deveria acontecer dado o crescimento sem limite de H
                continue
            rho_out[i, j] = brentq(f, 0.0, hi, args=(rhs_ij, varpi2_ij),
                                    xtol=1e-6 * max(rho_c, 1.0), rtol=1e-12)
    return rho_out


def hachisu_scf(rho, r, theta, rho_c, rotation=None, poloidal=None, toroidal=None,
                 mu_e=2.0, lmax=16, tol=1e-8, max_iter=200, verbose=False,
                 track_virial=False):
    """SCF parametrizado por (rho_c, termos) -- ver nota de projeto no topo
    do modulo.

    rho: chute inicial (nr, ntheta) — ver initial_guess()
    rho_c: densidade central alvo (g/cm^3), fixa H(r=0) = H(EOS(rho_c))
    rotation: scf.terms.rotation.Rotation ou None (sem rotacao)
    poloidal: scf.terms.poloidal.Poloidal ou None (sem campo poloidal)
    toroidal: scf.terms.toroidal_sc.ToroidalSC ou None (sem campo toroidal
        autoconsistente) -- MUTUAMENTE EXCLUSIVO com poloidal (D6/escopo)
    mu_e: peso molecular medio por eletron (Y_e = 1/mu_e)
    track_virial: se True, calcula o erro virial (diagnostics.virial_error,
        generalizado com 2T) a cada iteracao e retorna em "ve_history".

    Retorna dict com rho, Phi, u (0 se poloidal=None, para compatibilidade
    com leitores antigos), H, C, iterations, converged, history,
    ve_history, e os objetos rotation/poloidal/toroidal (ja' atualizados
    ao estado final -- usados pelo dashboard/diagnosticos para extrair
    Br/Bth/Bphi/T/E_pol/E_tor sem recalcular nada).
    """
    if poloidal is not None and toroidal is not None:
        raise ValueError(
            "poloidal e toroidal autoconsistente sao mutuamente exclusivos "
            "(fora de escopo: campo misto autoconsistente, D6) — use "
            "toroidal.impose_toroidal() depois da convergencia para o twisted torus"
        )

    x_c = (rho_c / B_of_mu_e(mu_e)) ** (1.0 / 3.0)
    H_c = enthalpy(x_c, mu_e)

    history = []
    ve_history = []
    converged = False
    it = 0
    C = 0.0
    H = np.zeros_like(rho)
    Phi = np.zeros_like(rho)

    for it in range(max_iter):
        Phi = solve_poisson(rho, r, theta, lmax=lmax)

        if poloidal is not None:
            poloidal.update(rho, r, theta, lmax=lmax)

        explicit_potential = np.zeros_like(rho)
        if rotation is not None:
            explicit_potential = explicit_potential + rotation.potential(r, theta)
        if poloidal is not None:
            explicit_potential = explicit_potential + poloidal.potential(r, theta)

        # ancora no centro: todos os termos se anulam em r=0 (varpi=0),
        # entao C reduz sempre a H_c + Phi_c (ver nota de projeto no topo)
        C = H_c + Phi[0, 0] - explicit_potential[0, 0]

        RHS = C - Phi + explicit_potential

        if toroidal is None:
            H = RHS
            rho_new = density_of_enthalpy(H, mu_e)
        else:
            rho_new = _solve_rho_implicit(RHS, toroidal, r, theta, mu_e, rho_c)
            varpi2 = (r[:, None] * np.sin(theta)[None, :]) ** 2
            H = RHS - toroidal.M_tor(rho_new * varpi2)

        rel_delta = np.max(np.abs(rho_new - rho)) / rho_c
        history.append(rel_delta)

        if track_virial:
            import diagnostics as _diag
            ve_dict = _diag.virial_error_terms(rho_new, Phi, H, r, theta, mu_e,
                                                rotation=rotation, poloidal=poloidal,
                                                toroidal=toroidal)
            ve_history.append(ve_dict["VE"])

        rho = rho_new

        if verbose:
            print(f"iter {it:3d}  rel_delta={rel_delta:.3e}")

        if rel_delta < tol:
            converged = True
            break

    u = poloidal.u if poloidal is not None and poloidal.u is not None else np.zeros_like(rho)

    return {"rho": rho, "Phi": Phi, "u": u, "H": H, "C": C,
            "iterations": it + 1, "converged": converged, "history": history,
            "ve_history": ve_history,
            "rotation": rotation, "poloidal": poloidal, "toroidal": toroidal}
