"""V-R3, V-R4 (prompt de rotacao) + checagem de sinal para o ramo toroidal
autoconsistente (scf/terms/toroidal_sc.py).

ACHADO IMPORTANTE sobre V-R4 (documentado aqui e reportado ao usuario, que
confirmou e localizou a causa exata -- ver scf/terms/toroidal_sc.py para a
versao completa): o prompt de fisica afirma a identidade

    int rho grad(M_tor).r dV = int B_phi^2/(8 pi) dV

mas a derivacao rigorosa (via divergencia do tensor de tensoes de Maxwell
-- ver test_maxwell_stress_sign_derivation abaixo) e a verificacao
numerica (converge limpo com a malha: 2.27% -> 1.00% -> 0.36% -> 0.14% em
nr=65->257) mostram que o sinal correto e' NEGATIVO:

    int rho grad(M_tor).r dV = - int B_phi^2/(8 pi) dV

CAUSA (nao e' rho aparecer dentro de M_tor -- essa e' a causa de a
inversao virar busca de raiz, uma questao de estrutura do algoritmo, sem
relacao com sinal): o prompt escreveu o Bernoulli mestre como
"H + Phi - C_rot - M_pol + M_tor = C", ou seja M_pol com sinal menos e
M_tor com sinal MAIS. Tomando o gradiente e comparando com a equacao de
momento, disso sai F_L,tor/rho = -grad(M_tor) (o oposto do lado poloidal,
F_L,pol/rho = +grad(M_pol) -- veja que o sinal trocado entre os dois vem
so' de como cada um foi escrito na equacao, um com + um com -). A
identidade do tensor de Maxwell int r.F_L dV = +int B^2/(8pi) dV e' GERAL
e nao depende de nenhuma convencao de Bernoulli. Substituindo
F_L,tor=-rho*grad(M_tor) nela e' que produz o sinal negativo. Ou seja: se
o prompt tivesse escrito "-M_tor" em vez de "+M_tor" no Bernoulli mestre,
a identidade do virial sairia com "+", nao "-" -- o prompt era
internamente inconsistente (uma convencao de sinal numa linha, a
identidade do virial escrita como se M_tor fosse potencial de forca na
outra). O codigo (com "+M_tor" no Bernoulli, sinal negativo no virial)
esta' certo E consistente consigo mesmo. O teste abaixo usa o sinal
CORRIGIDO (negativo), documentado aqui para quem for comparar com o
prompt original.
"""

import os
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from scf import hachisu_scf, initial_guess, total_mass
from terms.toroidal_sc import ToroidalSC
import diagnostics as diag

M_SUN = 1.989e33


def test_maxwell_stress_sign_derivation():
    """Verificacao simbolica (SymPy) de f_L = -rho*grad(Psi) para
    B_phi = K rho^m varpi^(2m-1), m=1,2,3/2 -- confirma a derivacao do
    docstring de toroidal_sc.py (residuo exatamente 0). Reproduzida aqui
    como teste automatizado, nao so' verificacao interativa."""
    r, theta = sp.symbols("r theta", positive=True, real=True)
    K = sp.symbols("K", positive=True, real=True)
    rho = sp.Function("rho", positive=True)(r, theta)

    for m_val in [1, 2, sp.Rational(3, 2)]:
        m = sp.nsimplify(m_val)
        varpi = r * sp.sin(theta)
        B_phi = K * rho**m * varpi ** (2 * m - 1)

        curlB_r = sp.diff(B_phi * sp.sin(theta), theta) / (r * sp.sin(theta))
        curlB_th = -sp.diff(r * B_phi, r) / r

        fL_r = sp.simplify((curlB_th * B_phi - 0) / (4 * sp.pi))
        fL_th = sp.simplify((0 - curlB_r * B_phi) / (4 * sp.pi))

        s = rho * varpi**2
        Psi = (m * K**2 / (4 * sp.pi * (2 * m - 1))) * s ** (2 * m - 1)
        gradPsi_r = sp.diff(Psi, r)
        gradPsi_th = sp.diff(Psi, theta) / r

        residual_r = sp.simplify(fL_r - (-rho * gradPsi_r))
        residual_th = sp.simplify(fL_th - (-rho * gradPsi_th))
        assert residual_r == 0, f"m={m}: f_L,r != -rho*dPsi/dr (residual {residual_r})"
        assert residual_th == 0, f"m={m}: f_L,theta != -rho*(1/r)dPsi/dtheta (residual {residual_th})"


def _run_toroidal(rho_c=1e9, R_guess=3.0e8, K=1.2e-3, m=1.0, nr=97, ntheta=49,
                   domain_factor=1.6, tol=1e-7, max_iter=200):
    r = np.linspace(0, domain_factor * R_guess, nr)
    theta = np.linspace(0, np.pi, ntheta)
    rho0 = initial_guess(r, theta, rho_c, R_guess)
    toroidal = ToroidalSC(K=K, m=m)
    result = hachisu_scf(rho0, r, theta, rho_c, toroidal=toroidal, lmax=16, tol=tol, max_iter=max_iter,
                          track_virial=True)
    return result, r, theta, toroidal


def test_sign_check_mass_increases_with_field():
    """Checagem fisica de sinal exigida pelo prompt: ligar o campo
    toroidal a rho_c fixo tem que AUMENTAR a massa (mesma logica usada
    para confirmar o sinal do termo poloidal originalmente)."""
    rho_c, R_guess = 1e9, 3.0e8
    result0, r, theta, _ = _run_toroidal(rho_c, R_guess, K=0.0)
    resultK, _, _, _ = _run_toroidal(rho_c, R_guess, K=1.2e-3)
    assert result0["converged"] and resultK["converged"]
    M0 = total_mass(result0["rho"], r, theta)
    MK = total_mass(resultK["rho"], r, theta)
    print(f"M(K=0)={M0/M_SUN:.4f} Msun   M(K=1.2e-3)={MK/M_SUN:.4f} Msun")
    assert MK > M0, "turning on the toroidal field at fixed rho_c should increase the mass -- sign check (see module docstring)"


def test_v_r3_prolate_deformation_virial_closes():
    result, r, theta, toroidal = _run_toroidal()
    assert result["converged"]
    rho, Phi, H = result["rho"], result["Phi"], result["H"]

    R_eq, R_pol = diag.equatorial_polar_radii(H, r, theta)
    print(f"R_eq={R_eq/1e5:.1f} km  R_pol={R_pol/1e5:.1f} km  R_pol/R_eq={R_pol/R_eq:.5f}")
    assert R_pol > R_eq, "pure toroidal field should produce PROLATE deformation (R_pol > R_eq)"

    ve = diag.virial_error_terms(rho, Phi, H, r, theta, 2.0, toroidal=toroidal)
    print(f"VE={ve['VE']:.4e}  E_tor={ve['E_tor']:.4e}  W={ve['W']:.4e}")
    assert ve["VE"] < 1e-3, f"VE={ve['VE']:.3e} above the V3 threshold -- virial should close for a converged equilibrium"


def test_v_r4_magnetic_virial_identity():
    """int rho grad(M_tor).r dV == -int B_phi^2/(8 pi) dV (sinal NEGATIVO
    -- ver docstring do modulo). Verifica convergencia com a malha para
    confirmar que e' uma identidade real (erro de discretizacao indo a
    zero), nao um numero que so' concorda por coincidencia numa resolucao."""
    rho_c, R_guess = 1e9, 3.0e8
    residuals = []
    for nr, ntheta in [(65, 33), (97, 49), (161, 81)]:
        result, r, theta, toroidal = _run_toroidal(rho_c, R_guess, nr=nr, ntheta=ntheta)
        assert result["converged"]
        rho = result["rho"]

        varpi = r[:, None] * np.sin(theta)[None, :]
        s = rho * varpi**2
        dM_ds = toroidal.dM_tor_ds(s)
        drho_dr = np.gradient(rho, r, axis=0)
        ds_dr = varpi**2 * drho_dr + rho * 2 * r[:, None] * np.sin(theta)[None, :] ** 2
        dM_dr = dM_ds * ds_dr
        LHS = diag.volume_integral(rho * r[:, None] * dM_dr, r, theta)

        Bphi = toroidal.B_phi(rho, varpi)
        E_tor = diag.volume_integral(Bphi**2 / (8 * np.pi), r, theta)

        rel = abs(LHS - (-E_tor)) / E_tor
        residuals.append(rel)
        print(f"nr={nr:4d} ntheta={ntheta:3d}  LHS={LHS:.6e}  -E_tor={-E_tor:.6e}  rel_diff={rel:.4e}")

    # o residuo tem que encolher com a malha (assinatura de identidade
    # real, nao coincidencia numa resolucao so')
    assert residuals[-1] < residuals[0], "residual should shrink with grid refinement"
    assert residuals[-1] < 5e-3, f"V-R4 residual {residuals[-1]:.2%} too large even at the finest tested mesh"


if __name__ == "__main__":
    test_maxwell_stress_sign_derivation()
    test_sign_check_mass_increases_with_field()
    test_v_r3_prolate_deformation_virial_closes()
    test_v_r4_magnetic_virial_identity()
    print("OK")
