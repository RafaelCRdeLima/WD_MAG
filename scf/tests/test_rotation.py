"""V-R6 (prompt de rotacao): o T que entra no virial bate com uma integral
direta e independente de T. V-R6 PASSA e e' solido.

V-R1 (rotacao rigida ~1.5 Msun, Hachisu et al. 2012 / ~1.534 Msun em
Boshkayev et al. 2013) NAO E' VALIDADO -- ver test_v_r1_sequence_termination
abaixo. Registro preciso do estado (nao "ponto de virada fisico", que foi
uma caracterizacao errada descartada apos revisao):

A sequencia com rotacao rigida em rho_c=1e12 termina NUMERICAMENTE em
q=R_pol/R_eq≈0.932, Omega_c≈26.47, com razao de perda de massa (mass_loss_
ratio) = 0.138. A causa NAO esta' resolvida e a terminacao NAO e' fisica.
Dois sinais contra terminacao por perda de massa real:

  (i)  mass_loss_ratio=0.138 no ultimo ponto valido, quando terminacao por
       perda de massa exige mass_loss_ratio -> 1 por construcao do proprio
       diagnostico (Omega^2(R_eq)*R_eq / gravidade efetiva -> 1). Parar em
       0.138 e' parar a 14% do caminho, nao no breakup.
  (ii) o modelo de Roche, aplicavel a configuracoes centralmente
       condensadas (anas brancas perto do limite sao essencialmente n=3),
       preve q=2/3≈0.667 na perda de massa real. O valor obtido (0.932) e'
       so' ~7% de achatamento, contra os ~33% esperados.

V-R1 portanto NAO valida contra a literatura (1.48-1.53 Msun) -- a
validacao do termo rotacional repousa em V-R2 (rotacao diferencial,
test_differential_rotation.py, 0.40% de erro contra Yoon & Langer 2005,
longe do breakup). O mecanismo em si esta' confirmado por dois caminhos
independentes: (a) dM/M0 escala proporcionalmente a T/|W| com coeficiente
estavel ~3.0 ao longo de todo o range testado (ver
test_v_r1_sequence_termination); (b) V-R6 abaixo.

Suspeitos para a terminacao prematura (nao investigados a fundo -- fora do
caminho critico do projeto, que quer rotacao diferencial, nao rigida):
  - resolucao radial na camada externa da malha (~38 pontos ao longo do
    raio na configuracao testada)
  - a propria parametrizacao por Omega_c: Hachisu evita impor Omega_c perto
    do fim da sequencia justamente por isso, parametrizando por razao axial
    e resolvendo Omega^2 a partir dela. Um teste direto (bisseccao expondo
    razao axial como alvo, resolvendo Omega_c por fora) mostrou que NENHUM
    alvo de razao axial abaixo de ~0.932 e' alcancavel -- a bisseccao
    satura exatamente no mesmo ponto, nao converge para o alvo pedido.

Bug real encontrado e corrigido durante esta investigacao (nao a causa da
terminacao prematura, mas afetava a PRECISAO de toda leitura de raio do
projeto): diagnostics.equatorial_polar_radii() (via surface_radius())
recebia rho em vez de H. rho e' clipado exatamente a 0.0 alem da superficie
(eos.density_of_enthalpy), entao a "interpolacao linear" degenerava e
sempre devolvia um ponto de grade cru -- R_pol/R_eq saia em razoes de
inteiros pequenos (37/39, 38/39, ...), uma escada, nao um continuo. H e'
continuo e cruza zero de verdade entre pontos de grade; interpolar nele da'
a posicao real (~46 km de diferenca numa malha nr=65 -- um efeito grande).
Corrigido em diagnostics.py/toroidal.py (ver seus docstrings) -- a
bisseccao em razao axial passou a acertar alvos moderados EXATAMENTE
(0.980000000304478 para um alvo de 0.98) em vez de travar num degrau, mas
o ponto de saturacao em si sobreviveu ao conserto (so' mudou de
q≈0.9487 para q≈0.9323) -- ou seja, o bug era real mas nao e' a causa raiz
da terminacao prematura.

QUANDO RETOMAR: trocar o parametro de controle perto do fim da sequencia
(razao axial ou T/|W| como alvo, Omega_c resolvido por bisseccao externa —
exatamente a tecnica de Hachisu). Com equatorial_polar_radii() agora
interpolado corretamente, essa bisseccao tem uma chance real de funcionar
(antes, travava numa escada de malha). Nao implementado aqui -- fora do
caminho critico do projeto (rotacao diferencial, de onde vem o alvo de
2.2 Msun, ja' funciona e esta' validada por V-R2).
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from scf import hachisu_scf, initial_guess, total_mass
from terms.rotation import Rotation
import diagnostics as diag

M_SUN = 1.989e33


def test_v_r1_sequence_termination():
    """NAO valida V-R1 (ver docstring do modulo) -- documenta o estado
    numerico atual com asserts que capturam os fatos exatos, para que uma
    mudanca futura (ex: fixar a causa da terminacao prematura) seja
    visivel como uma mudanca nestes numeros, nao como silencio."""
    rho_c = 1e12
    R_guess = 5.576e7
    nr, ntheta = 129, 65
    r = np.linspace(0, 2.0 * R_guess, nr)
    theta = np.linspace(0, np.pi, ntheta)
    rho_seed = initial_guess(r, theta, rho_c, R_guess)

    omegas = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 26.0]
    masses = []
    T_over_W_list = []
    last_mlr = 0.0
    last_q = 1.0
    for Omega_c in omegas:
        rotation = Rotation(Omega_c=Omega_c, A=float("inf")) if Omega_c > 0 else None
        result = hachisu_scf(rho_seed, r, theta, rho_c, rotation=rotation,
                              lmax=16, tol=1e-8, max_iter=300)
        assert result["converged"], f"Omega_c={Omega_c} did not converge"
        rho_seed = result["rho"]
        rho, Phi, H = result["rho"], result["Phi"], result["H"]
        M = total_mass(rho, r, theta) / M_SUN
        masses.append(M)
        ve = diag.virial_error_terms(rho, Phi, H, r, theta, 2.0, rotation=rotation)
        T_over_W_list.append(ve["T"] / abs(ve["W"]))
        if rotation is not None:
            R_eq, R_pol = diag.equatorial_polar_radii(H, r, theta)
            last_mlr = diag.equatorial_mass_loss_ratio(Phi, rotation, r, theta, R_eq)
            last_q = R_pol / R_eq
        print(f"Omega_c={Omega_c:5.1f}  M/Msun={M:.4f}  T/|W|={T_over_W_list[-1]:.5f}  "
              f"q=R_pol/R_eq={last_q:.4f}  mass_loss_ratio={last_mlr:.4f}")

    # mechanism check (survives independently of the termination question):
    # dM/M0 tracks T/|W| with a stable proportionality coefficient -- this
    # is what actually validates the rotational term is doing real work
    M0 = masses[0]
    ratios = [(m - M0) / M0 / tw for m, tw in zip(masses[1:], T_over_W_list[1:])]
    print(f"(dM/M0)/(T/|W|) ratios: {[f'{x:.3f}' for x in ratios]}")
    assert max(ratios) - min(ratios) < 0.1, \
        "the dM/M0 vs T/|W| proportionality should stay stable if the term is working correctly"

    # the documented (not physical) termination point -- if these numbers
    # move, the module docstring's "suspects" discussion needs revisiting
    assert last_mlr < 0.3, "mass_loss_ratio at the last resolvable point should stay far from 1 (see docstring)"
    assert 0.9 < last_q < 0.95, "termination axis ratio drifted -- revisit the module docstring's numbers"


def test_v_r6_kinetic_energy_two_ways():
    """T que entra no virial (rotation.energy()['T']) tem que bater com
    uma integral DIRETA e independente (nao chama rotation.energy nem
    reusa nenhum codigo de terms/rotation.py alem do Omega(varpi) —
    escrita fresca aqui, exercitando a mesma formula por um caminho de
    codigo diferente, no espirito do V-R4)."""
    rho_c = 1e10
    R_guess = 1.6e8
    nr, ntheta = 97, 49
    r = np.linspace(0, 1.8 * R_guess, nr)
    theta = np.linspace(0, np.pi, ntheta)
    rho0 = initial_guess(r, theta, rho_c, R_guess)

    rotation = Rotation(Omega_c=3.0, A=float("inf"))
    result = hachisu_scf(rho0, r, theta, rho_c, rotation=rotation, lmax=16, tol=1e-8, max_iter=300)
    assert result["converged"]
    rho = result["rho"]

    T_from_term = rotation.energy(rho, r, theta)["T"]

    # integral direta, independente: T = (1/2) int rho Omega^2 varpi^2 dV,
    # Omega=Omega_c constante (rotacao rigida), calculada aqui sem chamar
    # rotation.energy()
    varpi = r[:, None] * np.sin(theta)[None, :]
    Omega_direct = np.full_like(varpi, rotation.Omega_c)
    integrand = rho * Omega_direct**2 * varpi**2
    T_direct = 0.5 * diag.volume_integral(integrand, r, theta)

    rel_diff = abs(T_from_term - T_direct) / abs(T_direct)
    print(f"T_from_term={T_from_term:.6e}  T_direct={T_direct:.6e}  rel_diff={rel_diff:.2e}")
    assert rel_diff < 1e-10


if __name__ == "__main__":
    test_v_r1_sequence_termination()
    test_v_r6_kinetic_energy_two_ways()
    print("OK")
