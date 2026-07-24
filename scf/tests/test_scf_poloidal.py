"""Regressao para o SCF poloidal (k0!=0). Ver nota de projeto em scf.py sobre
o bug de gradshafranov.py (corrigido) que inflava o campo por ~7000x.

Nao testa o limite da sequencia (isso muda com rho_c, R, resolucao — fica
para a Aba 2 do dashboard); so' garante que, no regime de campo fraco, o
SCF poloidal converge, reduz suavemente ao caso sem campo, e fecha o
virial (V3)."""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from scf import hachisu_scf, initial_guess, total_mass
from terms.poloidal import Poloidal
import diagnostics as diag

M_SUN = 1.989e33


def test_weak_field_close_to_nonmagnetic():
    rho_c = 1e9
    R_guess = 3.0e8
    nr, ntheta = 161, 65
    r = np.linspace(0, 1.3 * R_guess, nr)
    theta = np.linspace(0, np.pi, ntheta)
    rho0 = initial_guess(r, theta, rho_c, R_guess)

    result0 = hachisu_scf(rho0, r, theta, rho_c, lmax=16, tol=1e-8, max_iter=200)
    result_k = hachisu_scf(rho0, r, theta, rho_c, poloidal=Poloidal(k0=1e-14),
                            lmax=16, tol=1e-8, max_iter=200)

    assert result0["converged"] and result_k["converged"]

    M0 = total_mass(result0["rho"], r, theta)
    Mk = total_mass(result_k["rho"], r, theta)
    assert abs(Mk - M0) / M0 < 1e-3, "campo fraco nao deveria mudar a massa apreciavelmente"

    W = diag.gravitational_energy(result_k["rho"], result_k["Phi"], r, theta)
    Pi = diag.pressure_integral(result_k["H"], r, theta)
    Br, Bth = diag.poloidal_field(result_k["u"], r, theta)
    _, _, Emag = diag.magnetic_energies(Br, Bth, np.zeros_like(Br), r, theta)
    VE = abs(W + 3 * Pi + Emag) / abs(W)
    print(f"VE (k0=1e-14) = {VE:.3e}")
    assert VE < 1e-3, f"VE={VE} acima do V3 do plano"


if __name__ == "__main__":
    test_weak_field_close_to_nonmagnetic()
    print("OK")
