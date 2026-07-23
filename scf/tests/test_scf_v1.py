"""V1 (plano_wd_magnetizada.md secao 4): sem campo, relacao M-R deve reproduzir
Chandrasekhar, M_max -> 1.44 Msun para mu_e=2 (Y_e=0.5), erro < 1% no ramo
ultrarrelativistico (rho_c alto). Usa a formulacao estavel do SCF (fixa
rho_c, ver nota de projeto em scf.py)."""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from scf import hachisu_scf, initial_guess, total_mass

M_SUN = 1.989e33  # g


def run_spherical(rho_c, R_guess, nr=300, ntheta=17):
    rmax = 1.3 * R_guess
    r = np.linspace(0, rmax, nr)
    theta = np.linspace(0, np.pi, ntheta)

    rho0 = initial_guess(r, theta, rho_c, R_guess)
    result = hachisu_scf(rho0, r, theta, rho_c, lmax=0, tol=1e-8, max_iter=100)
    M = total_mass(result["rho"], r, theta)
    return M, result


def test_chandrasekhar_limit():
    # (rho_c, R aproximado p/ dimensionar a malha) -- R nao precisa ser exato,
    # so' grande o suficiente pra conter a estrela com folga
    cases = [
        (1e7, 8.0e8),
        (1e8, 5.0e8),
        (1e9, 3.0e8),
        (1e10, 1.6e8),
        (1e11, 8.0e7),
        (1e12, 4.0e7),
    ]
    masses = []
    for rho_c, R_guess in cases:
        M, result = run_spherical(rho_c, R_guess)
        assert result["converged"], f"SCF nao convergiu para rho_c={rho_c:.1e}"
        masses.append(M / M_SUN)
        print(f"rho_c = {rho_c:.1e} g/cm3  ->  M = {M / M_SUN:.5f} Msun  "
              f"({result['iterations']} iter)")

    assert all(m2 >= m1 - 1e-6 for m1, m2 in zip(masses, masses[1:])), (
        "massa deveria crescer monotonamente com rho_c no ramo relativistico"
    )

    m_final = masses[-1]
    rel_err = abs(m_final - 1.44) / 1.44
    print(f"erro relativo ao limite de Chandrasekhar (mu_e=2): {rel_err:.3%}")
    assert rel_err < 0.01, f"M={m_final:.4f} Msun, erro {rel_err:.3%} acima de 1%"


if __name__ == "__main__":
    test_chandrasekhar_limit()
    print("OK")
