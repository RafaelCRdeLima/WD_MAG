"""Valida poisson.solve_poisson contra o potencial de uma esfera uniforme (l=0 exato)."""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from poisson import solve_poisson, G_CONST


def test_uniform_sphere():
    R = 1.0e8       # cm
    rho0 = 1.0e6    # g/cm^3
    rmax = 4 * R

    # nr alto pq rho tem descontinuidade em r=R (esfera uniforme); o erro do
    # trapezio la' e' de 1a ordem em dr. Perfis estelares reais sao suaves
    # (H->0 continuo na superficie), entao a convergencia no SCF sera melhor.
    nr, ntheta = 801, 65
    r = np.linspace(0, rmax, nr)
    theta = np.linspace(0, np.pi, ntheta)

    rho = np.where(r[:, None] <= R, rho0, 0.0) * np.ones((1, ntheta))

    Phi = solve_poisson(rho, r, theta, lmax=4)

    Phi_exact = np.where(
        r <= R,
        -2 * np.pi * G_CONST * rho0 * (R**2 - r**2 / 3),
        -(4.0 / 3.0) * np.pi * G_CONST * rho0 * R**3 / np.where(r > 0, r, np.nan),
    )
    Phi_exact[0] = -2 * np.pi * G_CONST * rho0 * R**2  # r=0

    # compara ao longo do equador (theta = pi/2), ignorando o polo/centro
    j_eq = ntheta // 2
    rel_err = np.abs(Phi[1:, j_eq] - Phi_exact[1:]) / np.abs(Phi_exact[1:])
    max_err = np.nanmax(rel_err)
    print(f"erro relativo maximo (Poisson, esfera uniforme): {max_err:.3e}")
    assert max_err < 1.5e-2, f"erro relativo {max_err} acima da tolerancia"


if __name__ == "__main__":
    test_uniform_sphere()
    print("OK")
