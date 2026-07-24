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


def test_high_lmax_large_domain_uniform_sphere():
    """REGRESSAO (achado investigando VE em K alto no ramo toroidal
    autoconsistente, dominio grande -- ver docs/teoria.md Sec 6): a forma
    antiga de _inner_terms/_outer_terms formava r'^(l+2) e r^-(l+1)
    SEPARADAMENTE. Para r~1e8-1e9 cm (escala de ana branca) e l=48,
    r'^50 ~ 1e400-1e450 -- transborda float64 (max ~1.8e308) mesmo para
    um raio estelar comum, nao so' em malhas artificialmente grandes.
    l_max=48 nunca tinha sido exercitado neste projeto (default e' 16)
    ate' essa investigacao achar o NaN. Esfera uniforme e' l=0 puro
    (rho nao depende de theta), entao o resultado exato NAO deve mudar
    com l_max -- checa que o solver nao transborda/NaN em l_max=48 E que
    a solucao l=0 continua correta nessa configuracao."""
    R = 3.0e8
    rho0 = 1.0e6
    rmax = 4 * R
    nr, ntheta = 801, 33
    r = np.linspace(0, rmax, nr)
    theta = np.linspace(0, np.pi, ntheta)
    rho = np.where(r[:, None] <= R, rho0, 0.0) * np.ones((1, ntheta))

    Phi = solve_poisson(rho, r, theta, lmax=48)
    assert np.all(np.isfinite(Phi)), "Phi has non-finite values at lmax=48 on a stellar-scale domain"

    Phi_exact = np.where(
        r <= R,
        -2 * np.pi * G_CONST * rho0 * (R**2 - r**2 / 3),
        -(4.0 / 3.0) * np.pi * G_CONST * rho0 * R**3 / np.where(r > 0, r, np.nan),
    )
    Phi_exact[0] = -2 * np.pi * G_CONST * rho0 * R**2

    j_eq = ntheta // 2
    rel_err = np.abs(Phi[1:, j_eq] - Phi_exact[1:]) / np.abs(Phi_exact[1:])
    max_err = np.nanmax(rel_err)
    print(f"erro relativo maximo (Poisson, lmax=48, dominio estelar): {max_err:.3e}")
    assert max_err < 1.5e-2, f"erro relativo {max_err} acima da tolerancia"


if __name__ == "__main__":
    test_uniform_sphere()
    test_high_lmax_large_domain_uniform_sphere()
    print("OK")
