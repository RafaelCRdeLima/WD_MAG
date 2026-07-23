"""V1 (independente do SCF 2D): integracao radial direta da equilibrio
hidrostatico esferico (dH/dr=-Gm/r^2, dm/dr=4 pi r^2 rho(H)) para validar
EOS + gravidade contra o limite de Chandrasekhar, sem depender do solver
de Poisson por Legendre. Serve de referencia independente para semear e
checar o SCF 2D (scf.py)."""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from eos import enthalpy, x_of_enthalpy, density_of_enthalpy

G_CONST = 6.674e-8
M_SUN = 1.989e33


def integrate_star(rho_c, dr=1.0e4, r_max=1.0e10):
    x_c = (rho_c / 1.96e6) ** (1.0 / 3.0)
    H = enthalpy(x_c)
    m = 0.0
    r = dr / 2  # evita r=0 exato na primeira derivada
    rho = rho_c

    while H > 0 and r < r_max:
        dHdr = -G_CONST * m / r**2
        dmdr = 4 * np.pi * r**2 * rho

        H_new = H + dHdr * dr
        m_new = m + dmdr * dr
        r_new = r + dr

        rho = density_of_enthalpy(H_new)
        H, m, r = H_new, m_new, r_new

    return m, r  # M, R no ponto H=0


def test_chandrasekhar_limit_shooting():
    rho_cs = np.array([1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12])
    masses = []
    radii = []
    for rho_c in rho_cs:
        M, R = integrate_star(rho_c, dr=R_STEP(rho_c))
        masses.append(M / M_SUN)
        radii.append(R)
        print(f"rho_c = {rho_c:.1e} g/cm3  ->  M = {M / M_SUN:.5f} Msun  R = {R / 1e5:.1f} km")

    m_max = masses[-1]
    rel_err = abs(m_max - 1.44) / 1.44
    print(f"erro relativo ao limite de Chandrasekhar (mu_e=2): {rel_err:.3%}")

    # massa deve crescer monotonamente com rho_c (ramo estavel) e saturar perto de 1.44
    assert all(m2 >= m1 for m1, m2 in zip(masses, masses[1:]))
    assert rel_err < 0.01, f"M_max={m_max:.4f} Msun, erro {rel_err:.3%} acima de 1%"


def R_STEP(rho_c):
    # passo mais fino para estrelas mais compactas (gradientes mais acentuados)
    if rho_c < 1e8:
        return 2.0e5
    elif rho_c < 1e10:
        return 2.0e4
    else:
        return 2.0e3


if __name__ == "__main__":
    test_chandrasekhar_limit_shooting()
    print("OK")
