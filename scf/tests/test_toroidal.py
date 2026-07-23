"""Valida toroidal.py: a razao de energia alvo e' atingida, e B_phi e'
continuo (vai a zero) na borda do toro u=u_c."""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from scf import hachisu_scf, initial_guess
import diagnostics as diag
import toroidal as tor


def _converged_poloidal():
    rho_c = 1e9
    R_guess = 3.0e8
    nr, ntheta = 161, 65
    r = np.linspace(0, 1.3 * R_guess, nr)
    theta = np.linspace(0, np.pi, ntheta)
    rho0 = initial_guess(r, theta, rho_c, R_guess)
    result = hachisu_scf(rho0, r, theta, rho_c, k0=1e-13, lmax=16, tol=1e-8, max_iter=200)
    assert result["converged"]
    return result, r, theta


def test_energy_ratio_target():
    result, r, theta = _converged_poloidal()
    u, rho = result["u"], result["rho"]

    target = 0.3
    Bphi, zeta, u_c = tor.solve_zeta_for_energy_ratio(u, rho, r, theta, target, m_tor=1)
    Br, Bth = diag.poloidal_field(u, r, theta)
    ratio_energy, ratio_amp = tor.bt_bp_ratios(Br, Bth, Bphi, r, theta)

    print(f"u_c={u_c:.3e}  zeta={zeta:.3e}  Bt/Bp(energia)={ratio_energy:.4f}  "
          f"Bt/Bp(amplitude)={ratio_amp:.4f}")
    assert abs(ratio_energy - target) / target < 1e-6

    frac = tor.closed_torus_volume_fraction(u, rho, r, theta, u_c)
    print(f"fracao de volume do toro: {frac:.4e}")
    assert 0 < frac < 1


def test_continuity_at_boundary():
    result, r, theta = _converged_poloidal()
    u, rho = result["u"], result["rho"]
    Bphi, u_c = tor.impose_toroidal(u, rho, r, theta, zeta=1.0, m_tor=1)

    # logo dentro (u ligeiramente > u_c) B_phi deve ser pequeno (continuidade);
    # logo fora (u <= u_c) e' exatamente zero por construcao
    du = 1e-3 * (np.max(u) - u_c)
    close_to_boundary = (u > u_c) & (u < u_c + du)
    if np.any(close_to_boundary):
        assert np.max(np.abs(Bphi[close_to_boundary])) < 1e-2 * np.max(np.abs(Bphi))
    assert np.all(Bphi[u <= u_c] == 0.0)


if __name__ == "__main__":
    test_energy_ratio_target()
    test_continuity_at_boundary()
    print("OK")
