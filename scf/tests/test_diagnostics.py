"""Valida diagnostics.surface_dipolarity(): B_pole via valor direto
(regularizado no eixo por poloidal_field) concorda com uma extrapolacao
quadratica independente a partir de pontos fora do eixo, a razao
B_pole/B_eq converge para ~2 (quase-dipolo no regime de campo fraco), e a
sanidade B_surf_max < B_central se mantem — motivado por um bug real
onde B_r no eixo era zerado por engano (ver poloidal_field)."""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from scf import hachisu_scf, initial_guess
from terms.poloidal import Poloidal
import diagnostics as diag


def _run(ntheta, rho_c=1e9, R_guess=3.0e8, k0=1e-13, nr=161):
    r = np.linspace(0, 1.3 * R_guess, nr)
    theta = np.linspace(0, np.pi, ntheta)
    rho0 = initial_guess(r, theta, rho_c, R_guess)
    result = hachisu_scf(rho0, r, theta, rho_c, poloidal=Poloidal(k0=k0),
                          lmax=16, tol=1e-8, max_iter=200)
    assert result["converged"]
    return result, r, theta


def test_surface_dipolarity_consistency_and_convergence():
    diag_by_ntheta = {}
    for ntheta in (65, 129):
        result, r, theta = _run(ntheta)
        rho, u, H = result["rho"], result["u"], result["H"]
        Br, Bth = diag.poloidal_field(u, r, theta)
        Bpol = np.sqrt(Br ** 2 + Bth ** 2)
        d = diag.surface_dipolarity(Bpol, H, r, theta)
        B_central = np.sqrt(Br[1, 0] ** 2 + Bth[1, 0] ** 2)
        diag_by_ntheta[ntheta] = d

        print(f"ntheta={ntheta}  B_pole={d['B_pole']:.6e}  "
              f"B_pole_extrap={d['B_pole_extrapolated']:.6e}  "
              f"dipolarity={d['dipolarity']:.4f}  B_surf_max={d['B_surf_max']:.4e}  "
              f"B_central={B_central:.4e}")

        # direct (axis-regularized) vs independent quadratic extrapolation
        rel_diff = abs(d["B_pole"] - d["B_pole_extrapolated"]) / d["B_pole"]
        assert rel_diff < 5e-3, f"direct/extrapolated B_pole disagree by {rel_diff:.2%}"

        # sanity: field is stronger in the interior than at the surface
        # for this k0=const source (peaks well inside the star)
        assert d["B_surf_max"] < B_central

        # near-dipole regime (weak field): dipolarity should be close to
        # (not necessarily exactly) 2
        assert 1.9 < d["dipolarity"] < 2.05

    # convergence: dipolarity should move monotonically closer to 2 as
    # resolution increases (weak-field/near-dipole configuration)
    dip_65 = diag_by_ntheta[65]["dipolarity"]
    dip_129 = diag_by_ntheta[129]["dipolarity"]
    assert abs(dip_129 - 2.0) < abs(dip_65 - 2.0)


if __name__ == "__main__":
    test_surface_dipolarity_consistency_and_convergence()
    print("OK")
