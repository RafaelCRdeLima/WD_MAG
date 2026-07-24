"""V-R2 (prompt de rotacao): rotacao diferencial j-constante alcanca
~2.2 Msun (Yoon & Langer 2005), MUITO alem do limite de rotacao rigida
(~1.5 Msun, V-R1, test_rotation.py). Ao contrario de V-R1, este alvo e'
alcancado longe do breakup Kepleriano (mass_loss_ratio fica bem abaixo de
1 em toda a sequencia) — a rotacao diferencial sustenta massa extra
deixando o nucleo girar mais rapido que o envelope, sem precisar levar o
equador ao limite de ruptura. Isso torna V-R2 numericamente mais limpo que
V-R1 (sem o "penhasco" perto do ponto critico)."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from scf import hachisu_scf, initial_guess, total_mass
from terms.rotation import Rotation
import diagnostics as diag

M_SUN = 1.989e33


def test_v_r2_differential_rotation_reaches_2p2_msun():
    rho_c = 1e10
    R_guess = 1.6e8
    nr, ntheta = 97, 49
    r = np.linspace(0, 2.2 * R_guess, nr)
    theta = np.linspace(0, np.pi, ntheta)
    rho_seed = initial_guess(r, theta, rho_c, R_guess)

    A = 0.3 * R_guess  # decay scale of the j-constant profile, in units of R_eq (documented in terms/rotation.py)
    omegas = [0.0, 8.0, 16.0, 24.0, 32.0, 35.0, 36.0, 37.0]
    masses = []
    last_mlr = 0.0
    for Omega_c in omegas:
        rotation = Rotation(Omega_c=Omega_c, A=A) if Omega_c > 0 else None
        result = hachisu_scf(rho_seed, r, theta, rho_c, rotation=rotation,
                              lmax=16, tol=1e-8, max_iter=400)
        assert result["converged"], f"Omega_c={Omega_c} (differential, A={A:.3e}) did not converge"
        rho_seed = result["rho"]
        M = total_mass(result["rho"], r, theta) / M_SUN
        masses.append(M)
        if rotation is not None:
            R_eq, _ = diag.equatorial_polar_radii(result["H"], r, theta)
            last_mlr = diag.equatorial_mass_loss_ratio(result["Phi"], rotation, r, theta, R_eq)
        print(f"Omega_c={Omega_c:5.1f}  M/Msun={M:.4f}  mass_loss_ratio={last_mlr:.4f}")

    assert all(m2 >= m1 - 1e-9 for m1, m2 in zip(masses, masses[1:])), \
        "mass should increase monotonically with Omega_c"

    M_final = masses[-1]
    # far from breakup -- this is genuinely differential-rotation support,
    # not a near-critical artifact (contrast with V-R1's rigid-rotation cliff)
    assert last_mlr < 0.5, f"mass_loss_ratio={last_mlr:.3f} closer to breakup than expected for this A"

    rel_err = abs(M_final - 2.2) / 2.2
    print(f"M_final={M_final:.4f} Msun, rel. error to 2.2 Msun target (Yoon & Langer 2005): {rel_err:.2%}")
    assert rel_err < 0.10, (
        f"M={M_final:.4f} Msun, {rel_err:.2%} away from the ~2.2 Msun differential-rotation "
        "target -- more than the 10% tolerance"
    )

    # sanity: differential rotation should reach well beyond the rigid limit
    # (~1.5 Msun, V-R1) at a comparable rho_c
    assert M_final > 1.7, "differential rotation should clear the rigid-rotation ceiling comfortably"


if __name__ == "__main__":
    test_v_r2_differential_rotation_reaches_2p2_msun()
    print("OK")
