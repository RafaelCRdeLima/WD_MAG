"""One certified toroidal + rigidly rotating configuration.

Run:
    scf/.venv/bin/python3 papers/wd-toroidal/figures/solve_rotating.py [P_spin_s]

Same certified path as solve_cut.py (R1: no physics reimplemented here) --
continuation in K from 0, domain grown until frac_pol <= 0.2 -- with a
rigid-rotation term added. The spin period is the input; Omega_c follows
from it.

Reference point is the rho_c = 1e10 row of the paper's mass table
(K = 3e-3, no rotation, M = 2.0722 Msun), so the rotating result can be
read as a delta against a number already certified.

All four acceptance gates this project applies to a rotating
configuration are checked, and the script refuses to report a
configuration that fails any of them:
  VE < 1e-3                     virial balance
  frac_pol <= 0.2               the star fits in its domain
  T/|W| < 0.14                  secular non-axisymmetric instability
                                (Ostriker & Bodenheimer 1973)
  mass-loss ratio < 1           equatorial breakup
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
for p in (REPO / "scf", REPO / "dashboard"):
    sys.path.insert(0, str(p))

import diagnostics as diag        # noqa: E402
import scf as scf_mod             # noqa: E402
import units                      # noqa: E402
from seed import r_guess          # noqa: E402
from sweep_worker import _solve_toroidal_certified   # noqa: E402
from terms.rotation import Rotation                  # noqa: E402
from terms.toroidal_sc import ToroidalSC             # noqa: E402

RHO_C = 1.0e10
MU_E = 2.0
K_TOR = 3.0e-3
M_TOR = 1.0
P_SPIN_DEFAULT = 25.0    # s

OUT = HERE / "rotating_point.json"


def solve(p_spin):
    """p_spin = None solves the non-rotating reference, through the
    identical numerical path, so the two differ only by the rotation
    term and the mass difference is not partly a difference in setup."""
    omega_c = 0.0 if p_spin is None else 2.0 * np.pi / p_spin
    R_guess = r_guess(RHO_C)
    print(f"rho_c = {RHO_C:.3e}  K = {K_TOR:.6g}  "
          f"P = {'none' if p_spin is None else f'{p_spin:.4g} s'}  "
          f"Omega_c = {omega_c:.4f} rad/s")

    rotation = None if p_spin is None else Rotation(Omega_c=omega_c)
    result, r, theta, overflow = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=R_guess, K_tor=K_TOR, m_tor_sc=M_TOR,
        rotation=rotation, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=16,
        tol=1e-8, max_iter=200,
    )
    if result is None:
        raise SystemExit("SCF did not converge along the continuation path")

    rho, Phi, H = result["rho"], result["Phi"], result["H"]
    toroidal = ToroidalSC(K=K_TOR, m=M_TOR)
    ve = diag.virial_error_terms(rho, Phi, H, r, theta, MU_E,
                                 rotation=rotation, poloidal=None,
                                 toroidal=toroidal)
    R_eq, R_pol = diag.equatorial_polar_radii(H, r, theta)
    W = ve["W"]

    out = {
        "rho_c": RHO_C, "mu_e": MU_E, "K_tor": K_TOR, "m_tor": M_TOR,
        "P_spin_s": p_spin, "Omega_c_rad_s": omega_c, "rigid": True,
        "M_Msun": units.g_to_msun(scf_mod.total_mass(rho, r, theta)),
        "R_eq_km": units.cm_to_km(R_eq), "R_pol_km": units.cm_to_km(R_pol),
        "VE": ve["VE"],
        "T_over_W": ve["T"] / abs(W),
        "E_tor_over_W": ve["E_mag"] / abs(W),
        "B_tor_max_G": float(np.max(np.abs(ve["Bphi"]))),
        "mass_loss_ratio": diag.equatorial_mass_loss_ratio(
            Phi, rotation, r, theta, R_eq),
        "frac_pol": overflow["frac_pol"], "frac_eq": overflow["frac_eq"],
        "v_eq_equator_cm_s": omega_c * R_eq,
        "Nr": len(r), "Ntheta": len(theta), "iterations": result["iterations"],
    }
    print(json.dumps(out, indent=2))

    gates = {
        "VE < 1e-3": out["VE"] < 1e-3,
        "frac_pol <= 0.2": out["frac_pol"] <= 0.2,
        "T/|W| < 0.14": out["T_over_W"] < 0.14,
        "mass-loss ratio < 1": out["mass_loss_ratio"] < 1.0,
    }
    for name, ok in gates.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    if not all(gates.values()):
        raise SystemExit("configuration is not certified; not written")
    out["gates_passed"] = list(gates)
    return out


def main(p_spin):
    reference = solve(None)
    rotating = solve(p_spin)
    delta = (rotating["M_Msun"] - reference["M_Msun"]) / reference["M_Msun"]
    print(f"\nM: {reference['M_Msun']:.4f} -> {rotating['M_Msun']:.4f} Msun "
          f"({100 * delta:+.2f}%)")
    OUT.write_text(json.dumps(
        {"reference_no_rotation": reference, "rotating": rotating,
         "delta_M_fraction": delta}, indent=2) + "\n")
    print("wrote", OUT.name)


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else P_SPIN_DEFAULT)
