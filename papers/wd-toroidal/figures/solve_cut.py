"""Solve the configuration plotted in the meridional-cut figure and cache
its fields.

Run once:
    scf/.venv/bin/python3 papers/wd-toroidal/figures/solve_cut.py

The target is the M = 2 Msun crossing of the toroidal sequence at
rho_c = 1e9 g/cm^3 -- the point marked with a star in the M(B) figure --
so the two figures show the same star.

Physics is not reimplemented here (R1): the solve goes through
sweep_worker._solve_toroidal_certified(), the same certified path the
sweep uses (continuation from K=0, domain grown until frac_pol <= 0.2,
Delta r at the validated ratio), and the diagnostics through
diagnostics.virial_error_terms(). This script only chooses the target,
checks the result against the acceptance gates, and crops and stores the
fields so the figure can be redrawn without solving again.
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
from terms.toroidal_sc import ToroidalSC             # noqa: E402

RHO_C = 1.0e9
MU_E = 2.0
M_TOR = 1.0
# investigations/rho_c_1e9_M2_configurations.csv, the interpolated
# M = 2 Msun crossing of this sequence.
K_TOR = 3.245e-3

OUT = HERE / "cut_rhoc1e9.npz"


def main():
    R_guess = r_guess(RHO_C)
    print(f"rho_c = {RHO_C:.3e}  K = {K_TOR:.6g}  R_guess = {R_guess:.4e} cm")

    result, r, theta, overflow = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=R_guess, K_tor=K_TOR, m_tor_sc=M_TOR,
        rotation=None, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=16,
        tol=1e-8, max_iter=200,
    )
    if result is None:
        raise SystemExit("SCF did not converge along the continuation path")

    rho, Phi, H = result["rho"], result["Phi"], result["H"]
    toroidal = ToroidalSC(K=K_TOR, m=M_TOR)
    ve = diag.virial_error_terms(rho, Phi, H, r, theta, MU_E,
                                 rotation=None, poloidal=None,
                                 toroidal=toroidal)
    Bphi = ve["Bphi"]
    M_sun = units.g_to_msun(scf_mod.total_mass(rho, r, theta))
    R_eq, R_pol = diag.equatorial_polar_radii(H, r, theta)

    scalars = {
        "rho_c": RHO_C, "mu_e": MU_E, "K_tor": K_TOR, "m_tor": M_TOR,
        "M_Msun": M_sun,
        "VE": ve["VE"],
        "E_tor_over_W": ve["E_mag"] / abs(ve["W"]),
        "B_tor_max_G": float(np.max(np.abs(Bphi))),
        "R_eq_cm": float(R_eq), "R_pol_cm": float(R_pol),
        "frac_pol": overflow["frac_pol"], "frac_eq": overflow["frac_eq"],
        "Nr": len(r), "Ntheta": len(theta), "domain_cm": float(r[-1]),
        "iterations": result["iterations"],
    }
    print(json.dumps(scalars, indent=2))

    # The gates this project applies to any number it reports. Fail loudly
    # rather than quietly drawing an uncertified star.
    if scalars["VE"] >= 1e-3:
        raise SystemExit(f"VE = {scalars['VE']:.3e} is not certified")
    if scalars["frac_pol"] > 0.2:
        raise SystemExit(f"frac_pol = {scalars['frac_pol']:.3f} > 0.2")

    # Only the inner region is drawn; storing the whole 10-20 R_guess
    # domain would be mostly vacuum.
    keep = r <= 1.35 * max(R_eq, R_pol)
    np.savez_compressed(
        OUT, r=r[keep], theta=theta, rho=rho[keep], Bphi=Bphi[keep],
        scalars=json.dumps(scalars),
    )
    print("wrote", OUT.name, f"({OUT.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
