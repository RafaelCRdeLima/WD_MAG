"""What the 2 Msun configuration actually looks like inside.

Run:  scf/.venv/bin/python3 investigations/spindle_structure.py

R_pol/R_eq = 1.10 describes the surface where H = 0 and says the star is
mildly prolate. The mass distribution is nothing of the sort. Measured on the
converged SCF solution, the density ratio between the axis and the equator at
the SAME spherical radius is:

    r/R_eq    rho(axis)     rho(equator)   ratio
     0.10     8.93e8        2.15e8          4.2
     0.20     6.35e8        5.52e7         11.5
     0.50     1.02e8        4.36e6         23.4
     0.80     8.92e6        4.74e5         18.8

Along the equator the density falls by a factor 230 from centre to half
radius; along the axis it falls by 10. The configuration is a dense axial
column inside an evacuated equatorial region -- a spindle, not a flattened
star -- because B_phi = K rho varpi puts its magnetic pressure exactly where
rho and varpi are both large.

This reframes the Phase 1 collapses. Every configuration in the mass scan
shares this structure, since they share rho_c and the field law and differ
only in K, which is why the death time ignored E_tor/|W|: the character of
the object does not change across the scan, only its degree.

It also means R_pol/R_eq badly understates the deformation, and the mixed
paper quotes that ratio as the measure of it.
"""

import sys
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
for p in (REPO / "scf", REPO / "dashboard"):
    sys.path.insert(0, str(p))
warnings.filterwarnings("ignore")

import diagnostics as diag                       # noqa: E402
import scf as scf_mod                            # noqa: E402
import units                                     # noqa: E402
from seed import r_guess                         # noqa: E402
from sweep_worker import _solve_toroidal_certified   # noqa: E402

OUT = Path(__file__).resolve().parent / "spindle_structure.csv"
CASES = ((1.5155e-3, "1.50"), (2.3503e-3, "1.70"),
         (2.8274e-3, "1.85"), (3.245e-3, "2.01"))


def main():
    rows = []
    print("  M (Msun)  K_tor       r/R_eq  rho(axis)   rho(eq)     ratio")
    for K, label in CASES:
        res, r, th, _ = _solve_toroidal_certified(
            rho_c=1.0e9, R_guess=r_guess(1.0e9), K_tor=K, m_tor_sc=1.0,
            rotation=None, mu_e=2.0, Nr_base=129, Ntheta=129, lmax=16,
            tol=1e-8, max_iter=200)
        if res is None:
            print(f"  {label}: no convergence"); continue
        rho, H = res["rho"], res["H"]
        M = units.g_to_msun(scf_mod.total_mass(rho, r, th))
        R_eq, R_pol = diag.equatorial_polar_radii(H, r, th)
        ja, je = 0, int(np.argmin(np.abs(th - np.pi / 2)))
        for f in (0.2, 0.5):
            i = int(np.argmin(np.abs(r - f * R_eq)))
            a, e = rho[i, ja], rho[i, je]
            rows.append((M, K, R_pol / R_eq, f, a, e, a / max(e, 1e-30)))
            print(f"  {M:7.4f}  {K:.4e}  {f:5.2f}  {a:.3e}  {e:.3e}  "
                  f"{a / max(e, 1e-30):8.2f}")
        print(f"            R_pol/R_eq = {R_pol / R_eq:.4f}\n")

    with OUT.open("w") as f:
        f.write("# axis-to-equator density contrast at matched spherical "
                "radius, rho_c=1e9, mu_e=2\n")
        f.write("M_msun,K_tor,R_pol_over_R_eq,r_over_Req,rho_axis,rho_eq,"
                "ratio\n")
        for row in rows:
            f.write(",".join(f"{v:.6e}" for v in row) + "\n")
    print(f"wrote {OUT.name}")


if __name__ == "__main__":
    main()
