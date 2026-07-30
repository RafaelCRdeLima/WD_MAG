"""Can the magnetized SCF star be handed to Castro with div B = 0 exactly?

Run:  scf/.venv/bin/python3 investigations/vector_potential_export.py

This is the de-risking step for the export (Step 1): it answers, in Python,
whether the construction works, before anyone writes the C++ that does it
inside the Castro problem.

The requirement is not "write B into the initial condition" -- interpolating
B onto faces leaves a divergence that constrained transport will then carry
forever. It is to write a VECTOR POTENTIAL whose discrete curl on the
staggered mesh reproduces the field, which makes div B = 0 an identity of
the discretization rather than something to be checked and cleaned. The
existing seed machinery (seed_field.H) already works this way; this reuses
the idea for an SCF field instead of a random one.

The construction is analytic for an axisymmetric field:

  poloidal   u = varpi * A_phi          ->  A_phi = u / varpi
             (the flux function IS the potential, no work needed)

  toroidal   need poloidal A with (curl A)_phi = B_phi. Taking A_varpi = 0,
             (curl A)_phi = -dA_z/dvarpi, so
                 A_z(varpi, z) = - int_0^varpi B_phi(varpi', z) dvarpi'

so A = (A_phi e_phi) + (A_z e_z), and B = curl A recovers both components.

What is measured here: the discrete divergence of the reconstructed field,
and how much of the field's amplitude survives the mapping onto Castro's
Cartesian grid.
"""

import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (REPO / "scf", REPO / "dashboard"):
    sys.path.insert(0, str(p))

warnings.filterwarnings("ignore")

from gradshafranov import solve_gradshafranov            # noqa: E402
from scipy.interpolate import RegularGridInterpolator    # noqa: E402
from seed import r_guess                                 # noqa: E402
from sweep_worker import _solve_toroidal_certified       # noqa: E402
from terms.toroidal_sc import ToroidalSC                 # noqa: E402

RHO_C = 1.0e9
MU_E = 2.0
K_TOR = 3.245e-3      # the M = 2 Msun crossing
K0_POL = 1.0e-13      # the poloidal level that keeps the pair certified
N_CART = 64           # Castro's production resolution
HALF = 9.0e8          # cm; this star is inflated, R_pol = 5.7e8
N_MER = 513           # meridional grid for building A


def main():
    res, r, th, _ = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=K_TOR, m_tor_sc=1.0,
        rotation=None, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=16,
        tol=1e-8, max_iter=200)
    rho = res["rho"]
    omega2 = (r[:, None] * np.sin(th)[None, :]) ** 2
    u = solve_gradshafranov(-4.0 * np.pi * omega2 * rho * K0_POL, r, th, lmax=16)
    Bphi = ToroidalSC(K=K_TOR, m=1.0).B_phi(rho, np.sqrt(omega2))
    print(f"SCF field:  |B_phi|max = {np.abs(Bphi).max():.4e} G")

    # --- A on a meridional (varpi, z) grid -----------------------------
    vp = np.linspace(0.0, r[-1], N_MER)
    zz = np.linspace(-r[-1], r[-1], 2 * N_MER - 1)
    VP, ZZ = np.meshgrid(vp, zz, indexing="ij")
    RR = np.sqrt(VP**2 + ZZ**2)
    TT = np.arccos(np.clip(ZZ / np.maximum(RR, 1e-30), -1.0, 1.0))

    u_m = RegularGridInterpolator((r, th), u, bounds_error=False,
                                  fill_value=0.0)((RR, TT))
    bphi_m = RegularGridInterpolator((r, th), Bphi, bounds_error=False,
                                     fill_value=0.0)((RR, TT))
    A_phi = np.where(VP > 0, u_m / np.maximum(VP, 1e-30), 0.0)
    A_z = -np.concatenate(
        [np.zeros((1, VP.shape[1])),
         np.cumsum(0.5 * (bphi_m[1:] + bphi_m[:-1]) * np.diff(vp)[:, None],
                   axis=0)], axis=0)

    f_phi = RegularGridInterpolator((vp, zz), A_phi, bounds_error=False,
                                    fill_value=0.0)
    f_z = RegularGridInterpolator((vp, zz), A_z, bounds_error=False,
                                  fill_value=0.0)

    # --- A on the staggered Cartesian mesh, half-shift geometry --------
    dx = 2.0 * HALF / N_CART
    lo = -((N_CART + 1) / 2.0) * dx        # Sec 6.6: a cell CENTRE at r=0
    ctr = lo + dx * (np.arange(N_CART) + 0.5)
    fac = lo + dx * np.arange(N_CART + 1)

    def A_at(X, Y, Z, comp):
        vpg = np.sqrt(X**2 + Y**2)
        pts = np.stack([vpg.ravel(), Z.ravel()], axis=-1)
        ap = f_phi(pts).reshape(X.shape)
        az = f_z(pts).reshape(X.shape)
        inv = np.where(vpg > 0, 1.0 / np.maximum(vpg, 1e-30), 0.0)
        return {"x": -ap * Y * inv, "y": ap * X * inv, "z": az}[comp]

    # each component of A lives on the edges parallel to it
    Ax = A_at(*np.meshgrid(ctr, fac, fac, indexing="ij"), "x")
    Ay = A_at(*np.meshgrid(fac, ctr, fac, indexing="ij"), "y")
    Az = A_at(*np.meshgrid(fac, fac, ctr, indexing="ij"), "z")

    # B = curl A, landing on faces
    Bx = (np.diff(Az, axis=1) - np.diff(Ay, axis=2)) / dx
    By = (np.diff(Ax, axis=2) - np.diff(Az, axis=0)) / dx
    Bz = (np.diff(Ay, axis=0) - np.diff(Ax, axis=1)) / dx

    div = (np.diff(Bx, axis=0) + np.diff(By, axis=1) + np.diff(Bz, axis=2)) / dx
    bmax = max(np.abs(Bx).max(), np.abs(By).max(), np.abs(Bz).max())
    rel = np.abs(div).max() / (bmax / dx)

    print(f"\nreconstructed on {N_CART}^3, dx = {dx:.3e} cm:")
    print(f"  |div B| max, normalised by max|B|/dx = {rel:.3e}")
    print(f"  max|B| = {bmax:.4e} G, against {np.abs(Bphi).max():.4e} G in the SCF")
    print(f"  amplitude retained = {100 * bmax / np.abs(Bphi).max():.1f}%")
    print(f"  cells across the star (R_pol = 5.7e8 cm) = {2 * 5.7e8 / dx:.0f}")

    if rel > 1e-12:
        raise SystemExit(f"divergence too large: {rel:.3e}")
    print("\ndiv B is zero to machine precision: the curl construction works,")
    print("and the C++ side of the export is a transcription of it.")


if __name__ == "__main__":
    main()
