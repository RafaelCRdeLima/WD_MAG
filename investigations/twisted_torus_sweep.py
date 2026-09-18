"""The trade-off DIARIO 15 says was never drawn: how super-Chandrasekhar can a
star be and still hold a confined toroidal field?

twisted_torus_beta.py answered one point -- M = 2.0072 Msun at rho_c = 1e9 --
and found beta_min = 6e-3, set by the TOROIDAL, not the poloidal: identical to
three digits across three decades of poloidal amplitude. The reason is
structural. A twisted torus needs beta(psi) to vanish on every line that
reaches the vacuum, so the toroidal field lives only inside the closed-line
region; squeezing the energy that supports the mass into that fraction of the
volume is what breaks beta, and it breaks it before any poloidal is added.

That makes the question a search, with three knobs that were held fixed:

  rho_c, K_tor  -- how much mass the field must support (the real trade-off)
  m             -- Fujisawa localisation of the POLOIDAL, which sets the
                   closed region and therefore where the toroidal may live
  zeta          -- the shape of beta(u) inside that region

Per (rho_c, K_tor) this solves the star once, then for each (m, zeta):

  1. beta_min with the CONFINED TOROIDAL ALONE. If that is already below 1,
     no poloidal can rescue the point and the row says so -- this is the cheap
     gate, and DIARIO 15 predicts it fails at M = 2 Msun.
  2. where it passes, the largest poloidal amplitude that keeps beta_min >= 1,
     and the B_t/B_p that amplitude buys. Comparable is what Braithwaite &
     Spruit require; B_t/B_p near 1 with beta_min >= 1 is the target.

One SLURM array task per (rho_c, K_tor). Run one point locally with
  scf/.venv/bin/python3 investigations/twisted_torus_sweep.py --index 0
and the whole grid with cluster/scicom/job_tt_sweep.sh.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (REPO / "scf", REPO / "dashboard"):
    sys.path.insert(0, str(p))

warnings.filterwarnings("ignore")

import diagnostics as diag                       # noqa: E402
import eos                                       # noqa: E402
import scf as scf_mod                            # noqa: E402
import units                                     # noqa: E402
from gradshafranov import solve_gradshafranov    # noqa: E402
from seed import r_guess                         # noqa: E402
from sweep_worker import _solve_toroidal_certified   # noqa: E402
from terms.toroidal_sc import ToroidalSC         # noqa: E402

sys.path.insert(0, str(HERE))
from export_ct_model import confined_flux        # noqa: E402

MU_E = 2.0
M_TOR = 1.0
B_C = 4.414e13
LMAX = 16
RHO_STAR = 1.0e6           # same cut as export_ct_model.report()
RHO_AMBIENT = 2.0e4

# The grid. K_tor as a fraction of 3.245e-3, the M = 2 Msun crossing at 1e9.
# It reaches down to 0.05 because the shape knobs turned out to be nearly
# exhausted: over the whole (m, zeta) grid at K_frac = 1, beta_min moves
# only from 3.9e-3 to 1.2e-2, a factor of 3 against the 85 that are
# missing. E_tor goes as K_tor^2, so beta_min = 1 needs K_frac ~ 0.1 and
# the frontier has to be bracketed from below, not approached from 2 Msun.
RHO_C_LIST = (5.0e8, 1.0e9, 2.0e9, 3.0e9, 5.0e9)
K_FRAC_LIST = (0.05, 0.08, 0.11, 0.15, 0.20, 0.30, 0.45, 0.65, 1.00,
               1.20)
K_TOR_REF = 3.245e-3
M_LIST = (0.0, -0.5, -1.0, -1.5, -2.0)      # poloidal localisation
ZETA_LIST = (0.5, 1.0, 1.5, 2.0, 3.0)       # beta(u) shape

GRID = [(rc, kf) for rc in RHO_C_LIST for kf in K_FRAC_LIST]


def surface_flux(u, H, th):
    vals = []
    for j in range(len(th)):
        ins = H[:, j] > 0.0
        if ins.any():
            vals.append(u[np.flatnonzero(ins)[-1], j])
    return max(vals) if vals else 0.0


def beta_min_of(P_gas, star, Br, Bth, Bphi):
    P_mag = (Br**2 + Bth**2 + Bphi**2) / (8.0 * np.pi)
    return float((P_gas[star] / np.maximum(P_mag[star], 1e-300)).min())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", type=int, default=None,
                    help=f"point in the (rho_c, K_frac) grid, 0..{len(GRID)-1}")
    ap.add_argument("--rho-c", type=float, default=None,
                    help="central density, off-grid (needs --k-frac)")
    ap.add_argument("--k-frac", type=float, default=None,
                    help="K_tor as a fraction of the 2 Msun value, off-grid")
    ap.add_argument("--tag", default=None,
                    help="name for the off-grid output, e.g. hi0")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--nr", type=int, default=129)
    a = ap.parse_args()

    if a.rho_c is not None:
        # Off-grid point. The frontier rises with rho_c -- 1.3227 at 5e8 to
        # 1.4186 at 5e9 -- and neutronization for mu_e = 2 only sets in at
        # 1.94e10, so whether the frontier crosses 1.44 above the grid has to
        # be measured, not extrapolated.
        if a.k_frac is None:
            raise SystemExit("--rho-c needs --k-frac")
        rho_c, k_frac = a.rho_c, a.k_frac
        label = a.tag or f"{rho_c:.0e}_{k_frac:.2f}"
        out = a.out or HERE / f"tt_sweep_hi_{label}.csv"
    else:
        if a.index is None or not 0 <= a.index < len(GRID):
            raise SystemExit(f"--index must be in 0..{len(GRID)-1}")
        rho_c, k_frac = GRID[a.index]
        label = f"{a.index:03d}"
        out = a.out or HERE / f"tt_sweep_{a.index:03d}.csv"
    K_tor = k_frac * K_TOR_REF

    print(f"[{label}] rho_c = {rho_c:.3e}, K_tor = {K_tor:.6g} "
          f"({k_frac:.2f} of the 2 Msun value)", flush=True)

    res, r, th, _ = _solve_toroidal_certified(
        rho_c=rho_c, R_guess=r_guess(rho_c), K_tor=K_tor, m_tor_sc=M_TOR,
        rotation=None, mu_e=MU_E, Nr_base=a.nr, Ntheta=a.nr, lmax=LMAX,
        tol=1e-8, max_iter=200)
    if res is None:
        print(f"[{label}] the toroidal solve did not converge", flush=True)
        out.write_text(f"# index={label} rho_c={rho_c:.6e} "
                       f"K_tor={K_tor:.6e} DID NOT CONVERGE\n")
        return

    rho, Phi, H = res["rho"], res["Phi"], res["H"]
    M = units.g_to_msun(scf_mod.total_mass(rho, r, th))
    W = abs(diag.gravitational_energy(rho, Phi, r, th))
    varpi = r[:, None] * np.sin(th)[None, :]
    P_gas = eos.pressure(eos.x_of_enthalpy(np.maximum(H, 0.0), mu_e=MU_E))
    star = rho > RHO_STAR
    inside = H > 0.0
    zero = np.zeros_like(rho)

    Bphi_sc = ToroidalSC(K=K_tor, m=M_TOR).B_phi(rho, varpi)
    _, target, _ = diag.magnetic_energies(zero, zero, Bphi_sc, r, th)

    P_amb = eos.pressure(eos.x_of_enthalpy(eos.enthalpy(
        (RHO_AMBIENT / eos.B_of_mu_e(MU_E)) ** (1.0 / 3.0), mu_e=MU_E),
        mu_e=MU_E))
    B_amb_max = float(np.sqrt(8.0 * np.pi * P_amb))

    print(f"[{label}] M = {M:.4f} Msun, |W| = {W:.4e}, "
          f"E_tor/|W| = {target / W:.4f}", flush=True)

    rows = []
    for m in M_LIST:
        try:
            u_shape, _ = confined_flux(rho, r, th, varpi, m, lmax=LMAX)
        except SystemExit as e:
            print(f"[{label}]   m = {m}: {e}", flush=True)
            continue
        u_s = surface_flux(u_shape, H, th)
        closed = inside & (u_shape > u_s)
        f_vol = (diag.volume_integral(closed.astype(float), r, th)
                 / diag.volume_integral(inside.astype(float), r, th))
        Br1, Bth1 = diag.poloidal_field(u_shape, r, th)
        unit_pol = float(np.hypot(Br1, Bth1).max())

        for zeta in ZETA_LIST:
            u_norm = max(u_shape.max() - u_s, 1e-300)
            w = np.where(closed, (u_shape - u_s) / u_norm, 0.0)
            shape = np.where(varpi > 0,
                             np.power(w, zeta) / np.maximum(varpi, 1e-30), 0.0)
            _, E_unit, _ = diag.magnetic_energies(zero, zero, shape, r, th)
            if E_unit <= 0 or not np.isfinite(E_unit):
                continue
            Bphi = np.sqrt(target / E_unit) * shape
            Bphi_max = float(np.abs(Bphi).max())

            # GATE: the toroidal alone, no poloidal at all
            b_tor_only = beta_min_of(P_gas, star, zero, zero, Bphi)
            if b_tor_only < 1.0:
                rows.append((m, zeta, f_vol, M, target / W, b_tor_only,
                             np.nan, np.nan, np.nan, np.nan,
                             Bphi_max / B_C, "toroidal alone fails"))
                continue

            # it passes: how much poloidal fits under beta_min = 1?
            lo, hi = 1e6, 1e15
            for _ in range(50):
                mid = np.sqrt(lo * hi)
                Br, Bth = diag.poloidal_field(u_shape * (mid / unit_pol), r, th)
                if beta_min_of(P_gas, star, Br, Bth, Bphi) >= 1.0:
                    lo = mid
                else:
                    hi = mid
            Br, Bth = diag.poloidal_field(u_shape * (lo / unit_pol), r, th)
            E_pol, _, _ = diag.magnetic_energies(Br, Bth, zero, r, th)
            bpol = np.hypot(Br, Bth)
            bp_ext = diag.surface_dipolarity(bpol, H, r, th)["B_pole"]
            btot = float(np.sqrt(Br**2 + Bth**2 + Bphi**2).max())
            beta_amb = (B_amb_max / max(bp_ext, 1e-300)) ** 2
            verdict = ("TWISTED TORUS" if Bphi_max / lo < 10.0
                       else "toroidal-dominated")
            rows.append((m, zeta, f_vol, M, target / W, b_tor_only,
                         Bphi_max / lo, target / max(E_pol, 1e-99),
                         bp_ext, beta_amb, btot / B_C, verdict))
            print(f"[{label}]   m={m:5.2f} zeta={zeta:4.2f}  "
                  f"B_t/B_p={Bphi_max / lo:8.2f}  E_t/E_p="
                  f"{target / max(E_pol, 1e-99):9.2f}  {verdict}", flush=True)

    with out.open("w") as f:
        f.write(f"# index={label} rho_c={rho_c:.6e} K_tor={K_tor:.6e} "
                f"k_frac={k_frac:.2f} M={M:.6f} Msun W={W:.6e} "
                f"E_tor_over_W={target / W:.6f}\n")
        f.write(f"# beta_min over rho > {RHO_STAR:.0e}; ambient "
                f"{RHO_AMBIENT:.0e} holds {B_amb_max:.4e} G at beta=1\n")
        f.write("m,zeta,closed_vol_frac,M_msun,Etor_over_W,beta_min_tor_only,"
                "Bt_over_Bp,Etor_over_Epol,B_pole_G,beta_ambient,"
                "Btot_over_Bc,verdict\n")
        for row in rows:
            f.write(",".join(f"{v:.6e}" if isinstance(v, float) else str(v)
                             for v in row) + "\n")
    print(f"[{label}] wrote {out.name} ({len(rows)} rows)", flush=True)


if __name__ == "__main__":
    main()
