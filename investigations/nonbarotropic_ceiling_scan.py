"""Does dropping barotropy beat the barotropic ceiling, and at what cost?

Run:  scf/.venv/bin/python3 investigations/nonbarotropic_ceiling_scan.py [nproc]

THE QUESTION, SHARPENED

The barotropic poloidal branch saturates at E_pol/|W| = 0.018, measured in this
project. Dropping barotropy is supposed to lift that. The previous scan could
not tell: it ran k0 from 2e-13 to 1.8e-12 and every point came back rejected,
because all three composition gates it used were reading the method's own
numerical scatter.

nonbarotropic_noise_floor.py settled what that scatter is. Run with NO field,
where mu_e must be exactly 2 and uniform, the solver returns a shell-averaged
spread of 0.0203 -- 2%, from differentiating twice along P -> dP/dr -> rho ->
mu_e. Every earlier gate compared against numbers of that size. The same
calibration showed the signal does climb out of the floor, at 1.54x by
k0 = 5e-13 and 8.98x by 1e-12, exactly where the mass starts to rise.

So this scan starts where the previous one stopped, and states composition in
MULTIPLES OF THE FLOOR rather than in absolute spread.

WHAT IT MEASURES

For each configuration: the mass, E_pol/|W| against the barotropic ceiling of
0.018, the peak field against B_c, and the shell-averaged mu_e profile -- its
mean, its range, and its spread relative to the floor. A configuration is
counted as genuinely stratified only above STRATIFIED_MULTIPLE; below that its
composition structure is not resolved and the Ledoux question does not arise.

The answer sought is not a single number but a trade: how far above 0.018 the
poloidal support can be pushed before either the field exceeds B_c or the
composition the equilibrium demands stops being a white dwarf's.

Writes nonbarotropic_ceiling_scan.csv.
"""

import itertools
import os
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (REPO / "scf", REPO / "dashboard", REPO / "investigations"):
    sys.path.insert(0, str(p))
warnings.filterwarnings("ignore")

import diagnostics as diag                       # noqa: E402
import eos                                       # noqa: E402
import scf as scf_mod                            # noqa: E402
import units                                     # noqa: E402
from gradshafranov import solve_gradshafranov    # noqa: E402
from nonbarotropic_mass_scan import (grads, integrate_pressure,   # noqa: E402
                                     mu_e_from_P_rho)
from poisson import solve_poisson                # noqa: E402
from seed import r_guess                         # noqa: E402
from sweep_worker import _solve_toroidal_certified   # noqa: E402

LMAX, NR, NTH, N_OUTER = 16, 129, 129, 8
B_C = 4.414e13
BAROTROPIC_CEILING = 0.018      # measured: the poloidal SCF branch saturates here
NOISE_FLOOR = 0.0203            # measured at k0 = 0 by nonbarotropic_noise_floor.py
STRATIFIED_MULTIPLE = 3.0       # below this the composition structure is unresolved
MU_E_MIN, MU_E_MAX = 1.90, 2.20

# Starts where the previous scan stopped. B/B_c was 0.69 at k0 = 1.25e-12 on
# the field-free background, so there is room to roughly 1.8e-12 before the
# Landau limit -- but the star inflates as k0 rises, so the gate is applied to
# each configuration rather than assumed from the grid.
K0_LIST = tuple(np.geomspace(1.0e-12, 3.0e-12, 7))
ALPHA_LIST = (1.00, 1.15, 1.30)
RHO_C_LIST = (1.0e9, 3.0e9)

OUT = HERE / "nonbarotropic_ceiling_scan.csv"


def solve_one(args):
    k0, alpha, rho_c = args
    try:
        res, r, th, _ = _solve_toroidal_certified(
            rho_c=rho_c, R_guess=r_guess(rho_c), K_tor=0.0, m_tor_sc=1.0,
            rotation=None, mu_e=2.0, Nr_base=NR, Ntheta=NTH, lmax=LMAX,
            tol=1e-8, max_iter=200)
        if res is None:
            return dict(k0=k0, alpha=alpha, rho_c=rho_c, status="no background")
        rho0, Phi, H0 = res["rho"], res["Phi"], res["H"]
        M0 = units.g_to_msun(scf_mod.total_mass(rho0, r, th))
        varpi = r[:, None] * np.sin(th)[None, :]
        x0 = eos.x_of_enthalpy(np.maximum(H0, 0.0), 2.0)
        P_axis = alpha * np.where(H0[:, 0] > 0, eos.pressure(x0[:, 0]), 0.0)

        src = -4.0 * np.pi * varpi**2 * rho0 * k0
        u = solve_gradshafranov(src, r, th, lmax=LMAX)
        g = np.where(varpi > 0,
                     -src / (4 * np.pi * np.maximum(varpi, 1e-30)**2), 0.0)
        R_pol = diag.equatorial_polar_radii(H0, r, th)[1]
        seeds = np.linspace(r[1], 1.6 * R_pol, 260)

        rho = rho0
        for _ in range(N_OUTER):
            P = integrate_pressure(P_axis, Phi, g, u, r, th, seeds)
            P_r, _ = grads(P, r, th)
            Phi_r, _ = grads(Phi, r, th)
            u_r, _ = grads(u, r, th)
            safe = np.abs(Phi_r) > 1e-4 * np.abs(Phi_r).max()
            rn = np.where(safe, (g * u_r - P_r) / np.where(safe, Phi_r, 1.0), 0)
            for j in range(len(th)):
                k = np.flatnonzero(safe[:, j])
                if k.size:
                    rn[:k[0], j] = rn[k[0], j]
            rn = np.where(P > 0, np.maximum(rn, 0.0), 0.0)
            rho = 0.5 * rho + 0.5 * rn
            Phi = 0.6 * Phi + 0.4 * solve_poisson(rho, r, th, lmax=LMAX)

        M = units.g_to_msun(scf_mod.total_mass(rho, r, th))
        mu, _ = mu_e_from_P_rho(P, rho)
        inside = (rho > 1e-4 * rho_c) & (P > 0) & np.isfinite(mu)
        if inside.sum() < 500:
            return dict(k0=k0, alpha=alpha, rho_c=rho_c, status="empty", M=M)

        Rg = np.broadcast_to(r[:, None], mu.shape)
        edges = np.linspace(0.0, float(Rg[inside].max()), 25)
        prof = []
        for a, b in zip(edges[:-1], edges[1:]):
            m = inside & (Rg >= a) & (Rg < b)
            if m.sum() > 20:
                w = rho[m]
                prof.append(float(np.sum(mu[m] * w) / np.sum(w)))
        prof = np.array(prof)
        if prof.size < 5:
            return dict(k0=k0, alpha=alpha, rho_c=rho_c,
                        status="too few shells", M=M)

        spread = float((prof.max() - prof.min()) / prof.mean())
        over_floor = spread / NOISE_FLOOR
        stratified = over_floor >= STRATIFIED_MULTIPLE
        # Ledoux only where the structure is resolved; a profile at the floor
        # is homogeneous as far as this method can see, which is stable.
        if stratified:
            d = np.diff(prof)
            ledoux = float(np.maximum(d, 0).sum()
                           / max(prof.max() - prof.min(), 1e-12))
        else:
            ledoux = 0.0

        Br, Bth = diag.poloidal_field(u, r, th)
        E_pol, _, _ = diag.magnetic_energies(Br, Bth, np.zeros_like(rho), r, th)
        W = abs(diag.gravitational_energy(rho, Phi, r, th))
        B_max = float(np.hypot(Br, Bth).max())
        EpW = E_pol / max(W, 1.0)

        physical = (MU_E_MIN <= prof.min() and prof.max() <= MU_E_MAX)
        ok = physical and B_max < B_C and ledoux < 0.30
        return dict(k0=k0, alpha=alpha, rho_c=rho_c,
                    status="ok" if ok else "gate", M=M, M0=M0,
                    mu_mean=float(prof.mean()), mu_lo=float(prof.min()),
                    mu_hi=float(prof.max()), over_floor=over_floor,
                    ledoux=ledoux, EpW=EpW, beats_ceiling=EpW > BAROTROPIC_CEILING,
                    B_over_Bc=B_max / B_C)
    except Exception as exc:                       # noqa: BLE001
        return dict(k0=k0, alpha=alpha, rho_c=rho_c,
                    status=f"error: {type(exc).__name__}: {exc}")


def main():
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else max(1, os.cpu_count() - 1)
    grid = list(itertools.product(K0_LIST, ALPHA_LIST, RHO_C_LIST))
    print(f"ceiling scan: {len(grid)} points on {nproc} processes")
    print(f"noise floor {NOISE_FLOOR:.4f}; stratified above "
          f"{STRATIFIED_MULTIPLE}x it; barotropic ceiling "
          f"E_pol/|W| = {BAROTROPIC_CEILING}\n")
    with Pool(nproc) as pool:
        rows = pool.map(solve_one, grid)

    print("  rho_c    k0        alp  status  M      dM/M0   mu_e mean/range      "
          "xfloor  E_pol/|W|  beats  B/Bc")
    for d in sorted(rows, key=lambda x: (x["rho_c"], x["k0"], x["alpha"])):
        h = (f"  {d['rho_c']:.0e}  {d['k0']:.2e} {d['alpha']:.2f}  "
             f"{d['status']:6s}")
        if "mu_mean" in d:
            print(h + f" {d['M']:6.3f} {(d['M']-d['M0'])/d['M0']:+6.3f}  "
                      f"{d['mu_mean']:.3f} {d['mu_lo']:.3f}-{d['mu_hi']:.3f}  "
                      f"{d['over_floor']:6.2f}  {d['EpW']:9.5f}  "
                      f"{'YES' if d['beats_ceiling'] else ' no'}  "
                      f"{d['B_over_Bc']:5.2f}")
        else:
            print(h + f" {d.get('M', '')}")

    beat = [d for d in rows if d.get("beats_ceiling") and d.get("B_over_Bc", 9) < 1]
    if beat:
        b = max(beat, key=lambda d: d["EpW"])
        print(f"\nbeats the barotropic ceiling below B_c: E_pol/|W| = "
              f"{b['EpW']:.5f} at rho_c = {b['rho_c']:.0e}, k0 = {b['k0']:.2e}")
        print(f"  M = {b['M']:.4f} Msun ({100*(b['M']-b['M0'])/b['M0']:+.1f}% "
              f"over field-free), mu_e {b['mu_lo']:.3f}-{b['mu_hi']:.3f} "
              f"({b['over_floor']:.1f}x floor), status {b['status']}")
    else:
        print("\nnothing beats the barotropic ceiling below B_c")

    with OUT.open("w") as f:
        f.write(f"# non-barotropic ceiling scan; floor={NOISE_FLOOR}, "
                f"barotropic ceiling={BAROTROPIC_CEILING}, B_c={B_C:.4e}\n")
        f.write("k0,alpha,rho_c,status,M_msun,M0_msun,mu_mean,mu_lo,mu_hi,"
                "over_floor,ledoux,E_pol_over_W,beats_ceiling,B_over_Bc\n")
        for d in rows:
            f.write(",".join(str(d.get(k, "")) for k in
                             ("k0", "alpha", "rho_c", "status", "M", "M0",
                              "mu_mean", "mu_lo", "mu_hi", "over_floor",
                              "ledoux", "EpW", "beats_ceiling",
                              "B_over_Bc")) + "\n")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
