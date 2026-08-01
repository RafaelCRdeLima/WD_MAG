"""Can a BAROTROPIC rotating white dwarf reach 2 Msun?

This is the question that decides whether a Castro relaxation run is worth
building, and it has to be asked barotropically.

The non-barotropic models cannot be handed to Castro. The code's EOS is ztwd:
P = P(rho) with mu_e fixed at 2. The extra pressure a non-barotropic
construction uses to hold its mass is exactly what ztwd discards, so the star
arrives out of equilibrium by precisely the amount of support under test and
falls. That failure mode is not hypothetical -- this project already hit it,
building under ztwd and evolving under gamma_law, and the star collapsed from
t = 0 with a field-free control collapsing identically.

So the initial condition has to be in genuine equilibrium under the same EOS
that will evolve it. What it may keep is differential rotation, which is the
standard super-Chandrasekhar route and which we have never tested
dynamically, plus a field.

Rotation is swept first and alone, because rotation is what delivers the
mass: the barotropic field is bounded by the ceiling measured in the
companion paper, E_mag/|W| <= 0.018, so it cannot be the main support. If
rotation alone does not clear 2 Msun in equilibrium and below mass shedding,
there is nothing to evolve and the Castro run should not be built.

Run:  scf/.venv/bin/python3 investigations/rotating_barotropic_scan.py
"""

import itertools
import multiprocessing as mp
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for _p in (REPO / "scf", REPO / "dashboard"):
    sys.path.insert(0, str(_p))
warnings.filterwarnings("ignore")

import diagnostics as diag                                  # noqa: E402
import scf as scf_mod                                       # noqa: E402
import units                                                # noqa: E402
from seed import r_guess                                    # noqa: E402
from sweep_worker import _solve_toroidal_certified          # noqa: E402
from terms.rotation import Rotation                         # noqa: E402

NR, NTH, LMAX = 129, 129, 16
RHO_C = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0e9
SHED_GATE = 0.95

# Omega_c quoted against the Keplerian frequency of the non-rotating star, so
# the numbers mean something before the solve. A is the j-constant shear
# length in units of R_eq; inf is rigid rotation.
OMEGA_FRAC = (0.0, 0.5, 0.8, 1.0, 1.2, 1.5)
# A = inf (rigid) failed to converge at every Omega >= 0.5, so the
# near-rigid end is approached through large finite A instead.
A_FRAC = (1.0, 2.0, 3.0, 4.0)
K_TOR_LIST = (0.0,)

KEPLER_REF = 0.0
OUT = HERE / "rotating_barotropic_scan.csv"


def solve_one(args):
    om_frac, a_frac, k_tor, om_kep = args
    try:
        R0 = r_guess(RHO_C)
        A = a_frac * R0 if np.isfinite(a_frac) else float("inf")
        rot = Rotation(om_frac * om_kep, A) if om_frac > 0 else None

        res, r, th, _ = _solve_toroidal_certified(
            rho_c=RHO_C, R_guess=R0, K_tor=k_tor, m_tor_sc=1.0,
            rotation=rot, mu_e=2.0, Nr_base=NR, Ntheta=NTH, lmax=LMAX,
            tol=1e-8, max_iter=400)
        if res is None:
            return dict(om_frac=om_frac, a_frac=a_frac, k_tor=k_tor,
                        status="no convergence")

        rho, Phi = res["rho"], res["Phi"]
        M = units.g_to_msun(scf_mod.total_mass(rho, r, th))
        W = diag.gravitational_energy(rho, Phi, r, th)
        T = rot.energy(rho, r, th)["T"] if rot is not None else 0.0

        jeq = len(th) // 2
        kk = np.flatnonzero(rho[:, jeq] > 0)
        R_eq = float(r[kk[-1]]) if kk.size else 0.0
        if kk.size and rot is not None:
            om2 = float(np.atleast_1d(rot.Omega(np.array([R_eq])))[0]) ** 2
            gr = abs(float(np.gradient(Phi[:, jeq], r)[kk[-1]]))
            shed = om2 * R_eq / max(gr, 1e-30)
        else:
            shed = 0.0

        return dict(om_frac=om_frac, a_frac=a_frac, k_tor=k_tor, status="ok",
                    M=M, TW=T / max(abs(W), 1.0), shed=shed, R_eq=R_eq,
                    ve=float(res.get("ve", float("nan"))))
    except Exception as exc:
        return dict(om_frac=om_frac, a_frac=a_frac, k_tor=k_tor,
                    status=f"{type(exc).__name__}: {exc}")


def main():
    R0 = r_guess(RHO_C)
    res, r, th, _ = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=R0, K_tor=0.0, m_tor_sc=1.0, rotation=None,
        mu_e=2.0, Nr_base=NR, Ntheta=NTH, lmax=LMAX, tol=1e-8, max_iter=400)
    M0 = scf_mod.total_mass(res["rho"], r, th)
    R_eq0 = diag.equatorial_polar_radii(res["H"], r, th)[0]
    om_kep = float(np.sqrt(units.G_CONST * M0 / R_eq0 ** 3))
    print(f"background: M = {units.g_to_msun(M0):.4f} Msun, "
          f"R_eq = {R_eq0:.3e} cm, Omega_K = {om_kep:.4e} rad/s\n")

    grid = [a + (om_kep,) for a in
            itertools.product(OMEGA_FRAC, A_FRAC, K_TOR_LIST)
            if not (a[0] == 0.0 and a[1] != A_FRAC[0])]
    print(f"rotating barotropic scan: {len(grid)} points\n")
    print("  Om/Om_K   A/R     K_tor    status   M       T/|W|   shed    "
          "R_eq")

    with mp.Pool(min(11, mp.cpu_count())) as pool:
        rows = pool.map(solve_one, grid)

    for d in sorted(rows, key=lambda x: (x["om_frac"], x["a_frac"],
                                         x["k_tor"])):
        head = (f"  {d['om_frac']:5.2f}   {d['a_frac']:5.2f}  "
                f"{d['k_tor']:.1e}")
        if "M" in d:
            print(head + f"   ok     {d['M']:6.3f}  {d['TW']:6.4f}  "
                         f"{d['shed']:6.3f}  {d['R_eq']:.3e}")
        else:
            print(head + f"   {d['status']}")

    good = [d for d in rows if "M" in d and d["shed"] < SHED_GATE]
    if good:
        b = max(good, key=lambda d: d["M"])
        print(f"\nheaviest below the shedding limit: M = {b['M']:.4f} Msun "
              f"at Omega/Omega_K = {b['om_frac']}, A/R_eq = {b['a_frac']}, "
              f"K_tor = {b['k_tor']:.1e}")
        print(f"  T/|W| = {b['TW']:.4f}, shedding parameter "
              f"{b['shed']:.3f}, R_eq = {b['R_eq']:.3e} cm")
        over = [d for d in good if d["M"] >= 2.0]
        print(f"  models at or above 2 Msun: {len(over)} of {len(good)}")
        if not over:
            print("  -> rotation alone does not reach 2 Msun here; there is "
                  "nothing to hand to Castro")
    else:
        print("\nnothing below the shedding limit")

    with OUT.open("w") as f:
        f.write("om_frac,a_frac,k_tor,status,M_msun,T_over_W,shed,R_eq,ve\n")
        for d in rows:
            f.write(",".join(str(d.get(k, "")) for k in
                             ("om_frac", "a_frac", "k_tor", "status", "M",
                              "TW", "shed", "R_eq", "ve")) + "\n")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
