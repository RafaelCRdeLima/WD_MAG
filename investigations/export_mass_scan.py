"""Models for a dynamical-stability scan in mass.

Run:  scf/.venv/bin/python3 investigations/export_mass_scan.py

The 2 Msun configuration collapses in about 2.5 s of simulated time, and the
collapse does not scale with resolution -- so it is the configuration, not the
mesh. That is a negative result about one point. The useful question is where
the boundary is: at what mass does a magnetically supported white dwarf stop
being dynamically unstable?

That turns a failure into a number the collaboration can use. Jorge's scenario
needs the star to sit still for Myr while magnetic braking acts; a maximum mass
for dynamical survival is a hard constraint on it, and one that comes before
Tayler, before braking timescales, and before carbon ignition.

Why the boundary should exist at all. Under homologous contraction at fixed
mass and frozen flux, |W| goes as 1/R, E_mag goes as 1/R, and at Gamma = 4/3
the internal energy goes as 1/R too. All three scale identically, so magnetic
support buys mass capacity but no stability -- it cannot lift a marginal star
off the knife edge. What decides the sign is the deviation of Gamma from 4/3,
which grows as the central density falls below the fully relativistic regime.
Less magnetic support means a less inflated star at the same rho_c, hence a
stiffer effective index. So lower masses should be stabler, and somewhere
between 1.35 (field-free, measured stable) and 2.0 (measured unstable) the
behaviour turns over.

Each mass is reached by secant on K_tor at fixed rho_c = 1e9, then exported
exactly as investigations/export_phase1_models.py does, with the same gates.

Writes models/scan_M<mass>.txt and their manifests.
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

import diagnostics as diag                            # noqa: E402
import scf as scf_mod                                 # noqa: E402
import units                                          # noqa: E402
from axisym_model_writer import (to_meridional, vector_potential,  # noqa: E402
                                 verify_curl_on_cartesian,
                                 verify_meridional_curl, write_model)
from gradshafranov import solve_gradshafranov         # noqa: E402
from seed import r_guess                              # noqa: E402
from sweep_worker import _solve_toroidal_certified    # noqa: E402
from terms.toroidal_sc import ToroidalSC              # noqa: E402

RHO_C, MU_E, M_TOR, LMAX = 1.0e9, 2.0, 1.0, 16
N_MER, HALF_CM, CORNER = 385, 9.0e8, 1.7320508
B_POLE_TARGET = 1.0e9
K0_REF = 1.0e-13

# 1.35 is the field-free star, measured dynamically stable; 2.01 is measured
# unstable. These bracket the turnover.
TARGETS = (1.5, 1.7, 1.85)
K_AT_2MSUN = 3.245e-3
MASS_TOL = 3.0e-3
MAX_SOLVES = 6

DIV_GATE, CURL_GATE = 1.0e-12, 5.0e-2
OUTDIR = REPO / "models"


def solve_at(K):
    res, r, th, ov = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=K, m_tor_sc=M_TOR,
        rotation=None, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=LMAX,
        tol=1e-8, max_iter=200)
    if res is None:
        return None
    M = units.g_to_msun(scf_mod.total_mass(res["rho"], r, th))
    return dict(res=res, r=r, th=th, ov=ov, M=M, K=K)


def main():
    OUTDIR.mkdir(exist_ok=True)
    print(f"mass scan at rho_c = {RHO_C:.1e}, mu_e = {MU_E}\n")

    K_guess = K_AT_2MSUN
    for target in sorted(TARGETS, reverse=True):
        s0 = solve_at(K_guess)
        n = 1
        if s0 is None:
            print(f"M = {target}: first solve failed"); continue
        if abs(s0["M"] - target) / target > MASS_TOL:
            s1 = solve_at(K_guess * 0.7)
            n += 1
            while (s1 is not None and n < MAX_SOLVES
                   and abs(s1["M"] - target) / target > MASS_TOL):
                dM = s1["M"] - s0["M"]
                if abs(dM) < 1e-12:
                    break
                K2 = s1["K"] - (s1["M"] - target) * (s1["K"] - s0["K"]) / dM
                s0, s1 = s1, solve_at(float(np.clip(K2, 1e-6, 1.0)))
                n += 1
            s0 = s1 if s1 is not None else s0
        if s0 is None:
            print(f"M = {target}: no convergence"); continue

        res, r, th, ov = s0["res"], s0["r"], s0["th"], s0["ov"]
        rho, Phi, H = res["rho"], res["Phi"], res["H"]
        varpi = r[:, None] * np.sin(th)[None, :]
        W = abs(diag.gravitational_energy(rho, Phi, r, th))
        R_eq, R_pol = diag.equatorial_polar_radii(H, r, th)

        # calibrate k0 to the scenario's surface dipole, as in Phase 1
        u_ref = solve_gradshafranov(-4.0 * np.pi * varpi**2 * rho * K0_REF,
                                    r, th, lmax=LMAX)
        Br0, Bth0 = diag.poloidal_field(u_ref, r, th)
        bp_ref = diag.surface_dipolarity(np.hypot(Br0, Bth0), H, r, th)["B_pole"]
        k0 = K0_REF * (B_POLE_TARGET / bp_ref)

        u = solve_gradshafranov(-4.0 * np.pi * varpi**2 * rho * k0, r, th,
                                lmax=LMAX)
        Bphi = ToroidalSC(K=s0["K"], m=M_TOR).B_phi(rho, varpi)
        Br, Bth = diag.poloidal_field(u, r, th)
        E_pol, E_tor, _ = diag.magnetic_energies(Br, Bth, Bphi, r, th)

        rmax = 1.02 * CORNER * HALF_CM
        vp = np.linspace(0.0, rmax, N_MER)
        zz = np.linspace(-rmax, rmax, 2 * N_MER - 1)
        rho_m, u_m, bphi_m = to_meridional(r, th, (rho, u, Bphi), vp, zz)
        A_phi, A_z = vector_potential(vp, u_m, bphi_m)

        err_pol, err_tor = verify_meridional_curl(vp, zz, A_phi, A_z, u_m,
                                                  bphi_m)
        rel_div, b_max, _ = verify_curl_on_cartesian(vp, zz, A_phi, A_z,
                                                     half=HALF_CM, n_cart=64)
        retained = b_max / max(np.abs(Bphi).max(), 1e-300)
        ok = (rel_div < DIV_GATE and err_pol < CURL_GATE
              and err_tor < CURL_GATE and retained <= 1.02)

        tag = f"{s0['M']:.2f}".replace(".", "p")
        print(f"M = {s0['M']:.4f} Msun  ({n} solves)  K_tor = {s0['K']:.4e}")
        print(f"   E_tor/|W| = {E_tor / W:.4f}, max|B_phi| = "
              f"{np.abs(Bphi).max():.3e} G, R_pol/R_eq = {R_pol / R_eq:.4f}")
        print(f"   VE = {ov.get('frac_pol', float('nan')):.3f} frac_pol, "
              f"div B = {rel_div:.2e}, retained = {100 * retained:.1f}%  "
              f"{'ok' if ok else 'GATE FAILED'}")
        if not ok:
            print("   not written"); continue

        man = write_model(
            vp, zz, rho_m, A_phi, A_z, OUTDIR / f"scan_M{tag}.txt",
            dict(rho_c=RHO_C, mu_e=MU_E, K_tor=s0["K"], m_tor=M_TOR, k0=k0,
                 M_msun=s0["M"], R_eq_cm=R_eq, R_pol_cm=R_pol,
                 E_tor_over_W=E_tor / W, E_tor_over_Epol=E_tor / E_pol,
                 Bphi_max_G=float(np.abs(Bphi).max())),
            dict(curl_err_poloidal=err_pol, curl_err_toroidal=err_tor,
                 relative_divB_64cubed=rel_div,
                 amplitude_retained_64cubed=retained))
        print(f"   wrote models/{man['file']}\n")
        K_guess = s0["K"]


if __name__ == "__main__":
    main()
