"""Where does a braking 2 Msun WD cross the ignition and neutronization lines?

Run:  scf/.venv/bin/python3 investigations/braking_sequence.py

THE QUESTION

A magnetically supported ultramassive white dwarf loses angular momentum and
field to magnetic braking. As the support goes, the star contracts and rho_c
rises, until it crosses one of two lines:

    rho ~ 3e9  g/cm^3   pycnonuclear C+C ignition   -> thermonuclear, SN Ia
    rho ~ 9.6e9 g/cm^3  electron capture on 20Ne    -> collapse, AIC to NS
    rho ~ 1.94e10       electron capture on 16O     -> collapse

Which line it meets first decides the fate. This maps the track: at FIXED
mass, how much toroidal support each central density needs, and therefore
where in that window the star sits as the support is removed.

HOW

The self-consistent toroidal branch is the one that reaches the required
energy (the mixed barotropic branch tops out four orders of magnitude short
-- see barotropic_ceiling.py). For each rho_c we solve for the amplitude
K_tor that puts the star at M_TARGET, by secant, warm-started from the
previous rho_c so it takes about three solves per point.

Read the ordering off the result, not the fate: the three lines are ordered
in density, so a track with rho_c RISING meets them in a fixed order. What
the sequence actually decides is where the star STARTS, hence how far it has
to travel, and whether it is still a certified equilibrium when it arrives.

Writes braking_sequence.csv.
"""

import sys
import time
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (REPO / "scf", REPO / "dashboard"):
    sys.path.insert(0, str(p))

warnings.filterwarnings("ignore")

import diagnostics as diag                       # noqa: E402
import scf as scf_mod                            # noqa: E402
import units                                     # noqa: E402
from seed import r_guess                         # noqa: E402
from sweep_worker import _solve_toroidal_certified   # noqa: E402
from terms.toroidal_sc import ToroidalSC         # noqa: E402

M_TARGET = 2.0
MU_E = 2.0
M_TOR = 1.0
NR = 65                    # 62 s/solve; M within 0.4% of the Nr=129 value
TOL = 1.0e-6
LMAX = 16
MASS_TOL = 2.0e-3          # relative, on M_TARGET
MAX_SOLVES = 6

RHO_CS = np.geomspace(8.0e8, 2.0e10, 11)
K_START = 3.245e-3         # the M = 2 Msun crossing at rho_c = 1e9

B_C = 4.414e13
RHO_CC = 3.0e9             # pycnonuclear C+C
RHO_NE = 9.6e9             # electron capture on 20Ne
RHO_O = 1.94e10            # electron capture on 16O
VE_GATE, FRACPOL_GATE = 1.0e-3, 0.2

OUT = HERE / "braking_sequence.csv"


def solve_at(rho_c, K):
    res, r, th, ov = _solve_toroidal_certified(
        rho_c=rho_c, R_guess=r_guess(rho_c), K_tor=K, m_tor_sc=M_TOR,
        rotation=None, mu_e=MU_E, Nr_base=NR, Ntheta=NR, lmax=LMAX,
        tol=TOL, max_iter=200)
    if res is None:
        return None
    M = units.g_to_msun(scf_mod.total_mass(res["rho"], r, th))
    return dict(res=res, r=r, th=th, ov=ov, M=M, K=K)


def main():
    print(f"target M = {M_TARGET} Msun, mu_e = {MU_E}, Nr = {NR}")
    print(f"lines: C+C {RHO_CC:.1e}, 20Ne {RHO_NE:.1e}, 16O {RHO_O:.2e} "
          f"g/cm^3\n")
    print("  rho_c        K_tor      M        E_tor/|W|  max|B_t|    /B_c  "
          "  VE       frac_pol  gate  solves")

    rows = []
    K_guess = K_START
    for rho_c in RHO_CS:
        # secant on M(K), warm-started from the previous rho_c
        K0, K1 = K_guess, K_guess * 0.85
        s0 = solve_at(rho_c, K0)
        if s0 is None:
            print(f"  {rho_c:.3e}   no convergence at K = {K0:.4e}")
            continue
        n = 1
        best = s0
        if abs(s0["M"] - M_TARGET) / M_TARGET > MASS_TOL:
            s1 = solve_at(rho_c, K1)
            n += 1
            while (s1 is not None
                   and abs(s1["M"] - M_TARGET) / M_TARGET > MASS_TOL
                   and n < MAX_SOLVES):
                dM = s1["M"] - s0["M"]
                if abs(dM) < 1e-12:
                    break
                K2 = s1["K"] - (s1["M"] - M_TARGET) * (s1["K"] - s0["K"]) / dM
                K2 = float(np.clip(K2, 1e-6, 1.0))
                s0, s1 = s1, solve_at(rho_c, K2)
                n += 1
            best = s1 if s1 is not None else s0
        if best is None:
            continue

        res, r, th, ov = best["res"], best["r"], best["th"], best["ov"]
        rho, Phi, H = res["rho"], res["Phi"], res["H"]
        tor = ToroidalSC(K=best["K"], m=M_TOR)
        ve = diag.virial_error_terms(rho, Phi, H, r, th, MU_E, rotation=None,
                                     poloidal=None, toroidal=tor)
        W = abs(ve["W"])
        Bt = np.abs(ve["Bphi"]).max()
        fp = ov.get("frac_pol", np.nan)
        ok = (ve["VE"] < VE_GATE) and (fp <= FRACPOL_GATE)
        R_eq, R_pol = diag.equatorial_polar_radii(H, r, th)
        rows.append((rho_c, best["K"], best["M"], ve["E_mag"] / W, Bt,
                     Bt / B_C, ve["VE"], fp, ok, R_eq, R_pol, n))
        print(f"  {rho_c:.3e}  {best['K']:.4e}  {best['M']:.4f}  "
              f"{ve['E_mag'] / W:9.5f}  {Bt:.3e}  {Bt / B_C:5.2f}  "
              f"{ve['VE']:.2e}  {fp:7.3f}   {'y' if ok else 'N'}     {n}")
        K_guess = best["K"]

    with OUT.open("w") as f:
        f.write(f"# braking sequence at fixed M = {M_TARGET} Msun, mu_e={MU_E}"
                f", Nr={NR}, m_tor={M_TOR}\n")
        f.write(f"# lines: C+C={RHO_CC:.3e}, 20Ne={RHO_NE:.3e}, "
                f"16O={RHO_O:.3e} g/cm^3; B_c={B_C:.4e} G\n")
        f.write("rho_c,K_tor,M_msun,E_mag_over_W,Bphi_max_G,Bphi_over_Bc,"
                "VE,frac_pol,certified,R_eq_cm,R_pol_cm,solves\n")
        for row in rows:
            f.write(",".join(f"{v:.6e}" if isinstance(v, float) else str(v)
                             for v in row) + "\n")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"elapsed {time.time() - t0:.0f} s")
