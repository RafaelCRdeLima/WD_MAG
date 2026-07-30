"""Build the 1D hydrostatic white dwarf that FLASH will be given.

Run:  scf/.venv/bin/python3 flash_crosscheck/make_wd_model.py

EOS choice. FLASH ships Gamma, Helmholtz, Multigamma and Tabulated. For a
white dwarf the appropriate one is Helmholtz (Timmes & Swesty): degenerate
electrons at finite temperature, plus ions and radiation. FLASH has no
ztwd, so the star cannot be handed over as a ztwd solution and evolved
under Helmholtz without an EOS mismatch -- the failure mode that caused
the Step 3 collapse (docs/teoria.md Sec 6.4).

The resolution is that the mismatch is negligible if the star is cold
enough, and this script MEASURES that rather than assuming it. The
structure is integrated in hydrostatic equilibrium with the project's own
ztwd EOS (scf/eos.py -- reused, not reimplemented, and already validated
against the Chandrasekhar limit), at a uniform temperature low enough that
the ion and radiation contributions Helmholtz adds are far below the drift
we are trying to measure. The printed diagnostic P_ion/P_degenerate is the
number that justifies handing this model to Helmholtz.

Validation target: the wd_braithwaite production star reports, in its own
run log, R_star = 238681011.6 cm and M = 1.346499895 Msun at
rho_c = 988393849.5 g/cm^3. This integration has to reproduce those, or
the model is not the same star the Castro results were measured on.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "scf"))

from eos import B_of_mu_e, enthalpy, pressure  # noqa: E402

G = 6.67428e-8
M_SUN = 1.98892e33
K_BOLTZ = 1.380658e-16
M_U = 1.660539e-24
A_RAD = 7.5657e-15        # erg cm^-3 K^-4

RHO_C = 988393849.5       # g/cm^3, the production star
MU_E = 2.0
ABAR, ZBAR = 4.0, 2.0     # mu2.net: He4-like, mu_e = A/Z = 2
T_ISO = 1.0e7             # K, uniform; justified by the printed P_ion/P_deg

DOMAIN_HALF = 4.90e8      # cm, same domain as the Castro run
RHO_FLOOR = 1.0e4          # g/cm^3; matches Castro small_dens = 1e4 (inputs.evolve:83)
NPTS_OUT = 2048

OUT = HERE / "wd_model.dat"


def rho_of_P(P_target, rho_hi=1e12):
    """Invert the ztwd P(rho) by bisection (P is monotonic in rho)."""
    lo, hi = 0.0, rho_hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if P_ztwd(mid) < P_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def P_ztwd(rho):
    if rho <= 0.0:
        return 0.0
    x = (rho / B_of_mu_e(MU_E)) ** (1.0 / 3.0)
    return float(pressure(x))


def main():
    P_c = P_ztwd(RHO_C)
    print(f"rho_c = {RHO_C:.6e} g/cm^3  ->  P_c = {P_c:.6e} erg/cm^3")

    # dP/dr = -G M rho / r^2 ;  dM/dr = 4 pi r^2 rho ;  rho = rho(P)
    def rhs(r, y):
        P, M = y
        if P <= 0.0:
            return [0.0, 0.0]
        rho = rho_of_P(P)
        dP = 0.0 if r <= 0.0 else -G * M * rho / r**2
        return [dP, 4.0 * np.pi * r**2 * rho]

    def surface(r, y):
        return y[0] - 1e-8 * P_c
    surface.terminal = True
    surface.direction = -1

    r0 = 1.0e3
    y0 = [P_c, 4.0 / 3.0 * np.pi * r0**3 * RHO_C]
    sol = solve_ivp(rhs, (r0, 1.0e10), y0, events=surface, rtol=1e-10,
                    atol=[1e-6 * P_c, 1e20], dense_output=True, max_step=1e6)

    R = float(sol.t_events[0][0])
    M = float(sol.y_events[0][0][1])
    t_dyn = np.sqrt(R**3 / (G * M))
    print(f"\nintegrated: R = {R:.6e} cm, M = {M / M_SUN:.6f} Msun, "
          f"t_dyn = {t_dyn:.7f} s")

    # The production star's own numbers, from its run log.
    R_ref, M_ref, t_dyn_ref = 238681011.6, 1.346499895, 0.2758062098
    print(f"production star:  R = {R_ref:.6e} cm, M = {M_ref:.6f} Msun, "
          f"t_dyn = {t_dyn_ref:.7f} s")
    print(f"  agreement: R {100 * abs(R / R_ref - 1):.3f}%, "
          f"M {100 * abs(M / M_SUN / M_ref - 1):.3f}%, "
          f"t_dyn {100 * abs(t_dyn / t_dyn_ref - 1):.3f}%")

    # Sample onto a uniform radial grid covering the domain diagonal.
    r = np.linspace(0.0, DOMAIN_HALF * np.sqrt(3.0), NPTS_OUT)
    rho = np.full_like(r, RHO_FLOOR)
    inside = (r > 0) & (r < R)
    P_interp = np.maximum(sol.sol(r[inside])[0], 0.0)
    rho[inside] = [max(rho_of_P(p), RHO_FLOOR) for p in P_interp]
    rho[0] = RHO_C

    # Is T_ISO cold enough that Helmholtz's extra terms do not disturb HSE?
    P_deg = np.array([P_ztwd(d) for d in rho])
    P_ion = rho * K_BOLTZ * T_ISO / (ABAR * M_U)
    P_rad = A_RAD * T_ISO**4 / 3.0
    core = rho > 0.1 * RHO_C
    print(f"\nat T = {T_ISO:.1e} K, inside the core:")
    print(f"  P_ion/P_deg  max = {(P_ion[core] / P_deg[core]).max():.3e}")
    print(f"  P_rad/P_deg  max = {(P_rad / P_deg[core]).max():.3e}")
    print("  (both must sit far below the ~2e-2 drift being measured)")

    T = np.full_like(r, T_ISO)
    with OUT.open("w") as f:
        f.write("# 1D hydrostatic white dwarf, ztwd structure, isothermal T\n")
        f.write(f"# rho_c={RHO_C:.10e} R={R:.10e} M_msun={M / M_SUN:.10f} "
                f"t_dyn={t_dyn:.10f} T_iso={T_ISO:.10e}\n")
        f.write("# r density temperature\n")
        for ri, di, ti in zip(r, rho, T):
            f.write(f"{ri:.10e} {di:.10e} {ti:.10e}\n")
    print(f"\nwrote {OUT.name}: {NPTS_OUT} points to r = {r[-1]:.4e} cm")
    print(f"  rho_c/rho_mean = {RHO_C / (M / (4 / 3 * np.pi * R**3)):.1f}"
          "  (the real central condensation, unlike a single polytrope)")
    print(f"\n--- for flash.par ---\n  sim_tempIso = {T_ISO:.6e}"
          f"\n  stop_time   = {4 * t_dyn:.6f}   (4 t_dyn)")


if __name__ == "__main__":
    main()
