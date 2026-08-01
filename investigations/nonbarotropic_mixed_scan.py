"""Mixed poloidal-toroidal, non-barotropic: can each component fix what the
other breaks?

Run:  scf/.venv/bin/python3 investigations/nonbarotropic_ceiling_scan.py [nproc]

WHY MIXED, FROM THIS PROJECT'S OWN MEASUREMENTS

Measured at equal magnetic energy, in the core (r < 0.2 R) the radial Lorentz
force of a toroidal field is -0.80 of gravity and that of a poloidal field is
+0.42. The toroidal branch reaches E_tor/|W| = 0.203 and 2 Msun but COMPRESSES
the centre, and in 3D the axial column runs away in about three dynamical
times. The poloidal branch holds the centre but saturates -- at 0.018 in
barotropy, at 0.036 non-barotropically -- and reaches only about 1.7 Msun, and
a purely poloidal field is Markey-Tayler unstable anyway.

Each component fixes what the other breaks. Barotropy is what kept them apart:
in a barotropic mixed equilibrium the toroidal component is confined to the
closed-line region and caps at about 4% of the magnetic energy, measured here.
Dropping barotropy lifts exactly that restriction, and the previous scans used
the freedom for the poloidal component alone.

WHAT IS ADDED

beta(u) enters the same way it does in Grad-Shafranov: with u obtained by
inverting Delta* on a chosen source, Delta* u = src exactly, and

    g = -[src + beta beta'(u)] / (4 pi varpi^2),      B_phi = beta(u)/varpi

beta must vanish where field lines reach the surface, so the toroidal
component is confined to u > u_s. That confinement follows from the azimuthal
force balance and the vacuum exterior, not from barotropy, so it survives.

THE GATE THAT MATTERS

The core Lorentz force is computable at construction time, so a configuration
that would compress its centre can be rejected without spending a 3D run to
find out. It is reported as f_core, the mass-weighted radial Lorentz force
inside 0.2 R as a fraction of gravity: negative compresses, positive supports.

Writes nonbarotropic_mixed_scan.csv.
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

LMAX, NR, NTH = 16, 129, 129
# 32, not the 8 the first scan used. Measured on one point, the virial error
# runs 4.8e-3 at 8 iterations, 1e-5 at 16 and settles near 4e-4 by 32 -- so
# every point of the first scan was under-converged, by five times the gate at
# the interesting configurations. The physics of the trade did not depend on
# it, but no number from that scan should be quoted.
N_OUTER = 32
RETURN_FIELDS = False
VE_GATE = 1.0e-3                # the same gate used throughout this project
B_C = 4.414e13
BAROTROPIC_CEILING = 0.018      # measured: the poloidal SCF branch saturates here
NOISE_FLOOR = 0.0203            # measured at k0 = 0 by nonbarotropic_noise_floor.py
STRATIFIED_MULTIPLE = 3.0       # below this the composition structure is unresolved
MU_E_MIN, MU_E_MAX = 1.90, 2.20
# The strictly defensible interior composition of an ultramassive white dwarf:
# 2.00 is C, O, Ne and He; 2.15 is Fe. Below 2.00 requires hydrogen, which has
# burned. Reported alongside the looser gate so near-misses stay visible
# instead of being collapsed into a pass/fail.
MU_E_STRICT_MIN, MU_E_STRICT_MAX = 2.00, 2.15

# Starts where the previous scan stopped. B/B_c was 0.69 at k0 = 1.25e-12 on
# the field-free background, so there is room to roughly 1.8e-12 before the
# Landau limit -- but the star inflates as k0 rises, so the gate is applied to
# each configuration rather than assumed from the grid.
K0_LIST = tuple(np.geomspace(1.0e-12, 2.5e-12, 4))
# alpha multiplies the anchor pressure, and mu_e = rho / (B1 x^3) with x set by
# P, so RAISING alpha LOWERS mu_e -- roughly as alpha^(-3/4) for a relativistic
# gas. The previous grid ran alpha >= 1.00 and returned mean mu_e between 1.62
# and 1.95 at every one of its 42 points: the physically defensible region,
# mu_e near 2.0 to 2.15, was never sampled. It lies below 1, near
# alpha = 1.93/2.05 to the -4/3, about 0.92.
ALPHA_LIST = (0.70, 0.80, 0.90, 1.00)
RHO_C_LIST = (1.0e9, 3.0e9)
# Toroidal strength, as the energy ratio E_tor/E_pol asked of the pair.
# 0 reproduces the purely poloidal scan; the barotropic mixed ceiling
# corresponds to about 0.04, so everything at 1 and above is territory
# barotropy cannot reach at all.
# Parametrised by E_tor/|W| ABSOLUTE, not by E_tor/E_pol.
#
# The ratio form was wrong and the smoke test showed it: asking E_tor/E_pol =
# 16 of a poloidal field already at E_pol/|W| = 0.034 demands a magnetic energy
# larger than the gravitational one, E_mag/|W| = 1.17, and the virial error
# went 0.0004, 0.19, 1.96. Those are not stars. A ratio has no scale in it, so
# it amplifies without bound; the absolute fraction of |W| is the quantity with
# physical meaning, and 0.203 is where the self-consistent toroidal branch
# delivers 2 Msun.
# Restricted to the range where an equilibrium is actually attainable. A
# sweep at zeta = 2 gave virial errors 0.0004, 0.0004, 0.0079, 0.035, 0.042,
# 0.026, 0.010 at E_tor/|W| = 0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10, so the
# 1e-3 gate is met only up to about 0.01 and everything beyond is rejected on
# arrival. Sampling it would spend the scan mapping non-equilibria.
TOR_W_LIST = (0.0, 0.005, 0.01, 0.02)
# beta ~ w^zeta with zeta = 2 so that beta AND beta' both vanish at the
# separatrix. At zeta = 1 the derivative jumps there, which is a current
# sheet, and the residual of the force balance piled up in the equatorial
# plane at 0.4-0.6 R_eq -- exactly on the separatrix -- carrying the whole
# virial error. Raising it to 2 cut the error at E_tor/|W| = 0.10 from 0.029
# to 0.010. This is the same choice Ciolfi-Rezzolla and Lander-Jones make.
ZETA = 2.0

OUT = HERE / "nonbarotropic_mixed_scan.csv"


def solve_one(args):
    k0, alpha, rho_c, tor_w = args
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

        # Toroidal component, confined to the closed-line region. u_s is the
        # largest u on the stellar surface: lines with u <= u_s reach the
        # vacuum, where B_phi must vanish, so beta must vanish on them too.
        #
        # This has to be recomputed against the CURRENT star, not the
        # field-free background, and that is not a refinement. Freezing u_s at
        # rho0 left the toroidal field confined to a region defined by a
        # 1.74 Msun star while the solve converged on a 1.21 Msun one, so part
        # of the field sat outside the matter with nothing to confine it. The
        # virial error read 0.17 and stayed there from 32 outer iterations to
        # 256 -- converged, and converged on something that was not in
        # equilibrium. beta is renormalised the same way, against the current
        # |W|, so E_tor/|W| means what it says at the solution rather than at
        # the starting guess.
        def toroidal(rho_now, Phi_now):
            u_s = max((u[np.flatnonzero(rho_now[:, j] > 0)[-1], j]
                       for j in range(len(th))
                       if np.any(rho_now[:, j] > 0)), default=0.0)
            u_norm = max(u.max() - u_s, 1e-300)
            closed = (rho_now > 0.0) & (u > u_s)
            w = np.where(closed, (u - u_s) / u_norm, 0.0)
            shape = np.where(varpi > 0, np.power(w, ZETA)
                             / np.maximum(varpi, 1e-30), 0.0)
            zero = np.zeros_like(rho_now)
            _, E_unit, _ = diag.magnetic_energies(zero, zero, shape, r, th)
            W_now = abs(diag.gravitational_energy(rho_now, Phi_now, r, th))
            if tor_w > 0.0 and E_unit > 0:
                beta_0 = np.sqrt(tor_w * W_now / E_unit)
            else:
                beta_0 = 0.0
            beta = beta_0 * np.power(w, ZETA)
            dbeta = np.where(closed, beta_0 * ZETA
                             * np.power(w, ZETA - 1.0) / u_norm, 0.0)
            Bphi = np.where(varpi > 0, beta / np.maximum(varpi, 1e-30), 0.0)
            S = src + beta * dbeta
            g = np.where(varpi > 0,
                         -S / (4 * np.pi * np.maximum(varpi, 1e-30)**2), 0.0)
            return Bphi, g

        Bphi, g = toroidal(rho0, Phi)
        R_pol = diag.equatorial_polar_radii(H0, r, th)[1]
        seeds = np.linspace(r[1], 1.6 * R_pol, 260)

        rho = rho0
        for _ in range(N_OUTER):
            Bphi, g = toroidal(rho, Phi)
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
        E_pol, E_tor, E_mag = diag.magnetic_energies(Br, Bth, Bphi, r, th)
        W_signed = diag.gravitational_energy(rho, Phi, r, th)
        W = abs(W_signed)
        B_max = float(np.hypot(Br, Bth).max())   # poloidal peak, for reference
        EpW = E_pol / max(W, 1.0)

        # Virial error, computed from the solver's OWN pressure rather than
        # from an equation of state -- the whole point here is that P is not a
        # function of rho alone. A configuration that fails this is not an
        # equilibrium and nothing else about it means anything.
        VE = abs(W_signed + 3.0 * diag.volume_integral(P, r, th) + E_mag) / W

        # Core Lorentz force, mass weighted inside 0.2 R, as a fraction of
        # gravity. f_L = g grad u, so this carries BOTH components through g.
        # Cells with varpi below 5% of the radius are excluded: the 1/varpi^2
        # in g is finite there only in the limit and is numerically useless.
        u_r_full, _ = grads(u, r, th)
        r_out = float(Rg[inside].max())
        core = inside & (Rg < 0.2 * r_out) & (varpi > 0.05 * r_out)
        if core.sum() > 20:
            wc = rho[core]
            grav = rho * np.abs(grads(Phi, r, th)[0])
            f_core = float(np.sum((g * u_r_full)[core] * wc) / np.sum(wc)
                           / max(np.sum(grav[core] * wc) / np.sum(wc), 1e-30))
        else:
            f_core = float("nan")
        B_tot_max = float(np.sqrt(Br**2 + Bth**2 + Bphi**2).max())
        bt_bp = (float(np.abs(Bphi).max()
                       / max(np.hypot(Br, Bth).max(), 1e-30)))

        if RETURN_FIELDS:
            return dict(rho=rho, P=P, Phi=Phi, u=u, g=g, Bphi=Bphi,
                        Br=Br, Bth=Bth, r=r, th=th, varpi=varpi,
                        inside=inside, VE=VE, M=M, E_tor=E_tor, E_mag=E_mag,
                        W=W_signed)

        # Composition gate, POINTWISE and mass weighted.
        #
        # The shell-averaged version this replaces was too weak to be worth
        # anything. Averaging mu_e over 25 mass-weighted shells reported
        # 2.048-2.630 for a configuration whose pointwise maximum was 65.2,
        # and reported 2.002-2.122 -- inside the strict window -- for the
        # 1.386 Msun point whose pointwise maximum is 2.627. A star is not
        # made of shell averages; if some of its mass wants mu_e = 2.6 there
        # is no cold white dwarf composition that supplies it.
        #
        # Mass weighting is what handles the surface. The inversion is
        # ill-conditioned where P -> 0, and with no toroidal field the raw
        # pointwise maximum sits at 0.98 R for exactly that reason -- but
        # that shell holds almost no mass, so it cannot move a mass-weighted
        # percentile. No ad hoc radial cut is needed.
        wm = (rho * (r[:, None] ** 2) * np.sin(th)[None, :])[inside]
        mv = mu[inside]
        order = np.argsort(mv)
        mv, wm = mv[order], wm[order]
        cw = np.cumsum(wm) / max(wm.sum(), 1e-300)
        mu_p01 = float(np.interp(0.01, cw, mv))
        mu_p99 = float(np.interp(0.99, cw, mv))
        bad = (mv < MU_E_MIN) | (mv > MU_E_MAX)
        f_bad = float(wm[bad].sum() / max(wm.sum(), 1e-300))

        # at most 1% of the mass allowed outside the window, and the 1-99
        # percentile band itself has to lie inside it
        physical = (f_bad <= 0.01
                    and MU_E_MIN <= mu_p01 and mu_p99 <= MU_E_MAX)
        # pointwise as well -- the shell-averaged form left here after the
        # gate was rewritten is what reported 1.3859 Msun as strictly
        # physical while 4.9% of its mass sat outside the window
        strict = (f_bad <= 0.01 and MU_E_STRICT_MIN <= mu_p01
                  and mu_p99 <= MU_E_STRICT_MAX)
        # B_c applies to the TOTAL field, not the poloidal part alone.
        ok = (physical and B_tot_max < B_C and ledoux < 0.30
              and VE < VE_GATE and f_core > 0.0)
        return dict(k0=k0, alpha=alpha, rho_c=rho_c,
                    status="ok" if ok else "gate", M=M, M0=M0,
                    mu_mean=float(prof.mean()), mu_lo=float(prof.min()),
                    mu_hi=float(prof.max()), over_floor=over_floor,
                    ledoux=ledoux, EpW=EpW, VE=VE,
                    strict=strict, beats_ceiling=EpW > BAROTROPIC_CEILING,
                    mu_p01=mu_p01, mu_p99=mu_p99, f_bad=f_bad,
                    tor_w=tor_w, EtW=E_tor / max(W, 1.0),
                    EmW=E_mag / max(W, 1.0),
                    bt_bp=bt_bp, f_core=f_core, B_over_Bc=B_tot_max / B_C)
    except Exception as exc:                       # noqa: BLE001
        return dict(k0=k0, alpha=alpha, rho_c=rho_c,
                    status=f"error: {type(exc).__name__}: {exc}")


def main():
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else max(1, os.cpu_count() - 1)
    grid = list(itertools.product(K0_LIST, ALPHA_LIST, RHO_C_LIST,
                                 TOR_W_LIST))
    print(f"ceiling scan: {len(grid)} points on {nproc} processes")
    print(f"noise floor {NOISE_FLOOR:.4f}; stratified above "
          f"{STRATIFIED_MULTIPLE}x it; barotropic ceiling "
          f"E_pol/|W| = {BAROTROPIC_CEILING}\n")
    with Pool(nproc) as pool:
        rows = pool.map(solve_one, grid)

    print("  rho_c    k0        alp  E_t/W   status  M      mu_e range     "
          "E_mag/|W|  Bt/Bp   f_core  B/Bc     VE   strict")
    for d in sorted(rows, key=lambda x: (x["rho_c"], x["k0"], x["alpha"],
                                         x.get("tor_w", 0))):
        head = (f"  {d['rho_c']:.0e}  {d['k0']:.2e} {d['alpha']:.2f} "
                f"{d.get('tor_w', 0):5.2f}  {d['status']:6s}")
        if "mu_mean" in d:
            print(head + f" {d['M']:6.3f} {d['mu_p01']:.3f}-{d['mu_p99']:.3f} "
                         f"{d['f_bad']:6.4f} "
                         f"{d['EmW']:9.5f} {d['bt_bp']:6.2f}  "
                         f"{d['f_core']:+6.2f}  {d['B_over_Bc']:5.2f}  "
                         f"{d['VE']:.5f}  {'YES' if d['strict'] else '  .'}")
        else:
            print(head + f" {d.get('M', '')}")

    # what the whole exercise is for: mass, physical composition, a core that
    # is held rather than squeezed, and a field the equation of state can
    # still describe
    good = [d for d in rows if d.get("status") == "ok"]
    if good:
        b = max(good, key=lambda d: d["M"])
        print(f"\nheaviest configuration passing every gate:")
        print(f"  M = {b['M']:.4f} Msun at rho_c = {b['rho_c']:.0e}, "
              f"k0 = {b['k0']:.2e}, alpha = {b['alpha']:.2f}, "
              f"E_tor/E_pol = {b['tor_w']:.2f}")
        print(f"  Bt/Bp = {b['bt_bp']:.2f}, E_mag/|W| = {b['EmW']:.5f}, "
              f"mu_e {b['mu_lo']:.3f}-{b['mu_hi']:.3f}")
        print(f"  core force {b['f_core']:+.2f} of gravity (positive = held), "
              f"B/B_c = {b['B_over_Bc']:.2f}, VE = {b['VE']:.5f}")
    else:
        held = [d for d in rows if d.get("f_core", -1) > 0]
        print(f"\nnothing passes every gate; {len(held)} of {len(rows)} at "
              f"least hold their core (f_core > 0)")

    beat = [d for d in rows if d.get("beats_ceiling") and d.get("B_over_Bc", 9) < 1]
    if beat:
        b = max(beat, key=lambda d: d["EpW"])
        print(f"\nbeats the barotropic ceiling below B_c: E_pol/|W| = "
              f"{b['EpW']:.5f} at rho_c = {b['rho_c']:.0e}, k0 = {b['k0']:.2e}")
        print(f"  M = {b['M']:.4f} Msun ({100*(b['M']-b['M0'])/b['M0']:+.1f}% "
              f"over field-free), mu_e {b['mu_lo']:.3f}-{b['mu_hi']:.3f} "
              f"({b['over_floor']:.1f}x floor), VE = {b['VE']:.5f}")
    else:
        print("\nnothing beats the barotropic ceiling below B_c")

    phys = [d for d in rows
            if d.get("strict") and d.get("B_over_Bc", 9) < 1
            and d.get("VE", 9) < VE_GATE]
    if phys:
        q = max(phys, key=lambda d: d["M"])
        print(f"\nheaviest with STRICTLY physical composition (mu_e in "
              f"[{MU_E_STRICT_MIN}, {MU_E_STRICT_MAX}]), below B_c, certified:")
        print(f"  M = {q['M']:.4f} Msun at rho_c = {q['rho_c']:.0e}, "
              f"k0 = {q['k0']:.2e}, alpha = {q['alpha']:.2f}")
        print(f"  mu_e {q['mu_lo']:.3f}-{q['mu_hi']:.3f}, "
              f"E_pol/|W| = {q['EpW']:.5f} "
              f"({'beats' if q['beats_ceiling'] else 'below'} the barotropic "
              f"ceiling), B/B_c = {q['B_over_Bc']:.2f}, VE = {q['VE']:.5f}")
    else:
        print("\nno configuration has strictly physical composition below "
              "B_c -- the composition demand, not the Landau limit, closes "
              "the route")

    with OUT.open("w") as f:
        f.write(f"# non-barotropic ceiling scan; floor={NOISE_FLOOR}, "
                f"barotropic ceiling={BAROTROPIC_CEILING}, B_c={B_C:.4e}\n")
        f.write("k0,alpha,rho_c,status,M_msun,M0_msun,mu_mean,mu_lo,mu_hi,"
                "over_floor,ledoux,mu_p01,mu_p99,f_bad,E_pol_over_W,E_mag_over_W,VE,tor_w,E_tor_over_W,"
                "bt_bp,f_core,strict,beats_ceiling,B_over_Bc\n")
        for d in rows:
            f.write(",".join(str(d.get(k, "")) for k in
                             ("k0", "alpha", "rho_c", "status", "M", "M0",
                              "mu_mean", "mu_lo", "mu_hi", "over_floor",
                              "ledoux", "mu_p01", "mu_p99", "f_bad",
                              "EpW", "EmW", "VE", "tor_w", "EtW",
                              "bt_bp", "f_core", "strict", "beats_ceiling",
                              "B_over_Bc")) + "\n")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
