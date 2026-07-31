"""How much mass will a poloidal-dominated field hold in a stratified star?

Run:  scf/.venv/bin/python3 investigations/nonbarotropic_mass_scan.py [nproc]

WHY THIS AND NOT THE TOROIDAL BRANCH

The self-consistent toroidal branch reaches 2 Msun and does not survive: in 3D
MHD it runs away in central density after about three dynamical times, at a
rate insensitive to resolution and to E_tor/|W|. Measured at equal magnetic
energy, the reason is plain -- in the core the toroidal field's radial Lorentz
force is -0.80 of gravity while the poloidal field's is +0.42. The toroidal
field supports the star globally through the virial while COMPRESSING its
centre, and the axial column is what collapses.

So the field has to be poloidal-dominated. Barotropy forbids that: the
poloidal SCF branch saturates at E_pol/|W| = 0.018 and the barotropic mixed
ceiling is 4% of the magnetic energy in the toroidal component. Both were
measured here, not assumed.

WHAT MAKES IT NON-BAROTROPIC, AND WHAT CLOSES IT

Barotropy forces the Lorentz force divided by rho to be a gradient, which pins
the poloidal source to rho M'(u). Dropping it leaves one condition, that a
single-valued P exists, and along a surface of constant Phi gravity does no
work, so with f_L = g grad u the momentum equation collapses to

    dP = g du     along each equipotential,          rho = (g u_r - P_r)/Phi_r

That leaves a free function -- P on the polar axis, one value per equipotential
-- which is genuine physical freedom, since a real star's stratification comes
from its history. Closing it needs physics, not numerics.

The closure used here is composition. This project's equation of state already
carries mu_e (eos.B_of_mu_e), so rho = B(mu_e) x^3 with x fixed by P: at given
pressure the density scales with mu_e. A mu_e gradient is therefore exactly the
baroclinicity the equilibrium needs, and it is physical for a white dwarf,
which is stratified in composition. Given P and rho from above, mu_e follows by
INVERSION -- nothing is fitted. What the free function must then satisfy is
that the mu_e it implies be

  * in a physical range, MU_E_MIN to MU_E_MAX below, and
  * stably stratified, i.e. non-increasing outward (Ledoux): heavier material
    below lighter, or the star is convectively unstable and the stratification
    it relies on does not exist.

Those two conditions are the gate. Configurations failing them are reported
and not counted.

WHAT IS SCANNED

k0 sets the poloidal amplitude; alpha stretches the anchor pressure profile on
the axis, which is the free function's one parameter here. The answer sought is
the largest mass whose implied composition profile passes the gate.

Writes nonbarotropic_mass_scan.csv.
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
for p in (REPO / "scf", REPO / "dashboard"):
    sys.path.insert(0, str(p))
warnings.filterwarnings("ignore")

import diagnostics as diag                       # noqa: E402
import eos                                       # noqa: E402
import scf as scf_mod                            # noqa: E402
import units                                     # noqa: E402
from gradshafranov import solve_gradshafranov    # noqa: E402
from poisson import solve_poisson                # noqa: E402
from scipy.interpolate import RegularGridInterpolator as RGI   # noqa: E402
from seed import r_guess                         # noqa: E402
from sweep_worker import _solve_toroidal_certified   # noqa: E402

RHO_C = 1.0e9
MU_E_REF = 2.0
LMAX = 16
NR, NTH = 129, 129
N_MARCH = 401
N_OUTER = 8

# Physical range for mu_e: 2.0 is C, O, Ne (Z/A = 1/2); 2.15 is Fe; below 2 needs
# hydrogen or helium, which an ultramassive interior does not have.
MU_E_MIN, MU_E_MAX = 1.95, 2.20
# Landau critical field. The equation of state has no field dependence, so a
# configuration whose peak exceeds B_c is inconsistent with its own
# microphysics -- the same gate the toroidal work used.
B_C = 4.414e13
# Fraction of the total mu_e range that may run the wrong way (outward
# increase) in the shell-averaged profile before the stratification counts as
# convectively unstable. Some non-monotonicity is expected from the numerics;
# a profile that rises by a third of its own range is not stratified.
LEDOUX_TOL = 0.30
# Below this fractional spread in the shell-averaged mu_e, the star counts as
# chemically homogeneous and the Ledoux criterion has nothing to act on.
HOMOGENEOUS_SPREAD = 0.01

# B_pol scales linearly with k0, and the field-free star gives max|B_pol| =
# 2.45e12 G at k0 = 1e-13, so B_c is reached near k0 = 1.8e-12. The first
# smoke test at 3e-12 came out at 1.66 B_c; the grid stops below that.
K0_LIST = tuple(np.geomspace(2.0e-13, 1.8e-12, 7))
ALPHA_LIST = (1.00, 1.10, 1.25, 1.40)
# rho_c is the lever the first grid never touched. The field-free mass rises
# toward the Chandrasekhar value with it, and the neutronization threshold for
# mu_e = 2 sits at 1.94e10, so there is room.
RHO_C_LIST = (1.0e9, 3.0e9, 8.0e9, 1.5e10)

OUT = HERE / "nonbarotropic_mass_scan.csv"


def grads(F, r, th):
    return np.gradient(F, r, axis=0), np.gradient(F, th, axis=1)


def integrate_pressure(P_axis, Phi, g, u, r, th, r_seeds):
    """dP = g du along equipotentials, marching in theta from the pole."""
    Phi_r, Phi_th = grads(Phi, r, th)
    u_r, u_th = grads(u, r, th)
    I = lambda A: RGI((r, th), A, bounds_error=False, fill_value=None)
    f_pr, f_pt, f_ur, f_ut, f_g = I(Phi_r), I(Phi_th), I(u_r), I(u_th), I(g)

    th_m = np.linspace(th[0], th[-1], N_MARCH)
    dth = th_m[1] - th_m[0]
    R = np.empty((len(r_seeds), len(th_m)))
    P = np.empty_like(R)
    R[:, 0], P[:, 0] = r_seeds, np.interp(r_seeds, r, P_axis)

    def deriv(rr, tt):
        tt = np.full_like(rr, tt) if np.isscalar(tt) else tt
        pts = np.stack([np.clip(rr, r[0], r[-1]),
                        np.clip(tt, th[0], th[-1])], axis=-1)
        pr = f_pr(pts)
        bad = np.abs(pr) < 1e-30
        pr = np.where(bad, 1.0, pr)
        drdt = np.where(bad, 0.0, -f_pt(pts) / pr)
        dudt = f_ut(pts) + f_ur(pts) * drdt
        return drdt, np.where(bad, 0.0, f_g(pts) * dudt)

    for i in range(len(th_m) - 1):
        t0, r0, p0 = th_m[i], R[:, i], P[:, i]
        k1r, k1p = deriv(r0, t0)
        k2r, k2p = deriv(r0 + 0.5 * dth * k1r, t0 + 0.5 * dth)
        k3r, k3p = deriv(r0 + 0.5 * dth * k2r, t0 + 0.5 * dth)
        k4r, k4p = deriv(r0 + dth * k3r, t0 + dth)
        R[:, i + 1] = np.clip(r0 + (dth / 6) * (k1r + 2*k2r + 2*k3r + k4r),
                              r[0], r[-1])
        P[:, i + 1] = p0 + (dth / 6) * (k1p + 2*k2p + 2*k3p + k4p)

    out = np.zeros((len(r), len(th)))
    for j, tj in enumerate(th):
        i = int(np.argmin(np.abs(th_m - tj)))
        o = np.argsort(R[:, i])
        rr, pp = R[o, i], P[o, i]
        keep = np.concatenate([[True], np.diff(rr) > 0])
        out[:, j] = np.interp(r, rr[keep], pp[keep], left=pp[keep][0], right=0.0)
    return out


def mu_e_from_P_rho(P, rho):
    """Invert the equation of state for composition.

    P fixes x through the T = 0 electron pressure, which does not involve
    mu_e at all; rho = B(mu_e) x^3 then gives mu_e directly. So the
    composition is READ OFF the solution rather than assumed.
    """
    A = eos.A_CONST
    y = np.maximum(P, 0.0) / A
    # invert P(x)/A = x(2x^2-3) sqrt(x^2+1) + 3 asinh(x) by bisection
    lo, hi = np.zeros_like(y), np.full_like(y, 200.0)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        f = mid * (2 * mid**2 - 3) * np.sqrt(mid**2 + 1) + 3 * np.arcsinh(mid)
        hi = np.where(f > y, mid, hi)
        lo = np.where(f > y, lo, mid)
    x = 0.5 * (lo + hi)
    B1 = eos.B_of_mu_e(1.0)
    denom = B1 * np.maximum(x, 1e-30) ** 3
    return np.where(P > 0, rho / denom, np.nan), x


def solve_one(args):
    k0, alpha, RHO_C = args
    try:
        # anchor: the field-free barotropic star supplies P on the axis
        res, r, th, _ = _solve_toroidal_certified(
            rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=0.0, m_tor_sc=1.0,
            rotation=None, mu_e=MU_E_REF, Nr_base=NR, Ntheta=NTH, lmax=LMAX,
            tol=1e-8, max_iter=200)
        if res is None:
            return dict(k0=k0, alpha=alpha, rho_c=RHO_C, status="no background")
        rho0, Phi, H0 = res["rho"], res["Phi"], res["H"]
        varpi = r[:, None] * np.sin(th)[None, :]
        x0 = eos.x_of_enthalpy(np.maximum(H0, 0.0), MU_E_REF)
        P_axis = alpha * np.where(H0[:, 0] > 0, eos.pressure(x0[:, 0]), 0.0)

        # poloidal field only: it is the component that pushes outward in the core
        src = -4.0 * np.pi * varpi**2 * rho0 * k0
        u = solve_gradshafranov(src, r, th, lmax=LMAX)
        g = np.where(varpi > 0, -src / (4 * np.pi * np.maximum(varpi, 1e-30)**2), 0.0)

        R_pol = diag.equatorial_polar_radii(H0, r, th)[1]
        seeds = np.linspace(r[1], 1.6 * R_pol, 260)

        rho = rho0
        for _ in range(N_OUTER):
            P = integrate_pressure(P_axis, Phi, g, u, r, th, seeds)
            P_r, _ = grads(P, r, th)
            Phi_r, _ = grads(Phi, r, th)
            u_r, _ = grads(u, r, th)
            safe = np.abs(Phi_r) > 1e-4 * np.abs(Phi_r).max()
            rho_new = np.where(safe, (g * u_r - P_r) / np.where(safe, Phi_r, 1.0), 0.0)
            # the centre: Phi_r vanishes there, so take rho from the innermost
            # safe shell rather than dividing by it
            for j in range(len(th)):
                k = np.flatnonzero(safe[:, j])
                if k.size:
                    rho_new[:k[0], j] = rho_new[k[0], j]
            rho_new = np.where(P > 0, np.maximum(rho_new, 0.0), 0.0)
            rho = 0.5 * rho + 0.5 * rho_new
            Phi = 0.6 * Phi + 0.4 * solve_poisson(rho, r, th, lmax=LMAX)

        M = units.g_to_msun(scf_mod.total_mass(rho, r, th))
        mu, x = mu_e_from_P_rho(P, rho)
        inside = (rho > 1e-4 * RHO_C) & (P > 0) & np.isfinite(mu)
        if inside.sum() < 500:
            return dict(k0=k0, alpha=alpha, rho_c=RHO_C, status="empty", M=M)

        # Both composition gates read the SHELL-AVERAGED profile, not the
        # pointwise field.
        #
        # mu_e comes from rho, which is a derivative of P, so pointwise it
        # carries that noise: at one grid point the 1-99 percentile range was
        # [1.955, 2.105], 7%, while the mass-weighted shell profile ran 1.995
        # to 2.006, 0.5%. Judging composition on percentiles of the pointwise
        # field therefore rejected configurations for numerical scatter. Both
        # earlier versions of these gates did exactly that.
        nsh = 24
        Rg = np.broadcast_to(r[:, None], mu.shape)
        r_star = float(Rg[inside].max())
        edges = np.linspace(0.0, r_star, nsh + 1)
        prof, mass_sh = [], []
        for a, b in zip(edges[:-1], edges[1:]):
            m = inside & (Rg >= a) & (Rg < b)
            if m.sum() > 20:
                w = rho[m]
                prof.append(float(np.sum(mu[m] * w) / np.sum(w)))
                mass_sh.append(float(w.sum()))
        prof, mass_sh = np.array(prof), np.array(mass_sh)
        if prof.size < 5:
            return dict(k0=k0, alpha=alpha, rho_c=RHO_C, status="too few shells", M=M)

        mu_lo, mu_hi = float(prof.min()), float(prof.max())
        spread = (mu_hi - mu_lo) / max(mu_lo, 1e-12)

        # Ledoux only bites where there IS a gradient. A star whose mean
        # composition is uniform to a fraction of a percent is homogeneous,
        # which is the commonest and most stable case, not an unstable one --
        # penalising it for lacking stratification was the previous gate's
        # error. Below HOMOGENEOUS_SPREAD the criterion does not apply; above
        # it, the profile must not rise outward by more than LEDOUX_TOL of its
        # own range.
        if spread < HOMOGENEOUS_SPREAD:
            frac_unstable = 0.0
        else:
            d = np.diff(prof)
            frac_unstable = float(np.maximum(d, 0).sum()
                                  / max(mu_hi - mu_lo, 1e-12))

        Br, Bth = diag.poloidal_field(u, r, th)
        E_pol, _, _ = diag.magnetic_energies(Br, Bth, np.zeros_like(rho), r, th)
        W = abs(diag.gravitational_energy(rho, Phi, r, th))

        B_pol_max = float(np.hypot(Br, Bth).max())
        ok = (MU_E_MIN <= mu_lo <= MU_E_MAX and MU_E_MIN <= mu_hi <= MU_E_MAX
              and frac_unstable < LEDOUX_TOL and B_pol_max < B_C)
        return dict(k0=k0, alpha=alpha, status="ok" if ok else "gate",
                    rho_c=RHO_C, M=M, mu_lo=mu_lo, mu_hi=mu_hi,
                    frac_unstable=frac_unstable, mu_spread=spread,
                    E_pol_over_W=E_pol / max(W, 1),
                    Bpol_max=B_pol_max, B_over_Bc=B_pol_max / B_C)
    except Exception as exc:                       # noqa: BLE001
        # the message, not just the type: a run that reported only
        # "AttributeError" 42 times cost a full round trip to the cluster
        return dict(k0=k0, alpha=alpha, rho_c=RHO_C,
                    status=f"error: {type(exc).__name__}: {exc}")


def main():
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else max(1, os.cpu_count() - 1)
    grid = list(itertools.product(K0_LIST, ALPHA_LIST, RHO_C_LIST))
    print(f"non-barotropic mass scan: {len(grid)} points on {nproc} processes")
    print(f"mu_e gate [{MU_E_MIN}, {MU_E_MAX}], Ledoux tolerance "
          f"{LEDOUX_TOL}\n")

    with Pool(nproc) as pool:
        rows = pool.map(solve_one, grid)

    print("  rho_c      k0          alpha  status  M (Msun)   mu_e range      "
          "spread  unstab  E_pol/|W|  B/B_c")
    for d in rows:
        head = (f"  {d.get('rho_c', 0):.1e}  {d['k0']:.3e}  {d['alpha']:.2f}  "
                f"{d['status']:6s}")
        if "M" in d and "mu_lo" in d:
            print(head
                  + f" {d['M']:8.4f}   "
                    f"{d['mu_lo']:.4f}-{d['mu_hi']:.4f}  "
                    f"{d.get('mu_spread', float('nan')):6.4f}  "
                    f"{d.get('frac_unstable', float('nan')):5.2f}  "
                    f"{d.get('E_pol_over_W', float('nan')):9.5f}  "
                    f"{d.get('B_over_Bc', float('nan')):5.2f}")
        else:
            print(head + f"  {d.get('M', '')}")

    good = [d for d in rows if d.get("status") == "ok"]
    if good:
        best = max(good, key=lambda d: d["M"])
        print(f"\nlargest mass passing the gate: {best['M']:.4f} Msun "
              f"at k0 = {best['k0']:.3e}, alpha = {best['alpha']:.2f}")
        print(f"  mu_e in [{best['mu_lo']:.3f}, {best['mu_hi']:.3f}], "
              f"E_pol/|W| = {best['E_pol_over_W']:.5f}, "
              f"max|B_pol| = {best['Bpol_max']:.3e} G")
    else:
        print("\nno configuration passed the gate -- see the mu_e ranges above")

    with OUT.open("w") as f:
        f.write(f"# non-barotropic poloidal scan, rho_c={RHO_C:.3e}, "
                f"mu_e gate [{MU_E_MIN},{MU_E_MAX}]\n")
        f.write("k0,alpha,rho_c,status,M_msun,mu_lo,mu_hi,frac_unstable,"
                "E_pol_over_W,Bpol_max_G,B_over_Bc,mu_spread\n")
        for d in rows:
            f.write(f"{d['k0']:.6e},{d['alpha']:.3f},{d.get('rho_c','')},{d['status']},"
                    f"{d.get('M', '')},{d.get('mu_lo', '')},"
                    f"{d.get('mu_hi', '')},{d.get('frac_unstable', '')},"
                    f"{d.get('E_pol_over_W', '')},{d.get('Bpol_max', '')},"
                    f"{d.get('B_over_Bc', '')},{d.get('mu_spread', '')}\n")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
