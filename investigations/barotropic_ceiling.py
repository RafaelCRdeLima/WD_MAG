"""How much toroidal field fits in a BAROTROPIC twisted torus? Measured here.

Run:  scf/.venv/bin/python3 investigations/barotropic_ceiling.py

This is the baseline for the non-barotropic route. The claim the mixed-field
paper rests on -- that a self-consistent barotropic mixed equilibrium admits
only a few percent of toroidal energy (docs/teoria.md Sec 1.9, note G4) -- is
taken from the literature. Before building a solver to beat that ceiling, it
has to be measured in this code, on this star.

The equilibrium. For an axisymmetric magnetostatic star the phi-component of
the Lorentz force must vanish, which forces

    varpi B_phi = beta(u)                                            (1)

with u the flux function. This holds whether or not the star is barotropic:
it says omega*B_phi is constant along poloidal field lines. Since B_phi must
vanish outside the star, beta must vanish on every line that reaches the
surface, so the toroidal field is confined to the closed-line region
u > u_s. That confinement is NOT what barotropy costs -- it survives the
generalisation.

What barotropy costs is the source. With (1/rho) grad P = grad h the whole
left side is a gradient, so the Lorentz force divided by rho must be one
too, and since it is parallel to grad u it must equal M'(u) grad u. That
pins the poloidal source to the form rho M'(u) and gives

    Delta* u = -4 pi varpi^2 rho f(u) - beta(u) beta'(u)             (2)

which is what this script solves, with f(u) = k0 constant and

    beta(u) = lambda * (u - u_s)^zeta   for u > u_s,   0 otherwise.

Sweeping lambda upward asks how much toroidal energy the closed region can
hold before the equilibrium stops existing. The answer is the ceiling.

The star is field-free at rho_c = 1e9, mu_e = 2 -- the standard setting for
this family, and deliberately NOT the 2 Msun toroidal star, whose structure
was set by a different toroidal law (B_phi = K rho varpi, which is not of
the form beta(u)/varpi and is only an equilibrium because a purely toroidal
field has no poloidal lines to constrain it).

Writes barotropic_ceiling.csv.
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

import diagnostics as diag                       # noqa: E402
import scf as scf_mod                            # noqa: E402
import units                                     # noqa: E402
from gradshafranov import solve_gradshafranov    # noqa: E402
from seed import r_guess                         # noqa: E402
from sweep_worker import _solve_toroidal_certified   # noqa: E402

RHO_C = 1.0e9
MU_E = 2.0
K0 = 1.0e-13          # the poloidal source level used by the mixed paper
ZETAS = (1.0, 1.1, 2.0)
# beta is parametrised by the toroidal field it asks for, not by a bare
# coefficient: u ~ 1e27 here, so a lambda of order unity would demand
# B_phi ~ 1e19 G and the iteration simply explodes. beta_0 = B_REF * varpi_ref
# makes the sweep a sweep in gauss.
B_REFS = (0.0, 1e10, 3e10, 1e11, 3e11, 1e12, 3e12, 1e13, 3e13, 1e14)
LMAX = 16
OUT = HERE / "barotropic_ceiling.csv"


def surface_flux(u, H, r, th):
    """u_s: the largest u on the stellar surface -- the last closed line.

    Lines with u <= u_s reach the surface, so beta must vanish there and the
    toroidal field lives only where u > u_s.
    """
    vals = []
    for j in range(len(th)):
        inside = H[:, j] > 0.0
        if not inside.any():
            continue
        vals.append(u[np.flatnonzero(inside)[-1], j])
    return max(vals) if vals else 0.0


def beta_of(u, u_s, H, beta_0, u_norm, zeta):
    """beta(u) = beta_0 * ((u - u_s)/u_norm)^zeta, zero outside the closed
    region. u_norm is FIXED from the beta_0 = 0 solution, so the functional
    form does not drift as u evolves during the iteration."""
    w = np.where((u > u_s) & (H > 0.0), (u - u_s) / u_norm, 0.0)
    return beta_0 * np.power(w, zeta), w


def solve_mixed(rho, H, r, th, k0, beta_0, u_norm, zeta, tol=1e-8,
                max_iter=600, relax=0.15):
    """Iterate Eq. (2) to convergence. Returns (u, u_s, iters, converged)."""
    omega2 = (r[:, None] * np.sin(th)[None, :]) ** 2
    src_pol = -4.0 * np.pi * omega2 * rho * k0
    u = solve_gradshafranov(src_pol, r, th, lmax=LMAX)
    if beta_0 == 0.0:
        return u, surface_flux(u, H, r, th), 0, True
    u0_scale = np.abs(u).max()

    for it in range(1, max_iter + 1):
        u_s = surface_flux(u, H, r, th)
        _, w = beta_of(u, u_s, H, beta_0, u_norm, zeta)
        # beta beta' = beta_0^2 zeta w^(2 zeta - 1) / u_norm
        bb = (beta_0 ** 2) * zeta * np.power(w, 2.0 * zeta - 1.0) / u_norm
        u_new = solve_gradshafranov(src_pol - bb, r, th, lmax=LMAX)
        if not np.isfinite(u_new).all():
            return u, u_s, it, False
        delta = np.abs(u_new - u).max() / max(np.abs(u).max(), 1e-300)
        u = (1.0 - relax) * u + relax * u_new
        if np.abs(u).max() > 1.0e3 * u0_scale:      # running away
            return u, surface_flux(u, H, r, th), it, False
        if delta < tol:
            return u, surface_flux(u, H, r, th), it, True
    return u, surface_flux(u, H, r, th), max_iter, False


def main():
    res, r, th, _ = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=0.0, m_tor_sc=1.0,
        rotation=None, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=LMAX,
        tol=1e-8, max_iter=200)
    if res is None:
        raise SystemExit("the field-free solve did not converge")
    rho, Phi, H = res["rho"], res["Phi"], res["H"]
    M = units.g_to_msun(scf_mod.total_mass(rho, r, th))
    W = abs(diag.gravitational_energy(rho, Phi, r, th))
    R_eq, R_pol = diag.equatorial_polar_radii(H, r, th)
    print(f"field-free background: M = {M:.4f} Msun, R_eq = {R_eq:.3e} cm, "
          f"|W| = {W:.4e} erg\n")

    varpi = r[:, None] * np.sin(th)[None, :]

    # reference scales, from the purely poloidal solution
    omega2 = varpi ** 2
    u0 = solve_gradshafranov(-4.0 * np.pi * omega2 * rho * K0, r, th, lmax=LMAX)
    u_s0 = surface_flux(u0, H, r, th)
    u_norm = u0.max() - u_s0
    varpi_ref = 0.5 * R_eq
    Br0, Bth0 = diag.poloidal_field(u0, r, th)
    print(f"poloidal reference: u_max = {u0.max():.4e}, u_s = {u_s0:.4e}, "
          f"u_norm = {u_norm:.4e}")
    print(f"                    max|B_pol| = {np.hypot(Br0, Bth0).max():.4e} G, "
          f"varpi_ref = {varpi_ref:.3e} cm")
    closed = ((u0 > u_s0) & (H > 0.0))
    vol_closed = diag.volume_integral(closed.astype(float), r, th)
    vol_star = diag.volume_integral((H > 0.0).astype(float), r, th)
    print(f"                    closed-line region = "
          f"{100 * vol_closed / vol_star:.1f}% of the stellar volume\n")

    rows = []
    for zeta in ZETAS:
        print(f"zeta = {zeta}")
        print("   B_ref (G)  iters  conv   E_tor/E_mag   E_tor/|W|   "
              "E_pol/|W|   max|B_t| (G)   VE")
        for B_ref in B_REFS:
            beta_0 = B_ref * varpi_ref
            u, u_s, it, ok = solve_mixed(rho, H, r, th, K0, beta_0, u_norm,
                                         zeta)
            Br, Bth = diag.poloidal_field(u, r, th)
            beta, _ = beta_of(u, u_s, H, beta_0, u_norm, zeta)
            Bphi = np.where(varpi > 0, beta / np.maximum(varpi, 1e-30), 0.0)
            E_pol, E_tor, E_mag = diag.magnetic_energies(Br, Bth, Bphi, r, th)
            VE = diag.virial_error(rho, Phi, H, Br, Bth, Bphi, r, th, MU_E)[0]
            frac = E_tor / E_mag if E_mag > 0 else 0.0
            rows.append((zeta, B_ref, it, ok, frac, E_tor / W, E_pol / W,
                         np.abs(Bphi).max(), VE))
            print(f"   {B_ref:9.1e}  {it:5d}  {'y' if ok else 'N':^4s}  "
                  f"{frac:11.4f}   {E_tor / W:9.5f}   {E_pol / W:9.5f}   "
                  f"{np.abs(Bphi).max():.3e}   {VE:.2e}")
        print()

    best = max((row for row in rows if row[3]), key=lambda t: t[4], default=None)
    if best is not None:
        print(f"highest converged toroidal fraction: E_tor/E_mag = "
              f"{best[4]:.4f} at zeta = {best[0]}, B_ref = {best[1]:.2e} G")

    with OUT.open("w") as f:
        f.write(f"# barotropic twisted torus on a field-free star, "
                f"rho_c={RHO_C:.3e}, mu_e={MU_E}, k0={K0:.3e}\n")
        f.write(f"# background M={M:.4f} Msun, |W|={W:.6e} erg\n")
        f.write("zeta,B_ref_G,iters,converged,E_tor_over_Emag,E_tor_over_W,"
                "E_pol_over_W,Bphi_max_G,VE\n")
        for z, lam, it, ok, fr, etw, epw, bm, ve in rows:
            f.write(f"{z},{lam:.6e},{it},{ok},{fr:.6e},{etw:.6e},{epw:.6e},"
                    f"{bm:.6e},{ve:.6e}\n")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
