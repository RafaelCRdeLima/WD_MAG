"""A non-barotropic mixed equilibrium: solve for the matter, not the field.

Run:  scf/.venv/bin/python3 investigations/nonbarotropic_equilibrium.py

THE FORMULATION

Axisymmetric magnetostatic equilibrium, no flows. The field is written
through the flux function u,

    B_pol = (1/varpi) grad u x e_phi,     B_phi = beta(u)/varpi

and the Lorentz force is then everywhere parallel to grad u,

    f_L = g grad u,    g = -[Delta* u + beta beta'(u)] / (4 pi varpi^2).

Barotropy would force g = rho M'(u). We drop that and keep everything else.
What remains is

    grad P = f_L - rho grad Phi                                       (1)

for which a single-valued P exists iff the curl vanishes. In spherical
coordinates that is one scalar condition,

    rho_r Phi_theta - rho_theta Phi_r  =  g_r u_theta - g_theta u_r == Q  (2)

and this is the whole content of the non-barotropic freedom. Equation (2)
is a linear first-order PDE for rho whose characteristics are the
equipotentials of Phi: along Phi = const, with theta as the parameter,

    dr/dtheta      = -Phi_theta / Phi_r
    drho/dtheta    = -Q / Phi_r                                        (3)

So rho is NOT free -- it is determined on each equipotential once its value
at one point of that equipotential is chosen. The free function is rho on
the polar axis, one value per equipotential: the barotropic backbone. The
baroclinic part is what Eq. (3) adds on top of it.

WHY THIS FORM AND NOT ANOTHER

g contains Delta* u, and computing that by finite differences and then
taking another derivative for Q would put third derivatives of a numerical
field into the source. Instead u is obtained by INVERTING Delta* on a
chosen source, so Delta* u is known exactly -- it is the source itself --
and Q needs only first derivatives. That is also the honest statement of
what non-barotropy buys: the poloidal source is prescribed rather than tied
to the density.

STATUS

Validated: run with --nofield, seeded from a field-free star, the scheme
returns that star's mass to 0.07% (1.3474 against 1.3464 Msun) and holds it
over the iteration. The characteristic marching and the free-boundary cut
are therefore doing what they should on the global invariant.

NOT yet working, both stated rather than buried:

  1. rho = (g u_r - P_r)/Phi_r is singular at the centre, where Phi_r -> 0.
     About 0.06% of cells come out negative there, with excursions of order
     100 rho_max. They carry negligible volume -- the mass is right to
     0.07% -- but they make any LOCAL diagnostic untrustworthy.

  2. Consequently the stratification diagnostic below, the spread of P over
     surfaces of constant rho, reads 0.25-0.27 even in the --nofield limit
     where it must be exactly zero. It is measuring the noise of (1), not
     the physics, and no number from it should be quoted yet.

AND A RESULT THAT IS NOT A BUG

The free function -- P on the polar axis, one value per equipotential -- is
genuine physical freedom, not a gauge choice to be fixed later. A
non-barotropic star's stratification is set by its thermal and compositional
history, which the equilibrium equations do not know. So for one prescribed
field there is a one-function FAMILY of equilibria with different masses,
and the mass is not predicted by the field alone. Seeding that function from
the 2 Msun toroidal star while prescribing a different (confined) field is
what gives 5.8 Msun below: not an error, a different member of the family.
Closing the system needs a physical stratification condition -- convective
stability is the obvious candidate -- not a better solver.

Writes nonbarotropic_equilibrium.csv.
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
from poisson import solve_poisson                # noqa: E402
from scipy.interpolate import RegularGridInterpolator as RGI   # noqa: E402
from seed import r_guess                         # noqa: E402
from sweep_worker import _solve_toroidal_certified   # noqa: E402
from terms.toroidal_sc import ToroidalSC         # noqa: E402

RHO_C, MU_E = 1.0e9, 2.0
K_TOR, M_TOR = 3.245e-3, 1.0     # the M = 2 Msun crossing, for the target
K0 = 1.0e-13                     # poloidal source amplitude
ZETA = 1.0
LMAX = 16
N_OUTER = 25                     # rho <-> Phi iterations
RELAX = 0.25                     # without it the fixed point oscillates 4<->8 Msun
N_THETA_MARCH = 721              # theta steps along each equipotential
OUT = HERE / "nonbarotropic_equilibrium.csv"

# Round-trip test. With the field off, dP = g du = 0 along equipotentials, so
# P is constant on them and rho = -P_r/Phi_r is plain hydrostatic balance: the
# scheme must return the star it was seeded from. This has to seed from a
# FIELD-FREE star -- seeding the 2 Msun star's pressure profile and then
# removing the field that inflates it is not a null test, it is a different
# physical question, and it was giving 6.6 Msun for that reason.
FIELD_OFF = "--nofield" in sys.argv


def surface_flux(u, H_or_rho, th):
    vals = []
    for j in range(len(th)):
        inside = H_or_rho[:, j] > 0.0
        if inside.any():
            vals.append(u[np.flatnonzero(inside)[-1], j])
    return max(vals) if vals else 0.0


def grads(F, r, th):
    """(dF/dr, dF/dtheta) on the (r, theta) grid."""
    return np.gradient(F, r, axis=0), np.gradient(F, th, axis=1)


def transport_pressure(P_axis, Phi, g, u, r, th, r_seeds):
    """Integrate dP = g du along equipotentials, from the pole outward.

    Along a surface of constant Phi gravity does no work, so the momentum
    equation collapses to dP/dtheta|_Phi = g du/dtheta|_Phi. Carrying P
    rather than rho is what makes the free boundary tractable: P is the
    quantity that has to vanish at the surface, and rho follows from
    rho = (g u_r - P_r)/Phi_r afterwards.

    The free function is P on the polar axis, one value per equipotential.
    """
    Phi_r, Phi_th = grads(Phi, r, th)
    u_r, u_th = grads(u, r, th)
    I = lambda A: RGI((r, th), A, bounds_error=False, fill_value=None)
    f_Phr, f_Phth, f_ur, f_uth, f_g = (I(Phi_r), I(Phi_th), I(u_r), I(u_th),
                                       I(g))

    th_m = np.linspace(th[0], th[-1], N_THETA_MARCH)
    dth = th_m[1] - th_m[0]
    R = np.empty((len(r_seeds), len(th_m)))
    Pm = np.empty_like(R)
    R[:, 0] = r_seeds
    Pm[:, 0] = np.interp(r_seeds, r, P_axis)

    def deriv(rr, tt):
        tt = np.full_like(rr, tt) if np.isscalar(tt) else tt
        pts = np.stack([np.clip(rr, r[0], r[-1]),
                        np.clip(tt, th[0], th[-1])], axis=-1)
        phr = f_Phr(pts)
        bad = np.abs(phr) < 1e-30
        phr = np.where(bad, 1.0, phr)
        drdt = np.where(bad, 0.0, -f_Phth(pts) / phr)
        dudt = f_uth(pts) + f_ur(pts) * drdt        # du/dtheta along Phi
        return drdt, np.where(bad, 0.0, f_g(pts) * dudt)

    for i in range(len(th_m) - 1):
        t0, r0, p0 = th_m[i], R[:, i], Pm[:, i]
        k1r, k1p = deriv(r0, t0)
        k2r, k2p = deriv(r0 + 0.5 * dth * k1r, t0 + 0.5 * dth)
        k3r, k3p = deriv(r0 + 0.5 * dth * k2r, t0 + 0.5 * dth)
        k4r, k4p = deriv(r0 + dth * k3r, t0 + dth)
        R[:, i + 1] = np.clip(r0 + (dth / 6.0) * (k1r + 2*k2r + 2*k3r + k4r),
                              r[0], r[-1])
        Pm[:, i + 1] = p0 + (dth / 6.0) * (k1p + 2*k2p + 2*k3p + k4p)

    P = np.zeros((len(r), len(th)))
    for j, tj in enumerate(th):
        i = np.argmin(np.abs(th_m - tj))
        order = np.argsort(R[:, i])
        rr, pp = R[order, i], Pm[order, i]
        keep = np.concatenate([[True], np.diff(rr) > 0])
        P[:, j] = np.interp(r, rr[keep], pp[keep], left=pp[keep][0], right=0.0)
    return P


def main():
    # --- the target: the certified 2 Msun toroidal star -----------------
    res, r, th, _ = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C),
        K_tor=(0.0 if FIELD_OFF else K_TOR), m_tor_sc=M_TOR,
        rotation=None, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=LMAX,
        tol=1e-8, max_iter=200)
    if res is None:
        raise SystemExit("the reference solve did not converge")
    rho0, Phi0, H0 = res["rho"], res["Phi"], res["H"]
    M0 = units.g_to_msun(scf_mod.total_mass(rho0, r, th))
    W0 = abs(diag.gravitational_energy(rho0, Phi0, r, th))
    varpi = r[:, None] * np.sin(th)[None, :]
    Bphi_sc = ToroidalSC(K=K_TOR, m=M_TOR).B_phi(rho0, varpi)
    _, E_tor_target, _ = diag.magnetic_energies(
        np.zeros_like(rho0), np.zeros_like(rho0), Bphi_sc, r, th)
    print(f"reference 2 Msun star: M = {M0:.4f} Msun, |W| = {W0:.4e} erg")
    print(f"  target E_tor = {E_tor_target:.4e} erg = "
          f"{E_tor_target / W0:.4f} |W|\n")

    # --- the prescribed field -------------------------------------------
    src = -4.0 * np.pi * varpi ** 2 * rho0 * K0     # == Delta* u, exactly
    u = solve_gradshafranov(src, r, th, lmax=LMAX)
    u_s = surface_flux(u, rho0, th)
    u_norm = max(u.max() - u_s, 1e-300)
    closed = (rho0 > 0.0) & (u > u_s)
    w = np.where(closed, (u - u_s) / u_norm, 0.0)

    shape = np.where(varpi > 0, np.power(w, ZETA) / np.maximum(varpi, 1e-30), 0)
    _, E_unit, _ = diag.magnetic_energies(np.zeros_like(rho0),
                                          np.zeros_like(rho0), shape, r, th)
    beta_0 = np.sqrt(E_tor_target / E_unit)
    beta = beta_0 * np.power(w, ZETA)
    dbeta = np.where(closed, beta_0 * ZETA * np.power(w, ZETA - 1.0) / u_norm,
                     0.0)
    Bphi = np.where(varpi > 0, beta / np.maximum(varpi, 1e-30), 0.0)
    Br, Bth = diag.poloidal_field(u, r, th)
    E_pol, E_tor, E_mag = diag.magnetic_energies(Br, Bth, Bphi, r, th)
    print(f"prescribed field: max|B_phi| = {np.abs(Bphi).max():.4e} G, "
          f"max|B_pol| = {np.hypot(Br, Bth).max():.4e} G")
    print(f"  E_tor/|W| = {E_tor / W0:.4f}, E_pol/|W| = {E_pol / W0:.5f}, "
          f"E_tor/E_pol = {E_tor / E_pol:.1f}\n")

    if FIELD_OFF:
        g = np.zeros_like(rho0)
        u = np.zeros_like(rho0)
        E_pol = E_tor = 0.0
        print("  FIELD OFF: this must reproduce the barotropic star\n")
    # g = -[Delta* u + beta beta'] / (4 pi varpi^2); Delta* u IS src
    S = src + beta * dbeta
    g_full = np.where(varpi > 0, -S / (4.0 * np.pi * np.maximum(varpi, 1e-30) ** 2),
                 0.0)
    g = g_full if not FIELD_OFF else np.zeros_like(rho0)
    u = u if not FIELD_OFF else np.zeros_like(rho0)
    g_r, g_th = grads(g, r, th)
    u_r, u_th = grads(u, r, th)
    Q = g_r * u_th - g_th * u_r

    # --- iterate P, rho <-> Phi -----------------------------------------
    import eos                                            # noqa: E402
    R_pol_ref = diag.equatorial_polar_radii(H0, r, th)[1]
    r_seeds = np.linspace(r[1], R_pol_ref, 220)
    x0 = eos.x_of_enthalpy(np.maximum(H0, 0.0), MU_E)
    P_ref = np.where(H0 > 0, eos.pressure(x0), 0.0)
    P_axis = P_ref[:, 0]

    Phi = Phi0.copy()
    print("  it    M (Msun)   min P       min rho     frac rho<0   "
          "P-rho spread   |drho|/rho0max")
    hist = []
    for it in range(1, N_OUTER + 1):
        P = transport_pressure(P_axis, Phi, g, u, r, th, r_seeds)
        P_r, _ = grads(P, r, th)
        Phi_r, _ = grads(Phi, r, th)
        u_r, _ = grads(u, r, th)
        safe = np.abs(Phi_r) > 1e-30
        rho_new = np.where(safe, (g * u_r - P_r) / np.where(safe, Phi_r, 1.0),
                           0.0)
        # the surface: outward from the centre along each ray, the star ends
        # at the first sign change of P. Without this the baroclinic term
        # keeps making matter outside and the mass runs away.
        for j in range(len(th)):
            bad_j = np.flatnonzero(P[:, j] <= 0.0)
            if bad_j.size:
                rho_new[bad_j[0]:, j] = 0.0
                P[bad_j[0]:, j] = 0.0
        rho_new = np.where(P > 0.0, rho_new, 0.0)
        neg = float((rho_new < 0).sum()) / rho_new.size
        rho_pos = np.maximum(rho_new, 0.0)
        M = units.g_to_msun(scf_mod.total_mass(rho_pos, r, th))

        # how non-barotropic: spread of P within bins of rho
        spread = np.nan
        m = rho_pos > 1e-3 * rho_pos.max()
        if m.sum() > 100:
            lr = np.log10(rho_pos[m])
            bins = np.linspace(lr.min(), lr.max(), 25)
            idx = np.digitize(lr, bins)
            rel = []
            for b in range(1, len(bins)):
                sel = P[m][idx == b]
                sel = sel[sel > 0]
                if sel.size > 5:
                    rel.append((sel.max() - sel.min()) / np.median(sel))
            spread = float(np.median(rel)) if rel else np.nan

        drho = (np.abs(rho_pos - rho0).max() / max(rho0.max(), 1e-300))
        print(f"  {it:3d}   {M:8.4f}   {P.min():.3e}   {rho_new.min():.3e}   "
              f"{100*neg:8.2f}%   {spread:.3e}   {drho:.3e}")
        hist.append((it, M, P.min(), rho_new.min(), neg, spread))
        Phi = ((1.0 - RELAX) * Phi
               + RELAX * solve_poisson(rho_pos, r, th, lmax=LMAX))

    W = abs(diag.gravitational_energy(rho_pos, Phi, r, th))
    print(f"\nfinal: M = {M:.4f} Msun (target {M0:.4f}), |W| = {W:.4e} erg")
    if not FIELD_OFF:
        print(f"       E_tor/|W| = {E_tor / W:.4f}, E_pol/|W| = {E_pol / W:.5f}")

    with OUT.open("w") as f:
        f.write(f"# non-barotropic mixed equilibrium, k0={K0:.3e}, "
                f"zeta={ZETA}, beta_0={beta_0:.6e}\n")
        f.write(f"# reference M={M0:.4f} Msun, |W|={W0:.6e}; "
                f"E_tor/|W|={E_tor / W0:.4f}, E_pol/|W|={E_pol / W0:.6f}\n")
        f.write("iter,M_msun,P_min,rho_min,frac_rho_negative,P_rho_spread\n")
        for row in hist:
            f.write(f"{row[0]},{row[1]:.6f},{row[2]:.6e},{row[3]:.6e},"
                    f"{row[4]:.6e},{row[5]:.6e}\n")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
