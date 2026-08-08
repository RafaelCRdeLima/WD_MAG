"""Campaign CT -- confined-torus current: a poloidal field born inside the star.

Run:  scf/.venv/bin/python3 investigations/export_ct_model.py --scan
      scf/.venv/bin/python3 investigations/export_ct_model.py --frac 0.3 --n 1 --bpol 1e13

WHY CT EXISTS, and why TT failed.

Campaign TT raised the poloidal amplitude of the SAME field we have been
evolving -- a vacuum dipole, produced by a Grad-Shafranov source with f(u)
constant. At B_pole = 3e12 G the star went to negative density at t = 0.06 s,
because that field extends outside the star and the ambient cannot hold it:
plasma beta was 0.088 in the envelope and 1.3e-4 in the ambient. The limit that
bites is not the Landau field, it is beta, and it bites at B ~ 3e10 G in the
ambient. A twisted torus is unreachable that way.

The fix is geometric, and it is standard practice for neutron stars although
apparently unused for white dwarfs. In the Grad-Shafranov equation the
azimuthal current is a free function of the flux. Choosing it CONSTANT -- what
we have been doing -- spreads current out to the surface and gives a vacuum
dipole outside. Choosing it to vanish below a threshold flux confines the
current to an interior torus, enlarges the closed-field-line region, and in the
limit leaves the poloidal field entirely inside the star. Ciolfi & Rezzolla
(2013, arXiv:1306.2803) use exactly this to reach toroidal-dominated twisted
tori for neutron stars; our problem is the mirror of theirs, and the lever is
the same. See also Pili & Bucciantini (2014, arXiv:1401.4308), whose XNS is the
code the family this configuration comes from was built with.

WHAT IT BUYS BEYOND THE TWISTED TORUS.

lambda_MRI = 2 pi v_Az / Omega scales with the VERTICAL field, so it is
adjustable by construction and not only by refinement. With B_z ~ 1e13 G
confined inside, lambda_MRI reaches ~11 cells at 256^3, against 0.004 cells
today. Stated honestly: that clears the six-cells-per-wavelength bar for the
LINEAR phase and misses the Q ~ 15-20 that MRI turbulence is normally held to.
So CT can let us WATCH the MRI grow and measure its linear rate; it cannot
support a claim about saturation or MRI-driven transport. Raising B_z further
improves Q and shrinks the number of wavelengths across the star at the same
time -- at 1.8e13 G, Q ~ 20 but lambda/R ~ 0.48, two wavelengths -- so the
useful window is narrow.

HOW THE SOLVE WORKS.

Making f depend on u makes the equation nonlinear, so this iterates: solve,
rebuild the source from the new flux, solve again, under-relaxed. The threshold
is set as a FRACTION of the running flux maximum, which makes the whole system
scale-invariant: the converged shape depends only on (frac, n) and the
amplitude is then a single linear rescaling. That is why the shape is converged
once and the amplitude fixed afterwards, rather than iterating on both.

PREDICTION, REGISTERED BEFORE THE FIRST RUN: E_tor/E_pol between 1 and 10, peak
field under B_c, beta > 1 everywhere inside, and Q between 8 and 15. If Q lands
below 6 the MRI will not appear even linearly, and CT reduces to the
twisted-torus campaign -- still worth running, for a different reason.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for _p in (REPO / "scf", REPO / "dashboard"):
    sys.path.insert(0, str(_p))
warnings.filterwarnings("ignore")

import diagnostics as diag                            # noqa: E402
import eos                                            # noqa: E402
import scf as scf_mod                                 # noqa: E402
import units                                          # noqa: E402
from gradshafranov import solve_gradshafranov         # noqa: E402
from seed import r_guess                              # noqa: E402
from sweep_worker import _solve_toroidal_certified    # noqa: E402
from terms.rotation import Rotation                   # noqa: E402
from terms.toroidal_sc import ToroidalSC              # noqa: E402

# Identical to the evolved configuration, so CT differs from it in the field
# geometry and nothing else.
RHO_C, MU_E = 3.0e9, 2.0
LMAX = 16
OMEGA_FRAC, A_FRAC = 1.5, 1.0
K_TOR, M_TOR = 5.0e-4, 1.0
B_C = 4.414e13
DX192, DX256 = 9.375e6, 7.031e6


def confined_flux(rho, r, th, varpi, m, eps=1.0e-6, lmax=LMAX,
                  iters=300, relax=0.3, tol=1e-9, verbose=False):
    """Solve Delta* u = -4 pi varpi^2 rho mu(u) with a NEGATIVE-POWER mu.

        mu(u) = (u + eps)^m

    This is the prescription of Fujisawa, Yoshida & Eriguchi (2012,
    arXiv:1204.5830), who reach interior fields two orders of magnitude above
    the surface with it. The localisation comes from m, not from any cutoff:
    with m < -1 the magnetic potential grows without bound as the flux falls
    towards the axis, so the poloidal field lines crowd inwards; with m > -1
    they spread out; m = 0 is the constant-mu case every earlier work used --
    and the one we had been using, which is why our field was a vacuum dipole.

    A FIRST ATTEMPT AT THIS USED THE OPPOSITE SIGN OF THE IDEA: a threshold
    that switched the current OFF where the flux was low. It gave an interior
    to exterior field ratio of only 3, against the ~100 needed, because a
    localised single-signed current still carries a dipole moment and still
    produces an exterior dipole. Concentrating the current is not the same as
    cutting it off, and the literature concentrates it.

    The threshold form does belong in this problem, but on the OTHER arbitrary
    function -- the one that carries the toroidal field, which must vanish
    outside the last closed line or the toroidal field leaks into the vacuum.
    Our toroidal field comes from ToroidalSC with B_phi proportional to rho,
    which already vanishes at the surface, so that side is covered.

    Because mu depends on u, the equation is nonlinear and this iterates,
    under-relaxed. u is returned normalised to unit maximum: the shape depends
    on m alone and the amplitude is a separate linear factor.
    """
    base = -4.0 * np.pi * varpi ** 2 * rho
    u = solve_gradshafranov(base, r, th, lmax=lmax)   # the mu = const seed
    umax = np.abs(u).max()
    if umax <= 0:
        raise SystemExit("the seed solve returned a null flux")
    u = u / umax

    err = np.inf
    for it in range(iters):
        mu = (np.maximum(u, 0.0) + eps) ** m
        new = solve_gradshafranov(base * mu, r, th, lmax=lmax)
        nmax = np.abs(new).max()
        if not np.isfinite(nmax) or nmax <= 0:
            raise SystemExit(f"solve diverged at m = {m}")
        new = new / nmax
        err = np.abs(new - u).max()
        u = (1.0 - relax) * u + relax * new
        u = u / np.abs(u).max()
        if err < tol:
            if verbose:
                print(f"      convergiu em {it + 1} iteracoes")
            break
    else:
        print(f"      AVISO: nao convergiu em {iters}; err = {err:.2e}")
    return u, err


def report(rho, r, th, varpi, H, W, u_shape, amp, Bphi, label=""):
    """Diagnostics for one amplitude of one converged shape."""
    u = u_shape * amp
    Br, Bth = diag.poloidal_field(u, r, th)
    E_pol, E_tor, E_mag = diag.magnetic_energies(Br, Bth, Bphi, r, th)
    bpol = np.hypot(Br, Bth)
    bpol_max = float(bpol.max())
    btot = float(np.sqrt(Br ** 2 + Bth ** 2 + Bphi ** 2).max())
    bp_ext = diag.surface_dipolarity(bpol, H, r, th)["B_pole"]

    # B_z in the meridional plane; the MRI wavelength is built from it
    Bz = Br * np.cos(th)[None, :] - Bth * np.sin(th)[None, :]
    inside = rho > 1.0e7
    Bz_typ = float(np.sqrt((Bz[inside] ** 2 * rho[inside]).sum()
                           / rho[inside].sum()))       # mass-weighted rms

    # plasma beta, from the code's own ztwd pressure
    x = eos.x_of_enthalpy(np.maximum(H, 0.0), mu_e=MU_E)
    P = eos.pressure(x)
    Pmag = (Br ** 2 + Bth ** 2 + Bphi ** 2) / (8.0 * np.pi)
    star = rho > 1.0e6
    beta_min = float((P[star] / np.maximum(Pmag[star], 1e-300)).min())

    return dict(E_tor_over_E_pol=E_tor / max(E_pol, 1e-99),
                E_pol_over_E_mag=E_pol / E_mag,
                bpol_max=bpol_max, btot_over_Bc=btot / B_C,
                bp_ext=bp_ext, ratio_int_ext=bpol_max / max(bp_ext, 1e-300),
                Bz_typ=Bz_typ, beta_min=beta_min,
                E_pol_over_W=E_pol / W)


def mri(Bz, rho_typ, Omega):
    vA = Bz / np.sqrt(4.0 * np.pi * rho_typ)
    return 2.0 * np.pi * vA / Omega


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", action="store_true",
                    help="sweep the localisation exponent m")
    ap.add_argument("-m", type=float, default=-1.5,
                    help="localisation exponent; m < -1 concentrates inwards")
    ap.add_argument("--bpol", type=float, default=1.0e13,
                    help="target interior peak |B_pol| in G")
    a = ap.parse_args()

    ref, r0, th0, _ = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=0.0, m_tor_sc=M_TOR,
        rotation=None, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=LMAX,
        tol=1e-8, max_iter=400)
    if ref is None:
        raise SystemExit("the non-rotating reference did not converge")
    M_ref = scf_mod.total_mass(ref["rho"], r0, th0)
    R_ref = diag.equatorial_polar_radii(ref["H"], r0, th0)[0]
    om_kep = float(np.sqrt(units.G_CONST * M_ref / R_ref ** 3))
    rot = Rotation(OMEGA_FRAC * om_kep, A_FRAC * R_ref)

    res, r, th, _ = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=K_TOR, m_tor_sc=M_TOR,
        rotation=rot, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=LMAX,
        tol=1e-8, max_iter=400)
    if res is None:
        raise SystemExit("the rotating solve did not converge")

    rho, Phi, H = res["rho"], res["Phi"], res["H"]
    varpi = r[:, None] * np.sin(th)[None, :]
    W = abs(diag.gravitational_energy(rho, Phi, r, th))
    R_eq, R_pol = diag.equatorial_polar_radii(H, r, th)
    Bphi = ToroidalSC(K=K_TOR, m=M_TOR).B_phi(rho, varpi)
    rho_typ = float(rho[rho > 1e7].mean())
    Om_c = float(rot.Omega_c)
    print(f"star: M = {units.g_to_msun(scf_mod.total_mass(rho, r, th)):.4f} Msun, "
          f"R_eq = {R_eq:.3e} cm, Omega_c = {Om_c:.3f} rad/s")
    print(f"      max|B_phi| = {np.abs(Bphi).max():.3e} G, "
          f"rho tipico no interior = {rho_typ:.2e}\n")

    shapes = ([0.0, -0.5, -1.0, -1.3, -1.5, -1.8, -2.0, -2.5]
              if a.scan else [a.m])

    # For each shape, find the LARGEST amplitude with beta_min >= 1 -- that is
    # the frontier the star can actually hold -- and report what it buys.
    hdr = (f"{'m':>6s} | {'B_pol pico':>11s} {'Et/Ep':>9s} {'|B|/Bc':>7s} "
           f"{'B_in/B_ex':>10s} {'beta_min':>9s} {'lam/dx256':>9s} {'lam/R':>6s}")
    print(hdr); print("-" * len(hdr))
    for m in shapes:
        try:
            u_shape, _ = confined_flux(rho, r, th, varpi, m)
        except SystemExit as e:
            print(f"{m:6.2f} | {e}")
            continue
        Br1, Bth1 = diag.poloidal_field(u_shape, r, th)
        unit = float(np.hypot(Br1, Bth1).max())
        lo, hi = 1e8, 1e15            # bisect on the interior peak |B_pol|
        for _ in range(60):
            mid = np.sqrt(lo * hi)
            if report(rho, r, th, varpi, H, W, u_shape, mid / unit,
                      Bphi)["beta_min"] >= 1.0:
                lo = mid
            else:
                hi = mid
        d = report(rho, r, th, varpi, H, W, u_shape, lo / unit, Bphi)
        L = mri(d["Bz_typ"], rho_typ, Om_c)
        print(f"{m:6.2f} | {lo:11.3e} {d['E_tor_over_E_pol']:9.3g} "
              f"{d['btot_over_Bc']:7.3f} {d['ratio_int_ext']:10.1f} "
              f"{d['beta_min']:9.3f} {L / DX256:9.2f} {L / R_eq:6.3f}")

    print("\nB_int/B_ext is the point of CT: the interior peak over the exterior")
    print("dipole. The vacuum-dipole construction gives ~1; large values mean")
    print("the field is confined and the ambient never sees it.")
    print("beta_min is over cells with rho > 1e6. Below 1 the field wins and")
    print("the star is pushed apart -- that is what killed TT at t = 0.06 s.")


if __name__ == "__main__":
    main()
