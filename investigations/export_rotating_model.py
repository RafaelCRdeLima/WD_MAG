"""Export a rotating, magnetized, BAROTROPIC white dwarf for Castro.

The target is the run this whole line of work has been heading towards: a
2 Msun white dwarf with differential rotation and a toroidal-dominated
interior beside a weak exterior dipole --- the collaboration's specification
--- evolved in full MHD to see whether it survives.

Why barotropic. Castro's EOS here is ztwd, P = P(rho) with mu_e fixed at 2.
A non-barotropic model's extra pressure is exactly what ztwd discards, so it
would arrive out of equilibrium by precisely the support under test. This
project has already made that mistake once, building under ztwd and evolving
under gamma_law; the star collapsed from t = 0 and a field-free control
collapsed identically. So the initial condition is built barotropically, in
genuine equilibrium under the same EOS that will evolve it.

Why it is worth running anyway. The 2 Msun purely toroidal configuration of
the companion paper collapsed in about 3.5 dynamical times, but it was not
rotating. Adding differential rotation adds a mechanism that was absent:
a differentially rotating magnetized star is MRI-unstable, so the shear that
supplies the mass may destroy itself. Neither outcome is known in advance,
which is what makes the run worth its time.

Rotation cannot ride in on the vector potential, and Castro's own rotation
support is a rotating frame at constant Omega, which a j-constant law is not.
So it is exported as an initial velocity field (model format v2) and left to
evolve.

Run:  scf/.venv/bin/python3 investigations/export_rotating_model.py
"""

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
import scf as scf_mod                                 # noqa: E402
import units                                          # noqa: E402
from axisym_model_writer import (to_meridional, vector_potential,  # noqa: E402
                                 verify_curl_on_cartesian,
                                 verify_meridional_curl, write_model)
from gradshafranov import solve_gradshafranov         # noqa: E402
from seed import r_guess                              # noqa: E402
from sweep_worker import _solve_toroidal_certified    # noqa: E402
from terms.rotation import Rotation                   # noqa: E402
from terms.toroidal_sc import ToroidalSC              # noqa: E402

RHO_C, MU_E = 3.0e9, 2.0
LMAX = 16
N_MER = 385
HALF_CM = 9.0e8
CORNER = 1.7320508

# Rotation: the heaviest point of investigations/rotating_barotropic_scan.py
# that stayed well below mass shedding. Omega_c is quoted against the
# Keplerian frequency of the non-rotating star at the same central density.
OMEGA_FRAC, A_FRAC = 1.5, 1.0

# Field: the collaboration's specification -- a toroidal-dominated interior
# with a weak exterior dipole. k0 is calibrated to the surface dipole, which
# is the observable, and K_TOR is set from a sweep: at 5e-4 the star reaches
# 2.005 Msun with max|B| = 0.73 B_c, while 1e-3 reaches 2.269 Msun but at
# 1.24 B_c, outside the range where a zero-temperature unquantised equation
# of state is valid.
#
# Note what this specification implies. A 1e9 G dipole beside a 3e13 G
# toroidal field is B_t/B_p ~ 4e3: the field is toroidal to within a part in
# a thousand, and a purely toroidal field is Tayler-unstable -- it is the
# configuration that collapsed in about 3.5 dynamical times in the companion
# paper. What is different here, and the reason the run is worth its time, is
# that this star rotates, and the collapsing one did not.
K0_REF = 1.0e-13
B_POLE_TARGET = 1.0e9
K_TOR, M_TOR = 5.0e-4, 1.0

DIV_GATE = 1.0e-12
CURL_GATE = 5.0e-2
SHED_GATE = 0.95
OUTDIR = REPO / "models"


def build(rho, r, th, k0, varpi):
    u = solve_gradshafranov(-4.0 * np.pi * varpi ** 2 * rho * k0, r, th,
                            lmax=LMAX)
    Bphi = ToroidalSC(K=K_TOR, m=M_TOR).B_phi(rho, varpi)
    Br, Bth = diag.poloidal_field(u, r, th)
    return u, Bphi, Br, Bth


def main():
    OUTDIR.mkdir(exist_ok=True)

    # the non-rotating star at the same rho_c, only to set Omega_K
    ref, r0, th0, _ = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=0.0, m_tor_sc=M_TOR,
        rotation=None, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=LMAX,
        tol=1e-8, max_iter=400)
    if ref is None:
        raise SystemExit("the non-rotating reference did not converge")
    M_ref = scf_mod.total_mass(ref["rho"], r0, th0)
    R_ref = diag.equatorial_polar_radii(ref["H"], r0, th0)[0]
    om_kep = float(np.sqrt(units.G_CONST * M_ref / R_ref ** 3))
    print(f"reference: M = {units.g_to_msun(M_ref):.4f} Msun, "
          f"R_eq = {R_ref:.4e} cm, Omega_K = {om_kep:.4e} rad/s")

    rot = Rotation(OMEGA_FRAC * om_kep, A_FRAC * R_ref)
    res, r, th, ov = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=K_TOR, m_tor_sc=M_TOR,
        rotation=rot, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=LMAX,
        tol=1e-8, max_iter=400)
    if res is None:
        raise SystemExit("the rotating solve did not converge")

    rho, Phi, H = res["rho"], res["Phi"], res["H"]
    M = units.g_to_msun(scf_mod.total_mass(rho, r, th))
    W = abs(diag.gravitational_energy(rho, Phi, r, th))
    T = rot.energy(rho, r, th)["T"]
    varpi = r[:, None] * np.sin(th)[None, :]
    R_eq, R_pol = diag.equatorial_polar_radii(H, r, th)

    jeq = len(th) // 2
    kk = np.flatnonzero(rho[:, jeq] > 0)
    om_eq = float(np.atleast_1d(rot.Omega(np.array([R_eq])))[0])
    grav_eq = abs(float(np.gradient(Phi[:, jeq], r)[kk[-1]]))
    shed = om_eq ** 2 * R_eq / max(grav_eq, 1e-30)
    print(f"rotating:  M = {M:.4f} Msun, R_eq = {R_eq:.4e}, "
          f"R_pol = {R_pol:.4e} cm")
    print(f"           T/|W| = {T / W:.4f}, shedding {shed:.3f} "
          f"(gate {SHED_GATE})")
    if shed >= SHED_GATE:
        raise SystemExit("the configuration is shedding mass")

    # calibrate k0 for the target surface dipole: B_pole is linear in k0
    _, _, Br0, Bth0 = build(rho, r, th, K0_REF, varpi)
    bp_ref = diag.surface_dipolarity(np.hypot(Br0, Bth0), H, r, th)["B_pole"]
    k0 = K0_REF * (B_POLE_TARGET / bp_ref)
    print(f"calibration: B_pole = {bp_ref:.4e} G at k0 = {K0_REF:.3e}"
          f"  ->  k0 = {k0:.4e} for {B_POLE_TARGET:.1e} G\n")

    u, Bphi, Br, Bth = build(rho, r, th, k0, varpi)
    E_pol, E_tor, E_mag = diag.magnetic_energies(Br, Bth, Bphi, r, th)
    bp = diag.surface_dipolarity(np.hypot(Br, Bth), H, r, th)["B_pole"]
    amp = np.abs(Bphi).max() / max(np.hypot(Br, Bth).max(), 1e-300)

    rmax = 1.02 * CORNER * HALF_CM
    vp = np.linspace(0.0, rmax, N_MER)
    zz = np.linspace(-rmax, rmax, 2 * N_MER - 1)
    rho_m, u_m, bphi_m = to_meridional(r, th, (rho, u, Bphi), vp, zz)
    A_phi, A_z = vector_potential(vp, u_m, bphi_m)

    # v_phi = Omega(varpi) varpi depends on varpi alone, so it is evaluated
    # directly on the meridional grid rather than interpolated, and is cut
    # off outside the star: the ambient floor must not be spun.
    v_phi = (np.atleast_1d(rot.Omega(vp))[:, None] * vp[:, None]
             * np.ones((1, len(zz))))
    v_phi = np.where(rho_m > 0.0, v_phi, 0.0)

    err_pol, err_tor = verify_meridional_curl(vp, zz, A_phi, A_z, u_m, bphi_m)
    rel_div, b_max, dx = verify_curl_on_cartesian(
        vp, zz, A_phi, A_z, half=HALF_CM, n_cart=64)

    b_tot_max = float(np.sqrt(Br**2 + Bth**2 + Bphi**2).max())
    print(f"field: B_pole = {bp:.4e} G, E_pol/|W| = {E_pol / W:.3e}, "
          f"E_tor/E_mag = {E_tor / E_mag:.6f}")
    print(f"       max|B_phi| = {np.abs(Bphi).max():.4e} G, "
          f"max|B|/B_c = {b_tot_max / 4.414e13:.3f}")
    print(f"       B_t/B_p = {amp:.4g} (amplitude)")
    print(f"       curl A vs B: poloidal {err_pol:.3e}, toroidal "
          f"{err_tor:.3e}   (gate {CURL_GATE:.0e})")
    print(f"       div B on 64^3: {rel_div:.3e}   (gate {DIV_GATE:.0e}), "
          f"amplitude retained "
          f"{100 * b_max / np.abs(Bphi).max():.1f}%")
    print(f"       max v_phi = {v_phi.max():.4e} cm/s")

    retained = b_max / np.abs(Bphi).max()
    if retained > 1.02:
        raise SystemExit(f"reconstruction gained amplitude "
                         f"({100 * retained:.1f}%) -- model grid too small")
    if not (rel_div < DIV_GATE):
        raise SystemExit("divergence gate failed")
    if not (err_pol < CURL_GATE and err_tor < CURL_GATE):
        raise SystemExit("curl gate failed")

    params = dict(rho_c=RHO_C, mu_e=MU_E, K_tor=K_TOR, m_tor=M_TOR, k0=k0,
                  M_msun=M, R_eq_cm=R_eq, R_pol_cm=R_pol, B_pole_G=bp,
                  E_pol_over_W=E_pol / W, E_tor_over_Emag=E_tor / E_mag,
                  Bphi_max_G=float(np.abs(Bphi).max()),
                  B_total_max_over_Bc=b_tot_max / 4.414e13,
                  Bt_over_Bp_amplitude=amp, omega_frac=OMEGA_FRAC,
                  A_over_Req=A_FRAC, Omega_c=rot.Omega_c, A_cm=rot.A,
                  T_over_W=T / W, shedding=shed,
                  v_phi_max_cms=float(v_phi.max()))
    checks = dict(curl_err_poloidal=err_pol, curl_err_toroidal=err_tor,
                  relative_divB_64cubed=rel_div,
                  amplitude_retained_64cubed=retained)
    man = write_model(vp, zz, rho_m, A_phi, A_z,
                      OUTDIR / "rotating_mixed.txt", params, checks,
                      v_phi=v_phi)
    print(f"\nwrote models/{man['file']} ({man['n_varpi']}x{man['n_z']}), "
          f"format {man['format']}")


if __name__ == "__main__":
    main()
