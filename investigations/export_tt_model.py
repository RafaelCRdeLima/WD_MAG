"""Campaign TT -- export a rotating star with a genuinely mixed field.

Run:  scf/.venv/bin/python3 investigations/export_tt_model.py --scan
      scf/.venv/bin/python3 investigations/export_tt_model.py --bpole 1.1e13

This is export_rotating_model.py with B_POLE_TARGET moved to the command line
and one claim dropped. The star, the rotation law, the toroidal field and the
whole export path are identical, so the two configurations differ in exactly
one parameter and their model files are comparable line for line.

WHY. The configuration evolved so far is toroidal to one part in 10^9 by
energy, which is the textbook Tayler-unstable case, and it behaved like one:
an m=1 mode destroyed the ordered field in a few Alfven times. Braithwaite-type
stable configurations are not toroidal-dominated -- they are twisted tori with
comparable poloidal and toroidal energy, and that geometry is the one the
equilibrium literature on super-Chandrasekhar white dwarfs never constructs.
Whether it survives is an IDEAL MHD question, so the answer does not depend on
the resistivity, which is the one thing this problem cannot get right: Rm ~ 6
on the finest grid affordable, against ~10^18 physically.

WHAT IS GIVEN UP. The poloidal field is imposed on a converged toroidal plus
rotation equilibrium, not solved alongside it, so the pair leaves virial
balance by its own poloidal energy. papers/wd-toroidal-poloidal measures the
cost on a non-rotating star and finds the virial error tracking E_pol/|W|.

The output is therefore NOT an equilibrium and must not be presented as one.
It is an initial condition to be relaxed, which is what that draft prescribes
for this range, and the accommodation phase in the inputs file is the
relaxation. Expect a larger initial transient than the toroidal-dominated run,
which already breathed by a factor 2.5 in volume.

--scan reports the trade-off and writes nothing. The star is solved once and
reused across the scan, so the extra points are nearly free.
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

# Unchanged from export_rotating_model.py.
RHO_C, MU_E = 3.0e9, 2.0
LMAX = 16
N_MER = 385
HALF_CM = 9.0e8
CORNER = 1.7320508
OMEGA_FRAC, A_FRAC = 1.5, 1.0
K0_REF = 1.0e-13
K_TOR, M_TOR = 5.0e-4, 1.0
DIV_GATE = 1.0e-12
CURL_GATE = 5.0e-2
SHED_GATE = 0.95
B_C = 4.414e13
OUTDIR = REPO / "models"

# 1.0e9 is the configuration already evolved. The rest walk towards a twisted
# torus; E_pol goes as B_pole^2, so the energy ratio falls as the square.
SCAN = [1.0e9, 1.0e12, 1.5e12, 2.0e12, 2.5e12, 3.0e12, 4.0e12, 5.0e12]


def build(rho, r, th, k0, varpi):
    u = solve_gradshafranov(-4.0 * np.pi * varpi ** 2 * rho * k0, r, th,
                            lmax=LMAX)
    Bphi = ToroidalSC(K=K_TOR, m=M_TOR).B_phi(rho, varpi)
    Br, Bth = diag.poloidal_field(u, r, th)
    return u, Bphi, Br, Bth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--bpole", type=float, default=None)
    a = ap.parse_args()
    if not a.scan and a.bpole is None:
        ap.error("give --scan or --bpole")

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

    rot = Rotation(OMEGA_FRAC * om_kep, A_FRAC * R_ref)
    res, r, th, _ = _solve_toroidal_certified(
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
    print(f"star: M = {M:.4f} Msun, R_eq = {R_eq:.4e} cm, "
          f"R_pol/R_eq = {R_pol/R_eq:.4f}")
    print(f"      T/|W| = {T/W:.4f}, shedding {shed:.3f} (gate {SHED_GATE})")
    if shed >= SHED_GATE:
        raise SystemExit("the configuration is shedding mass")

    # B_pole is linear in k0, and it is the SURFACE DIPOLE from
    # diag.surface_dipolarity, not the maximum of |B_r|. Getting that wrong is
    # what broke the first version of this script.
    _, _, Br0, Bth0 = build(rho, r, th, K0_REF, varpi)
    bp_ref = diag.surface_dipolarity(np.hypot(Br0, Bth0), H, r, th)["B_pole"]
    print(f"      calibration: B_pole = {bp_ref:.4e} G at k0 = {K0_REF:.3e}\n")

    targets = SCAN if a.scan else [a.bpole]
    print(f"{'B_pole (G)':>12s} {'E_tor/E_pol':>12s} {'B_t/B_p':>10s} "
          f"{'max|B|/B_c':>11s} {'E_pol/|W|':>11s}")
    keep = None
    for tgt in targets:
        k0 = K0_REF * (tgt / bp_ref)
        u, Bphi, Br, Bth = build(rho, r, th, k0, varpi)
        E_pol, E_tor, E_mag = diag.magnetic_energies(Br, Bth, Bphi, r, th)
        amp = np.abs(Bphi).max() / max(np.hypot(Br, Bth).max(), 1e-300)
        btot = float(np.sqrt(Br ** 2 + Bth ** 2 + Bphi ** 2).max())
        print(f"{tgt:12.2e} {E_tor/max(E_pol,1e-99):12.4g} {amp:10.4g} "
              f"{btot/B_C:11.3f} {E_pol/W:11.3e}")
        keep = (tgt, k0, u, Bphi, Br, Bth, E_pol, E_tor, E_mag, amp, btot)

    if a.scan:
        print("\nmax|B|/B_c is the binding constraint here, not the virial")
        print("error: the EOS is a zero-temperature unquantised one and is")
        print("not valid above the Landau field. Pick the most balanced")
        print("ratio that stays under 1, then relax it.")
        return

    tgt, k0, u, Bphi, Br, Bth, E_pol, E_tor, E_mag, amp, btot = keep
    bp = diag.surface_dipolarity(np.hypot(Br, Bth), H, r, th)["B_pole"]

    rmax = 1.02 * CORNER * HALF_CM
    vp = np.linspace(0.0, rmax, N_MER)
    zz = np.linspace(-rmax, rmax, 2 * N_MER - 1)
    rho_m, u_m, bphi_m = to_meridional(r, th, (rho, u, Bphi), vp, zz)
    A_phi, A_z = vector_potential(vp, u_m, bphi_m)

    # v_phi tapered across the sponge's own density bracket, exactly as in the
    # evolved model: a hard cut at the surface put 7.5e8 cm/s into one cell
    # beside a static ambient and killed the first run at t = 2.59 s.
    RHO_SPIN_LO, RHO_SPIN_HI = 1.0e4, 1.0e6
    t = np.clip((np.log10(np.maximum(rho_m, RHO_SPIN_LO))
                 - np.log10(RHO_SPIN_LO))
                / (np.log10(RHO_SPIN_HI) - np.log10(RHO_SPIN_LO)), 0.0, 1.0)
    v_phi = (np.atleast_1d(rot.Omega(vp))[:, None] * vp[:, None]
             * np.ones((1, len(zz))))
    v_phi = v_phi * (t * t * (3.0 - 2.0 * t))

    err_pol, err_tor = verify_meridional_curl(vp, zz, A_phi, A_z, u_m, bphi_m)
    rel_div, b_max, dx = verify_curl_on_cartesian(
        vp, zz, A_phi, A_z, half=HALF_CM, n_cart=64)

    print(f"\nfield: B_pole = {bp:.4e} G, E_tor/E_pol = {E_tor/E_pol:.4g}, "
          f"B_t/B_p = {amp:.4g}")
    print(f"       max|B_phi| = {np.abs(Bphi).max():.4e} G, "
          f"max|B|/B_c = {btot/B_C:.3f}")
    print(f"       curl A vs B: poloidal {err_pol:.3e}, toroidal {err_tor:.3e}"
          f"   (gate {CURL_GATE:.0e})")
    print(f"       div B on 64^3: {rel_div:.3e}   (gate {DIV_GATE:.0e}), "
          f"amplitude retained {100*b_max/np.abs(Bphi).max():.1f}%")
    print(f"       max v_phi = {v_phi.max():.4e} cm/s")
    print("       NOT an equilibrium -- see the module docstring. Relax it.")

    retained = b_max / np.abs(Bphi).max()
    if retained > 1.02:
        raise SystemExit(f"reconstruction gained amplitude "
                         f"({100*retained:.1f}%) -- model grid too small")
    if not (rel_div < DIV_GATE):
        raise SystemExit("divergence gate failed")
    if not (err_pol < CURL_GATE and err_tor < CURL_GATE):
        raise SystemExit("curl gate failed")

    params = dict(rho_c=RHO_C, mu_e=MU_E, K_tor=K_TOR, m_tor=M_TOR, k0=k0,
                  M_msun=M, R_eq_cm=R_eq, R_pol_cm=R_pol, B_pole_G=bp,
                  E_pol_over_W=E_pol / W, E_tor_over_Emag=E_tor / E_mag,
                  E_tor_over_E_pol=E_tor / E_pol,
                  Bphi_max_G=float(np.abs(Bphi).max()),
                  B_total_max_over_Bc=btot / B_C,
                  Bt_over_Bp_amplitude=amp, omega_frac=OMEGA_FRAC,
                  A_over_Req=A_FRAC, Omega_c=rot.Omega_c, A_cm=rot.A,
                  T_over_W=T / W, shedding=shed,
                  v_phi_max_cms=float(v_phi.max()),
                  is_equilibrium=False,
                  note="poloidal field imposed, not solved; relax before use")
    checks = dict(curl_err_poloidal=err_pol, curl_err_toroidal=err_tor,
                  relative_divB_64cubed=rel_div,
                  amplitude_retained_64cubed=retained)
    man = write_model(vp, zz, rho_m, A_phi, A_z,
                      OUTDIR / "mixed_tt.txt", params, checks, v_phi=v_phi)
    print(f"\nwrote models/{man['file']} ({man['n_varpi']}x{man['n_z']}), "
          f"format {man['format']}")


if __name__ == "__main__":
    main()
