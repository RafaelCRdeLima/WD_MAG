"""Export the TWISTED TORUS configuration of DIARIO 20 for Castro.

The companion export, export_rotating_model.py, writes the configuration the
equilibrium literature specifies: B_t/B_p ~ 4e3, which is toroidal to within a
part in a thousand and Tayler-unstable. This one writes the configuration the
parameter sweep found instead -- comparable poloidal and toroidal energy, the
Braithwaite & Spruit branch -- so the two can be evolved side by side.

The point of the pair is that only one variable differs in kind: both are
rotating, barotropic, built under the same ztwd the run uses, at the same
central density. What changes is the field geometry.

Which toroidal. The sweep replaced the self-consistent toroidal by a confined
profile of equal energy to compute beta, and a confined toroidal is NOT
exportable: it is a second imposed field on top of the imposed poloidal, and
the TT campaign died at t = 0.06 s on exactly that. So this exports the
SELF-CONSISTENT ToroidalSC field the star was actually solved with, and only
the poloidal is imposed -- the same arrangement as rotating_mixed.txt, which
ran 82 s. Checked beforehand: with the self-consistent toroidal the star holds
the field easily (beta_min = 20.8 with no poloidal), and there is a poloidal
amplitude where beta_min = 1 lands at E_tor/E_pol ~ 7.5, inside the comparable
band. That amplitude is what this calibrates to.

KNOWN BLOCKER, and it is not fixed here. The poloidal that comparable energies
require carries an exterior dipole near 6e11 G, and the 2e4 g/cm^3 ambient
holds 3.4e10 -- beta_ambient ~ 3e-3, within a factor of what killed the ML
campaign at t = 0.221 s. DIARIO 21 has the arithmetic and why neither a denser
ambient nor stronger confinement escapes it. This model is exported to be
PROBED in the box we have, not to survive it: the question the probe asks is
whether the initial transient already differs from the literature's
configuration, which is cheap to answer and decides whether the expensive
setup work is worth doing.

Run:  scf/.venv/bin/python3 investigations/export_twisted_torus_model.py
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

RHO_C, MU_E = 1.0e9, 2.0
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
BETA_TARGET = 1.0        # calibrate the poloidal to beta_min = 1
B_POLE_CONTROL = 1.0e9   # --control: the literature's specified dipole
RHO_STAR = 1.0e6         # same cut as export_ct_model.report()
# 0.11 of the 2 Msun crossing: the K_frac the sweep's best
# four-gate point sits at (DIARIO 20).
K_TOR, M_TOR = 0.11 * 3.245e-3, 1.0

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
    # The CONTROL is the same star with the literature's weak dipole, not
    # rotating_mixed.txt: that model is a different star (rho_c = 3e9, another
    # mass, another K_tor), so comparing against it would mix field geometry
    # with stellar structure. Here only the poloidal amplitude differs, by a
    # factor near 700, and everything else is identical by construction.
    #
    # Read the asymmetry honestly: the control's 1e9 G dipole sits fine in the
    # 2e4 ambient, while the twisted torus at 7.2e11 G does not. If the torus
    # dies early and the control does not, that is the ambient of DIARIO 21
    # and not physics. What the probe actually tests is whether that
    # prediction is right, cheaply, before weeks of boundary work.
    ap = argparse.ArgumentParser()
    ap.add_argument("--control", action="store_true",
                    help="calibrate to the literature dipole instead of beta")
    a = ap.parse_args()
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

    # Calibrate k0 so beta_min = 1: the largest poloidal the star holds.
    # The companion script targets a surface dipole instead, because there the
    # dipole is the observable being specified; here the constraint is that the
    # star can carry the field at all, and the dipole is whatever follows.
    import eos as _eos
    P_gas = _eos.pressure(_eos.x_of_enthalpy(np.maximum(H, 0.0), mu_e=MU_E))
    star = rho > RHO_STAR

    def beta_min_at(k0_try):
        _, Bphi_t, Br_t, Bth_t = build(rho, r, th, k0_try, varpi)
        P_mag = (Br_t**2 + Bth_t**2 + Bphi_t**2) / (8.0 * np.pi)
        return float((P_gas[star] / np.maximum(P_mag[star], 1e-300)).min())

    b_nofield = beta_min_at(0.0)
    print(f"toroidal alone: beta_min = {b_nofield:.3e}")

    if a.control:
        # B_pole is linear in k0, so one solve calibrates it
        _, _, Br0, Bth0 = build(rho, r, th, K0_REF, varpi)
        bp_ref = diag.surface_dipolarity(np.hypot(Br0, Bth0), H, r,
                                         th)["B_pole"]
        k0 = K0_REF * (B_POLE_CONTROL / bp_ref)
        print(f"CONTROL: k0 = {k0:.4e} for B_pole = {B_POLE_CONTROL:.1e} G, "
              f"beta_min = {beta_min_at(k0):.3e}\n")
    else:
        if b_nofield < BETA_TARGET:
            raise SystemExit("the self-consistent toroidal alone already "
                             "fails beta -- no poloidal amplitude rescues it")
        lo, hi = 1.0e-16, 1.0e-9
        for _ in range(60):
            mid = np.sqrt(lo * hi)
            if beta_min_at(mid) >= BETA_TARGET:
                lo = mid
            else:
                hi = mid
        k0 = lo
        print(f"calibration: k0 = {k0:.4e} gives beta_min = "
              f"{beta_min_at(k0):.4f}\n")

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
    # directly on the meridional grid rather than interpolated.
    #
    # It is TAPERED to zero across the density transition, not cut at the
    # stellar surface. Cutting it was the first version, and the first run
    # died of it: a hard cut puts the full 7.5e8 cm/s into one cell beside a
    # static ambient, and Castro aborted at t = 2.59 s with
    # "Invalid density = 871 at index 75, 77, 84" -- a cell at
    # (varpi/R_eq)^2 + (z/R_pol)^2 = 1.02, which is to say exactly on the
    # stellar surface. Mass was conserved to 6e-5 and the central density was
    # unchanged, so the star was healthy; what failed was the shear layer the
    # export had manufactured.
    #
    # The taper is a smoothstep in log10(rho) between the sponge's own
    # bracketing densities, so the rotation dies out over the same two
    # decades the sponge acts on, and both the value and its first derivative
    # vanish at each end.
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
                      OUTDIR / ("twisted_torus_ctl.txt" if a.control
                                else "twisted_torus.txt"),
                      params, checks,
                      v_phi=v_phi)
    print(f"\nwrote models/{man['file']} ({man['n_varpi']}x{man['n_z']}), "
          f"format {man['format']}")


if __name__ == "__main__":
    main()
