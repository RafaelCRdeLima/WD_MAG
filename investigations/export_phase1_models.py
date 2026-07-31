"""Phase 1 export: the models whose stability has to be tested.

Run:  scf/.venv/bin/python3 investigations/export_phase1_models.py

Phase 1 of docs/PLAN_ULTRAMASSIVE_SNIA.md is a go/no-go: if the configuration
that supports 2 Msun does not survive an Alfven time, everything downstream is
optimising something that does not exist. Its exit criterion asks for a
measured survival time at two values of B_t/B_p -- one of them the value the
scenario's 1e9 G surface dipole implies.

That value is the worry. E_pol scales as k0^2 at fixed E_tor, and the measured
pair is E_t/E_p = 884 at B_pole = 5.9e10 G, so a 1e9 G dipole implies
E_t/E_p ~ 3e6: the configuration most exposed to the Tayler m = 1 instability.
The second model is the balanced one confinement makes available, roughly
B_t/B_p = 1 in amplitude, as the control. If the first is unstable and the
second is not, that brackets the answer and tells the collaboration what the
dipole constraint actually costs.

Both models are written as vector potential on a meridional grid, so the
consuming problem takes a discrete curl and gets div B = 0 as an identity.
Nothing is written unless the two verifications pass: that the curl of the
written A returns the field it was built from, and that its discrete
divergence on a Cartesian mesh is zero to machine precision.

Writes models/phase1_*.txt (gitignored -- large derived data) and their
manifests (committed -- small, and they carry the check results).
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

import diagnostics as diag                            # noqa: E402
import scf as scf_mod                                 # noqa: E402
import units                                          # noqa: E402
from axisym_model_writer import (to_meridional, vector_potential,  # noqa: E402
                                 verify_curl_on_cartesian,
                                 verify_meridional_curl, write_model)
from gradshafranov import solve_gradshafranov         # noqa: E402
from seed import r_guess                              # noqa: E402
from sweep_worker import _solve_toroidal_certified    # noqa: E402
from terms.toroidal_sc import ToroidalSC              # noqa: E402

RHO_C, MU_E = 1.0e9, 2.0
K_TOR, M_TOR = 3.245e-3, 1.0        # the certified M = 2 Msun crossing
LMAX = 16
N_MER = 385                          # meridional resolution of the model file
HALF_CM = 9.0e8                      # simulation box half-width (Castro runs)
CORNER = 1.7320508                   # sqrt(3): the box DIAGONAL, not its face
K0_REF = 1.0e-13                     # reference poloidal amplitude
B_POLE_TARGET = 1.0e9                # the scenario's surface dipole
K0_BALANCED = 1.0e-12                # measured: B_t/B_p ~ 1 in amplitude

DIV_GATE = 1.0e-12                   # normalised |div B|
CURL_GATE = 5.0e-2                   # relative error of curl A against B

OUTDIR = REPO / "models"


def build(rho, r, th, k0, varpi):
    u = solve_gradshafranov(-4.0 * np.pi * varpi ** 2 * rho * k0, r, th,
                            lmax=LMAX)
    Bphi = ToroidalSC(K=K_TOR, m=M_TOR).B_phi(rho, varpi)
    Br, Bth = diag.poloidal_field(u, r, th)
    return u, Bphi, Br, Bth


def main():
    OUTDIR.mkdir(exist_ok=True)
    res, r, th, ov = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=K_TOR, m_tor_sc=M_TOR,
        rotation=None, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=LMAX,
        tol=1e-8, max_iter=200)
    if res is None:
        raise SystemExit("the 2 Msun solve did not converge")
    rho, Phi, H = res["rho"], res["Phi"], res["H"]
    M = units.g_to_msun(scf_mod.total_mass(rho, r, th))
    W = abs(diag.gravitational_energy(rho, Phi, r, th))
    varpi = r[:, None] * np.sin(th)[None, :]
    R_eq, R_pol = diag.equatorial_polar_radii(H, r, th)
    print(f"background: M = {M:.4f} Msun, R_eq = {R_eq:.4e}, "
          f"R_pol = {R_pol:.4e} cm, frac_pol = {ov['frac_pol']:.3f}")

    # calibrate k0 for the scenario dipole: B_pole is linear in k0
    _, _, Br0, Bth0 = build(rho, r, th, K0_REF, varpi)
    bp_ref = diag.surface_dipolarity(np.hypot(Br0, Bth0), H, r, th)["B_pole"]
    k0_scenario = K0_REF * (B_POLE_TARGET / bp_ref)
    print(f"calibration: B_pole = {bp_ref:.4e} G at k0 = {K0_REF:.3e}"
          f"  ->  k0 = {k0_scenario:.4e} for {B_POLE_TARGET:.1e} G\n")

    cases = [("scenario", k0_scenario), ("balanced", K0_BALANCED)]
    # The meridional grid must cover the CORNERS of the Cartesian box, not
    # its faces: a point at (HALF, HALF, HALF) sits at radius sqrt(3)*HALF.
    # Sizing it to the face instead leaves A undefined in the corners, where
    # it falls to the interpolator fill value and the discrete curl of that
    # step produces a field larger than the star's -- which shows up as an
    # "amplitude retained" above 100%, the check that caught it.
    rmax = 1.02 * CORNER * HALF_CM
    vp = np.linspace(0.0, rmax, N_MER)
    zz = np.linspace(-rmax, rmax, 2 * N_MER - 1)
    print(f"model grid to {rmax:.4e} cm covers the corners of a "
          f"+-{HALF_CM:.2e} cm box\n")

    for name, k0 in cases:
        u, Bphi, Br, Bth = build(rho, r, th, k0, varpi)
        E_pol, E_tor, E_mag = diag.magnetic_energies(Br, Bth, Bphi, r, th)
        bp = diag.surface_dipolarity(np.hypot(Br, Bth), H, r, th)["B_pole"]
        amp = np.abs(Bphi).max() / max(np.hypot(Br, Bth).max(), 1e-300)

        rho_m, u_m, bphi_m = to_meridional(r, th, (rho, u, Bphi), vp, zz)
        A_phi, A_z = vector_potential(vp, u_m, bphi_m)

        err_pol, err_tor = verify_meridional_curl(vp, zz, A_phi, A_z, u_m,
                                                  bphi_m)
        rel_div, b_max, dx = verify_curl_on_cartesian(
            vp, zz, A_phi, A_z, half=HALF_CM, n_cart=64)

        print(f"[{name}] k0 = {k0:.4e}")
        print(f"   B_pole = {bp:.4e} G, E_tor/E_pol = {E_tor / E_pol:.4g}, "
              f"B_t/B_p = {amp:.3g} (amplitude)")
        print(f"   E_tor/|W| = {E_tor / W:.4f}, max|B_phi| = "
              f"{np.abs(Bphi).max():.4e} G")
        print(f"   curl A vs B: poloidal {err_pol:.3e}, toroidal "
              f"{err_tor:.3e}   (gate {CURL_GATE:.0e})")
        print(f"   div B on 64^3: {rel_div:.3e}   (gate {DIV_GATE:.0e}), "
              f"amplitude retained {100 * b_max / np.abs(Bphi).max():.1f}%")

        retained = b_max / np.abs(Bphi).max()
        if retained > 1.02:
            raise SystemExit(
                f"[{name}] reconstruction gained amplitude ({100*retained:.1f}%)"
                " -- the model grid does not cover the Cartesian box")
        if not (rel_div < DIV_GATE):
            raise SystemExit(f"[{name}] divergence gate failed")
        if not (err_pol < CURL_GATE and err_tor < CURL_GATE):
            raise SystemExit(f"[{name}] curl gate failed")

        params = dict(rho_c=RHO_C, mu_e=MU_E, K_tor=K_TOR, m_tor=M_TOR,
                      k0=k0, M_msun=M, R_eq_cm=R_eq, R_pol_cm=R_pol,
                      B_pole_G=bp, E_tor_over_Epol=E_tor / E_pol,
                      E_tor_over_W=E_tor / W, Bt_over_Bp_amplitude=amp,
                      Bphi_max_G=float(np.abs(Bphi).max()))
        checks = dict(curl_err_poloidal=err_pol, curl_err_toroidal=err_tor,
                      relative_divB_64cubed=rel_div,
                      amplitude_retained_64cubed=b_max / np.abs(Bphi).max())
        man = write_model(vp, zz, rho_m, A_phi, A_z,
                          OUTDIR / f"phase1_{name}.txt", params, checks)
        print(f"   wrote models/{man['file']} "
              f"({man['n_varpi']}x{man['n_z']})\n")


if __name__ == "__main__":
    main()
