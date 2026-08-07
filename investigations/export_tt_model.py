"""Campaign TT -- export a rotating star with a genuinely mixed field.

Run:  scf/.venv/bin/python3 investigations/export_tt_model.py --scan
      scf/.venv/bin/python3 investigations/export_tt_model.py --bpole 1.5e12

This is export_rotating_model.py with one parameter moved and one claim
dropped. The star, the rotation law and the toroidal field are identical to
the run already on disk; only the poloidal amplitude changes.

WHY, in one paragraph. The configuration evolved so far is toroidal to one
part in 10^7 by energy, which is the textbook Tayler-unstable case, and it
behaved like one: an m=1 mode destroyed the ordered field in a few Alfven
times. Braithwaite-type stable configurations are not toroidal-dominated --
they are twisted tori with comparable poloidal and toroidal energy. That
geometry is the one the equilibrium literature on super-Chandrasekhar white
dwarfs never constructs, so whether it survives is open, and it is an IDEAL
MHD question: stability does not depend on the resistivity, which is the one
thing this problem cannot get right (Rm ~ 6 on the finest grid we can afford,
against ~10^18 physically).

WHAT IS GIVEN UP, stated plainly. The poloidal field here is imposed on a
converged toroidal+rotation equilibrium, not solved self-consistently, so it
leaves the pair out of virial balance by its own energy. papers/wd-toroidal-
poloidal measures the cost: the virial error tracks E_pol/|W| almost one for
one, and the pair leaves the certified band once the exterior dipole passes
~10^11 G. At B_t/B_p = 8.8 the error is 2.3e-2, twenty times the threshold.

So the output of this script is NOT an equilibrium and must not be presented
as one. It is an initial condition to be relaxed, which is what that same
draft prescribes for this range. The accommodation/damping phase in the
inputs file is the relaxation, and the state the star settles into -- not the
state written here -- is what the run is about. Expect a larger initial
transient than the toroidal-dominated run, which already breathed by a factor
2.5 in volume.

--scan reports the trade-off and writes nothing. Pick a point from it, then
run again with --bpole to write the model.
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

# Everything below is copied unchanged from export_rotating_model.py so the
# two configurations differ in exactly one parameter.
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
OUTDIR = REPO / "models"

# The scan. 1.0e9 is the configuration already evolved; the rest raise the
# poloidal amplitude towards a twisted torus. E_pol goes as B_pole^2, so the
# ratio falls as the square.
SCAN = [1.0e9, 1.0e11, 5.0e11, 1.0e12, 1.5e12, 3.0e12]


def build_poloidal(rho, r, th, k0, varpi):
    u = solve_gradshafranov(-4.0 * np.pi * varpi ** 2 * rho * k0, r, th,
                            lmax=LMAX)
    Br, Bth = diag.poloidal_field(u, r, th)
    return u, Br, Bth


def solve_star():
    """The rotating, toroidal equilibrium. Identical to the evolved run."""
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
    return res, r, th, rot


def energies(rho, r, th, Bphi, Br, Bth):
    varpi = r[:, None] * np.sin(th)[None, :]
    w = 2.0 * np.pi * varpi
    dr = np.gradient(r)[:, None]
    dth = np.gradient(th)[None, :]
    dV = w * r[:, None] * dr * dth
    E_tor = float((Bphi ** 2 / (8 * np.pi) * dV).sum())
    E_pol = float(((Br ** 2 + Bth ** 2) / (8 * np.pi) * dV).sum())
    return E_tor, E_pol


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", action="store_true",
                    help="report the trade-off and write nothing")
    ap.add_argument("--bpole", type=float, default=None,
                    help="surface dipole target in G; writes the model")
    a = ap.parse_args()
    if not a.scan and a.bpole is None:
        ap.error("give --scan or --bpole")

    res, r, th, rot = solve_star()
    rho, Phi, H = res["rho"], res["Phi"], res["H"]
    varpi = r[:, None] * np.sin(th)[None, :]
    M = units.g_to_msun(scf_mod.total_mass(rho, r, th))
    W = abs(diag.gravitational_energy(rho, Phi, r, th))
    R_eq, R_pol = diag.equatorial_polar_radii(H, r, th)
    Bphi = ToroidalSC(K=K_TOR, m=M_TOR).B_phi(rho, varpi)
    print(f"star: M = {M:.4f} Msun, R_eq = {R_eq:.4e} cm, "
          f"R_pol/R_eq = {R_pol/R_eq:.4f}, |W| = {W:.4e} erg")

    # B_pole is linear in k0, so one reference solve calibrates all of them.
    u_ref, Br_ref, Bth_ref = build_poloidal(rho, r, th, K0_REF, varpi)
    pole_ref = float(abs(Br_ref[:, 0]).max())
    print(f"calibration: k0 = {K0_REF:.3e} gives B_pole = {pole_ref:.4e} G\n")

    targets = SCAN if a.scan else [a.bpole]
    print(f"{'B_pole (G)':>12s} {'E_tor/E_pol':>12s} {'max|Bpol| (G)':>14s} "
          f"{'E_pol/|W|':>11s} {'virial est.':>12s}")
    chosen = None
    for tgt in targets:
        k0 = K0_REF * tgt / pole_ref
        u, Br, Bth = build_poloidal(rho, r, th, k0, varpi)
        E_tor, E_pol = energies(rho, r, th, Bphi, Br, Bth)
        bpol_max = float(np.hypot(Br, Bth).max())
        # the draft measures the virial error tracking E_pol/|W| one for one
        print(f"{tgt:12.2e} {E_tor/max(E_pol,1e-99):12.3e} {bpol_max:14.4e} "
              f"{E_pol/W:11.3e} {E_pol/W:12.3e}")
        chosen = (tgt, k0, u, Br, Bth, E_tor, E_pol)

    if a.scan:
        print("\nThe last column is the expected virial error: the draft in")
        print("papers/wd-toroidal-poloidal finds it tracks E_pol/|W| almost")
        print("one for one. The certified band is ~1e-3; a twisted torus is")
        print("not inside it, which is the point -- pick the ratio you want")
        print("to evolve, not the one that certifies, and relax it.")
        return

    tgt, k0, u, Br, Bth, E_tor, E_pol = chosen
    print(f"\nwriting the model at B_pole = {tgt:.3e} G, "
          f"E_tor/E_pol = {E_tor/E_pol:.3g}")
    print("NOT an equilibrium -- see the module docstring. Relax it.")

    OUTDIR.mkdir(exist_ok=True)
    out = OUTDIR / "mixed_tt.txt"
    # the same export path the evolved model went through, so the two are
    # byte-for-byte comparable in format
    A_phi = vector_potential(u, r, th)
    verify_meridional_curl(A_phi, Br, Bth, r, th, gate=CURL_GATE)
    mer = to_meridional(rho, Phi, H, A_phi, Bphi, rot, r, th,
                        n=N_MER, half_cm=HALF_CM, corner=CORNER)
    verify_curl_on_cartesian(mer, gate=CURL_GATE, div_gate=DIV_GATE)
    write_model(out, mer, rho_c=RHO_C, mu_e=MU_E)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
