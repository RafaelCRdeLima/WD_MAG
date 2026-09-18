"""The gate confinement_cost.py never opened: does the comparable-split row hold?

confinement_cost.py reports peak fields against B_c and calls a row "ok" when
both stay under it. B_c is the Landau limit -- the field at which the
field-free EOS stops describing its own star. It says nothing about whether
the star can hold the field, which is beta = P_gas / P_mag.

The k0 = 1e-12 row of confinement_cost.csv has B_t/B_p = 0.97 and
E_tor/E_pol = 8.8: a genuine twisted torus in a star whose 2.007 Msun the
toroidal field supports (E_tor/|W| = 0.203, against 1.346 Msun field-free).
DIARIO 14 flagged it as a candidate, not a counterexample, precisely because
beta was never measured on it.

This measures it, with the same definition export_ct_model.py uses so the
numbers compare to the 6.6 table: beta_min over rho > 1e6 g/cm^3.

Two separate questions, and they have different answers:
  (A) beta_min inside the star   -- can this configuration exist at all?
  (B) beta in the Castro ambient -- can we evolve it in the box we have?

(B) is a setup question, not a stellar one: the ambient is an artefact of
filling the domain, and DIARIO 6.8 already named the escapes (vacuum boundary,
or a denser ambient in a smaller box).

Run:  scf/.venv/bin/python3 investigations/twisted_torus_beta.py
Writes twisted_torus_beta.csv.
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
import eos                                       # noqa: E402
import scf as scf_mod                            # noqa: E402
import units                                     # noqa: E402
from gradshafranov import solve_gradshafranov    # noqa: E402
from seed import r_guess                         # noqa: E402
from sweep_worker import _solve_toroidal_certified   # noqa: E402
from terms.toroidal_sc import ToroidalSC         # noqa: E402

RHO_C, MU_E = 1.0e9, 2.0
K_TOR, M_TOR = 3.245e-3, 1.0        # the M = 2 Msun crossing
K0_LIST = (1.0e-13, 3.0e-13, 1.0e-12, 3.0e-12, 1.0e-11)
ZETA = 1.0
B_C = 4.414e13
LMAX = 16
RHO_STAR = 1.0e6                    # same cut as export_ct_model.report()
RHO_AMBIENT = 2.0e4                 # what the Castro box is filled with
OUT = HERE / "twisted_torus_beta.csv"


def surface_flux(u, H, th):
    vals = []
    for j in range(len(th)):
        inside = H[:, j] > 0.0
        if inside.any():
            vals.append(u[np.flatnonzero(inside)[-1], j])
    return max(vals) if vals else 0.0


def main():
    res, r, th, _ = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=K_TOR, m_tor_sc=M_TOR,
        rotation=None, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=LMAX,
        tol=1e-8, max_iter=200)
    if res is None:
        raise SystemExit("the 2 Msun toroidal solve did not converge")
    rho, Phi, H = res["rho"], res["Phi"], res["H"]
    M = units.g_to_msun(scf_mod.total_mass(rho, r, th))
    W = abs(diag.gravitational_energy(rho, Phi, r, th))

    varpi = r[:, None] * np.sin(th)[None, :]
    P_gas = eos.pressure(eos.x_of_enthalpy(np.maximum(H, 0.0), mu_e=MU_E))
    star = rho > RHO_STAR
    inside = H > 0.0
    vol_star = diag.volume_integral(inside.astype(float), r, th)

    Bphi_sc = ToroidalSC(K=K_TOR, m=M_TOR).B_phi(rho, varpi)
    _, target, _ = diag.magnetic_energies(np.zeros_like(rho),
                                          np.zeros_like(rho), Bphi_sc, r, th)

    # the ambient the Castro box is filled with holds this much field at beta=1
    P_amb = eos.pressure(eos.x_of_enthalpy(
        eos.enthalpy((RHO_AMBIENT / eos.B_of_mu_e(MU_E)) ** (1.0 / 3.0),
                     mu_e=MU_E), mu_e=MU_E))
    B_amb_max = np.sqrt(8.0 * np.pi * P_amb)

    print(f"reference: M = {M:.4f} Msun, |W| = {W:.4e} erg, "
          f"E_tor/|W| = {target / W:.4f}")
    print(f"ambient rho = {RHO_AMBIENT:.1e} g/cm^3 holds at most "
          f"{B_amb_max:.3e} G at beta = 1\n")
    print("   k0      B_t/B_p  E_t/E_p   beta_min   rho@bmin   r/R_eq   "
          "f(beta<1)  B_pole      beta_amb   verdict")

    rows = []
    for k0 in K0_LIST:
        u = solve_gradshafranov(-4.0 * np.pi * varpi**2 * rho * k0, r, th,
                                lmax=LMAX)
        u_s = surface_flux(u, H, th)
        Br, Bth = diag.poloidal_field(u, r, th)
        dip = diag.surface_dipolarity(np.hypot(Br, Bth), H, r, th)
        closed = inside & (u > u_s)

        u_norm = max(u.max() - u_s, 1e-300)
        w = np.where(closed, (u - u_s) / u_norm, 0.0)
        shape = np.where(varpi > 0, np.power(w, ZETA) / np.maximum(varpi, 1e-30),
                         0.0)
        _, E_unit, _ = diag.magnetic_energies(np.zeros_like(rho),
                                              np.zeros_like(rho), shape, r, th)
        if E_unit <= 0 or not np.isfinite(E_unit):
            continue
        Bphi = np.sqrt(target / E_unit) * shape

        P_mag = (Br**2 + Bth**2 + Bphi**2) / (8.0 * np.pi)
        beta = P_gas / np.maximum(P_mag, 1e-300)
        bmin = float(beta[star].min())
        k = np.unravel_index(np.argmin(np.where(star, beta, np.inf)),
                             beta.shape)
        R_eq = diag.equatorial_polar_radii(H, r, th)[0]

        # mass fraction of the star that sits below beta = 1
        dm = rho * (beta < 1.0) * star
        f_lowbeta = (diag.volume_integral(dm, r, th)
                     / diag.volume_integral(rho * star, r, th))

        Bpol_max = float(np.hypot(Br, Bth).max())
        Bphi_max = float(np.abs(Bphi).max())
        beta_amb = (B_amb_max / max(dip["B_pole"], 1e-300)) ** 2

        verdict = "HOLDS" if bmin >= 1.0 else "beta < 1"
        rows.append((k0, Bphi_max / Bpol_max, target / diag.magnetic_energies(
            Br, Bth, np.zeros_like(rho), r, th)[0], bmin, rho[k], r[k[0]] / R_eq,
            f_lowbeta, dip["B_pole"], beta_amb, verdict))
        print(f"  {k0:.0e}  {Bphi_max / Bpol_max:7.2f} "
              f"{rows[-1][2]:8.2f}  {bmin:9.3e}  {rho[k]:9.2e}  "
              f"{r[k[0]] / R_eq:6.3f}  {f_lowbeta:8.4f}  "
              f"{dip['B_pole']:.3e}  {beta_amb:9.2e}  {verdict}")

    with OUT.open("w") as f:
        f.write(f"# 2 Msun toroidal-supported star, rho_c={RHO_C:.3e}, "
                f"M={M:.4f} Msun, |W|={W:.6e} erg\n")
        f.write(f"# beta_min over rho > {RHO_STAR:.0e} (same cut as "
                f"export_ct_model.report); ambient {RHO_AMBIENT:.0e} holds "
                f"{B_amb_max:.4e} G at beta=1\n")
        f.write("k0,Bt_over_Bp,Etor_over_Epol,beta_min,rho_at_beta_min,"
                "r_over_Req_at_beta_min,mass_frac_beta_below_1,B_pole_G,"
                "beta_ambient,verdict\n")
        for row in rows:
            f.write(",".join(f"{v:.6e}" if isinstance(v, float) else str(v)
                             for v in row) + "\n")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
