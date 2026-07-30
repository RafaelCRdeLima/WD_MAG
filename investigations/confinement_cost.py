"""Can a CONFINED toroidal field still support 2 Msun? The Landau limit says.

Run:  scf/.venv/bin/python3 investigations/confinement_cost.py

Two questions decide whether the non-barotropic route is worth building:
does it give B_tor > B_pol, and does it keep the mass high? The first is
about the source, which non-barotropy frees. The second is about the
confinement, which it does NOT free -- varpi B_phi = beta(u) and beta must
vanish on lines reaching the vacuum, so the toroidal field lives only where
u > u_s, whatever the equation of state.

That is the whole difficulty, and it is quantitative. The certified 2 Msun
configuration carries E_tor/|W| = 0.203 with B_phi = K rho varpi, which is
non-zero EVERYWHERE the density is. A mixed equilibrium must fit the same
energy into the closed-line region alone. Less volume at the same energy
means a larger peak field -- and that configuration already peaks at
4.28e13 G against a Landau critical field B_c = 4.414e13 G, so there is
almost no headroom before the field-independent equation of state stops
describing its own star.

A second squeeze was expected and does NOT happen, which is worth
recording. The guess was that a stronger exterior dipole means more flux
escaping, a larger u_s, a smaller closed region, and so a direct
competition between the dipole the paper wants and the toroidal field the
mass needs. Measured, the closed region is 37.2% of the stellar volume at
every poloidal amplitude tried, across three decades. With f(u) = k0
constant the Grad-Shafranov source is linear in k0, so u and u_s scale
together and the geometry is invariant. The competition is real but it
runs through the Landau limit, not through the volume.

This measures, for a range of poloidal amplitudes: the exterior dipole,
the closed-region volume, the peak toroidal field a confined beta(u) needs
to reproduce E_tor/|W| = 0.203, and the peak poloidal field, both against
B_c.

Writes confinement_cost.csv.
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
from terms.toroidal_sc import ToroidalSC         # noqa: E402

RHO_C, MU_E = 1.0e9, 2.0
K_TOR, M_TOR = 3.245e-3, 1.0        # the M = 2 Msun crossing
K0_LIST = (1.0e-13, 3.0e-13, 1.0e-12, 3.0e-12, 1.0e-11)
ZETA = 1.0                          # the most favourable, from the ceiling run
B_C = 4.414e13                      # Landau critical field, m_e^2 c^3 / e hbar
LMAX = 16
OUT = HERE / "confinement_cost.csv"


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
    Bphi_sc = ToroidalSC(K=K_TOR, m=M_TOR).B_phi(rho, varpi)
    _, E_tor_sc, _ = diag.magnetic_energies(np.zeros_like(rho),
                                            np.zeros_like(rho), Bphi_sc, r, th)
    target = E_tor_sc
    inside = H > 0.0
    vol_star = diag.volume_integral(inside.astype(float), r, th)

    print(f"reference: M = {M:.4f} Msun, |W| = {W:.4e} erg")
    print(f"  self-consistent toroidal branch, B_phi = K rho varpi:")
    print(f"    E_tor/|W| = {E_tor_sc / W:.4f}, max|B_phi| = "
          f"{np.abs(Bphi_sc).max():.4e} G = {np.abs(Bphi_sc).max() / B_C:.3f} B_c")
    print(f"    it fills the whole star, not a closed region\n")
    print("A confined beta(u) must put that same energy in less volume:\n")
    print("   k0        B_pole (G)  closed  max|B_t| need  /B_c   max|B_p|   "
          "/B_c   B_t/B_p  E_t/E_p")

    rows = []
    for k0 in K0_LIST:
        u = solve_gradshafranov(-4.0 * np.pi * varpi**2 * rho * k0, r, th,
                                lmax=LMAX)
        u_s = surface_flux(u, H, th)
        Br, Bth = diag.poloidal_field(u, r, th)
        dip = diag.surface_dipolarity(np.hypot(Br, Bth), H, r, th)
        closed = inside & (u > u_s)
        f_vol = diag.volume_integral(closed.astype(float), r, th) / vol_star

        # beta_0 = 1 shape, then scale: E_tor is quadratic in beta_0
        u_norm = max(u.max() - u_s, 1e-300)
        w = np.where(closed, (u - u_s) / u_norm, 0.0)
        shape = np.where(varpi > 0, np.power(w, ZETA) / np.maximum(varpi, 1e-30),
                         0.0)
        _, E_unit, _ = diag.magnetic_energies(np.zeros_like(rho),
                                              np.zeros_like(rho), shape, r, th)
        if E_unit <= 0 or not np.isfinite(E_unit):
            print(f"   {k0:.0e}   -- no closed region --")
            continue
        beta_0 = np.sqrt(target / E_unit)
        Bphi_max = beta_0 * np.abs(shape).max()
        ratio = Bphi_max / B_C
        Bpol_max = np.hypot(Br, Bth).max()
        E_pol, _, _ = diag.magnetic_energies(Br, Bth, np.zeros_like(rho),
                                             r, th)
        # both peaks must stay under B_c for the field-free EOS to hold
        verdict = ("ok" if max(ratio, Bpol_max / B_C) < 1.0
                   else "ABOVE B_c")
        rows.append((k0, dip["B_pole"], f_vol, Bphi_max, ratio, Bpol_max,
                     Bpol_max / B_C, Bphi_max / Bpol_max, target / E_pol,
                     verdict))
        print(f"   {k0:.0e}   {dip['B_pole']:.3e}  {100*f_vol:5.1f}%  "
              f"{Bphi_max:.3e}   {ratio:5.2f}  {Bpol_max:.3e}  "
              f"{Bpol_max / B_C:6.2f}  {Bphi_max / Bpol_max:7.2f}  "
              f"{target / E_pol:7.2f}  {verdict}")

    print("\n  E_tor/|W| = 0.203 held fixed in every row: this asks what peak")
    print("  field the confinement costs, not whether the mass is reachable.")

    with OUT.open("w") as f:
        f.write(f"# 2 Msun star, rho_c={RHO_C:.3e}, K_tor={K_TOR:.6g}, "
                f"M={M:.4f} Msun, |W|={W:.6e} erg\n")
        f.write(f"# target E_tor = {target:.6e} erg = {target / W:.4f} |W|, "
                f"zeta={ZETA}, B_c={B_C:.4e} G\n")
        f.write("k0,B_pole_G,closed_volume_fraction,Bphi_max_needed_G,"
                "Bphi_over_Bc,Bpol_max_G,Bpol_over_Bc,Bt_over_Bp,"
                "Etor_over_Epol,verdict\n")
        for k0, bp, fv, bm, ra, bpm, bpr, amp, en, vd in rows:
            f.write(f"{k0:.6e},{bp:.6e},{fv:.6e},{bm:.6e},{ra:.6e},"
                    f"{bpm:.6e},{bpr:.6e},{amp:.6e},{en:.6e},{vd}\n")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
