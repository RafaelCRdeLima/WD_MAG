"""The confinement gain along the braking sequence, recomputed per point.

Run:  scf/.venv/bin/python3 investigations/confinement_gain_sequence.py

braking_sequence.py showed that a space-filling B_phi = K rho varpi holding
2 Msun needs a peak field above the Landau critical field B_c across the
whole ignition/collapse window -- 1.98 B_c at the C+C line, 6.73 at 16O. A
field-independent equation of state is then inconsistent with its own star
exactly where the nuclear physics happens.

confinement_cost.py measured, at rho_c = 1e9 only, that a toroidal field
confined to the closed-poloidal-line region reaches the SAME E_tor/|W| at
1/2.71 of the peak field, because peak field at fixed energy is set by how
peaked the profile is, not only by how much volume it has. Carried as a
constant that put the C+C point at 0.73 B_c -- below threshold.

That constant is the weak link: the closed region is bounded by the last
field line that stays inside the star, and the star's shape changes along
the sequence. So the gain is recomputed here at every rho_c, from the
converged configuration at that point.

This decides which equation of state the collaboration needs. If the
confined ignition point stays below B_c, the calculation lives in a regime
where a field-independent EOS is defensible and the priority is finite
temperature plus Coulomb (Skye or PC). If it does not, a magnetised EOS with
Landau quantisation and P_perp != P_par is required instead.

Scoping, as in confinement_cost.py: the confined field is imposed on the
density that the self-consistent space-filling branch produced, and the
question asked is what peak field the same energy costs. Whether that
configuration is itself an equilibrium is the non-barotropic problem.

Writes confinement_gain_sequence.csv.
"""

import sys
import time
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

MU_E, M_TOR = 2.0, 1.0
NR, TOL, LMAX = 65, 1.0e-6, 16
K0 = 1.0e-13        # the closed-region geometry is independent of k0 (the
                    # Grad-Shafranov source is linear in it), so any value
                    # gives the same gain; this one matches earlier work
ZETA = 1.0
B_C = 4.414e13
RHO_CC, RHO_NE, RHO_O = 3.0e9, 9.6e9, 1.94e10

SRC = HERE / "braking_sequence.csv"
OUT = HERE / "confinement_gain_sequence.csv"


def surface_flux(u, rho, th):
    vals = []
    for j in range(len(th)):
        inside = rho[:, j] > 0.0
        if inside.any():
            vals.append(u[np.flatnonzero(inside)[-1], j])
    return max(vals) if vals else 0.0


def main():
    seq = []
    with SRC.open() as f:
        for line in f:
            if line.startswith("#") or line.startswith("rho_c"):
                continue
            p = line.strip().split(",")
            seq.append((float(p[0]), float(p[1])))       # rho_c, K_tor

    print(f"recomputing the confinement gain at {len(seq)} points, "
          f"Nr = {NR}\n")
    print("  rho_c       closed vol   B_t space-fill  /B_c   B_t confined  "
          "/B_c   gain   B_pole (G)")

    rows = []
    for rho_c, K in seq:
        res, r, th, _ = _solve_toroidal_certified(
            rho_c=rho_c, R_guess=r_guess(rho_c), K_tor=K, m_tor_sc=M_TOR,
            rotation=None, mu_e=MU_E, Nr_base=NR, Ntheta=NR, lmax=LMAX,
            tol=TOL, max_iter=200)
        if res is None:
            print(f"  {rho_c:.3e}   no convergence")
            continue
        rho, H = res["rho"], res["H"]
        M = units.g_to_msun(scf_mod.total_mass(rho, r, th))
        varpi = r[:, None] * np.sin(th)[None, :]
        inside = rho > 0.0

        Bphi_sc = ToroidalSC(K=K, m=M_TOR).B_phi(rho, varpi)
        _, E_tor_sc, _ = diag.magnetic_energies(
            np.zeros_like(rho), np.zeros_like(rho), Bphi_sc, r, th)
        Bt_sc = np.abs(Bphi_sc).max()

        u = solve_gradshafranov(-4.0 * np.pi * varpi**2 * rho * K0, r, th,
                                lmax=LMAX)
        u_s = surface_flux(u, rho, th)
        closed = inside & (u > u_s)
        f_vol = (diag.volume_integral(closed.astype(float), r, th)
                 / diag.volume_integral(inside.astype(float), r, th))
        u_norm = max(u.max() - u_s, 1e-300)
        w = np.where(closed, (u - u_s) / u_norm, 0.0)
        shape = np.where(varpi > 0,
                         np.power(w, ZETA) / np.maximum(varpi, 1e-30), 0.0)
        _, E_unit, _ = diag.magnetic_energies(
            np.zeros_like(rho), np.zeros_like(rho), shape, r, th)
        if not np.isfinite(E_unit) or E_unit <= 0:
            print(f"  {rho_c:.3e}   no closed region")
            continue
        beta_0 = np.sqrt(E_tor_sc / E_unit)
        Bt_conf = beta_0 * np.abs(shape).max()
        gain = Bt_sc / Bt_conf

        Br, Bth = diag.poloidal_field(u, r, th)
        B_pole = diag.surface_dipolarity(np.hypot(Br, Bth), H, r, th)["B_pole"]

        rows.append((rho_c, M, f_vol, Bt_sc, Bt_sc / B_C, Bt_conf,
                     Bt_conf / B_C, gain, B_pole))
        print(f"  {rho_c:.3e}   {100*f_vol:7.2f}%   {Bt_sc:.3e}   "
              f"{Bt_sc/B_C:5.2f}   {Bt_conf:.3e}   {Bt_conf/B_C:5.2f}  "
              f"{gain:5.2f}   {B_pole:.3e}")

    a = np.array([[r[0], r[6], r[7]] for r in rows])          # rho, /Bc, gain
    print("\n  interpolated at the three lines:")
    for x, name in ((RHO_CC, "C+C  "), (RHO_NE, "20Ne "), (RHO_O, "16O  ")):
        bc = np.interp(np.log10(x), np.log10(a[:, 0]), a[:, 1])
        gn = np.interp(np.log10(x), np.log10(a[:, 0]), a[:, 2])
        flag = "BELOW B_c" if bc < 1.0 else "above B_c"
        print(f"    {name} ({x:.2e}): confined {bc:.2f} B_c, "
              f"gain {gn:.2f}   {flag}")

    with OUT.open("w") as f:
        f.write(f"# confinement gain along the 2 Msun braking sequence, "
                f"mu_e={MU_E}, Nr={NR}, zeta={ZETA}, k0={K0:.3e}\n")
        f.write(f"# B_c={B_C:.4e} G; lines C+C={RHO_CC:.3e}, "
                f"20Ne={RHO_NE:.3e}, 16O={RHO_O:.3e}\n")
        f.write("rho_c,M_msun,closed_volume_fraction,Bt_spacefilling_G,"
                "Bt_sf_over_Bc,Bt_confined_G,Bt_conf_over_Bc,gain,B_pole_G\n")
        for row in rows:
            f.write(",".join(f"{v:.6e}" for v in row) + "\n")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"elapsed {time.time() - t0:.0f} s")
