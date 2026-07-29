"""Build the 1D hydrostatic polytrope both codes will be given.

Run:  scf/.venv/bin/python3 flash_crosscheck/make_polytrope.py

Why a polytrope and not the real ztwd star: FLASH ships Gamma,
Helmholtz, Multigamma and Tabulated EOS units -- no ztwd. Handing Castro
a ztwd star and FLASH a Helmholtz star would mean any difference in
rho_c(t) could be the equation of state rather than the scheme, which is
exactly the confound that produced the Step 3 collapse (docs/teoria.md
Sec 6.4). A gamma-law polytrope is available in BOTH codes (Castro's
Microphysics gamma_law, FLASH's Eos/Gamma), so the EOS drops out of the
comparison and what is left is the question we actually want answered:
does an Eulerian Godunov scheme hold a self-gravitating hydrostatic star,
or does rho_c drift at the percent level per dynamical time?

n = 3/2 (gamma = 5/3), not n = 3. An n = 3 polytrope is marginally
stable, so a drifting rho_c could be genuine near-neutral physics rather
than a numerical artifact -- the one confound this test cannot afford.
n = 3/2 is unambiguously stable, so any drift measured is the scheme.

The star is scaled to the SAME central density and radius as the
wd_braithwaite background star, so its dynamical time matches and the
curves are directly comparable to the Castro reference in the paper's
Fig. 1 (as context, not as the primary test -- that one is Castro
polytrope vs FLASH polytrope).

Writes model_polytrope.dat in the single-header-line format Castro's
model_parser reads and scf/castro_model_writer.py already emits
("# r density temperature pressure X"), so the Castro side needs no new
reader; the FLASH side reads the same file.
"""

from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent

G = 6.67428e-8          # cm^3 g^-1 s^-2, the value Castro's own constants use
K_BOLTZ = 1.380658e-16  # erg/K
M_U = 1.660539e-24      # g

# Matched to the wd_braithwaite background star (docs/teoria.md Sec 6.9).
RHO_C = 9.883938495e8   # g/cm^3
R_STAR = 2.36e8         # cm

N_POLY = 1.5            # gamma = 1 + 1/n = 5/3
GAMMA = 1.0 + 1.0 / N_POLY

# Composition: mu2.net's single species, A = 4, Z = 2 (He4-like), so the
# gamma-law temperature is defined with mu = abar/(1 + zbar) for full
# ionization. Both codes are given the same A and Z.
ABAR, ZBAR = 4.0, 2.0
MU = ABAR / (1.0 + ZBAR)

NPTS = 4096
DOMAIN_HALF = 4.90e8    # cm, the Castro run's domain half-width
RHO_FLOOR = 1.0e-4      # g/cm^3, ambient outside the star


def lane_emden(n):
    """theta'' + (2/xi) theta' + theta^n = 0, to the first zero."""
    def rhs(xi, y):
        theta, dtheta = y
        # series start avoids the 1/xi singularity at the origin
        lap = -max(theta, 0.0) ** n - (2.0 / xi) * dtheta if xi > 0 else 0.0
        return [dtheta, lap]

    def hit_surface(xi, y):
        return y[0]
    hit_surface.terminal = True
    hit_surface.direction = -1

    xi0 = 1e-8
    y0 = [1.0 - xi0**2 / 6.0, -xi0 / 3.0]
    sol = solve_ivp(rhs, (xi0, 20.0), y0, events=hit_surface,
                    rtol=1e-12, atol=1e-14, dense_output=True, max_step=0.01)
    xi1 = float(sol.t_events[0][0])
    dtheta1 = float(sol.y_events[0][0][1])
    return sol, xi1, dtheta1


def main():
    sol, xi1, dtheta1 = lane_emden(N_POLY)
    print(f"Lane-Emden n={N_POLY}: xi1 = {xi1:.6f}, "
          f"-xi1^2 theta'(xi1) = {-xi1**2 * dtheta1:.6f}")

    a = R_STAR / xi1                       # radial scale length
    # a^2 = (n+1) K rho_c^(1/n - 1) / (4 pi G)
    K = 4.0 * np.pi * G * a**2 * RHO_C ** (1.0 - 1.0 / N_POLY) / (N_POLY + 1.0)
    M = 4.0 * np.pi * a**3 * RHO_C * (-xi1**2 * dtheta1)
    t_dyn = np.sqrt(R_STAR**3 / (G * M))   # same formula as seed_field.H:85
    print(f"K = {K:.6e} (cgs, P = K rho^{GAMMA:.4f})")
    print(f"M = {M:.6e} g = {M / 1.98892e33:.4f} Msun")
    print(f"R = {R_STAR:.4e} cm ; t_dyn = {t_dyn:.6f} s")
    print(f"mean rho = {M / (4 / 3 * np.pi * R_STAR**3):.4e} g/cm^3 ; "
          f"rho_c/rho_mean = {RHO_C / (M / (4 / 3 * np.pi * R_STAR**3)):.2f}")

    # Sample onto a uniform radial grid covering the whole Castro domain,
    # so neither code has to extrapolate at its outer boundary.
    r = np.linspace(0.0, DOMAIN_HALF * np.sqrt(3.0), NPTS)
    rho = np.full_like(r, RHO_FLOOR)
    interior = r < R_STAR
    xi = r[interior] / a
    theta = np.clip(sol.sol(np.maximum(xi, 1e-8))[0], 0.0, None)
    rho[interior] = np.maximum(RHO_C * theta ** N_POLY, RHO_FLOOR)

    P = K * rho ** GAMMA
    T = P * MU * M_U / (rho * K_BOLTZ)     # gamma-law, consistent with P, rho

    hdr = "# r density temperature pressure X"
    out = HERE / "model_polytrope.dat"
    with out.open("w") as f:
        f.write(hdr + "\n")
        for ri, di, ti, pi in zip(r, rho, T, P):
            f.write(f"{ri:.10e} {di:.10e} {ti:.10e} {pi:.10e} 1.0000000000e+00\n")
    print(f"\nwrote {out.name}: {NPTS} points to r = {r[-1]:.4e} cm")
    print(f"  T(centre) = {T[0]:.4e} K, T(surface) = {T[interior][-1]:.4e} K")
    print(f"  P(centre) = {P[0]:.4e} erg/cm^3")

    # Both codes need these as runtime parameters; print them so the two
    # inputs files cannot silently disagree with the model.
    print("\n--- put these in BOTH inputs files ---")
    print(f"  gamma        = {GAMMA:.10f}")
    print(f"  A (abar)     = {ABAR}")
    print(f"  Z (zbar)     = {ZBAR}")
    print(f"  rho_ambient  = {RHO_FLOOR:.4e}")
    print(f"  t_dyn        = {t_dyn:.10f} s")
    print(f"  stop_time    = {4.0 * t_dyn:.6f} s   (4 t_dyn)")


if __name__ == "__main__":
    main()
