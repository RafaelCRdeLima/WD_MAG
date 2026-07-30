"""A 2 Msun star with a mixed field: what the SCF route can and cannot give.

Run:  scf/.venv/bin/python3 investigations/mixed_2msun.py

Why this configuration and not another. The mass cannot come from the
poloidal branch: measured, that branch stops at E_pol/|W| ~ 0.018 (+3.3%
in mass) with VE already at 8e-4, and neither the Grad-Shafranov ratio
recursion nor continuation in k0 moves it. It cannot come from a
self-consistent mixed equilibrium either: in barotropy beta = beta(u)
confines the toroidal field to the closed-line region and the achievable
ratio is a few percent (docs/teoria.md Sec 1.9, note G4). So the mass has
to come from the self-consistent toroidal branch, which reaches
E_tor/|W| = 0.20 and 2 Msun, and the poloidal component can only be added
on top of it.

That addition is not self-consistent -- hachisu_scf() refuses poloidal and
toroidal together (D6) -- so the poloidal field is solved on the converged
toroidal density and imposed. This script measures what that costs, which
is the number the route turns on: the virial error of the combined
configuration comes out equal to the imposed E_pol/|W| to within a few
percent, because the imposed field's energy IS the unbalanced term in the
virial. There is no free lunch in it, and the size of the bill is what
decides whether "impose and relax" is a method or a wish.

Writes mixed_2msun.csv.
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

RHO_C = 1.0e9          # below the mu_e = 2 neutronization threshold
MU_E = 2.0
K_TOR = 3.245e-3       # the interpolated M = 2 Msun crossing
M_TOR = 1.0
K0_LIST = (1.0e-13, 3.0e-13, 6.0e-13, 1.0e-12)
VE_GATE = 1.0e-3

OUT = HERE / "mixed_2msun.csv"


def main():
    res, r, th, overflow = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=K_TOR, m_tor_sc=M_TOR,
        rotation=None, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=16,
        tol=1e-8, max_iter=200)
    if res is None:
        raise SystemExit("the toroidal solve did not converge")

    rho, Phi, H = res["rho"], res["Phi"], res["H"]
    ve = diag.virial_error_terms(rho, Phi, H, r, th, MU_E, rotation=None,
                                 poloidal=None, toroidal=ToroidalSC(K=K_TOR, m=M_TOR))
    W = abs(ve["W"])
    Bphi = ve["Bphi"]
    M = units.g_to_msun(scf_mod.total_mass(rho, r, th))

    print("toroidal branch, self-consistent:")
    print(f"  M = {M:.4f} Msun   VE = {ve['VE']:.2e}   frac_pol = {overflow['frac_pol']:.3f}")
    print(f"  E_tor/|W| = {ve['E_mag'] / W:.4f}   B_tor,max = {np.abs(Bphi).max():.3e} G")

    omega2 = (r[:, None] * np.sin(th)[None, :]) ** 2
    rows = []
    print("\npoloidal imposed on that density (not self-consistent):")
    print("  k0        E_pol/|W|   E_pol/E_mag   B_pole (G)    VE(total)  certified")
    for k0 in K0_LIST:
        u = solve_gradshafranov(-4.0 * np.pi * omega2 * rho * k0, r, th, lmax=16)
        Br, Bth = diag.poloidal_field(u, r, th)
        E_pol, E_tor, E_mag = diag.magnetic_energies(Br, Bth, Bphi, r, th)
        VE_tot = diag.virial_error(rho, Phi, H, Br, Bth, Bphi, r, th, MU_E)[0]
        dip = diag.surface_dipolarity(np.sqrt(Br**2 + Bth**2), H, r, th)
        ok = VE_tot < VE_GATE
        rows.append((k0, E_pol / W, E_pol / E_mag, dip["B_pole"], VE_tot, ok))
        print(f"  {k0:.0e}   {E_pol / W:9.5f}   {E_pol / E_mag:10.4f}   "
              f"{dip['B_pole']:.3e}   {VE_tot:.2e}   {'yes' if ok else 'NO'}")

    print("\n  VE(total) tracks E_pol/|W| -- the imposed energy is the")
    print("  unbalanced term, so imposing x costs a virial error of about x.")

    with OUT.open("w") as f:
        f.write("# 2 Msun toroidal star (rho_c=%.3e, K=%.6g, m=%g) with an\n"
                % (RHO_C, K_TOR, M_TOR))
        f.write("# imposed poloidal field. M=%.4f Msun, VE_toroidal=%.3e,\n"
                % (M, ve["VE"]))
        f.write("# E_tor/|W|=%.4f, B_tor_max=%.4e G\n" % (ve["E_mag"] / W,
                                                          np.abs(Bphi).max()))
        f.write("k0,E_pol_over_W,E_pol_over_Emag,B_pole_G,VE_total,certified\n")
        for k0, epw, epm, bp, vt, ok in rows:
            f.write(f"{k0:.6e},{epw:.6e},{epm:.6e},{bp:.6e},{vt:.6e},{ok}\n")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
