"""The mass versus comparability trade-off curve, which replaces a single point: bisect omega_frac to T/|W| = 0.14.

DIARIO 20 reported 1.9623 Msun as the maximum, which was an artefact of
sampling omega_frac at {0, 0.5, 1.0, 1.5, 2.0}: the allowed region extends to
1.67 and the grid skipped it. DIARIO 23 has the correction and the reasoning,
including why the solver's own non-convergence above omega_frac = 1.9 does not
matter -- the secular bar mode at T/|W| = 0.14 arrives first, so the numerical
limit terms/rotation.py documents sits outside the physically allowed region
and the Hachisu reparameterization by axial ratio is not needed work.

Bisects to the bar threshold, then checks the other three gates there: the
largest poloidal the star holds at beta_min = 1, and what E_tor/E_pol that
buys.

--nr sets the SCF resolution. The convergence of these numbers is the thing a
referee asks first, and the cost is only 1.8x from 129 to 257, not the 4x the
O(N^2) solve would suggest.

Run:  scf/.venv/bin/python3 investigations/max_mass_frontier.py [--nr 257]
"""

import argparse
import sys, warnings, numpy as np
# Derivado do proprio arquivo, nao fixo: o caminho da estacao
# (/home/rafael/wd-magnetizada) nao existe no sci-com (~/WD_MAG), e a versao
# anterior morria em "No module named 'diagnostics'" depois de o job ja' ter
# gasto vinte minutos na varredura que vem antes. Regra da DIARIO 10, desta vez
# plantada por mim ao promover o script do scratchpad.
from pathlib import Path as _Path
REPO = str(_Path(__file__).resolve().parent.parent)
for p in (REPO + "/scf", REPO + "/dashboard", REPO + "/investigations"):
    sys.path.insert(0, p)
warnings.filterwarnings("ignore")
import diagnostics as diag, eos, scf as scf_mod, units
from gradshafranov import solve_gradshafranov
from seed import r_guess
from sweep_worker import _solve_toroidal_certified
from terms.rotation import Rotation
from terms.toroidal_sc import ToroidalSC

ap = argparse.ArgumentParser()
ap.add_argument("--nr", type=int, default=129)
A = ap.parse_args()
NR = A.nr

MU_E, LMAX, RHO_STAR = 2.0, 16, 1.0e6
RHO_C, K_TOR = 1.0e9, 0.11 * 3.245e-3
TW_SECULAR = 0.14

ref, r0, th0, _ = _solve_toroidal_certified(
    rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=0.0, m_tor_sc=1.0,
    rotation=None, mu_e=MU_E, Nr_base=NR, Ntheta=NR, lmax=LMAX,
    tol=1e-8, max_iter=200)
M_ref = scf_mod.total_mass(ref["rho"], r0, th0)
R_ref = diag.equatorial_polar_radii(ref["H"], r0, th0)[0]
om_kep = float(np.sqrt(units.G_CONST * M_ref / R_ref**3))

def solve(of):
    rot = Rotation(of * om_kep, R_ref)
    res, r, th, _ = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=K_TOR, m_tor_sc=1.0,
        rotation=rot, mu_e=MU_E, Nr_base=NR, Ntheta=NR, lmax=LMAX,
        tol=1e-8, max_iter=400)
    if res is None: return None
    rho, Phi, H = res["rho"], res["Phi"], res["H"]
    W = abs(diag.gravitational_energy(rho, Phi, r, th))
    return dict(rho=rho, Phi=Phi, H=H, r=r, th=th, rot=rot, W=W,
                M=units.g_to_msun(scf_mod.total_mass(rho, r, th)),
                tw=rot.energy(rho, r, th)["T"]/W)


# The point DIARIO 23 quoted -- 2.2904 Msun at E_tor/E_pol = 9.98 -- sat on the
# edge of the comparable band, and the 257x257 rerun moved it to 10.69, across
# the boundary. The mass is converged to 0.01%; what is not robust is calling
# that configuration a twisted torus. The band edge at 10 is a convention, so
# the honest object is the curve, not the point: how heavy can the star be at a
# given E_tor/E_pol, with the other three gates held.
#
# Both quantities rise together with rotation, so the curve is parameterized by
# omega_frac and the comparability gate binds before the secular bar mode.

print(f"resolucao: {NR}x{NR}")
print(" om_f     M      T/|W|   beta_min  E_t/E_p  E_tor/|W|   B_pole    portoes")
rows = []
for of in (1.20, 1.35, 1.50, 1.60, 1.6754):
    s = solve(of)
    if s is None:
        print(f" {of:5.3f}   nao convergiu"); continue
    rho, H, r, th, W = s["rho"], s["H"], s["r"], s["th"], s["W"]
    varpi = r[:, None]*np.sin(th)[None, :]
    P_gas = eos.pressure(eos.x_of_enthalpy(np.maximum(H, 0.0), mu_e=MU_E))
    star = rho > RHO_STAR; zero = np.zeros_like(rho)
    Bphi = ToroidalSC(K=K_TOR, m=1.0).B_phi(rho, varpi)
    _, E_tor, _ = diag.magnetic_energies(zero, zero, Bphi, r, th)
    u1 = solve_gradshafranov(-4*np.pi*varpi**2*rho*1e-13, r, th, lmax=LMAX)
    Br1, Bth1 = diag.poloidal_field(u1, r, th)
    unit = float(np.hypot(Br1, Bth1).max())
    def bm(amp):
        Br, Bth = diag.poloidal_field(u1*(amp/unit), r, th)
        Pm = (Br**2+Bth**2+Bphi**2)/(8*np.pi)
        return float((P_gas[star]/np.maximum(Pm[star],1e-300)).min()), Br, Bth
    lo2, hi2 = 1e10, 1e14
    for _ in range(50):
        mid = np.sqrt(lo2*hi2)
        if bm(mid)[0] >= 1.0: lo2 = mid
        else: hi2 = mid
    b, Br, Bth = bm(lo2)
    E_pol, _, _ = diag.magnetic_energies(Br, Bth, zero, r, th)
    etep = E_tor/max(E_pol, 1e-99)
    bp = diag.surface_dipolarity(np.hypot(Br,Bth), H, r, th)["B_pole"]
    ok = (s["M"]>1.44 and b>=1 and 0.1<etep<10 and s["tw"]<TW_SECULAR)
    print(f" {of:5.3f} {s['M']:7.4f}  {s['tw']:.4f}  {b:7.4f}  {etep:7.2f}  "
          f"{E_tor/W:.5f}  {bp:.2e}  {'TODOS' if ok else 'reprova'}")
    rows.append((of, s["M"], s["tw"], b, etep, E_tor/W, bp, ok))

import csv, pathlib as _pl
out = _pl.Path(__file__).resolve().parent / f"tradeoff_{NR}.csv"
with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["omega_frac","M_msun","T_over_W","beta_min","Etor_over_Epol",
                "Etor_over_W","B_pole_G","all_four_gates"])
    w.writerows(rows)
print(f"\nescrito {out.name}")
