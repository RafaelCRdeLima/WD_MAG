"""What does this solver produce when the answer is known to be nothing?

Run:  scf/.venv/bin/python3 investigations/nonbarotropic_noise_floor.py

Three composition gates were built and rebuilt on this scan and all three
rejected everything. Each time the cause turned out to be the measurement
rather than the configuration. The step that should have come first is this
one: run the solver where the answer is known exactly, and see what it gives.

With k0 = 0 the field vanishes, so g = 0 and dP = g du = 0. Pressure is then
constant on equipotentials, rho = -P_r/Phi_r is plain hydrostatic balance, and
the star is the barotropic field-free model the pressure was anchored to. The
composition must come back EXACTLY uniform at mu_e = 2, and the mass must equal
the background mass.

Whatever spread in mu_e appears instead is the floor: the noise the chain
P -> dP/dr -> rho -> mu_e introduces. Below that floor no statement about
stratification means anything, and every gate built so far has been comparing
signal against numbers of that size without knowing it.

The nonzero-field rows then say where, if anywhere, the composition signal
climbs out of it.
"""

import sys
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
for p in (REPO / "scf", REPO / "dashboard", REPO / "investigations"):
    sys.path.insert(0, str(p))
warnings.filterwarnings("ignore")

import diagnostics as diag                       # noqa: E402
import eos                                       # noqa: E402
import nonbarotropic_mass_scan as N              # noqa: E402
import scf as scf_mod                            # noqa: E402
import units                                     # noqa: E402
from gradshafranov import solve_gradshafranov    # noqa: E402
from poisson import solve_poisson                # noqa: E402
from seed import r_guess                         # noqa: E402
from sweep_worker import _solve_toroidal_certified   # noqa: E402

RHO_C = 1.0e9


def run(k0, alpha=1.0):
    res, r, th, _ = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=0.0, m_tor_sc=1.0,
        rotation=None, mu_e=2.0, Nr_base=N.NR, Ntheta=N.NTH, lmax=N.LMAX,
        tol=1e-8, max_iter=200)
    rho0, Phi, H0 = res["rho"], res["Phi"], res["H"]
    M0 = units.g_to_msun(scf_mod.total_mass(rho0, r, th))
    varpi = r[:, None] * np.sin(th)[None, :]
    x0 = eos.x_of_enthalpy(np.maximum(H0, 0.0), 2.0)
    P_axis = alpha * np.where(H0[:, 0] > 0, eos.pressure(x0[:, 0]), 0.0)

    src = -4.0 * np.pi * varpi**2 * rho0 * k0
    u = solve_gradshafranov(src, r, th, lmax=N.LMAX)
    g = np.where(varpi > 0,
                 -src / (4 * np.pi * np.maximum(varpi, 1e-30)**2), 0.0)
    R_pol = diag.equatorial_polar_radii(H0, r, th)[1]
    seeds = np.linspace(r[1], 1.6 * R_pol, 260)

    rho = rho0
    for _ in range(N.N_OUTER):
        P = N.integrate_pressure(P_axis, Phi, g, u, r, th, seeds)
        P_r, _ = N.grads(P, r, th)
        Phi_r, _ = N.grads(Phi, r, th)
        u_r, _ = N.grads(u, r, th)
        safe = np.abs(Phi_r) > 1e-4 * np.abs(Phi_r).max()
        rn = np.where(safe, (g * u_r - P_r) / np.where(safe, Phi_r, 1.0), 0.0)
        for j in range(len(th)):
            k = np.flatnonzero(safe[:, j])
            if k.size:
                rn[:k[0], j] = rn[k[0], j]
        rn = np.where(P > 0, np.maximum(rn, 0.0), 0.0)
        rho = 0.5 * rho + 0.5 * rn
        Phi = 0.6 * Phi + 0.4 * solve_poisson(rho, r, th, lmax=N.LMAX)

    M = units.g_to_msun(scf_mod.total_mass(rho, r, th))
    mu, _ = N.mu_e_from_P_rho(P, rho)
    inside = (rho > 1e-4 * RHO_C) & (P > 0) & np.isfinite(mu)
    Rg = np.broadcast_to(r[:, None], mu.shape)
    edges = np.linspace(0.0, float(Rg[inside].max()), 25)
    prof = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = inside & (Rg >= a) & (Rg < b)
        if m.sum() > 20:
            w = rho[m]
            prof.append(float(np.sum(mu[m] * w) / np.sum(w)))
    prof = np.array(prof)
    Br, Bth = diag.poloidal_field(u, r, th)
    E_pol, _, _ = diag.magnetic_energies(Br, Bth, np.zeros_like(rho), r, th)
    W = abs(diag.gravitational_energy(rho, Phi, r, th))
    return dict(M=M, M0=M0, prof=prof, EpW=E_pol / max(W, 1.0))


def main():
    print(f"noise floor of the non-barotropic solver, rho_c = {RHO_C:.1e}\n")
    print("  k0         M (Msun)  dM/M0     mu_e mean   mu_e spread   E_pol/|W|")
    base = None
    for k0 in (0.0, 5.0e-14, 1.0e-13, 2.0e-13, 5.0e-13, 1.0e-12):
        d = run(k0)
        pr = d["prof"]
        spread = (pr.max() - pr.min()) / pr.mean()
        if base is None:
            base = spread
        print(f"  {k0:.1e}   {d['M']:7.4f}  {(d['M']-d['M0'])/d['M0']:+8.5f}  "
              f"{pr.mean():9.5f}   {spread:10.5f}    {d['EpW']:.6f}"
              + ("   <- FLOOR (no field)" if k0 == 0.0 else
                 f"   {spread/base:5.2f}x floor"))
    print("\n  A configuration whose mu_e spread is not several times the")
    print("  k0 = 0 value is not stratified -- it is the method's own scatter.")


if __name__ == "__main__":
    main()
