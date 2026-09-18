"""Field AND rotation: is there a super-Chandrasekhar star holding a twisted torus?

twisted_torus_sweep.py answered the field-only question and the answer was no:
E_tor/|W| tops out near 0.0101 whatever the shape, the field buys +0.017 Msun,
and the frontier converges to M_Ch from below without crossing it (DIARIO 19).

But our own 2.005 Msun model is rotational -- T/|W| = 0.099 against
E_tor/|W| = 0.011 -- so rotation clearly does pass the limit. The open question
is whether a twisted torus SURVIVES in a star that rotation has already carried
past Chandrasekhar, and rotation cuts both ways here: centrifugal support
lowers the gas pressure that has to hold the field, so beta_min could get worse
even as the mass gets better. That is measured here, not assumed.

A point counts as an answer only if it passes FOUR gates at once:

  M > 1.44 Msun          super-Chandrasekhar at all
  beta_min >= 1          the star can hold its own confined toroidal
  0.1 < E_t/E_p < 10     comparable energies, the Braithwaite & Spruit branch
  T/|W| < 0.14           below the secular bar mode, and not shedding

The fourth is the one the field-only sweep never needed and the one that makes
this honest: a star held up by rotation it cannot keep is not a result. The
dynamical bar threshold is 0.27; 0.14 is the secular one, and DIARIO 1 records
the production model at 0.0993, below both.

Rotation follows the production model's law, Omega(varpi) = Omega_c A^2/(A^2 +
varpi^2) with A = R_eq of the non-rotating reference, and Omega_c quoted as a
fraction of that reference's Keplerian frequency -- omega_frac = 1.5 is the
production value.

Run one point with --index; the whole grid with cluster/scicom/job_tt_rot.sh.
"""

import argparse
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
from terms.rotation import Rotation              # noqa: E402
from terms.toroidal_sc import ToroidalSC         # noqa: E402

sys.path.insert(0, str(HERE))
from export_ct_model import confined_flux        # noqa: E402

MU_E = 2.0
M_TOR = 1.0
B_C = 4.414e13
LMAX = 16
RHO_STAR = 1.0e6
RHO_AMBIENT = 2.0e4
M_CH = 1.44                 # mu_e = 2
TW_SECULAR = 0.14           # secular bar mode; dynamical is 0.27
SHED_GATE = 0.9             # same gate export_rotating_model.py uses

K_TOR_REF = 3.245e-3
RHO_C_LIST = (1.0e9, 3.0e9, 5.0e9)
OMEGA_FRAC_LIST = (0.0, 0.5, 1.0, 1.5, 2.0)   # 0 reproduces the field-only run
K_FRAC_LIST = (0.05, 0.11, 0.20, 0.45, 1.00)
A_FRAC = 1.0
M_LIST = (0.0, -1.0, -2.0)
ZETA_LIST = (0.5, 1.5, 3.0)

GRID = [(rc, of, kf) for rc in RHO_C_LIST
        for of in OMEGA_FRAC_LIST for kf in K_FRAC_LIST]


def surface_flux(u, H, th):
    vals = []
    for j in range(len(th)):
        ins = H[:, j] > 0.0
        if ins.any():
            vals.append(u[np.flatnonzero(ins)[-1], j])
    return max(vals) if vals else 0.0


def beta_min_of(P_gas, star, Br, Bth, Bphi):
    P_mag = (Br**2 + Bth**2 + Bphi**2) / (8.0 * np.pi)
    return float((P_gas[star] / np.maximum(P_mag[star], 1e-300)).min())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", type=int, required=True,
                    help=f"point in the grid, 0..{len(GRID)-1}")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--nr", type=int, default=129)
    a = ap.parse_args()
    if not 0 <= a.index < len(GRID):
        raise SystemExit(f"--index must be in 0..{len(GRID)-1}")
    rho_c, omega_frac, k_frac = GRID[a.index]
    K_tor = k_frac * K_TOR_REF
    out = a.out or HERE / f"tt_rot_{a.index:03d}.csv"
    tag = f"{a.index:03d}"

    print(f"[{tag}] rho_c = {rho_c:.2e}, omega_frac = {omega_frac:.2f}, "
          f"K_frac = {k_frac:.2f}", flush=True)

    # The non-rotating reference sets Omega_K and the length scale A, exactly
    # as export_rotating_model.py does, so omega_frac means the same thing
    # here as in the production model's manifest.
    ref, r0, th0, _ = _solve_toroidal_certified(
        rho_c=rho_c, R_guess=r_guess(rho_c), K_tor=0.0, m_tor_sc=M_TOR,
        rotation=None, mu_e=MU_E, Nr_base=a.nr, Ntheta=a.nr, lmax=LMAX,
        tol=1e-8, max_iter=200)
    if ref is None:
        raise SystemExit(f"[{tag}] the non-rotating reference did not converge")
    M_ref = scf_mod.total_mass(ref["rho"], r0, th0)
    R_ref = diag.equatorial_polar_radii(ref["H"], r0, th0)[0]
    om_kep = float(np.sqrt(units.G_CONST * M_ref / R_ref**3))
    rot = (Rotation(omega_frac * om_kep, A_FRAC * R_ref)
           if omega_frac > 0 else None)

    res, r, th, _ = _solve_toroidal_certified(
        rho_c=rho_c, R_guess=r_guess(rho_c), K_tor=K_tor, m_tor_sc=M_TOR,
        rotation=rot, mu_e=MU_E, Nr_base=a.nr, Ntheta=a.nr, lmax=LMAX,
        tol=1e-8, max_iter=400)
    if res is None:
        print(f"[{tag}] the rotating solve did not converge", flush=True)
        out.write_text(f"# index={tag} rho_c={rho_c:.6e} "
                       f"omega_frac={omega_frac:.4f} k_frac={k_frac:.4f} "
                       f"DID NOT CONVERGE\n")
        return

    rho, Phi, H = res["rho"], res["Phi"], res["H"]
    M = units.g_to_msun(scf_mod.total_mass(rho, r, th))
    W = abs(diag.gravitational_energy(rho, Phi, r, th))
    T_rot = rot.energy(rho, r, th)["T"] if rot is not None else 0.0
    t_over_w = T_rot / W
    R_eq, R_pol = diag.equatorial_polar_radii(H, r, th)

    # mass shedding, the same estimate export_rotating_model.py gates on
    if rot is not None:
        jeq = len(th) // 2
        kk = np.flatnonzero(rho[:, jeq] > 0)
        om_eq = float(np.atleast_1d(rot.Omega(np.array([R_eq])))[0])
        grav_eq = abs(float(np.gradient(Phi[:, jeq], r)[kk[-1]]))
        shed = om_eq**2 * R_eq / max(grav_eq, 1e-30)
    else:
        shed = 0.0

    varpi = r[:, None] * np.sin(th)[None, :]
    P_gas = eos.pressure(eos.x_of_enthalpy(np.maximum(H, 0.0), mu_e=MU_E))
    star = rho > RHO_STAR
    inside = H > 0.0
    zero = np.zeros_like(rho)
    Bphi_sc = ToroidalSC(K=K_tor, m=M_TOR).B_phi(rho, varpi)
    _, target, _ = diag.magnetic_energies(zero, zero, Bphi_sc, r, th)

    P_amb = eos.pressure(eos.x_of_enthalpy(eos.enthalpy(
        (RHO_AMBIENT / eos.B_of_mu_e(MU_E)) ** (1.0 / 3.0), mu_e=MU_E),
        mu_e=MU_E))
    B_amb_max = float(np.sqrt(8.0 * np.pi * P_amb))

    print(f"[{tag}] M = {M:.4f} Msun, T/|W| = {t_over_w:.4f}, "
          f"shedding = {shed:.3f}, E_tor/|W| = {target / W:.5f}", flush=True)

    rows = []
    for m in M_LIST:
        try:
            u_shape, _ = confined_flux(rho, r, th, varpi, m, lmax=LMAX)
        except SystemExit as e:
            print(f"[{tag}]   m = {m}: {e}", flush=True)
            continue
        u_s = surface_flux(u_shape, H, th)
        closed = inside & (u_shape > u_s)
        Br1, Bth1 = diag.poloidal_field(u_shape, r, th)
        unit_pol = float(np.hypot(Br1, Bth1).max())

        for zeta in ZETA_LIST:
            u_norm = max(u_shape.max() - u_s, 1e-300)
            w = np.where(closed, (u_shape - u_s) / u_norm, 0.0)
            shape = np.where(varpi > 0,
                             np.power(w, zeta) / np.maximum(varpi, 1e-30), 0.0)
            _, E_unit, _ = diag.magnetic_energies(zero, zero, shape, r, th)
            if E_unit <= 0 or not np.isfinite(E_unit):
                continue
            Bphi = np.sqrt(target / E_unit) * shape if target > 0 else zero
            b_tor_only = beta_min_of(P_gas, star, zero, zero, Bphi)

            if b_tor_only < 1.0:
                rows.append((m, zeta, M, t_over_w, shed, target / W,
                             b_tor_only, np.nan, np.nan, np.nan, np.nan,
                             "toroidal alone fails"))
                continue

            lo, hi = 1e6, 1e15
            for _ in range(50):
                mid = np.sqrt(lo * hi)
                Br, Bth = diag.poloidal_field(u_shape * (mid / unit_pol), r, th)
                if beta_min_of(P_gas, star, Br, Bth, Bphi) >= 1.0:
                    lo = mid
                else:
                    hi = mid
            Br, Bth = diag.poloidal_field(u_shape * (lo / unit_pol), r, th)
            E_pol, _, _ = diag.magnetic_energies(Br, Bth, zero, r, th)
            etep = target / max(E_pol, 1e-99)
            bp_ext = diag.surface_dipolarity(np.hypot(Br, Bth), H, r, th)["B_pole"]
            btot = float(np.sqrt(Br**2 + Bth**2 + Bphi**2).max())
            beta_amb = (B_amb_max / max(bp_ext, 1e-300)) ** 2

            # the four gates, together
            ok = (M > M_CH and 0.1 < etep < 10.0
                  and t_over_w < TW_SECULAR and shed < SHED_GATE)
            verdict = "ALL FOUR GATES" if ok else (
                "bar-unstable" if t_over_w >= TW_SECULAR else
                "shedding" if shed >= SHED_GATE else
                "sub-Chandrasekhar" if M <= M_CH else "split not comparable")
            rows.append((m, zeta, M, t_over_w, shed, target / W, b_tor_only,
                         float(np.abs(Bphi).max()) / lo, etep, bp_ext,
                         beta_amb, verdict))
            print(f"[{tag}]   m={m:5.2f} zeta={zeta:4.2f}  M={M:.4f}  "
                  f"T/W={t_over_w:.4f}  E_t/E_p={etep:8.2f}  {verdict}",
                  flush=True)

    with out.open("w") as f:
        f.write(f"# index={tag} rho_c={rho_c:.6e} omega_frac={omega_frac:.4f} "
                f"k_frac={k_frac:.4f} M={M:.6f} T_over_W={t_over_w:.6f} "
                f"shedding={shed:.6f} W={W:.6e}\n")
        f.write(f"# gates: M > {M_CH}, beta_min >= 1, 0.1 < E_t/E_p < 10, "
                f"T/|W| < {TW_SECULAR}, shedding < {SHED_GATE}\n")
        f.write("m,zeta,M_msun,T_over_W,shedding,Etor_over_W,"
                "beta_min_tor_only,Bt_over_Bp,Etor_over_Epol,B_pole_G,"
                "beta_ambient,verdict\n")
        for row in rows:
            f.write(",".join(f"{v:.6e}" if isinstance(v, float) else str(v)
                             for v in row) + "\n")
    print(f"[{tag}] wrote {out.name} ({len(rows)} rows)", flush=True)


if __name__ == "__main__":
    main()
