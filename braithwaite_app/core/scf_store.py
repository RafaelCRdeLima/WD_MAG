"""Thin import shim to the dashboard's own store.py -- reused, not
reimplemented (per the project rule: no physics or persistence logic
duplicated between the dashboard and this app). This module only adds
`sys.path` wiring so `braithwaite_app/` can import `dashboard/store.py`
and `scf/castro_model_writer.py` without copying them.
"""

import sys

from core.paths import DASHBOARD_DIR, SCF_DIR, REPO_ROOT

for _p in (str(DASHBOARD_DIR), str(SCF_DIR), str(REPO_ROOT)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np  # noqa: E402

import store  # noqa: E402  (dashboard/store.py)
import castro_model_writer  # noqa: E402  (scf/castro_model_writer.py)
import scf as scf_mod  # noqa: E402  (scf/scf.py -- the real solver, same one Tab 1 uses)
import diagnostics as diag  # noqa: E402  (scf/diagnostics.py)
import seed as seed_mod  # noqa: E402  (dashboard/seed.py -- r_guess)
from eos import neutronization_threshold_rho_c  # noqa: E402  (scf/eos.py)

# Same defaults Tab 1's UI ships with (dashboard/pages/1_equilibrium.py:
# Nr=129, Ntheta=129, lmax=16, tol=1e-6, max_iter=200) -- verified this
# session to reproduce the actual science star's model.dat to 1 ULP in an
# unused placeholder column, exact on every density value (129/129 lines).
SCF_GRID_NR = 129
SCF_GRID_NTHETA = 129
SCF_LMAX = 16
SCF_TOL = 1e-6
SCF_MAX_ITER = 200

# Matches castro.small_dens=1e4 in wd_braithwaite's inputs (see
# CASTRO_CORE_PATCHES.md / this session's "FLUFF HANDLING" notes) -- the
# model's own density floor must sit at or below the Castro floor or init
# aborts on floored exterior points.
MODEL_DENSITY_FLOOR_FRACTION = 1.0e-4


def list_background_stars():
    """Field-free, non-rotating equilibria from the SCF registry -- the
    same filter Tab 5's Streamlit skeleton used (param_k0 == param_K_tor
    == param_Omega_c == 0), reused here rather than redefined so the two
    UIs can never disagree about which runs qualify.
    """
    idx = store.load_index()
    if idx.empty:
        return idx

    def _col(df, name):
        return df[name] if name in df.columns else 0.0

    k0 = _col(idx, "param_k0")
    K = _col(idx, "param_K_tor")
    Om = _col(idx, "param_Omega_c")
    mask = (k0 == 0.0) & (K == 0.0) & (Om == 0.0)
    return idx[mask].fillna(0.0)


def load_background_star(run_hash: str) -> dict:
    return store.load_run(run_hash)


def converge_field_free_star(rho_c: float, mu_e: float = 2.0) -> dict:
    """Runs the REAL SCF solver (scf.hachisu_scf, the same function Tab 1
    calls) for a field-free, non-rotating star. Genuinely fast for this
    case -- no k0/rotation self-consistency coupling to slow the Picard
    iteration -- measured 16-25ms this session (10 and 8 iterations to
    tol=1e-6). Called synchronously; does not need to be backgrounded.
    """
    R_guess = seed_mod.r_guess(rho_c)
    r = np.linspace(0, 1.3 * R_guess, SCF_GRID_NR)
    theta = np.linspace(0, np.pi, SCF_GRID_NTHETA)
    rho0 = scf_mod.initial_guess(r, theta, rho_c, R_guess)
    result = scf_mod.hachisu_scf(
        rho0, r, theta, rho_c, mu_e=mu_e,
        lmax=SCF_LMAX, tol=SCF_TOL, max_iter=SCF_MAX_ITER,
    )
    return {"result": result, "r": r, "theta": theta, "R_guess": R_guess}


def compute_ve(scf_out: dict, mu_e: float) -> float:
    """Virial error for a field-free, non-rotating result -- via
    diag.virial_error, the single source of truth for the residual
    formula (not reimplemented; see diagnostics.py's own docstring on
    why duplicating this once already cost the project a bug)."""
    result = scf_out["result"]
    rho, Phi, H = result["rho"], result["Phi"], result["H"]
    r, theta = scf_out["r"], scf_out["theta"]
    zeros = np.zeros_like(rho)
    VE, _W, _Pi, _E_mag, _T2 = diag.virial_error(
        rho, Phi, H, zeros, zeros, zeros, r, theta, mu_e, T=0.0
    )
    return float(VE)


def build_model_dat(rho_c: float, mu_e: float, out_path) -> dict:
    """Converges the real SCF equilibrium and writes model.dat + its
    manifest via castro_model_writer.write_model_file (reused, not
    reimplemented). Raises RuntimeError on non-convergence -- never
    falls back to reusing whatever model.dat happened to already be on
    disk (that silent-reuse gap is exactly what this function replaces).

    Returns the write_model_file manifest, extended with "VE" and
    "scf_params_hash" (the latter for core/persistence.py's star cache
    key, via the same store.params_hash() mechanism Tab 1/2 use).
    """
    scf_out = converge_field_free_star(rho_c, mu_e)
    result = scf_out["result"]
    if not result["converged"]:
        raise RuntimeError(
            f"SCF did not converge for rho_c={rho_c:.3e} g/cm^3, mu_e={mu_e} "
            f"(iterations={result['iterations']}, last delta_rho/rho_c="
            f"{result['history'][-1]:.3e}) -- refusing to write model.dat "
            "or reuse a stale one."
        )

    params = {"rho_c": rho_c, "mu_e": mu_e, "k0": 0.0, "K_tor": 0.0, "Omega_c": 0.0}
    manifest = castro_model_writer.write_model_file(
        scf_out["r"], scf_out["theta"], result["rho"], params, out_path,
        run_hash="braithwaite_app", git_commit=store.git_commit_hash(),
        density_floor_fraction=MODEL_DENSITY_FLOOR_FRACTION,
    )
    manifest["VE"] = compute_ve(scf_out, mu_e)
    scf_params = {
        "rho_c": rho_c, "mu_e": mu_e, "Nr": SCF_GRID_NR, "Ntheta": SCF_GRID_NTHETA,
        "lmax": SCF_LMAX, "tol": SCF_TOL, "max_iter": SCF_MAX_ITER,
    }
    manifest["scf_params_hash"] = store.params_hash(scf_params)
    return manifest


def neutronization_check(rho_c: float, mu_e: float = 2.0) -> dict:
    """Same gate Tab 3's export uses (R5), reused rather than
    reimplemented: rho_c must stay below the inverse beta-decay
    (neutronization) threshold for this mu_e (Boshkayev et al. 2013,
    ApJ 762, 117) -- above it, the constant-mu_e cold EOS this whole
    project uses no longer describes a real white dwarf, and nothing
    downstream (including a Braithwaite relaxation study) is physical.
    No option to force it, same as Tab 3.
    """
    threshold = neutronization_threshold_rho_c(mu_e)
    return {
        "threshold": threshold,
        "ok": rho_c < threshold,
    }
