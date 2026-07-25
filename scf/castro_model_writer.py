"""Writes a saved, field-free, non-rotating SCF equilibrium as a 1D radial
initial model file in Castro's Util/model_parser plain-text format (the
"new", single-header-line style: `# r density temperature pressure X`),
so a Castro problem can read it with read_model_file() exactly like
Castro's own StarGrav/toy_convect problems do.

Pure functions only (numpy arrays + plain dicts in, text file out) --
no dependency on dashboard/store.py, per R1 (physics lives in scf/,
dashboard orchestrates). The CLI at the bottom of this file is the only
place that reaches into dashboard/store.py, and only when run directly.

The written temperature/pressure columns are PLACEHOLDERS (uniform, low
density-independent T), not the SCF's actual cold-degenerate-EOS values:
the Castro problem this feeds (Exec/science/wd_braithwaite) uses
gamma_law, a different EOS entirely, so matching the SCF's real pressure
profile into it would be physically inconsistent regardless of what's
written here. Density -- the only quantity the Braithwaite field-seed
generator's confinement actually depends on -- is the SCF's real,
converged profile, unmodified. This is a deliberate, documented scope
limit for Step 2 (see the Braithwaite plan: this step delivers the
static B-field initial condition, not a Castro-native reproduction of
the SCF's stellar structure).
"""

import json
from pathlib import Path

import numpy as np

# R*T with R = k_B/m_u in erg/(g*K), mu=1 -- only used to fill the (unused
# by the consuming problem) pressure column with a value of the right
# order of magnitude, so model_parser's format detector/warnings behave
# like they do for every other Castro model file.
_R_GAS_CGS = 8.31446e7
PLACEHOLDER_TEMPERATURE_K = 1.0e6


def check_field_free_non_rotating(params: dict, run_hash: str = "<run>") -> None:
    """Raises ValueError if params carries a nonzero field or rotation
    parameter. The Braithwaite background star (plan Step 1) must have
    neither -- Step 2 seeds the field itself."""
    for key in ("k0", "K_tor", "Omega_c"):
        val = params.get(key, 0.0) or 0.0
        if val:
            raise ValueError(
                f"run {run_hash} is not field-free/non-rotating ({key}={val}) "
                "-- the Braithwaite background star must have no imposed "
                "field and no rotation (Step 1 of the plan)."
            )


def check_spherical_symmetry(rho: np.ndarray, rho_c: float, rel_tol: float = 1e-3) -> None:
    """Raises ValueError if rho(r, theta) varies with theta beyond
    rel_tol * rho_c. Field-free + non-rotating implies spherical symmetry
    -- this verifies it instead of assuming it, since reducing to a 1D
    profile silently drops any theta-dependence there might be.

    rel_tol defaults to 1e-3, not machine precision: a converged SCF
    solution carries its own residual (finite lmax truncation, grid
    discretization) of a few x 1e-4 relative even for a field-free,
    non-rotating star -- this check is catching genuine asymmetry (a
    rotating or toroidal-field star, O(1) relative), not solver noise.
    """
    spread = rho.max(axis=1) - rho.min(axis=1)
    tol = rel_tol * max(rho_c, 1.0)
    bad = spread > tol
    if np.any(bad):
        raise ValueError(
            f"not spherically symmetric to tolerance {tol:.3e} "
            f"(max spread over theta = {spread.max():.3e}) -- cannot "
            "reduce to a 1D model file."
        )


def build_model_lines(r: np.ndarray, rho_r: np.ndarray, rho_c: float,
                       density_floor_fraction: float = 1e-10) -> list[str]:
    """Returns the lines of a Castro model_parser-format text file for a
    single species "X" (matching Exec/mhd_tests/*/gammalaw.net)."""
    rho_floored = np.maximum(rho_r, rho_c * density_floor_fraction)
    lines = ["# r density temperature pressure X"]
    for ri, rhoi in zip(r, rho_floored):
        p_nominal = rhoi * _R_GAS_CGS * PLACEHOLDER_TEMPERATURE_K
        lines.append(
            f"{ri:.10e} {rhoi:.10e} {PLACEHOLDER_TEMPERATURE_K:.6e} "
            f"{p_nominal:.10e} 1.0"
        )
    return lines


def write_model_file(r: np.ndarray, theta: np.ndarray, rho: np.ndarray,
                      params: dict, out_path, run_hash: str = "<run>",
                      git_commit: str = "no-git",
                      density_floor_fraction: float = 1e-10) -> dict:
    """Validates and writes `out_path`, plus a `<out_path>.manifest.json`
    sidecar with provenance. Returns the manifest dict."""
    check_field_free_non_rotating(params, run_hash)
    rho_c = params["rho_c"]
    check_spherical_symmetry(rho, rho_c)

    # mean over theta, not a single column -- more robust against the
    # small per-column solver noise check_spherical_symmetry tolerates.
    rho_r = rho.mean(axis=1)
    lines = build_model_lines(r, rho_r, rho_c, density_floor_fraction)

    out_path = Path(out_path)
    out_path.write_text("\n".join(lines) + "\n")

    surface_mask = rho_r > rho_c * 1e-6
    R_star_cm = float(r[surface_mask][-1]) if np.any(surface_mask) else float(r[-1])

    manifest = {
        "source_run_hash": run_hash,
        "git_commit_scf": git_commit,
        "rho_c_gcm3": float(rho_c),
        "R_star_cm": R_star_cm,
        "n_points": int(len(r)),
        "placeholder_temperature_K": PLACEHOLDER_TEMPERATURE_K,
        "note": (
            "temperature/pressure columns are gamma_law-EOS placeholders, "
            "not the SCF's cold-degenerate values -- density is the real, "
            "converged SCF profile. See module docstring."
        ),
    }
    (out_path.parent / f"{out_path.name}.manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True)
    )
    return manifest


if __name__ == "__main__":
    import argparse
    import subprocess
    import sys

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_hash", help="hash of a saved, field-free, non-rotating run")
    ap.add_argument("out_path", help="path to write the Castro model file to")
    ap.add_argument("--runs-dir", default=None)
    args = ap.parse_args()

    _repo_root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(_repo_root / "dashboard"))
    import store  # noqa: E402

    run = store.load_run(args.run_hash, runs_dir=args.runs_dir)

    def _git_hash(path):
        try:
            out = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"],
                                  capture_output=True, text=True, timeout=5)
            return out.stdout.strip() if out.returncode == 0 else "no-git"
        except Exception:
            return "no-git"

    manifest = write_model_file(
        run["fields"]["r"], run["fields"]["theta"], run["fields"]["rho"],
        run["params"], args.out_path, run_hash=args.run_hash,
        git_commit=_git_hash(_repo_root),
    )
    print(json.dumps(manifest, indent=2))
