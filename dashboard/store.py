"""Persistence and provenance for the dashboard (rules R2, R3).

Each run becomes a directory runs/<hash>/:
    params.json    -- complete input parameters
    scalars.json   -- derived scalars (M, R_eq, VE, E_mag/|W|, ...)
    fields.npz     -- rho, Phi, u, H, Bphi on the (r,theta) grid, r, theta
    manifest.json  -- git hash, dependency versions, timestamp

Plus an index runs/index.csv so Tab 4 (run registry) can load quickly
without opening every directory.

R3: this module only PERSISTS what the dashboard computed via scf.*. It
never launches Castro runs nor decides physics.
"""

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent  # wd-magnetizada/
DEFAULT_RUNS_DIR = REPO_ROOT / "dashboard" / "runs"

# Bumped when the SET OF SCALARS a run carries changes shape (not for
# ordinary bug fixes to existing scalars). Rotation + self-consistent
# toroidal added T, T/|W|, Omega_c, mass_loss_ratio, rotation_period_s to
# scalars.json -- runs saved before this bump lack those columns. This was
# a documented gap (docs/teoria.md Sec 7: a real regression happened once
# during development, silently, before this check existed) -- run_exists()
# now treats a schema mismatch as a cache miss instead of silently
# returning a run with missing columns.
#
# v2 -> v3: Tab 2 (sweep) grid extended beyond (rho_c, k0) to rotation
# (rigid/differential) and the self-consistent toroidal field (K, m) --
# sweep_worker.run_one() now also returns B_tor,max, T, T/|W|, the two
# Bt/Bp ratios, equatorial mass-loss ratio, and rho_c_valid (item 3's
# neutronization gate flag). Runs cached under v2 lack every one of these.
SCHEMA_VERSION = 3


def params_hash(params: dict) -> str:
    """Stable hash of the parameters — cache key and run directory name."""
    payload = json.dumps(params, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()[:12]


def git_commit_hash(path=REPO_ROOT) -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:
        pass
    return "no-git"


def git_dirty(path=REPO_ROOT) -> bool:
    try:
        out = subprocess.run(
            ["git", "-C", str(path), "status", "--porcelain"],
            capture_output=True, text=True, timeout=5,
        )
        return bool(out.stdout.strip())
    except Exception:
        return True


def dependency_versions() -> dict:
    versions = {"python": sys.version.split()[0]}
    for mod in ("numpy", "scipy", "streamlit", "plotly", "h5py"):
        try:
            versions[mod] = __import__(mod).__version__
        except Exception:
            versions[mod] = "not installed"
    return versions


def _runs_dir(runs_dir=None) -> Path:
    d = Path(runs_dir) if runs_dir else DEFAULT_RUNS_DIR
    d.mkdir(parents=True, exist_ok=True)
    return d


def run_exists(params: dict, runs_dir=None) -> str | None:
    """Returns the hash if a run with these parameters AND the current
    SCHEMA_VERSION already exists, else None. A schema mismatch is a
    cache miss (forces recompute) rather than silently returning a run
    whose scalars.json is missing newer columns (see SCHEMA_VERSION)."""
    h = params_hash(params)
    run_dir = _runs_dir(runs_dir) / h
    if not run_dir.exists():
        return None
    try:
        with open(run_dir / "manifest.json") as f:
            manifest = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return None
    if manifest.get("schema_version") != SCHEMA_VERSION:
        return None
    return h


def save_run(params: dict, scalars: dict, fields: dict, runs_dir=None) -> str:
    """Saves a complete run. fields: dict[str, np.ndarray]. Returns the hash."""
    h = params_hash(params)
    run_dir = _runs_dir(runs_dir) / h
    run_dir.mkdir(parents=True, exist_ok=True)

    with open(run_dir / "params.json", "w") as f:
        json.dump(params, f, indent=2, sort_keys=True)
    with open(run_dir / "scalars.json", "w") as f:
        json.dump(scalars, f, indent=2, sort_keys=True, default=float)

    np.savez_compressed(run_dir / "fields.npz", **fields)

    manifest = {
        "hash": h,
        "schema_version": SCHEMA_VERSION,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit_hash(),
        "git_dirty": git_dirty(),
        "dependencies": dependency_versions(),
    }
    with open(run_dir / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)

    _append_index(h, params, scalars, manifest, runs_dir)
    return h


def load_run(run_hash: str, runs_dir=None) -> dict:
    run_dir = _runs_dir(runs_dir) / run_hash
    with open(run_dir / "params.json") as f:
        params = json.load(f)
    with open(run_dir / "scalars.json") as f:
        scalars = json.load(f)
    with open(run_dir / "manifest.json") as f:
        manifest = json.load(f)
    fields = dict(np.load(run_dir / "fields.npz"))
    return {"hash": run_hash, "params": params, "scalars": scalars,
            "manifest": manifest, "fields": fields}


def _index_path(runs_dir=None) -> Path:
    return _runs_dir(runs_dir) / "index.csv"


def _append_index(h, params, scalars, manifest, runs_dir=None):
    row = {"hash": h, "timestamp": manifest["timestamp"],
           "schema_version": manifest.get("schema_version"),
           "git_commit": manifest["git_commit"], "reference": False}
    row.update({f"param_{k}": v for k, v in params.items()})
    row.update(scalars)

    index_path = _index_path(runs_dir)
    if index_path.exists():
        df = pd.read_csv(index_path)
        df = df[df["hash"] != h]  # replace if it already exists (re-run)
        df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    else:
        df = pd.DataFrame([row])
    df.to_csv(index_path, index=False)


def load_index(runs_dir=None) -> pd.DataFrame:
    index_path = _index_path(runs_dir)
    if not index_path.exists():
        return pd.DataFrame()
    return pd.read_csv(index_path)


def mark_reference(run_hash: str, is_reference: bool = True, runs_dir=None):
    index_path = _index_path(runs_dir)
    if not index_path.exists():
        return
    df = pd.read_csv(index_path)
    df.loc[df["hash"] == run_hash, "reference"] = is_reference
    df.to_csv(index_path, index=False)
