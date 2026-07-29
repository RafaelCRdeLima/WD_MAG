"""The shared, protected, versioned store -- star rows and seed rows
together, one file, never deleted by any Castro-run cleanup (lives under
core/paths.RESULTS_DIR, outside castro/Exec/science/wd_braithwaite/).
This is the fix for what actually went wrong this session: 7 of 10
seeds' numbers lived only in conversation text after their plotfiles
were deleted. Plotfiles are disposable; rows in this file are not.

Cache-key discipline mirrors dashboard/store.py's SCHEMA_VERSION
precedent exactly: STAR_CACHE_SCHEMA_VERSION bumps whenever the star
construction PHYSICS changes (not on ordinary bug fixes), and rows
written under an older version are treated as cache misses, not
silently reused -- store.py already learned this lesson once (v3->v4
wiped 86 runs after the Poisson solver was fixed).
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from core.paths import DIPOLE_DIAGNOSTIC_CSV, PRERELAX_DIAGNOSTIC_CSV, RESULTS_CSV
from core.scf_store import store  # dashboard/store.py, reused for params_hash + git_commit_hash

STAR_CACHE_SCHEMA_VERSION = 1

_COLUMNS = [
    "row_type", "schema_version", "timestamp", "star_cache_key",
    # star rows
    "rho_c", "mu_e", "resolution", "scf_params_hash", "git_commit_scf",
    "git_commit_wd_braithwaite", "window_valid", "window_lo", "window_hi",
    "window_reason", "VE", "gravity_patch_executable_hash",
    # seed rows
    "seed", "e_mag_over_w_target", "t_ttdyn_measured", "E_mag_over_W",
    "E_tor_over_Emag", "divB_interior_max", "plotfile_path",
]


def _load_df() -> pd.DataFrame:
    if not RESULTS_CSV.exists():
        return pd.DataFrame(columns=_COLUMNS)
    return pd.read_csv(RESULTS_CSV)


def _append_row(row: dict):
    df = _load_df()
    row = {**{c: None for c in _COLUMNS}, **row}
    df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)


def star_cache_key(rho_c: float, mu_e: float, resolution: int,
                    scf_params_hash: str = "", git_commit_scf: str = "",
                    git_commit_wd_braithwaite: str = "") -> str:
    """Everything that defines the star -- not just (rho_c, resolution).
    Reuses store.params_hash(), the same mechanism Tab 1/2 already use,
    rather than inventing a new hashing scheme.
    """
    payload = {
        "rho_c": rho_c, "mu_e": mu_e, "resolution": resolution,
        "scf_params_hash": scf_params_hash, "git_commit_scf": git_commit_scf,
        "git_commit_wd_braithwaite": git_commit_wd_braithwaite,
        "star_cache_schema_version": STAR_CACHE_SCHEMA_VERSION,
    }
    return store.params_hash(payload)


def save_star_result(cache_key: str, rho_c: float, mu_e: float, resolution: int,
                      window_result: dict, VE: float = float("nan"),
                      scf_params_hash: str = "", git_commit_scf: str = "",
                      git_commit_wd_braithwaite: str = "",
                      gravity_patch_executable_hash: str = "") -> None:
    row = {
        "row_type": "star",
        "schema_version": STAR_CACHE_SCHEMA_VERSION,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "star_cache_key": cache_key,
        "rho_c": rho_c, "mu_e": mu_e, "resolution": resolution,
        "scf_params_hash": scf_params_hash, "git_commit_scf": git_commit_scf,
        "git_commit_wd_braithwaite": git_commit_wd_braithwaite,
        "window_valid": window_result["valid"],
        "window_lo": window_result["window"][0] if window_result["valid"] else None,
        "window_hi": window_result["window"][1] if window_result["valid"] else None,
        "window_reason": window_result.get("reason"),
        "VE": VE,
        "gravity_patch_executable_hash": gravity_patch_executable_hash,
    }
    _append_row(row)


def load_star_result(cache_key: str) -> dict | None:
    """Cache lookup -- returns None (cache miss) if absent OR if the
    stored row's schema_version doesn't match STAR_CACHE_SCHEMA_VERSION,
    same "mismatch = miss" convention as store.py's run_exists()."""
    df = _load_df()
    hits = df[(df["row_type"] == "star") & (df["star_cache_key"] == cache_key)]
    if hits.empty:
        return None
    hits = hits[hits["schema_version"] == STAR_CACHE_SCHEMA_VERSION]
    if hits.empty:
        return None
    return hits.iloc[-1].to_dict()


def save_seed_result(star_cache_key_value: str, seed: int, resolution: int,
                      e_mag_over_w_target: float, measurement: dict,
                      plotfile_path: str) -> None:
    row = {
        "row_type": "seed",
        "schema_version": STAR_CACHE_SCHEMA_VERSION,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "star_cache_key": star_cache_key_value,
        "resolution": resolution,
        "seed": seed,
        "e_mag_over_w_target": e_mag_over_w_target,
        "t_ttdyn_measured": measurement["t_ttdyn"],
        "E_mag_over_W": measurement["E_mag_over_W"],
        "E_tor_over_Emag": measurement["E_tor_over_Emag"],
        "divB_interior_max": measurement["divB_interior_max"],
        "plotfile_path": plotfile_path,
    }
    _append_row(row)


def load_all_results() -> pd.DataFrame:
    return _load_df()


def load_seeds_for_star(star_cache_key_value: str) -> pd.DataFrame:
    df = _load_df()
    return df[(df["row_type"] == "seed") & (df["star_cache_key"] == star_cache_key_value)]


# --- Dipole diagnostic: a SEPARATE file/schema from the rows above, on purpose. ---
#
# Every number here (core/dipole.py's m(r) measurement) was taken while
# the background star was already outside its own validity window
# (core/star_builder.py::find_measurement_window) -- confirmed on two
# independent seeds (42ext, 123ext) sharing one background star:
# X_2pct = 0.034-0.036 t/t_dyn, while the field itself only reaches its
# own quasi-steady band by t_field_relax~0.4 t/t_dyn. The star leaves
# validity roughly 10x faster than the field needs just to relax, let
# alone however much longer an exterior dipole would need to organize
# (m(r) was still not r-independent, and still visibly growing at the
# largest radius tested, at t/t_dyn~4 in both runs -- see docs/teoria.md
# for the full writeup). This is a LIMITATION OF THE MEASUREMENT SETUP
# (same well-balancing gap as everywhere else in this project, Sec 6.8),
# not a physical finding that no exterior dipole exists -- the question
# was never actually tested on a star that stayed physical long enough
# to answer it. That is exactly why this does NOT live in RESULTS_CSV
# next to the (in-window, trustworthy) E_tor/E_mag numbers: mixing a
# structurally-caveated number into the same table as validated ones is
# how a caveat gets silently lost the next time someone reads the file.
DIPOLE_DIAGNOSTIC_SCHEMA_VERSION = 1
DIPOLE_CAVEAT = (
    "Measured with the background star already OUTSIDE its validity window "
    "(t/t_dyn > X_2pct). Real, non-noise field (verified: absolute B_r on the "
    "sampling shell is structured and well above the numerical floor, and "
    "exterior |B| never exceeds interior |B| at any tested radius -- see "
    "session notes). NOT a measurement of the relaxed field's true exterior "
    "structure: m(r) is not r-independent and was still growing at the "
    "largest radius tested. Tool limitation (Castro's MHD path has no "
    "well-balancing, docs/teoria.md Sec 6.8), not evidence against a "
    "physical exterior dipole."
)

_DIPOLE_COLUMNS = [
    "schema_version", "timestamp", "run_id", "seed", "plotfile",
    "t_s", "t_ttdyn", "rho_c_ic", "window_valid", "window_x2pct_ttdyn",
    "t_field_relax_ttdyn", "r_star_cm",
    "m_mag_1.2R_G_cm3", "m_mag_1.5R_G_cm3", "m_mag_2.0R_G_cm3",
    "b_pole_implied_1.2R_G", "ratio_1.2R_over_1.5R", "ratio_1.5R_over_2.0R",
    "r_independent", "caveat",
]


def _load_dipole_df() -> pd.DataFrame:
    if not DIPOLE_DIAGNOSTIC_CSV.exists():
        return pd.DataFrame(columns=_DIPOLE_COLUMNS)
    return pd.read_csv(DIPOLE_DIAGNOSTIC_CSV)


def save_dipole_diagnostic_row(
    run_id: str, seed: int, plotfile: str, t_s: float, t_ttdyn: float,
    rho_c_ic: float, window_valid: bool, window_x2pct_ttdyn: float,
    t_field_relax_ttdyn: float, r_star_cm: float,
    m_mag_by_radius: dict, r_independent: bool,
) -> None:
    """One row per (run, timestep). `m_mag_by_radius` keys are the tested
    radii as multiples of R_star, e.g. {1.2: ..., 1.5: ..., 2.0: ...} in
    G*cm^3 (core/dipole.py's output units) -- caller passes exactly what
    dipole_moment_from_grid() returned, this function does no physics of
    its own, only ratios for the r-independence columns.
    """
    m12 = m_mag_by_radius.get(1.2)
    m15 = m_mag_by_radius.get(1.5)
    m20 = m_mag_by_radius.get(2.0)
    row = {
        "schema_version": DIPOLE_DIAGNOSTIC_SCHEMA_VERSION,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "run_id": run_id, "seed": seed, "plotfile": plotfile,
        "t_s": t_s, "t_ttdyn": t_ttdyn, "rho_c_ic": rho_c_ic,
        "window_valid": window_valid, "window_x2pct_ttdyn": window_x2pct_ttdyn,
        "t_field_relax_ttdyn": t_field_relax_ttdyn, "r_star_cm": r_star_cm,
        "m_mag_1.2R_G_cm3": m12, "m_mag_1.5R_G_cm3": m15, "m_mag_2.0R_G_cm3": m20,
        "b_pole_implied_1.2R_G": (2.0 * m12 / r_star_cm**3) if m12 and r_star_cm else None,
        "ratio_1.2R_over_1.5R": (m12 / m15) if m12 and m15 else None,
        "ratio_1.5R_over_2.0R": (m15 / m20) if m15 and m20 else None,
        "r_independent": r_independent,
        "caveat": DIPOLE_CAVEAT,
    }
    df = _load_dipole_df()
    row = {**{c: None for c in _DIPOLE_COLUMNS}, **row}
    df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    DIPOLE_DIAGNOSTIC_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(DIPOLE_DIAGNOSTIC_CSV, index=False)


def load_dipole_diagnostics() -> pd.DataFrame:
    return _load_dipole_df()


# --- Pre-relaxation diagnostic: same "separate file, caveat baked into every row" ---
# discipline as the dipole diagnostic above, for a different reason. The
# measurement itself (rho_c(t) during a hydro pre-relaxation attempt with
# castro.use_pslope=1 + castro.grav_source_type=4) is solid: near-flat
# from t/t_dyn~5.7 to ~19, then a monotonic, accelerating, non-oscillating
# rise to >1,000,000% of rho_c_ic by t/t_dyn~35, starting almost exactly
# at the damping ramp-off window (18-20). What is NOT yet settled is WHY:
#
#   (A) damping was masking a pure numerical instability in the
#       pslope+grav_source_type=4 combination -- it would collapse at
#       ANY damping ramp-off time, even one early enough that rho_c is
#       still near rho_c_ic.
#   (B) the drift itself pushed rho_c toward a genuine physical
#       threshold, and damping was suppressing a real collapse that
#       only becomes inevitable once rho_c has climbed far enough --
#       ramp-off timing wouldn't matter, rho_c level would.
#
# Both predict the same qualitative shape (flat, then runaway right
# after damping ends) from THIS ONE RUN alone -- indistinguishable
# without a second run that moves the ramp-off time independently of
# where rho_c happens to be. See the early-ramp-off test (a second row
# of runs in this same file, run_id containing "earlyramp") for the
# discriminating measurement. Do not write "well-balancing doesn't work"
# (hypothesis A) into docs/teoria.md Sec 6.8 until that test's rows are
# also in this file and confirm it.
PRERELAX_DIAGNOSTIC_SCHEMA_VERSION = 1
PRERELAX_CAUSE_STATUS_UNCONFIRMED = "unconfirmed -- pending early-ramp-off discriminating test"
PRERELAX_CAVEAT = (
    "Real rho_c(t) measurement (castro.use_pslope=1 + grav_source_type=4, hydro-only "
    "build, field-free). The COLLAPSE ITSELF is measured fact: near-flat for ~14 "
    "t_dyn, then monotonic accelerating runaway starting at damping ramp-off. The "
    "CAUSE is not yet decided between (A) pure numerical instability masked by "
    "damping [would collapse at ANY ramp-off time] vs (B) drift pushing rho_c to a "
    "genuine physical threshold [ramp-off timing wouldn't matter, rho_c level would] "
    "-- ambiguous from a single run's ramp-off timing. See this module's docstring "
    "above the save function. Do not cite as 'well-balancing failed' until the "
    "early-ramp-off discriminating test's rows are also present."
)

# Set once the early-ramp-off test (run_id "prerelax_v2_earlyramp") came
# back: damping ramped off at t/t_dyn~6, with rho_c still at -7.6% dev
# (i.e. NOT elevated toward any threshold -- if anything below rho_c_ic)
# -- and the same runaway began within ~1.5 t_dyn of ramp-off completion
# (dev crosses from -6.0% at t/t_dyn=7.02 to +79% by 9.06). rho_c level
# cannot be the trigger here since it never climbed before the collapse
# started; ramp-off timing alone predicts the onset in both this run and
# the original (ramp at t/t_dyn~20 -> collapse then; ramp at ~6 ->
# collapse then). Hypothesis (A) confirmed, (B) ruled out.
PRERELAX_CAUSE_STATUS_CONFIRMED_A = (
    "CONFIRMED (A): damping was masking a pure numerical instability in "
    "castro.use_pslope=1 + grav_source_type=4 -- collapse onset tracks damping "
    "ramp-off time, not rho_c level. Discriminating test: prerelax_v2_earlyramp "
    "(ramp-off at t/t_dyn~6, rho_c still at -7.6%, collapse begins within ~1.5 "
    "t_dyn of ramp-off) vs. prerelax_v1_extended_damping (ramp-off at t/t_dyn~20, "
    "rho_c at +9.5%, collapse begins then). Hypothesis (B) (drift-to-threshold) "
    "ruled out: rho_c had NOT climbed before the second run's collapse started."
)

_PRERELAX_COLUMNS = [
    "schema_version", "timestamp", "run_id", "t_s", "t_ttdyn",
    "rho_c", "rho_c_ic", "dev_pct",
    "damping_end_time_ttdyn", "damping_ramp_start_ttdyn",
    "cause_status", "caveat",
]


def _load_prerelax_df() -> pd.DataFrame:
    if not PRERELAX_DIAGNOSTIC_CSV.exists():
        return pd.DataFrame(columns=_PRERELAX_COLUMNS)
    return pd.read_csv(PRERELAX_DIAGNOSTIC_CSV)


def save_prerelax_diagnostic_series(
    run_id: str, rho_c_series: list, rho_c_ic: float, t_dyn_s: float,
    damping_end_time_ttdyn: float, damping_ramp_start_ttdyn: float,
    cause_status: str = PRERELAX_CAUSE_STATUS_UNCONFIRMED,
    caveat: str = PRERELAX_CAVEAT,
) -> None:
    """One row per (run_id, timestep). `rho_c_series` is exactly what
    core/star_builder.py::parse_rho_c_log() returns -- [(t_ttdyn, rho_c), ...]
    -- this function does no physics of its own, only packages it with
    the caveat and provenance needed to interpret it responsibly later.
    `cause_status`/`caveat` default to the unconfirmed pair; pass the
    CONFIRMED_A constants once (and only once) the discriminating test
    result is in hand.
    """
    timestamp = datetime.now(timezone.utc).isoformat()
    rows = []
    for t_ttdyn, rho_c in rho_c_series:
        rows.append({
            "schema_version": PRERELAX_DIAGNOSTIC_SCHEMA_VERSION,
            "timestamp": timestamp, "run_id": run_id,
            "t_s": t_ttdyn * t_dyn_s, "t_ttdyn": t_ttdyn,
            "rho_c": rho_c, "rho_c_ic": rho_c_ic,
            "dev_pct": 100.0 * (rho_c - rho_c_ic) / rho_c_ic,
            "damping_end_time_ttdyn": damping_end_time_ttdyn,
            "damping_ramp_start_ttdyn": damping_ramp_start_ttdyn,
            "cause_status": cause_status,
            "caveat": caveat,
        })
    df = _load_prerelax_df()
    new_df = pd.DataFrame(rows, columns=_PRERELAX_COLUMNS)
    df = pd.concat([df, new_df], ignore_index=True)
    PRERELAX_DIAGNOSTIC_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(PRERELAX_DIAGNOSTIC_CSV, index=False)


def mark_prerelax_cause_confirmed_a() -> int:
    """Retroactively updates every existing row's cause_status/caveat to
    the CONFIRMED_A finding -- called once, after the early-ramp-off
    test's own rows are already saved, so the file never has "unconfirmed"
    sitting next to the very data that confirmed it. Returns the number
    of rows updated.
    """
    df = _load_prerelax_df()
    if df.empty:
        return 0
    n = int((df["cause_status"] == PRERELAX_CAUSE_STATUS_UNCONFIRMED).sum())
    df.loc[df["cause_status"] == PRERELAX_CAUSE_STATUS_UNCONFIRMED, "cause_status"] = PRERELAX_CAUSE_STATUS_CONFIRMED_A
    df.loc[df["caveat"] == PRERELAX_CAVEAT, "caveat"] = PRERELAX_CAUSE_STATUS_CONFIRMED_A
    df.to_csv(PRERELAX_DIAGNOSTIC_CSV, index=False)
    return n


def load_prerelax_diagnostics() -> pd.DataFrame:
    return _load_prerelax_df()
