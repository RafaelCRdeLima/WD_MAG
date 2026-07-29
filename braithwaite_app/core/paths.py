"""Shared path constants for the Braithwaite desktop app.

Kept in one place because the persisted-results directory (data/) is the
project's source of truth for extracted numbers -- it must never sit
under a path that any Castro-run cleanup touches (see the module
docstring in core/persistence.py for why this matters: 7 of 10 seeds
were lost this session when plotfiles were deleted before their numbers
were durably saved anywhere).
"""

from pathlib import Path

APP_DIR = Path(__file__).resolve().parent.parent          # braithwaite_app/
REPO_ROOT = APP_DIR.parent                                 # wd-magnetizada/
DASHBOARD_DIR = REPO_ROOT / "dashboard"
SCF_DIR = REPO_ROOT / "scf"
WD_BRAITHWAITE_DIR = REPO_ROOT / "castro" / "Exec" / "science" / "wd_braithwaite"

RESULTS_DIR = APP_DIR / "data"                              # persisted numbers (source of truth)
RUN_LOGS_DIR = APP_DIR / "data" / "run_logs"
RESULTS_CSV = RESULTS_DIR / "results.csv"
# Deliberately NOT in RESULTS_CSV / _COLUMNS: every row here was measured
# outside its star's validity window (see core/persistence.py's
# save_dipole_diagnostic_row docstring) -- a separate file keeps that
# caveat structural, not just a comment someone could miss.
DIPOLE_DIAGNOSTIC_CSV = RESULTS_DIR / "dipole_diagnostic.csv"
# Same reasoning: a rho_c(t) trace from a hydro pre-relaxation attempt --
# real measurements, but the CAUSE of what they show (damping masking a
# numerical instability vs. drift-triggered physical collapse) is not
# yet decided (see core/persistence.py's save_prerelax_diagnostic_series
# docstring). Separate file/schema so an unconfirmed interpretation can
# never silently read as settled fact.
PRERELAX_DIAGNOSTIC_CSV = RESULTS_DIR / "prerelax_diagnostic.csv"

for _d in (RESULTS_DIR, RUN_LOGS_DIR):
    _d.mkdir(parents=True, exist_ok=True)
