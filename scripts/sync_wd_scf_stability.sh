#!/usr/bin/env bash
# Syncs the SCF-field stability Castro problem (Phase 1) between its live build location
# (castro/Exec/science/wd_scf_stability/ -- inside the gitignored castro/
# clone, where it actually compiles and runs) and its versioned mirror
# (castro_problems/wd_scf_stability/ -- tracked in this repo).
#
# castro/ cannot be selectively un-ignored (git: "it is not possible to
# re-include a file if a parent directory of that file is excluded" --
# confirmed directly, see .gitignore's note next to the castro/ line).
# This mirror + sync script is the documented workaround.
#
# Only source files are synced -- never build/run artifacts
# (tmp_build_dir/, *.ex, *.o, plt*, chk*, Backtrace*, the *_diag.out
# logs) or the model file generated per-run by
# scf/castro_model_writer.py (model.dat, model.dat.manifest.json).
#
# Usage:
#   scripts/sync_wd_scf_stability.sh save     # live -> mirror (before committing)
#   scripts/sync_wd_scf_stability.sh restore  # mirror -> live (after a fresh clone/reset)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LIVE="$REPO_ROOT/castro/Exec/science/wd_scf_stability"
MIRROR="$REPO_ROOT/castro_problems/wd_scf_stability"

SOURCE_FILES=(
    GNUmakefile
    Make.package
    _prob_params
    inputs.ic_check
    inputs.evolve
    inputs.control_nofield
    inputs.prod96
    mu2.net
    problem_initialize.H
    problem_initialize_state_data.H
    problem_initialize_mhd_data.H
    problem_source.H
    problem_restart.H
    scf_model.H
    scf_model_data.H
    scf_model_data.cpp
    scf_setup.H
    Problem_Derive.H
    Problem_Derive.cpp
    Problem_Derives.H
)

usage() {
    echo "usage: $0 {save|restore}" >&2
    exit 1
}

[ $# -eq 1 ] || usage

case "$1" in
    save)
        [ -d "$LIVE" ] || { echo "error: $LIVE does not exist -- nothing to save" >&2; exit 1; }
        mkdir -p "$MIRROR"
        for f in "${SOURCE_FILES[@]}"; do
            if [ -f "$LIVE/$f" ]; then
                cp "$LIVE/$f" "$MIRROR/$f"
                echo "saved: $f"
            else
                echo "skipped (not present in live dir): $f"
            fi
        done
        ;;
    restore)
        [ -d "$MIRROR" ] || { echo "error: $MIRROR does not exist -- nothing to restore" >&2; exit 1; }
        mkdir -p "$LIVE"
        for f in "${SOURCE_FILES[@]}"; do
            if [ -f "$MIRROR/$f" ]; then
                cp "$MIRROR/$f" "$LIVE/$f"
                echo "restored: $f"
            else
                echo "skipped (not present in mirror): $f"
            fi
        done
        echo "note: castro/'s own submodules (external/amrex, external/Microphysics)" \
             "still need 'git submodule update --init' if this is a fresh clone."
        ;;
    *)
        usage
        ;;
esac
