#!/usr/bin/env bash
# Build the wd_scf_stability problem on CENAPAD-SP (Lovelace).
#
# Run this ON THE CLUSTER, from a directory with plenty of quota. It clones
# Castro at the exact version this project uses, pins the submodules to the
# exact commits, re-applies the two core patches this problem needs, and drops
# our problem in.
#
# The patches are NOT optional. Without Gravity.cpp, the half-shift geometry
# that puts a cell centre at r = 0 makes Castro divide by zero and abort with
# NaN density on the first step. Without Castro_io.cpp, the ztwd EOS fails to
# build ('network_rp has not been declared') because ztwd has no runtime
# parameters and so never pulls in extern_parameters.H transitively.
set -euo pipefail

CASTRO_TAG=26.07
AMREX_SHA=6ff3cc07b020281ba189de79815a0e1ccd57b446
MICRO_SHA=b489615a42cf9dc3f3ec12d377f3c5b48fc45d7f

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${1:-$PWD/wd-mag}"
mkdir -p "$ROOT" && cd "$ROOT"

if [ ! -d Castro ]; then
    git clone https://github.com/AMReX-Astro/Castro.git
fi
cd Castro
git fetch --tags
git checkout "$CASTRO_TAG"

git submodule update --init external/amrex external/Microphysics
git -C external/amrex checkout "$AMREX_SHA"
git -C external/Microphysics checkout "$MICRO_SHA"

# idempotent: skip if already applied
if git apply --check "$HERE/castro_core.patch" 2>/dev/null; then
    git apply "$HERE/castro_core.patch"
    echo "core patches applied"
else
    echo "core patches already applied (or do not fit) -- checking"
    git apply --reverse --check "$HERE/castro_core.patch" 2>/dev/null \
        && echo "  already applied, fine" \
        || { echo "  ERROR: patch neither applies nor is applied"; exit 1; }
fi

mkdir -p Exec/science/wd_scf_stability
cp "$HERE"/../../castro_problems/wd_scf_stability/* Exec/science/wd_scf_stability/

echo
echo "Castro is at $ROOT/Castro"
echo "Problem at   $ROOT/Castro/Exec/science/wd_scf_stability"
echo
echo "Still needed there: the model files. From your laptop:"
echo "  scp -P 31459 models/phase1_*.txt \\"
echo "      rcrlima@cenapad.unicamp.br:$ROOT/Castro/Exec/science/wd_scf_stability/"
