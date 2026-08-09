#!/usr/bin/env bash
# Pm scan: MRI transport against magnetic Prandtl number, net vertical flux.
#
# Why q = 1.5 and not our star's q. Keplerian shear is where Lesur & Longaretti
# (2007) and Fromang et al. (2007) published alpha(Pm), so this run both
# measures the trend we need AND checks our use of the code against numbers
# someone else got. Scanning our own q = 0-2 comes after that check passes;
# doing it first would leave any disagreement unattributable.
#
# Pm = reynolds_magnetic / reynolds. Re is held at 1000 and Rm varied, so the
# viscous scale is fixed across the scan and only the resistive one moves.
#
# CAVEAT, stated before the numbers exist: at 64^3 the top of the scan is
# marginal. Rm = 16000 puts the resistive scale near the grid, and the
# published high-Pm runs use 128^3. Pm = 8 and 16 here are provisional until
# one of them is repeated at 128^3.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNOOPY="$ROOT/src/snoopy-v6.0-official/snoopy"
CFG="$ROOT/src/snoopy-v6.0-official/src/problem/mri/snoopy.cfg"

RE=1000
PM_LIST=(1 2 4 8 16)
T_FINAL=200.0          # ~32 orbits; saturation by ~15, average over the rest
THREADS=4
CONCURRENT=3           # 3 x 4 = 12 cores

mkdir -p "$ROOT/runs/pm_scan"

launch() {
    local pm=$1
    local rm_val=$(( RE * pm ))
    local dir="$ROOT/runs/pm_scan/pm$(printf '%02d' "$pm")"
    mkdir -p "$dir/data"
    sed -e "s/reynolds = 1000.0;/reynolds = ${RE}.0;/" \
        -e "s/reynolds_magnetic = 1000.0;/reynolds_magnetic = ${rm_val}.0;/" \
        -e "s/t_final = 620.0;/t_final = ${T_FINAL};/" \
        -e "s/snapshot_step = 1.0;/snapshot_step = 1000.0;/" \
        -e "s/dump_step = 1.0;/dump_step = 50.0;/" \
        -e "s/timevar_step = 0.1;/timevar_step = 0.1;/" \
        "$CFG" > "$dir/snoopy.cfg"
    ( cd "$dir" && OMP_NUM_THREADS=$THREADS "$SNOOPY" >run.log 2>&1 \
        && echo "done Pm=$pm" || echo "FAILED Pm=$pm" ) &
}

echo "Pm scan: Re=$RE, Pm=${PM_LIST[*]}, t_final=$T_FINAL, ${THREADS}t x ${CONCURRENT}"
n=0
for pm in "${PM_LIST[@]}"; do
    launch "$pm"
    n=$((n + 1))
    if [ $((n % CONCURRENT)) -eq 0 ]; then wait; fi
done
wait
echo "scan complete"
