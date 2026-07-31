#!/usr/bin/env bash
# Follow the Phase 1 production run without typing qstat over and over.
#
#   bash ~/WD_MAG/cluster/cenapad/watch_job.sh
#
# Prints one line per check: queue state, elapsed steps, simulated time, and
# the central density relative to its initial value -- which is the quantity
# that decides whether the measurement window is still open.
set -u

DIR="${1:-$HOME/wd-mag/Castro/Exec/science/wd_scf_stability}"
LOG="$DIR/run_prod96.log"
EVERY="${2:-120}"

while :; do
    STATE="$(qstat -u "$USER" 2>/dev/null | awk '/wdscf96/ {print $10}')"
    if [ -f "$LOG" ]; then
        STEPS=$(grep -c '^STEP' "$LOG" 2>/dev/null || echo 0)
        LAST=$(grep '^STEP' "$LOG" 2>/dev/null | tail -1)
        RHO=$(grep 'MAXIMUM DENSITY' "$LOG" 2>/dev/null | tail -1 | awk '{print $NF}')
        RHO0=$(grep 'MAXIMUM DENSITY' "$LOG" 2>/dev/null | head -1 | awk '{print $NF}')
        REL=$(awk -v a="$RHO" -v b="$RHO0" 'BEGIN{if(b>0) printf "%.4f", a/b; else print "?"}')
        printf '%s  state=%-3s steps=%-6s rho_c/rho_c0=%-8s  %s\n' \
               "$(date +%H:%M:%S)" "${STATE:-done}" "$STEPS" "$REL" "$LAST"
    else
        printf '%s  state=%-3s (waiting for the run to start)\n' \
               "$(date +%H:%M:%S)" "${STATE:-queued}"
    fi
    [ -z "$STATE" ] && [ -f "$LOG" ] && { echo "job finished"; break; }
    sleep "$EVERY"
done
