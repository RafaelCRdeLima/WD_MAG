#!/usr/bin/env bash
# q scan: MRI transport against shear rate, at fixed Pm = 4.
#
# Our star's rotation law is Komatsu j-constant, Omega = Omega_c A^2/(A^2+w^2),
# giving q = -dlnOmega/dlnw = 2w^2/(A^2+w^2), so q runs from 0 on the axis to 2
# far out. The Pm scan was run at Keplerian q = 1.5 to check the setup against
# published numbers; this covers the range the star actually spans.
#
# q enters twice, which is why this is the cheaper of the two remaining scans:
#   - it sets the transport we measure here;
#   - it sets MInIT's own coefficient, alpha^MRI = 1 - 4/q.
#
# q = 1.9 rather than 2.0: at q = 2 the flow is Rayleigh-marginal and the MRI
# criterion degenerates. 1.9 stays inside without being indistinguishable.
#
# q = 1.5 at Pm = 4 already exists in runs/pm_scan/pm04 and is not repeated;
# the analysis picks it up from there.
#
# Wavelength check, k_MRI = sqrt(1-(2-q)^2/4) Omega/v_Az with v_Az = bz0 = 0.1:
#   q=0.5  lambda=0.95   q=1.0  lambda=0.73   q=1.9  lambda=0.63
# all inside Lz = 1, so the fastest-growing mode fits at every q.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNOOPY="$ROOT/src/snoopy-v6.0-official/snoopy"
CFG="$ROOT/src/snoopy-v6.0-official/src/problem/mri/snoopy.cfg"

RE=1000
PM=4
Q_LIST=(0.5 1.0 1.9)
T_FINAL=200.0
THREADS=3          # 3 runs x 3 threads = 9 cores, under the 10-core cap
mkdir -p "$ROOT/runs/q_scan"

for q in "${Q_LIST[@]}"; do
    dir="$ROOT/runs/q_scan/q${q/./p}"
    mkdir -p "$dir/data"
    sed -e "s/reynolds = 1000.0;/reynolds = ${RE}.0;/" \
        -e "s/reynolds_magnetic = 1000.0;/reynolds_magnetic = $((RE * PM)).0;/" \
        -e "s/shear = 1.5;/shear = ${q};/" \
        -e "s/t_final = 620.0;/t_final = ${T_FINAL};/" \
        -e "s/snapshot_step = 1.0;/snapshot_step = 1000.0;/" \
        -e "s/dump_step = 1.0;/dump_step = 50.0;/" \
        "$CFG" > "$dir/snoopy.cfg"
    ( cd "$dir" && OMP_NUM_THREADS=$THREADS "$SNOOPY" >run.log 2>&1 \
        && echo "done q=$q" || echo "FAILED q=$q" ) &
done

echo "q scan: Re=$RE, Pm=$PM, q=${Q_LIST[*]}, ${THREADS}t each, all concurrent"
wait
echo "q scan complete"
