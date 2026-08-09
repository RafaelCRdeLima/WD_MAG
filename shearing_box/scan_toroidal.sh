#!/usr/bin/env bash
# Net toroidal flux alongside the vertical one, at our star's field ratio.
#
# ---------------------------------------------------------------------------
# WHAT THIS ANSWERS, AND WHAT IT DOES NOT
# ---------------------------------------------------------------------------
# It DOES answer: does a net toroidal field, at our star's B_tor/B_pol of 2 to
# 9.5, change the transport the VERTICAL-field MRI produces? That is a real
# question and MInIT has no term for it -- its k_MRI is built from v_Az alone.
#
# It does NOT answer: is the non-axisymmetric TOROIDAL-field MRI the dominant
# instability at our ratio? That was the original intent and this geometry
# cannot do it. The most unstable wavelength is lambda ~ 6.5 v_A, so
#
#     by0 = 0.20  ->  lambda_tor = 1.3   needs Lz >~ 2
#     by0 = 0.50  ->  lambda_tor = 3.3   needs Lz >~ 5
#     by0 = 0.95  ->  lambda_tor = 6.2   needs Lz >~ 9
#
# against Lz = 1 here. The azimuthal mode does not fit in the box at ANY of our
# ratios, so it simply cannot grow, and its absence from these runs will be a
# property of the domain rather than of the physics. Checked before launching
# rather than discovered in the output.
#
# The proper test needs a tall box and goes to the cluster. At Q = 16 instead of
# the present 42, dz = 0.041, and Lz = 2 with Lx = Ly = 4 gives 98 x 98 x 49,
# about 1.8x this cost -- affordable for the low end of the ratio range. The
# high end needs Lz ~ 9 and roughly 33x, which needs chaining.
#
# ---------------------------------------------------------------------------
# by0 = 0 is not repeated: runs/pm_scan/pm04 is exactly that, same Re, Pm and q.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNOOPY="$ROOT/src/snoopy-v6.0-official/snoopy"
CFG="$ROOT/src/snoopy-v6.0-official/src/problem/mri/snoopy.cfg"

RE=1000
PM=4
BZ=0.1
BY_LIST=(0.2 0.5 0.95)      # ratios 2, 5, 9.5
T_FINAL=200.0
THREADS=3                    # 3 x 3 = 9 cores, under the 10-core cap
mkdir -p "$ROOT/runs/tor_scan"

# Timestep warning: incompressible MHD has no sound speed, so dt is set by the
# Alfven speed. The saturated turbulent field is b_rms ~ 0.88; adding by0 = 0.95
# in quadrature gives ~1.3, so the strongest run should take roughly 1.5x the
# 2.3 h a vertical-only box took. Watch it against the four-hour budget.
for by in "${BY_LIST[@]}"; do
    dir="$ROOT/runs/tor_scan/by${by/./p}"
    mkdir -p "$dir/data"
    sed -e "s/reynolds = 1000.0;/reynolds = ${RE}.0;/" \
        -e "s/reynolds_magnetic = 1000.0;/reynolds_magnetic = $((RE * PM)).0;/" \
        -e "s/t_final = 620.0;/t_final = ${T_FINAL};/" \
        -e "s/snapshot_step = 1.0;/snapshot_step = 1000.0;/" \
        -e "s/dump_step = 1.0;/dump_step = 25.0;/" \
        -e "s/bz0 = 0.1;/bz0 = ${BZ};/" \
        -e "s/by0 = 0.0;/by0 = ${by};/" \
        "$CFG" > "$dir/snoopy.cfg"
    ( cd "$dir" && OMP_NUM_THREADS=$THREADS "$SNOOPY" >run.log 2>&1 \
        && echo "done by0=$by" || echo "FAILED by0=$by" ) &
done

echo "toroidal scan: Re=$RE Pm=$PM bz0=$BZ by0=${BY_LIST[*]}, ${THREADS}t each"
wait
echo "toroidal scan complete"
