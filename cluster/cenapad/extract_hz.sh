#!/usr/bin/env bash
#PBS -N hzextract
#PBS -q testes
#PBS -l nodes=1:ppn=4
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -o hzextract.out
#
# Reduce the campaign HZ runs to two small CSVs that can actually come home.
#
# The PBS directives above are comments to bash, so this runs either way:
#
#     bash ~/WD_MAG/cluster/cenapad/extract_hz.sh      # on lovelace, directly
#     qsub  ~/WD_MAG/cluster/cenapad/extract_hz.sh     # if the sweep is long
#
# testes caps at 4 ncpus and 1 h (DIARIO/ONBOARDING, queue limits read off
# qstat -Qf) -- asking for the usual ppn=48 there fails with a message that
# does not name the offending resource. Four is what fits, and the sweep is
# I/O bound anyway.
#
# WHY THIS EXISTS
#
# A 192^3 plotfile is a few hundred MB and each run writes of order a hundred
# of them. Bringing them down over the two scp hops the separate frontend and
# lovelace filesystems force is not a transfer, it is an afternoon. Every
# earlier campaign reduced on the cluster and moved the CSV -- bt_bp_192.csv is
# 14 kB and answered the question the 70 plotfiles were written for.
#
# WHAT COMES OUT
#
#   hz_results/thermal_hz192.csv      fthermal over dir_hz192
#   hz_results/thermal_hz192ctl.csv   fthermal over dir_hz192ctl   <- the control
#   hz_results/field_hz192.csv        fbtbp   over dir_hz192  (E_mag, rho, radii)
#   hz_results/status.txt             how each run actually ended
#   hz_results.tgz                    all of the above, a few hundred kB
#
# The control is the point (DIARIO 10.1): field_scale = 0 gives the identical
# star with no field, so it is IN equilibrium. Quiet at ~1e7 K means the
# heating in dir_hz192 is the field's and is a measurement; the same heating
# without a field means it is numerics and campaign HZ has answered its own
# question in the negative.
set -u

# ----------------------------------------------------------------------------
# Paths. Overridable, because the layout on lovelace is not the layout on the
# workstation: the repo is ~/WD_MAG and Castro is a separate ~/wd-mag/Castro,
# with different capitalisation on both halves. Assuming otherwise is the
# phantom-directory bug that cost three queue windows (DIARIO 10).
# ----------------------------------------------------------------------------
EXEC="${EXEC:-$HOME/wd-mag/Castro/Exec/science/wd_scf_stability}"
AMREX="${AMREX:-$HOME/wd-mag/Castro/external/amrex}"
REPO="${REPO:-$HOME/WD_MAG}"
OUT="${OUT:-$EXEC/hz_results}"

echo "exec:  $EXEC"
echo "amrex: $AMREX"
echo "repo:  $REPO"
echo "out:   $OUT"

# ----------------------------------------------------------------------------
# Gates, before anything expensive. A gate costs seconds; a wasted hour in the
# queue costs an hour. VERIFY AT THE DESTINATION, never at the source.
# ----------------------------------------------------------------------------
fail=0
for d in "$EXEC" "$AMREX/Tools/Plotfile" "$REPO/tools"; do
    if [ ! -d "$d" ]; then echo "MISSING: $d"; fail=1; fi
done
for f in "$REPO/tools/fthermal.cpp" "$REPO/tools/fbtbp.cpp" "$REPO/tools/patch_plotfile_inf.sh"; do
    if [ ! -f "$f" ]; then echo "MISSING: $f  (git pull in $REPO?)"; fail=1; fi
done
[ "$fail" -eq 0 ] || { echo "gates failed, nothing done"; exit 1; }

mkdir -p "$OUT"

# ----------------------------------------------------------------------------
# 1. The headline, first and cheap: how did each run actually end?
#
# "The runs finished" has several meanings and they are not equivalent. The
# chain stops on reaching stop_time, on CHAIN_MAX submissions, on a window
# that took no steps, or on an abort -- and only the first is a finished run.
# The .out files name which. Reading this before the sweep means a truncated
# run is known before an hour is spent reducing it.
# ----------------------------------------------------------------------------
{
    echo "=== campaign HZ, run status ==="
    echo "extracted: $(date)"
    for tag in hz192 hz192ctl; do
        d="$EXEC/dir_$tag"
        echo
        echo "--- dir_$tag ---"
        if [ ! -d "$d" ]; then echo "  no such directory"; continue; fi
        log="$d/run_$tag.log"
        if [ -f "$log" ]; then
            # No "|| echo 0" on a grep -c: it already prints 0, and it exits 1
            # when it finds nothing, so the fallback fired ON TOP of the count
            # and the line came out as "aborts:   0" followed by a stray
            # "0 (cumulative over all windows)". Cosmetic, but this block is
            # read to decide whether a run is trustworthy and it should not
            # look like it is reporting two different numbers.
            echo "  steps:    $(grep -ac '^STEP' "$log" 2>/dev/null)"
            echo "  last:     $(grep -a '^STEP' "$log" 2>/dev/null | tail -1)"
            echo "  aborts:   $(grep -ac 'amrex::Abort' "$log" 2>/dev/null) (cumulative over all windows)"
            echo "  last abort: $(grep -a 'amrex::Abort' "$log" 2>/dev/null | tail -1)"
        else
            echo "  no $log"
        fi
        echo "  chain_count: $(cat "$d/chain_count" 2>/dev/null || echo '(none)')"
        echo "  plotfiles:   $(ls -d "$d"/plt????? 2>/dev/null | wc -l)"
        echo "  checkpoints: $(ls -d "$d"/chk????? 2>/dev/null | wc -l)"
        echo "  last chk:    $(ls -d "$d"/chk????? 2>/dev/null | sort | tail -1)"
        # Why the chain stopped, in its own words.
        echo "  chain says:"
        grep -ah '^chain:\|^target:' "$EXEC"/wd$tag.out "$EXEC"/wd${tag}.o* 2>/dev/null | tail -6 | sed 's/^/    /'
    done
} | tee "$OUT/status.txt"

# ----------------------------------------------------------------------------
# 2. Build the two diagnostics into AMReX's own Tools/Plotfile, which is where
#    finterior and fline were built (castro_problems/wd_braithwaite/
#    CASTRO_CORE_PATCHES.md). Serial: the sweep is I/O bound and Tools/Plotfile
#    defaults to USE_MPI=FALSE.
# ----------------------------------------------------------------------------
module purge 2>/dev/null || true
module load openmpi/5.0.8-gcc-15.2.0 2>/dev/null || true
if ! command -v python >/dev/null 2>&1 && command -v python3 >/dev/null 2>&1; then
    mkdir -p "$HOME/bin"; ln -sf "$(command -v python3)" "$HOME/bin/python"; export PATH="$HOME/bin:$PATH"
fi

cp "$REPO/tools/fthermal.cpp" "$REPO/tools/fbtbp.cpp" "$AMREX/Tools/Plotfile/"
( cd "$AMREX/Tools/Plotfile" && make programs=fthermal -j4 2>&1 | tail -5 )
( cd "$AMREX/Tools/Plotfile" && make programs=fbtbp   -j4 2>&1 | tail -5 )

# Named by dimension and compiler, so find rather than guess.
FTHERMAL="$(ls "$AMREX"/Tools/Plotfile/fthermal*.ex 2>/dev/null | head -1)"
FBTBP="$(ls "$AMREX"/Tools/Plotfile/fbtbp*.ex 2>/dev/null | head -1)"
echo "fthermal: ${FTHERMAL:-NOT BUILT}"
echo "fbtbp:    ${FBTBP:-NOT BUILT}"
[ -n "$FTHERMAL" ] || { echo "fthermal did not build; stopping before the sweep"; exit 1; }

# ----------------------------------------------------------------------------
# 3. Patch the plotfile headers.
#
# emag_density and etor_density are corrupt derives -- they read off the end of
# the data box -- and the overflow reaches the per-FAB min/max block of
# Level_*/Cell_H as an `inf`. AMReX's Real parser aborts on that token before
# reading any data, so every plotfile tool fails in the header. On the 192^3
# ztwd run this hit all 70 plotfiles. Only metadata is touched, Cell_D is left
# alone, and Cell_H.orig is kept.
#
# hz192 still writes those derives (inputs.hz192 line 138), so expect hits.
# ----------------------------------------------------------------------------
for tag in hz192 hz192ctl; do
    d="$EXEC/dir_$tag"
    [ -d "$d" ] || continue
    ( cd "$d" && bash "$REPO/tools/patch_plotfile_inf.sh" plt????? 2>&1 | tail -3 )
done

# ----------------------------------------------------------------------------
# 4. The sweeps. Whitespace tables out of the tools, converted to CSV here so
#    the comment lines survive as comments -- every earlier CSV in
#    investigations/ carries its provenance in a leading # block and
#    plot_*.py reads them with that assumption.
# ----------------------------------------------------------------------------
to_csv () {  # stdin: fixed-width table with # comments -> stdout: CSV
    # AMReX prints its own banner to STDOUT, not stderr:
    #
    #     Initializing AMReX (525c31011b65)...
    #     AMReX (525c31011b65) initialized
    #     AMReX (525c31011b65) finalized
    #
    # Those are not comments, so they came through as three data rows and a
    # sweep of ten plotfiles reported thirteen. It survived local testing only
    # because the test pipeline had a `grep -v "AMReX ("` on it that this
    # script did not, which is the whole reason a diagnostic has to be run
    # through the SAME path it will be used through.
    grep -v -e '^Initializing AMReX' -e '^AMReX (' \
    | awk '/^#/ { print; next }
           NF   { s=""; for (i=1;i<=NF;i++) s = s (i>1 ? "," : "") $i; print s }'
}

# A sweep is worth exactly as many rows as it read plotfiles. Anything less
# means it was cut short, and a short CSV is the dangerous failure here: it
# still parses, still plots, and quietly describes a window that is not the
# one the run covers.
#
# This exists because it happened. The first extraction was interrupted during
# the hz192ctl sweep and the CSVs came home at ZERO bytes -- awk block-buffers
# at 4 kB, so a file that is being written looks identical to a file that
# failed. They were packed and shipped anyway, and the emptiness was only
# noticed on the far side of two scp hops.
check_rows () {  # file, expected
    local f="$1" want="$2" got
    got=$(grep -vc '^#' "$f" 2>/dev/null || echo 0)
    if [ "$got" -eq "$want" ]; then
        echo "  -> $f ($got rows)"
    else
        echo "  -> $f (GOT $got ROWS, EXPECTED $want) -- INCOMPLETE"
        INCOMPLETE=$((INCOMPLETE + 1))
    fi
}
INCOMPLETE=0

for tag in hz192 hz192ctl; do
    d="$EXEC/dir_$tag"
    [ -d "$d" ] || { echo "skip $tag: no $d"; continue; }
    n=$(ls -d "$d"/plt????? 2>/dev/null | wc -l)
    [ "$n" -gt 0 ] || { echo "skip $tag: no plotfiles"; continue; }
    echo "sweeping $tag: $n plotfiles"

    ( cd "$d" && "$FTHERMAL" plt????? ) \
        | sed "1i # fthermal over dir_$tag, $n plotfiles. Campaign HZ thermal readout." \
        | to_csv > "$OUT/thermal_$tag.csv"
    check_rows "$OUT/thermal_$tag.csv" "$n"
done

# The field only where there is one. In the control it is identically zero by
# construction, and a column of zeros is not a measurement -- but E_mag in
# dir_hz192 is what the thermal budget has to be compared against, so it is
# needed here and not only in the with-field run's own history.
d="$EXEC/dir_hz192"
if [ -d "$d" ] && [ "$(ls -d "$d"/plt????? 2>/dev/null | wc -l)" -gt 0 ] && [ -n "$FBTBP" ]; then
    nf=$(ls -d "$d"/plt????? 2>/dev/null | wc -l)
    echo "sweeping hz192 field: $nf plotfiles"
    ( cd "$d" && "$FBTBP" plt????? ) \
        | sed '1i # fbtbp over dir_hz192. E_mag, mass, radii and rotation.' \
        | to_csv > "$OUT/field_hz192.csv"
    check_rows "$OUT/field_hz192.csv" "$nf"
fi

# ----------------------------------------------------------------------------
# 5. Package -- but only if there is nothing incomplete to package.
#
# Shipping a truncated CSV is worse than shipping nothing: it arrives looking
# like a result. Refusing here costs a re-run of a sweep that takes minutes;
# not refusing cost two scp hops and a conclusion drawn from a file that had
# no rows in it.
# ----------------------------------------------------------------------------
if [ "$INCOMPLETE" -gt 0 ]; then
    echo
    echo "REFUSING TO PACKAGE: $INCOMPLETE sweep(s) came out short."
    echo "The CSVs are in $OUT and can be inspected, but hz_results.tgz was NOT"
    echo "written, so a partial extraction cannot be mistaken for a finished one."
    echo "Re-run this script; the header patching and the builds are idempotent."
    ls -la "$OUT"
    exit 1
fi

tar czf "$EXEC/hz_results.tgz" -C "$(dirname "$OUT")" "$(basename "$OUT")"
echo
echo "done: $EXEC/hz_results.tgz  ($(du -h "$EXEC/hz_results.tgz" | cut -f1))"
ls -la "$OUT"
echo
echo "bring it home, two hops -- the frontend and lovelace filesystems are separate:"
echo "  lovelace\$  scp $EXEC/hz_results.tgz frontend:"
echo "  laptop\$    scp -P 31459 rcrlima@cenapad.unicamp.br:hz_results.tgz ~/"
