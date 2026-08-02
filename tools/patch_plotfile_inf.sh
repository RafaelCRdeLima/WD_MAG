#!/bin/bash
# Make plotfiles readable again by AMReX's own tools.
#
# The corrupt emag_density and etor_density derives overflow to inf, and the
# overflow is written into the per-FAB min/max block of Level_*/Cell_H. AMReX's
# Real parser aborts on that token -- "Read of Vector<Vector<Real>> failed" --
# before any data is read, so fextrema, fvolumesum and fbtbp cannot open the
# file at all. On the 192^3 run this hit all 70 plotfiles, including plt00000.
#
# The fix is to replace the inf/nan tokens with 0 in Cell_H. Only the min/max
# metadata is touched; every Cell_D data file is left alone, and nothing we
# integrate reads that metadata. On dir_rot96/plt00050 a single token had to be
# replaced, after which the density integrated to the correct mass.
#
# The original Cell_H is kept as Cell_H.orig, so this is reversible:
#     for f in plt*/Level_*/Cell_H.orig; do mv "$f" "${f%.orig}"; done
#
# Usage:  patch_plotfile_inf.sh plt?????

set -u

if [ "$#" -eq 0 ]; then
    echo "usage: $(basename "$0") plotfile [plotfile ...]" >&2
    exit 2
fi

patched=0
skipped=0

for plt in "$@"; do
    for hdr in "$plt"/Level_*/Cell_H; do
        [ -f "$hdr" ] || continue

        if ! grep -qiE 'inf|nan' "$hdr"; then
            skipped=$((skipped + 1))
            continue
        fi

        # Keep the first original only: re-running must not overwrite the
        # backup with an already-patched copy.
        [ -f "$hdr.orig" ] || cp -p "$hdr" "$hdr.orig"

        n=$(grep -oiE '[-+]?(inf|nan)' "$hdr" | wc -l)
        sed -i -E 's/[-+]?[iI][nN][fF]/0/g; s/[-+]?[nN][aA][nN]/0/g' "$hdr"
        echo "$plt: $n token(s) replaced"
        patched=$((patched + 1))
    done
done

echo "patched $patched header(s), $skipped already clean"
