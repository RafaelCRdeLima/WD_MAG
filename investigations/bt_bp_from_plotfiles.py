#!/usr/bin/env python3
"""Toroidal/poloidal magnetic energy split, rebuilt from B_x, B_y, B_z.

Why not read emag_density and etor_density straight from the plotfile: they
are corrupt. Both derives, and Castro's own Div_B, register several
face-centred components (Mag_Type_x/y/z) into one cell-centred FAB and then
read dat(i+1,...), which runs off the end of the data box. Measured on
dir_rot96/plt00000, before a single step:

    B_x, B_y      +-8.67e12       plausible, matches the model
    B_z           -4.1e7 .. 2.1e9 plausible
    Div_B         up to 3.1e146   should be ~0
    emag_density  up to 4.2e306   should be ~B^2/2 ~ 4e25
    etor_density  up to 3.3e44

The overflow reaches the plotfile metadata as `inf`, and AMReX's own
fextrema/fvolumesum then abort on the header before reading anything -- all
70 plotfiles of the 192^3 run are unreadable by them. yt parses `inf` as a
Python float and reads the files fine, which is why this is written in yt.
(The alternative, if yt is unavailable: sed the inf/nan tokens in
Level_0/Cell_H to 0 on a copy. One token per file was enough on the 96^3
run, and density then integrates to the correct mass.)

B_x, B_y and B_z themselves are written correctly at every time checked, so
the split is recoverable without re-running:

    B_tor   = (-y*B_x + x*B_y) / varpi
    B_pol^2 = B_x^2 + B_y^2 + B_z^2 - B_tor^2

This assumes the rotation axis is z, which it is for these models.

Restricted to cells with density above RHO_CUT so the ambient does not
dominate the volume integral. Override with the RHO_CUT environment
variable; the ambient is 2.0e4 g/cm^3.

    python3 bt_bp_from_plotfiles.py plt?????

Writes bt_bp.csv in the working directory.
"""
import os
import sys

import numpy as np
import yt

RHO_CUT = float(os.environ.get("RHO_CUT", 1.0e5))

yt.set_log_level("error")


def split(fn):
    ds = yt.load(fn)
    ad = ds.all_data()

    rho = np.asarray(ad["boxlib", "density"])
    sel = rho > RHO_CUT
    if not sel.any():
        raise RuntimeError(f"no cell above rho = {RHO_CUT:g}")

    bx = np.asarray(ad["boxlib", "B_x"])[sel]
    by = np.asarray(ad["boxlib", "B_y"])[sel]
    bz = np.asarray(ad["boxlib", "B_z"])[sel]
    x = np.asarray(ad["index", "x"])[sel]
    y = np.asarray(ad["index", "y"])[sel]
    dv = np.asarray(ad["index", "cell_volume"])[sel]

    varpi = np.hypot(x, y)
    safe = np.where(varpi > 0.0, varpi, 1.0)
    btor = np.where(varpi > 0.0, (-y * bx + x * by) / safe, 0.0)
    bpol2 = np.maximum(bx * bx + by * by + bz * bz - btor * btor, 0.0)

    e_tor = 0.5 * float(np.sum(btor * btor * dv))
    e_pol = 0.5 * float(np.sum(bpol2 * dv))
    return float(ds.current_time), e_tor, e_pol, float(rho[sel].max())


def main(paths):
    print(f"# RHO_CUT = {RHO_CUT:g} g/cm^3")
    print(f"# {'t':>9} {'E_tor':>12} {'E_pol':>12} {'Et/Ep':>12} {'Et/Emag':>10} {'rho_max':>11}")
    rows = []
    for fn in paths:
        try:
            t, e_tor, e_pol, rho_max = split(fn)
        except Exception as exc:                       # a truncated plotfile must not stop the sweep
            print(f"# {fn}: FAILED -- {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
            continue
        e_mag = e_tor + e_pol
        ratio = e_tor / e_pol if e_pol > 0 else float("inf")
        frac = e_tor / e_mag if e_mag > 0 else float("nan")
        print(f"{t:11.5f} {e_tor:12.5e} {e_pol:12.5e} {ratio:12.4e} {frac:10.7f} {rho_max:11.4e}",
              flush=True)
        rows.append((t, e_tor, e_pol, ratio, frac, rho_max))

    if not rows:
        return 1
    rows.sort()
    with open("bt_bp.csv", "w") as fh:
        fh.write("t,E_tor,E_pol,Et_over_Ep,Et_over_Emag,rho_max\n")
        for row in rows:
            fh.write(",".join(f"{v:.8e}" for v in row) + "\n")
    print(f"# {len(rows)} of {len(paths)} plotfiles -> bt_bp.csv", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
