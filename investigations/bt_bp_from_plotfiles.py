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
CASTRO_TO_GAUSS = np.sqrt(4.0 * np.pi)
B_CRIT = 4.4140e13   # the field above which the unquantised ztwd EOS is out of range

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

    # Peak strengths as well as energies: the model manifest reports
    # Bt_over_Bp_amplitude, a ratio of maxima, and that is what it has to be
    # compared against. The energy ratio weights the whole volume and is a
    # different number.
    #
    # Reported in GAUSS. Castro's MHD state is in Heaviside-Lorentz units,
    # B' = B/sqrt(4 pi) -- see problem_initialize_mhd_data.H. Reading the raw
    # state as gauss understates every field by 3.5449. The energies need no
    # conversion: 0.5*B'^2 is already B^2/(8 pi) in erg/cm^3.
    b_tor_max = CASTRO_TO_GAUSS * float(np.abs(btor).max())
    b_pol_max = CASTRO_TO_GAUSS * float(np.sqrt(bpol2).max())
    b_max = CASTRO_TO_GAUSS * float(np.sqrt(bx * bx + by * by + bz * bz).max())

    return (float(ds.current_time), e_tor, e_pol, float(rho[sel].max()),
            b_tor_max, b_pol_max, b_max)


def main(paths):
    print(f"# RHO_CUT = {RHO_CUT:g} g/cm^3")
    print(f"# energies in erg; peak fields in GAUSS (state is Heaviside-Lorentz,"
          f" converted by sqrt(4 pi)); B_c = {B_CRIT:.3e} G")
    print(f"# {'t':>9} {'E_tor':>12} {'E_pol':>12} {'Et/Ep':>12} {'Et/Emag':>10}"
          f" {'Btor_max_G':>12} {'Bpol_max_G':>12} {'B_max_G':>12} {'Bt/Bp_amp':>12}"
          f" {'B/B_c':>10} {'rho_max':>11}")
    rows = []
    for fn in paths:
        try:
            t, e_tor, e_pol, rho_max, b_tor_max, b_pol_max, b_max = split(fn)
        except Exception as exc:                       # a truncated plotfile must not stop the sweep
            print(f"# {fn}: FAILED -- {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
            continue
        e_mag = e_tor + e_pol
        ratio = e_tor / e_pol if e_pol > 0 else float("inf")
        frac = e_tor / e_mag if e_mag > 0 else float("nan")
        amp = b_tor_max / b_pol_max if b_pol_max > 0 else float("inf")
        bbc = b_max / B_CRIT
        print(f"{t:11.5f} {e_tor:12.5e} {e_pol:12.5e} {ratio:12.4e} {frac:10.7f}"
              f" {b_tor_max:12.5e} {b_pol_max:12.5e} {b_max:12.5e} {amp:12.5e}"
              f" {bbc:10.5f} {rho_max:11.4e}",
              flush=True)
        rows.append((t, e_tor, e_pol, ratio, frac, b_tor_max, b_pol_max, b_max, amp, bbc, rho_max))

    if not rows:
        return 1
    rows.sort()
    with open("bt_bp.csv", "w") as fh:
        fh.write("t,E_tor,E_pol,Et_over_Ep,Et_over_Emag,"
                 "Btor_max_G,Bpol_max_G,B_max_G,Bt_over_Bp_amp,B_over_Bc,rho_max\n")
        for row in rows:
            fh.write(",".join(f"{v:.8e}" for v in row) + "\n")
    print(f"# {len(rows)} of {len(paths)} plotfiles -> bt_bp.csv", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
