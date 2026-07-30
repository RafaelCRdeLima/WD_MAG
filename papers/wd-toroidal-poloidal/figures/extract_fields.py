"""Cache the 2D slices and 3D field lines of one relaxed field.

Run:
    scf/.venv/bin/python3 papers/wd-toroidal-poloidal/figures/extract_fields.py

Target: seed 20 at t/t_dyn = 0.827 (plt_seed20_caeade00070), a plotfile
that falls INSIDE the validity window and whose toroidal fraction,
0.3384, is the closest in the batch to the ensemble mean -- a
representative realization, not a hand-picked one.

Plotfile reading, the gauss conversion and the poloidal/toroidal split
come from braithwaite_app/core/field_reader.py (R1). Two properties of
that module matter here:

  * Bt and Bp are built from the CELL-CENTERED B_x/B_y/B_z, not from
    Castro's derived emag_density/etor_density, so they are free of the
    ghost-cell gap of docs/teoria.md Sec 6.5. We verified the two routes
    agree anyway: cell by cell inside the star the two toroidal
    fractions differ by ~1e-17 (machine precision), and the
    volume-integrated ratio comes out 0.3175 against 0.3177.
  * Bt is B_phi about the z axis. This star is non-rotating and the
    field is random, so z is a bookkeeping axis, not a physically
    preferred one -- the same caveat the app carries.

Field lines are traced by streamlines.py (a plain RK4 on the trilinearly
interpolated grid) rather than by yt's Streamlines, which integrates on
the AMR hierarchy this single-level data does not need and does not
return the per-vertex magnitudes the figure colors by.

The plotfiles live under the (gitignored) castro tree, so this cannot be
re-run from a fresh clone -- hence its output is committed.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "braithwaite_app"))

import streamlines as sl                            # noqa: E402
from core import field_reader as fr                 # noqa: E402

PLOTFILE = (REPO / "castro" / "Exec" / "science" / "wd_braithwaite"
            / "plt_seed20_caeade00070")
SEED = 20
T_DYN_S = 0.2758062098
RHO_FLOOR_FRAC = 1e-3      # the module's own stellar-surface convention
SLICE_RES = 256

# Seeding: strongest-field cells, thinned so the lines do not all start
# in the same structure, with a fixed RNG seed so the figure is
# reproducible. Same seed points for all three panels, so the panels
# differ only by which field is being followed.
N_LINES = 8
SEED_PCTILE = 85
MIN_SEED_SEP_CM = 1.05e8
RNG_SEED = 7
MAX_STEPS = 250            # per direction, for the 3D view: long enough to
                           # show structure, short enough that a single line
                           # can still be followed by eye
WANDER_STEPS = 4000        # one long line per field, for the wander curve

OUT = HERE / "fields3d.npz"

FIELDS = [("Bmag", "Bmag"), ("Bt", "Bt (toroidal, B_phi)"),
          ("Bp", "Bp (poloidal)")]


def main():
    store = {}

    # --- 2D: meridional (x-z) slice, i.e. normal to y -------------------
    # yt returns a y-normal slice as (x, z) on (axis0, axis1); the figure
    # transposes so that z is vertical.
    for key, spec_name in FIELDS:
        s = fr.load_slice(PLOTFILE, spec_name, "y", resolution=SLICE_RES)
        arr = np.asarray(s["array"])
        store[f"slice_{key}"] = np.abs(arr)
        if "slice_density" not in store:
            store["slice_density"] = np.asarray(s["density_array"])
            store["extent_cm"] = np.asarray(s["extent_cm"])
        print(f"  slice {key}: {arr.shape}, |.| max {np.abs(arr).max():.3e} G")

    # --- 3D: field lines of the total, poloidal and toroidal fields -----
    vg = fr.load_vector_grid(PLOTFILE)
    dims = vg["dims"]
    origin = np.asarray(vg["origin_cm"])
    spacing = np.asarray(vg["spacing_cm"])
    rho = vg["density"]
    inside = rho > rho.max() * RHO_FLOOR_FRAC

    total, pol, tor = sl.decompose(vg["B_x"], vg["B_y"], vg["B_z"],
                                   origin, spacing, dims)
    r_star_pre = fr.estimate_star_radius(rho, tuple(spacing))
    vectors = {"Bmag": total, "Bp": pol, "Bt": tor}

    mag_total = np.linalg.norm(total, axis=-1)
    candidates = np.argwhere(
        inside & (mag_total > np.percentile(mag_total[inside], SEED_PCTILE)))
    rng = np.random.default_rng(RNG_SEED)
    rng.shuffle(candidates)
    seeds = []
    for c in candidates:
        p = origin + spacing * (c + 0.5)
        if all(np.linalg.norm(p - q) > MIN_SEED_SEP_CM for q in seeds):
            seeds.append(p)
        if len(seeds) >= N_LINES:
            break
    store["seeds_cm"] = np.asarray(seeds)
    print(f"  {len(seeds)} seed points")

    for key in ("Bmag", "Bp", "Bt"):
        field = vectors[key]
        lines = sl.trace(field, seeds, origin, spacing, dims, inside,
                         max_steps=MAX_STEPS)
        mags = [sl.sample_along(field, ln, origin, spacing, dims)
                for ln in lines]
        # Ragged, so stored flat with an offset table.
        store[f"line_{key}_xyz"] = np.concatenate(lines)
        store[f"line_{key}_mag"] = np.concatenate(mags)
        store[f"line_{key}_off"] = np.cumsum(
            [0] + [len(ln) for ln in lines])
        print(f"  lines {key}: {len(lines)}, "
              f"{np.mean([len(x) for x in lines]):.0f} points each")

        # One long line, for the quantitative wander curve. The 3D view
        # cannot show both an individual line's topology and how far it
        # travels -- at the length where the wandering becomes visible the
        # panel is a solid ball. So the distance-from-seed against
        # arc-length is measured here and plotted in 2D instead.
        longline = sl.trace(field, [seeds[0]], origin, spacing, dims, inside,
                            max_steps=WANDER_STEPS)
        if longline:
            ln = longline[0]
            # The seed sits where trace() inserted it, which is len(back),
            # NOT the midpoint -- the two branches stop independently, so
            # slicing at len(ln)//2 lands on the wrong point whenever they
            # differ in length.
            k = int(np.argmin(np.linalg.norm(ln - np.asarray(seeds[0]), axis=1)))
            fwd = ln[k:]
            seg = np.linalg.norm(np.diff(fwd, axis=0), axis=1)
            store[f"wander_{key}_arc"] = np.concatenate([[0.0], np.cumsum(seg)])
            store[f"wander_{key}_dist"] = np.linalg.norm(fwd - fwd[0], axis=1)
            print(f"    long line: {len(fwd)} pts, arc = "
                  f"{store[f'wander_{key}_arc'][-1] / r_star_pre:.1f} R_star")

    # Volume statistics quoted in the text, computed here so the numbers
    # and the figure cannot drift apart.
    iso = float(np.percentile(mag_total[inside], 90))
    stats = {}
    for key in ("Bmag", "Bp", "Bt"):
        m = np.linalg.norm(vectors[key], axis=-1)
        stats[key] = {
            "max_G": float(m[inside].max()),
            "frac_above_iso": float(((m > iso) & inside).sum() / inside.sum()),
        }
    mp = np.linalg.norm(pol, axis=-1)[inside]
    mt = np.linalg.norm(tor, axis=-1)[inside]
    stats["median_Bp_over_Bt"] = float(np.median(mp / np.maximum(mt, 1e-30)))
    stats["frac_volume_Bt_gt_Bp"] = float((mt > mp).sum() / inside.size
                                          * inside.size / inside.sum())
    print(f"  volume fractions above {iso:.3e} G: "
          + ", ".join(f"{k} {100 * stats[k]['frac_above_iso']:.2f}%"
                      for k in ("Bmag", "Bp", "Bt")))
    print(f"  median |Bp|/|Bt| = {stats['median_Bp_over_Bt']:.2f}; "
          f"|Bt|>|Bp| in {100 * stats['frac_volume_Bt_gt_Bp']:.1f}% of volume")

    r_star = fr.estimate_star_radius(rho, tuple(spacing))
    store["meta"] = json.dumps({
        "plotfile": PLOTFILE.name, "seed": SEED,
        "time_s": vg["time_s"], "t_ttdyn": vg["time_s"] / T_DYN_S,
        "iso_gauss": iso, "stats": stats,
        "r_star_cm": float(r_star),
        "rho_peak": float(rho.max()), "rho_floor_frac": RHO_FLOOR_FRAC,
        "dims": [int(d) for d in dims],
        "n_lines": len(seeds), "max_steps": MAX_STEPS,
    })
    np.savez_compressed(OUT, **store)
    print("wrote", OUT.name, f"({OUT.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
