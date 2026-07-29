"""Cache the 2D slices and 3D isosurfaces of one relaxed field.

Run:
    scf/.venv/bin/python3 papers/wd-toroidal-poloidal/figures/extract_fields.py

Target: seed 20 at t/t_dyn = 0.827 (plt_seed20_caeade00070), a plotfile
that falls INSIDE the validity window and whose toroidal fraction,
0.3384, is the closest in the batch to the ensemble mean -- a
representative realization, not a hand-picked one.

Everything is read through braithwaite_app/core/field_reader.py (R1: the
plotfile reader, the gauss conversion and the poloidal/toroidal split
already live there). Two properties of that module matter here:

  * Bt and Bp are built from the CELL-CENTERED B_x/B_y/B_z, not from
    Castro's derived emag_density/etor_density, so they are free of the
    ghost-cell gap of docs/teoria.md Sec 6.5. We verified the two routes
    agree anyway: cell by cell inside the star the two toroidal
    fractions differ by ~1e-17 (machine precision), and the
    volume-integrated ratio comes out 0.3175 against 0.3177.
  * Bt is B_phi about the z axis. This star is non-rotating and the
    field is random, so z is a bookkeeping axis, not a physically
    preferred one -- the same caveat the app carries.

ds.force_periodicity() is required for yt's isosurface extraction, which
otherwise refuses to fill ghost cells on a non-periodic domain. It wraps
only the outermost cell layer, which is vacuum here; every isosurface
extracted lies well inside the star.

The plotfiles live under the (gitignored) castro tree, so this cannot be
re-run from a fresh clone -- hence its output is committed.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "braithwaite_app"))

import yt                                          # noqa: E402
from core import field_reader as fr                 # noqa: E402

PLOTFILE = (REPO / "castro" / "Exec" / "science" / "wd_braithwaite"
            / "plt_seed20_caeade00070")
SEED = 20
T_DYN_S = 0.2758062098
RHO_FLOOR_FRAC = 1e-3      # the module's own stellar-surface convention
SLICE_RES = 256

OUT = HERE / "fields3d.npz"

FIELDS = [("Bmag", "Bmag"), ("Bt", "Bt (toroidal, B_phi)"),
          ("Bp", "Bp (poloidal)")]


def main():
    store = {}

    # --- 2D: meridional (x-z) slice, i.e. normal to y -------------------
    for key, spec_name in FIELDS:
        s = fr.load_slice(PLOTFILE, spec_name, "y", resolution=SLICE_RES)
        arr = np.asarray(s["array"])
        store[f"slice_{key}"] = np.abs(arr)
        if "slice_density" not in store:
            store["slice_density"] = np.asarray(s["density_array"])
            store["interior_mask"] = np.asarray(s["interior_mask"])
            store["extent_cm"] = np.asarray(s["extent_cm"])
        print(f"  slice {key}: {arr.shape}, |.| max {np.abs(arr).max():.3e} G")

    # --- 3D: isosurfaces at one common threshold ------------------------
    ds = yt.load(str(PLOTFILE))
    fr._add_derived_fields(ds)
    ds.force_periodicity()
    dims = ds.domain_dimensions
    cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=dims)
    rho = np.asarray(cg["boxlib", "density"])
    inside = rho > rho.max() * RHO_FLOOR_FRAC
    bmag = np.asarray(cg["gas", "Bmag"])
    from units import castro_to_gauss               # noqa: E402
    bmag_g = castro_to_gauss(bmag)

    # One threshold for all three fields, so the enclosed volumes are
    # directly comparable: the 90th percentile of |B| inside the star.
    iso_g = float(np.percentile(bmag_g[inside], 90))
    iso_code = iso_g / (bmag_g.max() / bmag.max())   # back to code units
    print(f"  isosurface level: {iso_g:.3e} G (|B| 90th pct inside star)")

    ad = ds.all_data()
    for key, _ in FIELDS:
        field = ("gas", key) if key != "Bmag" else ("gas", "Bmag")
        if key == "Bt":
            # Bt is signed; the surface is |Bt| = iso, which yt cannot ask
            # for directly, so register the magnitude once here.
            def _abs_bt(field, data):
                return np.abs(np.asarray(data["gas", "Bt"]))
            ds.add_field(("gas", "absBt"), function=_abs_bt,
                         sampling_type="cell", units="auto",
                         dimensions="dimensionless", take_log=False)
            field = ("gas", "absBt")
        surf = ds.surface(ad, field, iso_code)
        tri = np.asarray(surf.triangles)
        store[f"tri_{key}"] = tri
        print(f"  isosurface {key}: {len(tri)} triangles")

    star = ds.surface(ad, ("boxlib", "density"), rho.max() * RHO_FLOOR_FRAC)
    store["tri_star"] = np.asarray(star.triangles)
    print(f"  isosurface star: {len(store['tri_star'])} triangles")

    store["meta"] = json.dumps({
        "plotfile": PLOTFILE.name, "seed": SEED,
        "time_s": float(ds.current_time.to("s").value),
        "t_ttdyn": float(ds.current_time.to("s").value) / T_DYN_S,
        "iso_gauss": iso_g,
        "bmag_max_gauss": float(bmag_g[inside].max()),
        "rho_peak": float(rho.max()), "rho_floor_frac": RHO_FLOOR_FRAC,
        "dims": [int(d) for d in dims],
        "domain_left_cm": [float(v) for v in ds.domain_left_edge.to("cm").value],
        "domain_right_cm": [float(v) for v in ds.domain_right_edge.to("cm").value],
    })
    np.savez_compressed(OUT, **store)
    print("wrote", OUT.name, f"({OUT.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
