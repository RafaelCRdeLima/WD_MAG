"""Spatial field visualization -- 2D slices and 3D vector data straight
from a plotfile, via `yt` (the standard AMReX/Castro plotfile reader; no
custom parser reimplemented here, same reuse-not-rebuild rule as
extraction.py's use of `finterior`). Everything here is read-only and
diagnostic-only -- it never feeds a physics number into
core/persistence.py (that stays on the finterior-verified path in
extraction.py); this module exists purely so a human can look at the
field's actual 2D/3D shape.

Field values are converted from Castro's internal B'=B/sqrt(4*pi)
convention to gauss (dashboard/units.py :: castro_to_gauss(), the same
single source of truth Tab 3 uses -- see docs/teoria.md Sec 1.10) so a
number shown here means the same gauss the rest of the project reports.

Poloidal/toroidal split: Castro's own derived fields `emag_density`
(total B^2/8pi) and `etor_density` (the B_phi-only piece, same
quantities `core/extraction.py` reads via finterior for the E_tor/E_mag
diagnostic) are enough to build a per-cell decomposition without
reimplementing the poloidal/toroidal projection ourselves --
epol_density = emag_density - etor_density, toroidal_frac =
etor_density / emag_density (0 = purely poloidal, 1 = purely toroidal
at that point).
"""

import re
import sys
from pathlib import Path

from core.paths import DASHBOARD_DIR

_dashboard_path = str(DASHBOARD_DIR)
if _dashboard_path not in sys.path:
    sys.path.insert(0, _dashboard_path)

import numpy as np  # noqa: E402
from units import castro_to_gauss  # noqa: E402  (dashboard/units.py)

import yt  # noqa: E402

yt.set_log_level(50)  # this app has its own status labels -- yt's own logging is noise here

_PLOTFILE_RE = re.compile(r"^plt_(.+?)(\d{5})$")

# One entry per selectable 2D-slice/3D-color field. `log`: whether the
# quantity is positive-definite and spans enough orders of magnitude to
# want a log color norm (energy densities, |B|, density) vs a signed or
# already-bounded quantity that must stay linear (B components,
# toroidal_frac) -- using log on a signed field would silently drop
# every negative value from the plot instead of erroring, so this is
# tracked explicitly per field rather than guessed from the data.
FIELD_SPECS = {
    "Bmag": {"yt": ("gas", "Bmag"), "unit": "G", "log": True, "diverging": False},
    "B_x": {"yt": ("boxlib", "B_x"), "unit": "G", "log": False, "diverging": True},
    "B_y": {"yt": ("boxlib", "B_y"), "unit": "G", "log": False, "diverging": True},
    "B_z": {"yt": ("boxlib", "B_z"), "unit": "G", "log": False, "diverging": True},
    "Bt (toroidal, B_phi)": {"yt": ("gas", "Bt"), "unit": "G", "log": False, "diverging": True},
    "Bp (poloidal)": {"yt": ("gas", "Bp"), "unit": "G", "log": True, "diverging": False},
    "density": {"yt": ("boxlib", "density"), "unit": "g/cm^3", "log": True, "diverging": False},
    "emag_density": {"yt": ("boxlib", "emag_density"), "unit": "erg/cm^3", "log": True, "diverging": False},
    "etor_density": {"yt": ("boxlib", "etor_density"), "unit": "erg/cm^3 (toroidal)", "log": True, "diverging": False},
    "epol_density": {"yt": ("gas", "epol_density"), "unit": "erg/cm^3 (poloidal)", "log": True, "diverging": False},
    "toroidal_frac": {"yt": ("gas", "toroidal_frac"), "unit": "E_tor/E_mag local", "log": False, "diverging": True, "vrange": (0.0, 1.0)},
}


def list_run_ids(run_dir: Path) -> list[str]:
    """Distinct run_ids with at least one real (non-backup) plotfile in
    run_dir, newest-activity-first is not attempted here -- just sorted
    for a stable dropdown order.
    """
    ids = set()
    for entry in Path(run_dir).iterdir():
        if not entry.is_dir() or ".old." in entry.name:
            continue
        m = _PLOTFILE_RE.match(entry.name)
        if m:
            ids.add(m.group(1))
    return sorted(ids)


def list_plotfiles_for_run(run_dir: Path, run_id: str) -> list[Path]:
    """Real plotfiles for one run_id, in step order. Excludes the
    `.old.<pid>` backup copies Castro leaves behind on a restart --
    without this filter, `\\d{5}$` can accidentally match the tail of a
    long pid suffix and produce a spurious run_id/step pair, since the
    regex only anchors to the end of the string, not to "no more
    characters at all."
    """
    out = []
    for entry in Path(run_dir).iterdir():
        if ".old." in entry.name:
            continue
        m = _PLOTFILE_RE.match(entry.name)
        if m and m.group(1) == run_id:
            out.append((int(m.group(2)), entry))
    return [p for _, p in sorted(out)]


def _add_derived_fields(ds):
    """Registers every derived field this module can show. Cheap and
    idempotent enough to call unconditionally per-load (yt.load() opens
    a fresh dataset object each time, so there's no stale-registration
    risk across calls)."""

    def _bmag(field, data):
        return (
            data["boxlib", "B_x"] ** 2
            + data["boxlib", "B_y"] ** 2
            + data["boxlib", "B_z"] ** 2
        ) ** 0.5

    ds.add_field(
        ("gas", "Bmag"), function=_bmag, sampling_type="cell",
        units="auto", dimensions="dimensionless", take_log=False,
    )

    def _epol(field, data):
        return data["boxlib", "emag_density"] - data["boxlib", "etor_density"]

    ds.add_field(
        ("gas", "epol_density"), function=_epol, sampling_type="cell",
        units="auto", dimensions="dimensionless", take_log=False,
    )

    def _tfrac(field, data):
        emag = data["boxlib", "emag_density"]
        etor = data["boxlib", "etor_density"]
        ratio = etor / emag
        ratio[~np.isfinite(ratio)] = 0.0  # emag~0 in the field-free exterior/vacuum -- not toroidal, not poloidal, just absent
        return ratio

    ds.add_field(
        ("gas", "toroidal_frac"), function=_tfrac, sampling_type="cell",
        units="auto", dimensions="dimensionless", take_log=False,
    )

    # Bt = B_phi (azimuthal component about the z axis -- same axis
    # convention as estimate_axis_radii's R_pol, not a claim of real
    # rotational symmetry for this random field), Bp = the rest of |B|.
    # Coordinates/components are stripped of units (.d) before the numpy
    # ops -- unyt's np.where refuses to mix a unyt_array with a bare
    # float literal (code_length vs dimensionless), so the "guard
    # against varpi=0 on the z-axis itself" step has to happen on plain
    # ndarrays.
    def _bt(field, data):
        x = data["index", "x"].to("cm").d
        y = data["index", "y"].to("cm").d
        bx = data["boxlib", "B_x"].d
        by = data["boxlib", "B_y"].d
        varpi = np.sqrt(x**2 + y**2)
        varpi_safe = np.where(varpi > 0, varpi, 1.0)
        bt = (-bx * y + by * x) / varpi_safe
        return np.where(varpi > 0, bt, 0.0)

    ds.add_field(
        ("gas", "Bt"), function=_bt, sampling_type="cell",
        units="auto", dimensions="dimensionless", take_log=False,
    )

    def _bp(field, data):
        bmag2 = (
            data["boxlib", "B_x"].d ** 2
            + data["boxlib", "B_y"].d ** 2
            + data["boxlib", "B_z"].d ** 2
        )
        # np.asarray, not .d: unlike the native "boxlib" fields above,
        # this reads back our OWN "gas","Bt" derived field, and yt does
        # not reliably hand that back as a unyt_array (depends on the
        # field-detection pass vs. the real evaluation pass) -- asarray
        # strips units if present and is a no-op if it's already a
        # plain ndarray, so it's correct either way.
        bt2 = np.asarray(data["gas", "Bt"]) ** 2
        return np.sqrt(np.clip(bmag2 - bt2, 0.0, None))

    ds.add_field(
        ("gas", "Bp"), function=_bp, sampling_type="cell",
        units="auto", dimensions="dimensionless", take_log=False,
    )


def load_slice(plotfile: Path, field: str, axis: str, resolution: int = 256) -> dict:
    """A 2D fixed-resolution-buffer slice through the domain center.
    Always also returns the density slice (same frb, negligible extra
    cost) so the caller can draw a stellar-surface contour regardless of
    which field is being displayed -- otherwise there is no way to tell,
    from the field plot alone, what part of the image is inside the star
    versus the surrounding field-free/vacuum-floor exterior.
    """
    ds = yt.load(str(plotfile))
    _add_derived_fields(ds)
    spec = FIELD_SPECS[field]

    slc = yt.SlicePlot(ds, axis, spec["yt"])
    width_cm = float((ds.domain_right_edge - ds.domain_left_edge)[0].to("cm").value)
    frb = slc.data_source.to_frb((width_cm, "cm"), resolution)
    arr = np.asarray(frb[spec["yt"]])
    if field.startswith("B"):
        arr = castro_to_gauss(arr)
    density_arr = np.asarray(frb[("boxlib", "density")])

    rho_peak = float(np.nanmax(density_arr))
    interior_mask = density_arr > rho_peak * 1e-3 if rho_peak > 0 else np.zeros_like(density_arr, dtype=bool)

    if field == "toroidal_frac":
        # emag_density ~ 0 in the field-free exterior means the ratio is
        # 0/0-noise there (numerically anything from 0 to 1, not a real
        # ratio) -- mask it out rather than show meaningless static.
        emag_arr = np.asarray(frb[("boxlib", "emag_density")])
        floor = emag_arr.max() * 1e-6
        arr = np.where(emag_arr > floor, arr, np.nan)
        interior_vrange = spec["vrange"]
    else:
        # Scale the color range to what the field actually does INSIDE
        # the star, not the full domain -- the field-free exterior spans
        # many more orders of magnitude down toward zero than the
        # interior does, and letting it set vmin turns the whole interior
        # into a single saturated color (the same "irrelevant region
        # contaminates the view" failure the density-boundary contour
        # already fixed for "where is the star", now fixed for "what
        # color range matters").
        interior_vals = arr[interior_mask]
        interior_vals = interior_vals[np.isfinite(interior_vals)]
        if interior_vals.size:
            if spec["diverging"]:
                vmax = float(np.percentile(np.abs(interior_vals), 99)) or 1.0
                interior_vrange = (-vmax, vmax)
            else:
                lo = float(np.percentile(interior_vals, 1))
                hi = float(np.percentile(interior_vals, 99))
                if lo <= 0:  # log norm needs a strictly positive floor
                    positive = interior_vals[interior_vals > 0]
                    lo = float(np.percentile(positive, 1)) if positive.size else hi * 1e-6
                if lo >= hi:
                    hi = lo * 10 + 1e-300
                interior_vrange = (lo, hi)
        else:
            interior_vrange = spec.get("vrange")

    return {
        "array": arr,
        "density_array": density_arr,
        "interior_mask": interior_mask,
        "interior_vrange": interior_vrange,
        "extent_cm": (-width_cm / 2, width_cm / 2, -width_cm / 2, width_cm / 2),
        "time_s": float(ds.current_time.to("s").value),
        "field": field,
        "axis": axis,
        "unit": spec["unit"],
        "log": spec["log"],
        "diverging": spec["diverging"],
    }


def load_vector_grid(plotfile: Path) -> dict:
    """The full B (gauss), density and toroidal-fraction grid at native
    plotfile resolution, as plain numpy arrays -- input for the 3D
    streamline/isosurface view (ui/field_view.py), which builds the
    pyvista mesh itself so this module stays yt-only / pyvista-free.
    """
    ds = yt.load(str(plotfile))
    dims = ds.domain_dimensions
    cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=dims)

    bx = castro_to_gauss(np.asarray(cg["boxlib", "B_x"]))
    by = castro_to_gauss(np.asarray(cg["boxlib", "B_y"]))
    bz = castro_to_gauss(np.asarray(cg["boxlib", "B_z"]))
    density = np.asarray(cg["boxlib", "density"])
    emag = np.asarray(cg["boxlib", "emag_density"])
    etor = np.asarray(cg["boxlib", "etor_density"])
    with np.errstate(divide="ignore", invalid="ignore"):
        toroidal_frac = np.where(emag > 0, etor / emag, 0.0)

    left = np.asarray(ds.domain_left_edge.to("cm").value)
    right = np.asarray(ds.domain_right_edge.to("cm").value)
    spacing = (right - left) / np.asarray(dims)

    return {
        "B_x": bx, "B_y": by, "B_z": bz, "density": density,
        "toroidal_frac": toroidal_frac,
        "dims": tuple(int(d) for d in dims),
        "origin_cm": tuple(left), "spacing_cm": tuple(spacing),
        "time_s": float(ds.current_time.to("s").value),
    }


def estimate_star_radius(density: np.ndarray, spacing_cm: tuple) -> float:
    """Equivalent-volume radius of the density>0.1%-of-peak cells -- the
    same boundary convention used everywhere else in this module (the 2D
    contour, the 3D density envelope), reduced to one scalar. Public
    because it's a shared reference length: ui/field_view.py uses it to
    place 3D streamline seeds inside the actual star (not the domain),
    and core/dipole.py uses it to choose measurement radii outside the
    star -- both need the exact same boundary definition or their
    results wouldn't be comparable to each other.
    """
    rho_max = float(density.max())
    if rho_max <= 0:
        return 0.0
    cell_volume = spacing_cm[0] * spacing_cm[1] * spacing_cm[2]
    n_inside = int(np.count_nonzero(density > rho_max * 1e-3))
    volume = n_inside * cell_volume
    return (3 * volume / (4 * np.pi)) ** (1 / 3)


def estimate_axis_radii(density: np.ndarray, dims: tuple, spacing_cm: tuple, origin_cm: tuple) -> dict:
    """R_eq and R_pol (cm) from the same 0.1%-of-peak density boundary
    used everywhere else in this module (estimate_star_radius, the 2D
    contour, the 3D envelope) -- but measured as a length along a
    direction, not folded into a single equivalent-volume scalar, so an
    oblate/prolate deformation actually shows up as R_pol != R_eq
    instead of cancelling out in a volume average. z is the reference
    axis (R_pol along +/-z, R_eq averaged over +/-x and +/-y in the
    z-through-center plane) -- the same axis convention `ui/field_view.py`
    already exposes as the default "corte" (2D slice) direction, not a
    claim that this random Braithwaite field has any real rotational
    symmetry about z. Each ray is walked outward from the domain's own
    center cell and the crossing is linearly interpolated between the
    last cell above threshold and the first below it, for sub-cell
    accuracy (the exact-zero-density-clip failure mode that forced
    H-based surface detection on the SCF side, docs/teoria.md Sec 1.13,
    does not apply here -- Castro's density floor is `castro.small_dens`,
    not an exact clip, so it varies continuously through the boundary).
    """
    rho_max = float(density.max())
    if rho_max <= 0:
        return {"r_eq_cm": 0.0, "r_pol_cm": 0.0}
    threshold = rho_max * 1e-3

    nx, ny, nz = dims
    ox, oy, oz = origin_cm
    sx, sy, sz = spacing_cm
    xs = ox + sx * (np.arange(nx) + 0.5)
    ys = oy + sy * (np.arange(ny) + 0.5)
    zs = oz + sz * (np.arange(nz) + 0.5)
    ix0 = int(np.argmin(np.abs(xs)))
    iy0 = int(np.argmin(np.abs(ys)))
    iz0 = int(np.argmin(np.abs(zs)))

    def ray_radius(profile: np.ndarray, coords: np.ndarray, i0: int) -> float:
        radii = []
        for step in (1, -1):
            i, prev_i = i0, i0
            while 0 <= i < len(profile) and profile[i] > threshold:
                prev_i = i
                i += step
            if 0 <= i < len(profile):
                d0, d1 = float(profile[prev_i]), float(profile[i])
                r0, r1 = float(coords[prev_i]), float(coords[i])
                frac = (threshold - d0) / (d1 - d0) if d1 != d0 else 0.0
                radii.append(abs(r0 + frac * (r1 - r0)))
            # else: this ray never dropped below threshold before the
            # domain edge -- excluded rather than reporting a domain-clipped radius
        return float(np.mean(radii)) if radii else float("nan")

    r_pol = ray_radius(density[ix0, iy0, :], zs, iz0)
    r_eq_x = ray_radius(density[:, iy0, iz0], xs, ix0)
    r_eq_y = ray_radius(density[ix0, :, iz0], ys, iy0)
    r_eq = float(np.nanmean([r_eq_x, r_eq_y]))

    return {"r_eq_cm": r_eq, "r_pol_cm": r_pol}


def estimate_bt_bp_max(b_x: np.ndarray, b_y: np.ndarray, b_z: np.ndarray,
                        dims: tuple, origin_cm: tuple, spacing_cm: tuple) -> dict:
    """Peak |Bt| and peak Bp (poloidal magnitude) over the whole grid, in
    gauss -- b_x/b_y/b_z are expected already in gauss (load_vector_grid's
    output, not the raw Castro-internal units). Same z-axis convention
    (Bt = B_phi about z) as _add_derived_fields' 2D "Bt"/"Bp" and
    estimate_axis_radii's R_pol -- this is the 3D-grid, numpy-only
    equivalent of that yt derived field, used here because
    load_vector_grid doesn't route through yt derived fields.
    """
    nx, ny, nz = dims
    ox, oy, _oz = origin_cm
    sx, sy, _sz = spacing_cm
    xs = ox + sx * (np.arange(nx) + 0.5)
    ys = oy + sy * (np.arange(ny) + 0.5)
    x, y = np.meshgrid(xs, ys, indexing="ij")
    varpi = np.sqrt(x**2 + y**2)[:, :, None]  # broadcast over z
    varpi_safe = np.where(varpi > 0, varpi, 1.0)
    bt = np.where(varpi > 0, (-b_x * y[:, :, None] + b_y * x[:, :, None]) / varpi_safe, 0.0)
    bp = np.sqrt(np.clip(b_x**2 + b_y**2 + b_z**2 - bt**2, 0.0, None))

    bt_max = float(np.abs(bt).max())
    bp_max = float(bp.max())
    return {
        "bt_max_G": bt_max, "bp_max_G": bp_max,
        "ratio": bt_max / bp_max if bp_max > 0 else float("nan"),
    }


def plotfile_time_s(plotfile: Path) -> float:
    """Just the simulation time (s) stored in a plotfile's header -- no
    field data is touched, so this is cheap enough to call once per
    plotfile when filtering/sorting a run's step list (unlike
    load_slice/load_vector_grid, which actually read the grid).
    """
    ds = yt.load(str(plotfile))
    return float(ds.current_time.to("s").value)


_LOG_T_DYN_RE = re.compile(r"t_dyn = ([\d.eE+-]+) s")


def read_run_t_dyn(run_dir: Path, run_id: str) -> float | None:
    """t_dyn (s) for a run, read straight from its own log -- every run
    prints this at startup (problem_initialize.H), field-bearing or not,
    since it's a property of the background star, not of the seeded
    field. Returns None if the log is missing or doesn't have the line
    yet (e.g. a run that crashed before finishing initialization).
    """
    log_path = Path(run_dir) / f"run_{run_id}.log"
    if not log_path.exists():
        return None
    m = _LOG_T_DYN_RE.search(log_path.read_text(errors="replace"))
    return float(m.group(1)) if m else None


def find_validity_window_for_run(run_dir: Path, run_id: str, t_dyn_s: float) -> tuple | None:
    """The [t_field_relax, X_2pct] window (t/t_dyn) for the background
    star behind `run_id`, computed from a FIELD-FREE reference run (any
    run_id starting with "star", sharing the same t_dyn -- i.e. the same
    background star) via core/star_builder.py's find_measurement_window().

    Deliberately never computed from `run_id`'s own rho_c(t) if it has a
    field: a field-bearing run's early density dip is the field's own
    dynamical push, not the star's structural drift, and using it would
    give a wrong window (this exact mistake was made and caught earlier
    in this project -- see docs/teoria.md's validity-window section).
    Returns None if no matching field-free reference run is found on
    disk, or its own window isn't valid.
    """
    from core.star_builder import find_measurement_window, parse_rho_c_log

    run_dir = Path(run_dir)
    for candidate in list_run_ids(run_dir):
        if not candidate.startswith("star"):
            continue
        cand_t_dyn = read_run_t_dyn(run_dir, candidate)
        if cand_t_dyn is None or abs(cand_t_dyn - t_dyn_s) > 1e-6 * t_dyn_s:
            continue
        series = parse_rho_c_log(run_dir / f"run_{candidate}.log", cand_t_dyn)
        if not series:
            continue
        window = find_measurement_window(series, series[0][1])
        if window["valid"]:
            return window["window"]
    return None
