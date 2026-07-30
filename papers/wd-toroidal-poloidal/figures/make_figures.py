"""Figures for the mixed poloidal-toroidal paper.

Run:  scf/.venv/bin/python3 papers/wd-toroidal-poloidal/figures/make_figures.py

Reads two committed data files and nothing else:

  series.npz   time series cached by extract_series.py (rho_c(t) for the
               field-free background star and for the three extended
               seeded runs; E_tor/E_mag(t) and E_mag/|W|(t) for the ten
               seeds of batch C)
  ../../../braithwaite_app/data/results.csv
               the persisted per-seed measurements, read at run time

Nothing is recomputed here -- every diagnostic came out of the app's own
interior-masked extraction path (see extract_series.py).
"""

import csv
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RESULTS = REPO / "braithwaite_app" / "data" / "results.csv"
SERIES = HERE / "series.npz"

COL_IN = 88.0 / 25.4      # A&A single column

C_FIELD = "#2a78d6"
C_SEED = "#52514e"
C_ALT = "#eb6834"
C_MUTED = "#898781"
C_GRID = "#e1e0d9"
C_BAND = "#f2f1ec"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8,
    "axes.labelsize": 8,
    "legend.fontsize": 6.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.frameon": False,
    "figure.dpi": 200,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})


def style_axes(ax):
    ax.grid(True, axis="y", color=C_GRID, linewidth=0.4, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_linewidth(0.6)
        ax.spines[side].set_color(C_MUTED)


def load():
    d = np.load(SERIES)
    return d, json.loads(str(d["meta"]))


def seed_rows():
    rows = [r for r in csv.DictReader(RESULTS.open())
            if r["row_type"] == "seed"]
    by_res = {}
    for r in rows:
        # One entry per (resolution, seed); later re-runs of a seed
        # reproduce it, so the first is representative.
        by_res.setdefault((r["resolution"], int(float(r["seed"]))), r)
    return by_res


# ----------------------------------------------------------------------
# Fig. 1 -- the measurement window, and the drift it does not bound
# ----------------------------------------------------------------------

def fig_window():
    d, meta = load()
    lo, x1, x2 = meta["t_field_relax"], meta["X_1pct"], meta["X_2pct"]
    rho_ic = meta["rho_c_ic"]

    fig, (ax, axb) = plt.subplots(
        2, 1, figsize=(COL_IN, COL_IN * 1.12),
        gridspec_kw=dict(height_ratios=[1.25, 1.0], hspace=0.42),
    )
    for a in (ax, axb):
        style_axes(a)
        a.axhline(0.0, color=C_MUTED, linewidth=0.5, zorder=1)

    def dev(key):
        return 100.0 * (d[f"{key}_rho_c"] / rho_ic - 1.0)

    # (a) the window itself, and what the seeded runs actually do in it
    for pct, alpha in ((1.0, 0.55), (2.0, 0.28)):
        ax.axhspan(-pct, pct, color="#dce9fc", alpha=alpha, zorder=0)
    ax.axvspan(lo, x2, color=C_BAND, zorder=0)
    ax.plot(d["star_t"], dev("star"), color=C_FIELD, linewidth=1.3,
            label="field-free (window calibration)", zorder=3)
    for i, run in enumerate(meta["ext_runs"]):
        ax.plot(d[f"{run}_t"], dev(run), color=C_SEED, linewidth=0.9,
                linestyle=[(0, (4, 2)), (0, (1, 1.5)), (0, (5, 1, 1, 1))][i],
                label="seeded runs" if i == 0 else None, zorder=2)
    for xv, lab in ((x1, r"$X_{1\%}$"), (x2, r"$X_{2\%}$")):
        ax.axvline(xv, color=C_MUTED, linewidth=0.6, linestyle=(0, (1, 2)),
                   zorder=1)
        ax.text(xv, 1.5, lab, fontsize=6.5, color=C_MUTED, ha="center",
                va="bottom")
    ax.text(0.5 * (lo + x2), 3.15, "window", fontsize=6.5, color=C_MUTED,
            ha="center", va="top")
    ax.set_xlim(0, 1.42)
    ax.set_ylim(-13.0, 3.6)
    ax.set_xlabel(r"$t/t_{\rm dyn}$", labelpad=1.5)
    ax.set_ylabel(r"$\rho_c$ deviation (%)")
    ax.legend(loc="lower left", handlelength=2.6, borderaxespad=0.3,
              bbox_to_anchor=(0.0, 0.02))

    # (b) the field-free star over its whole baseline: it never settles
    axb.axhspan(-2, 2, color="#dce9fc", alpha=0.28, zorder=0)
    axb.axvspan(lo, x2, color=C_BAND, zorder=0)
    axb.plot(d["star_t"], dev("star"), color=C_FIELD, linewidth=1.1,
             zorder=3)
    axb.set_xlim(0, 16.5)
    axb.set_ylim(-7.5, 9.5)
    axb.set_xlabel(r"$t/t_{\rm dyn}$", labelpad=1.5)
    axb.set_ylabel(r"$\rho_c$ deviation (%)")
    axb.text(1.6, 8.2, "field-free, same run, full baseline", fontsize=6.5,
             color=C_MUTED, ha="left", va="top")

    fig.savefig(HERE / "window.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------
# Fig. 2 -- what the relaxation does to the ratio, and to the energy
# ----------------------------------------------------------------------

def fig_relaxation():
    d, meta = load()
    lo, x2 = meta["t_field_relax"], meta["X_2pct"]
    seeds = meta["seeds_C"]
    rows = seed_rows()

    fig, (ax, axb) = plt.subplots(
        2, 1, figsize=(COL_IN, COL_IN * 1.10), sharex=True,
        gridspec_kw=dict(hspace=0.10),
    )
    for a in (ax, axb):
        style_axes(a)
        a.axvspan(lo, x2, color=C_BAND, zorder=0)

    for s in seeds:
        ax.plot(d[f"seed{s}_t"], d[f"seed{s}_ratio"], color=C_FIELD,
                linewidth=0.8, alpha=0.75, zorder=2)
        axb.plot(d[f"seed{s}_t"], d[f"seed{s}_emag"], color=C_FIELD,
                 linewidth=0.8, alpha=0.75, zorder=2)

    ax.axhline(1 / 3, color=C_ALT, linewidth=0.8, linestyle=(0, (4, 2)),
               zorder=3)
    # The 0.5 threshold has to be on the axis, or "well below 0.5" is a
    # claim the reader cannot check against the figure.
    ax.axhline(0.5, color=C_MUTED, linewidth=0.7, linestyle=(0, (4, 3)),
               zorder=3)
    ax.text(0.02, 0.507, r"$0.5$: equal poloidal and toroidal energy",
            fontsize=6.5, color=C_MUTED, va="bottom")
    ax.text(0.5 * (lo + x2), 0.478, "window", fontsize=6.5, color=C_MUTED,
            ha="center", va="top")
    ax.set_ylabel(r"$E_{\rm tor}/E_{\rm mag}$")
    ax.set_ylim(0.24, 0.545)

    # The one 128^3 run, at the time it was measured: the energy left in
    # the field at a given time depends on resolution, the ratio barely.
    r128 = rows.get(("128", 42))
    if r128:
        axb.plot([float(r128["t_ttdyn_measured"])],
                 [float(r128["E_mag_over_W"])], linestyle="none",
                 marker="D", markersize=4.5, color=C_ALT,
                 markeredgecolor="white", markeredgewidth=0.5,
                 zorder=5, label=r"seed 42, $128^3$")
        r64 = rows.get(("64", 42))
        if r64:
            axb.plot([float(r64["t_ttdyn_measured"])],
                     [float(r64["E_mag_over_W"])], linestyle="none",
                     marker="o", markersize=4.0, markerfacecolor="white",
                     markeredgecolor=C_ALT, markeredgewidth=0.9, zorder=5,
                     label=r"seed 42, $64^3$")
        axb.legend(loc="lower left", handlelength=1.2, borderaxespad=0.3,
                   ncol=2, columnspacing=1.0)
    axb.set_yscale("log")
    axb.set_ylabel(r"$E_{\rm mag}/|W|$")
    axb.set_xlabel(r"$t/t_{\rm dyn}$")
    axb.set_xlim(-0.03, 1.30)
    axb.set_ylim(1.2e-3, 7e-2)
    axb.set_yticks([2e-3, 5e-3, 1e-2, 2e-2, 5e-2])
    axb.set_yticklabels(["0.002", "0.005", "0.01", "0.02", "0.05"])
    axb.minorticks_off()
    axb.annotate(r"$3.2\times$", xy=(0.735, 6.6e-3), fontsize=6.5,
                 color=C_ALT, ha="right")
    axb.annotate("", xy=(0.803, 1.35e-2), xytext=(0.803, 4.2e-3),
                 arrowprops=dict(arrowstyle="<->", linewidth=0.6,
                                 color=C_ALT, shrinkA=0.5, shrinkB=0.5))

    fig.savefig(HERE / "relaxation.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------
# Fig. 3 -- the ratio does not move away from where the seed put it
# ----------------------------------------------------------------------

def fig_paired():
    d, meta = load()
    lo, x2 = meta["t_field_relax"], meta["X_2pct"]
    seeds = meta["seeds_C"]

    initial, windowed = [], []
    for s in seeds:
        t, ratio = d[f"seed{s}_t"], d[f"seed{s}_ratio"]
        m = (t >= lo) & (t <= x2)
        initial.append(float(ratio[0]))
        windowed.append(float(ratio[m].mean()))
    initial, windowed = np.array(initial), np.array(windowed)
    delta = windowed - initial
    se = delta.std(ddof=1) / np.sqrt(len(delta))

    fig, ax = plt.subplots(figsize=(COL_IN, COL_IN * 0.80))
    style_axes(ax)
    ax.grid(False)

    ax.axhline(1 / 3, color=C_ALT, linewidth=0.8, linestyle=(0, (4, 2)),
               zorder=1)
    ax.text(2.26, 1 / 3, r"$1/3$", fontsize=7, color=C_ALT, ha="left",
            va="center")

    for a, b, s in zip(initial, windowed, seeds):
        up = b >= a
        ax.plot([1, 2], [a, b], color=C_FIELD if up else C_SEED,
                linewidth=0.9, alpha=0.85, zorder=2)
        ax.plot([1, 2], [a, b], linestyle="none", marker="o",
                markersize=3.4, color=C_FIELD if up else C_SEED,
                markeredgecolor="white", markeredgewidth=0.5, zorder=3)

    for xpos, vals, lab in ((1, initial, "seeded"), (2, windowed, "window")):
        ax.plot([xpos - 0.13, xpos + 0.13], [vals.mean()] * 2,
                color="#0b0b0b", linewidth=1.6, zorder=4)
        ax.text(xpos, 0.259, f"mean {vals.mean():.3f}", fontsize=6.5,
                color="#0b0b0b", ha="center")

    ax.annotate(
        f"paired change ${delta.mean():+.3f} \\pm {se:.3f}$\n"
        r"($n=10$, indistinguishable from zero)",
        xy=(1.5, 0.435), fontsize=7, color="#0b0b0b", ha="center",
        va="center")

    ax.set_xlim(0.74, 2.42)
    ax.set_ylim(0.252, 0.455)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([r"$t = 0$", r"inside window"])
    ax.set_ylabel(r"$E_{\rm tor}/E_{\rm mag}$")
    ax.tick_params(top=False, right=False)

    fig.savefig(HERE / "paired.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------
# Figs. 4 and 5 -- the spatial structure of one relaxed field
# Slices and isosurfaces cached by extract_fields.py.
# ----------------------------------------------------------------------

FIELDS3D = HERE / "fields3d.npz"
PANELS = [("Bmag", r"$|\mathbf{B}|$"),
          ("Bt", r"$|B_{\rm t}|$  (toroidal)"),
          ("Bp", r"$|B_{\rm p}|$  (poloidal)")]

# Two-column A&A figure.
WIDE_IN = 180.0 / 25.4


def _load3d():
    if not FIELDS3D.exists():
        print(f"skipping field figures: run extract_fields.py first "
              f"({FIELDS3D.name})")
        return None, None
    d = np.load(FIELDS3D)
    return d, json.loads(str(d["meta"]))


def fig_slices():
    from matplotlib.colors import LinearSegmentedColormap, LogNorm

    d, meta = _load3d()
    if d is None:
        return

    ext = np.asarray(d["extent_cm"]) / 1e8
    rho = np.asarray(d["slice_density"]).T
    inside = rho > meta["rho_peak"] * meta["rho_floor_frac"]

    # One sequential ramp, one hue, monotone in lightness, and one shared
    # range across the three panels -- the panels are the same quantity
    # measured three ways, so a per-panel autoscale would erase exactly
    # the comparison the figure exists to make.
    cmap = LinearSegmentedColormap.from_list(
        "b", ["#f4f8fe", "#cde2fb", "#9ec5f4", "#5598e7", "#256abf",
              "#184f95", "#0d366b"])
    cmap.set_bad("#ffffff")
    total = np.asarray(d["slice_Bmag"]).T
    vmax = float(np.percentile(total[inside], 99.5))
    vmin = vmax / 3e2

    fig, axes = plt.subplots(1, 3, figsize=(WIDE_IN, WIDE_IN * 0.40))
    for ax, (key, title) in zip(axes, PANELS):
        arr = np.ma.masked_less_equal(np.asarray(d[f"slice_{key}"]).T, 0.0)
        mesh = ax.imshow(arr, origin="lower", extent=ext, cmap=cmap,
                         norm=LogNorm(vmin=vmin, vmax=vmax),
                         interpolation="bilinear", rasterized=True,
                         zorder=1)
        ax.contour(np.linspace(ext[0], ext[1], rho.shape[1]),
                   np.linspace(ext[2], ext[3], rho.shape[0]),
                   inside.astype(float), levels=[0.5], colors="#0b0b0b",
                   linewidths=0.7, zorder=3)
        ax.set_title(title, fontsize=8, pad=3)
        ax.set_xlabel(r"$x$  ($10^8$ cm)", labelpad=1.5)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=6.5, pad=1.5)
        for s in ax.spines.values():
            s.set_color(C_MUTED)
    axes[0].set_ylabel(r"$z$  ($10^8$ cm)", labelpad=1.0)
    for ax in axes[1:]:
        ax.set_yticklabels([])

    cb = fig.colorbar(mesh, ax=axes, fraction=0.020, pad=0.015)
    cb.set_label("field strength (G)", fontsize=7.5, labelpad=2)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=6.5, width=0.5, length=2, pad=1.5)

    fig.savefig(HERE / "slices.pdf")
    plt.close(fig)


def _lines(d, key):
    """The cached ragged polylines for one field, in units of 1e8 cm."""
    xyz = np.asarray(d[f"line_{key}_xyz"]) / 1e8
    mag = np.asarray(d[f"line_{key}_mag"])
    off = np.asarray(d[f"line_{key}_off"])
    return [(xyz[off[i]:off[i + 1]], mag[off[i]:off[i + 1]])
            for i in range(len(off) - 1)]


def fig_fieldlines():
    from matplotlib.colors import LinearSegmentedColormap, LogNorm
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    d, meta = _load3d()
    if d is None:
        return

    cmap = LinearSegmentedColormap.from_list(
        "b", ["#cde2fb", "#9ec5f4", "#5598e7", "#256abf", "#184f95",
              "#0d366b"])

    # One color range across the panels: the panels follow three
    # different fields but the quantity is the same, so a per-panel
    # autoscale would hide that the toroidal part is the weaker one.
    all_mag = np.concatenate([np.asarray(d[f"line_{k}_mag"])
                              for k, _ in PANELS])
    vmax = float(np.percentile(all_mag, 99.5))
    norm = LogNorm(vmin=vmax / 60.0, vmax=vmax)

    r_star = meta["r_star_cm"] / 1e8
    lim = 1.05 * r_star

    fig = plt.figure(figsize=(WIDE_IN, WIDE_IN * 0.38))
    # Explicit column for the colorbar: letting it steal space from the 3D
    # axes puts it on top of the wander panel.
    # Colourbar horizontal, under the three 3D panels. Any vertical
    # placement puts its tick labels and its own axis label into the wander
    # panel's y-axis label; a 3D axis's get_position() returns the whole
    # subplot slot, which is much wider than the visible box, so hand
    # placement from it lands inside the neighbour.
    gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 0.72],
                          height_ratios=[1, 0.075],
                          left=0.01, right=0.985, bottom=0.13, top=0.96,
                          wspace=0.30, hspace=0.05)
    axes = [fig.add_subplot(gs[0, i], projection="3d") for i in range(3)]
    cax = fig.add_subplot(gs[1, 0:3])
    axw = fig.add_subplot(gs[0:2, 3])

    # A faint wireframe at the stellar surface, for scale.
    u = np.linspace(0, 2 * np.pi, 37)
    v = np.linspace(0, np.pi, 19)
    sx = r_star * np.outer(np.cos(u), np.sin(v))
    sy = r_star * np.outer(np.sin(u), np.sin(v))
    sz = r_star * np.outer(np.ones_like(u), np.cos(v))

    for ax, (key, title) in zip(axes, PANELS):
        ax.plot_wireframe(sx, sy, sz, rstride=6, cstride=6,
                          color=C_GRID, linewidth=0.25, zorder=1)

        for xyz, mag in _lines(d, key):
            segs = np.stack([xyz[:-1], xyz[1:]], axis=1)
            lc = Line3DCollection(
                segs, cmap=cmap, norm=norm, linewidths=0.85,
                capstyle="round", zorder=3)
            lc.set_array(0.5 * (mag[:-1] + mag[1:]))
            ax.add_collection3d(lc)

            # Direction: a few arrowheads per line, tangent to it.
            n = len(xyz)
            for frac in (0.3, 0.62, 0.9):
                i = int(frac * (n - 2))
                p, t = xyz[i], xyz[i + 1] - xyz[i]
                if np.linalg.norm(t) == 0:
                    continue
                t = t / np.linalg.norm(t)
                ax.quiver(*p, *t, length=0.30 * r_star, normalize=True,
                          arrow_length_ratio=0.55, linewidth=0.7,
                          color=cmap(norm(mag[i])), zorder=4)

        ax.set_title(title, fontsize=8, pad=-4)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=18, azim=-56)
        ax.tick_params(labelsize=5.5, pad=-3.5)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.set_facecolor("white")
            axis.pane.set_edgecolor(C_GRID)
            axis.line.set_color(C_MUTED)
            axis._axinfo["grid"]["color"] = C_GRID
            axis._axinfo["grid"]["linewidth"] = 0.3
            axis.set_ticks([-2, 0, 2])
    axes[0].set_xlabel(r"$x$ ($10^8$ cm)", fontsize=6.5, labelpad=-8)
    axes[0].set_ylabel(r"$y$", fontsize=6.5, labelpad=-8)
    axes[0].set_zlabel(r"$z$", fontsize=6.5, labelpad=-8)

    # (d) how far a line actually travels. The 3D panels are truncated at a
    # fixed integration length -- their ends are a drawing choice, not
    # physics -- so the wandering is shown here instead, where it is legible.
    # Total field only. The poloidal line stays inside a region of ~0.3
    # R_star and its curve is a flat band; the toroidal tracer oscillates
    # between two neighbouring states rather than following a circle, so
    # neither curve says anything the reader can use, and plotting them
    # would only invite the eye to compare three things of which two are
    # artifacts.
    r_star = meta["r_star_cm"] / 1e8
    arc = np.asarray(d["wander_Bmag_arc"]) / 1e8 / r_star
    dist = np.asarray(d["wander_Bmag_dist"]) / 1e8 / r_star
    axw.plot(arc, dist, color=C_FIELD, linewidth=0.7, zorder=3)
    axw.axhline(1.0, color=C_MUTED, linewidth=0.7, linestyle=(0, (4, 3)),
                zorder=1)
    axw.text(0.5, 1.03, "stellar radius", fontsize=6, color=C_MUTED)
    style_axes(axw)
    axw.set_xlabel("arc length along the line " r"($R_\star$)", fontsize=7,
                   labelpad=1.5)
    axw.set_ylabel("distance from start " r"($R_\star$)", fontsize=7,
                   labelpad=1.5)
    axw.tick_params(labelsize=6, pad=1.5)
    axw.set_xlim(0, 70)
    axw.set_ylim(0, 1.55)
    axw.set_title(r"(d) one $|\mathbf{B}|$ line, followed far", fontsize=7,
                  pad=3, color=C_MUTED)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("field strength along the line (G)", fontsize=7,
                 labelpad=2)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=6.5, width=0.5, length=2, pad=1.5)

    fig.savefig(HERE / "fieldlines.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------
# Fig. 6 -- each component in the plane where it actually lives
# ----------------------------------------------------------------------

def fig_components():
    from matplotlib.colors import LinearSegmentedColormap, LogNorm, TwoSlopeNorm

    d, meta = _load3d()
    if d is None or "eq_Bt_signed" not in d:
        print("skipping components figure: re-run extract_fields.py")
        return

    rho_floor = meta["rho_peak"] * meta["rho_floor_frac"]

    # Equatorial slice: yt returns a z-normal slice as [y, x], so no
    # transpose (verified by rebuilding B_phi from sliced B_x, B_y).
    bt = np.asarray(d["eq_Bt_signed"]) / 1e13
    rho_eq = np.asarray(d["eq_density"])
    ext_eq = np.asarray(d["eq_extent_cm"]) / 1e8

    # Meridional slice: y-normal comes back as [x, z], so transpose to put
    # z on the vertical.
    bx = np.asarray(d["mer_B_x"]).T
    bz = np.asarray(d["mer_B_z"]).T
    bp = np.asarray(d["slice_Bp"]).T
    rho_mer = np.asarray(d["slice_density"]).T
    ext = np.asarray(d["extent_cm"]) / 1e8

    fig, (axe, axm) = plt.subplots(1, 2, figsize=(WIDE_IN * 0.78,
                                                  WIDE_IN * 0.36),
                                   gridspec_kw=dict(wspace=0.42))

    # --- (a) B_phi in the equatorial plane, signed -----------------------
    # Diverging: two hues with a neutral midpoint, symmetric limits, so the
    # zero contour reads as zero and the two senses of circulation are
    # equally weighted.
    div = LinearSegmentedColormap.from_list(
        "div", ["#0d366b", "#2a78d6", "#9ec5f4", "#efeee9", "#f6b79b",
                "#eb6834", "#8c3210"])
    vlim = float(np.percentile(np.abs(bt[rho_eq > rho_floor]), 99))
    im = axe.imshow(np.ma.masked_where(rho_eq <= rho_floor, bt),
                    origin="lower", extent=ext_eq, cmap=div,
                    norm=TwoSlopeNorm(vmin=-vlim, vcenter=0.0, vmax=vlim),
                    interpolation="bilinear", rasterized=True, zorder=1)
    gx = np.linspace(ext_eq[0], ext_eq[1], bt.shape[1])
    gy = np.linspace(ext_eq[2], ext_eq[3], bt.shape[0])
    # Zero contour only where there is a star: outside, B_phi is the
    # numerical floor and its sign changes are noise, which would draw a
    # maze over the whole panel.
    inside_eq = rho_eq > rho_floor
    axe.contour(gx, gy, np.ma.masked_where(~inside_eq, bt), levels=[0.0],
                colors="#0b0b0b", linewidths=0.35, alpha=0.55, zorder=2)
    axe.contour(gx, gy, (rho_eq > rho_floor).astype(float), levels=[0.5],
                colors="#0b0b0b", linewidths=0.7, zorder=3)
    # In-plane vectors of B_t: purely azimuthal, B_phi * e_phi with
    # e_phi = (-y, x)/varpi. Their geometry is fixed by construction --
    # concentric circles -- so what the arrows add is the sense of
    # circulation, read directly instead of decoded from the colour.
    X, Y = np.meshgrid(gx, gy)
    varpi = np.sqrt(X**2 + Y**2)
    safe = np.where(varpi > 0, varpi, 1.0)
    btx = np.where(inside_eq, bt * (-Y / safe), np.nan)
    bty = np.where(inside_eq, bt * (X / safe), np.nan)
    axe.streamplot(gx, gy, btx, bty, color="#0b0b0b", linewidth=0.4,
                   density=1.7, arrowsize=0.55, zorder=4)

    axe.set_title(r"(a) $B_\varphi$ in the equatorial plane", fontsize=7.5,
                  pad=3)
    axe.set_xlabel(r"$x$  ($10^8$ cm)", labelpad=1.5)
    axe.set_ylabel(r"$y$  ($10^8$ cm)", labelpad=1.0)
    cb = fig.colorbar(im, ax=axe, fraction=0.046, pad=0.02)
    cb.set_label(r"$B_\varphi$  ($10^{13}$ G)", fontsize=7)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=6, width=0.5, length=2, pad=1.5)

    # --- (b) B_p in the meridional plane, as streamlines -----------------
    cmap = LinearSegmentedColormap.from_list(
        "b", ["#f4f8fe", "#cde2fb", "#9ec5f4", "#5598e7", "#256abf",
              "#184f95", "#0d366b"])
    inside = rho_mer > rho_floor
    vmax = float(np.percentile(bp[inside], 99))
    im2 = axm.imshow(np.ma.masked_where(~inside | (bp <= 0.0), bp),
                     origin="lower",
                     extent=ext, cmap=cmap,
                     norm=LogNorm(vmin=vmax / 3e2, vmax=vmax),
                     interpolation="bilinear", rasterized=True, zorder=1)
    mx = np.linspace(ext[0], ext[1], bp.shape[1])
    mz = np.linspace(ext[2], ext[3], bp.shape[0])
    # Same reason: streamlines through the vacuum are streamlines of the
    # floor, and they swamp the structure that is actually there.
    bx_in = np.where(inside, bx, np.nan)
    bz_in = np.where(inside, bz, np.nan)
    axm.streamplot(mx, mz, bx_in, bz_in, color="#0b0b0b", linewidth=0.4,
                   density=1.7, arrowsize=0.55, zorder=3)
    axm.contour(mx, mz, inside.astype(float), levels=[0.5],
                colors="#0b0b0b", linewidths=0.7, zorder=4)
    axm.set_title(r"(b) $\mathbf{B}_{\rm p}$ in the meridional plane",
                  fontsize=7.5, pad=3)
    axm.set_xlabel(r"$x$  ($10^8$ cm)", labelpad=1.5)
    axm.set_ylabel(r"$z$  ($10^8$ cm)", labelpad=1.0)
    cb2 = fig.colorbar(im2, ax=axm, fraction=0.046, pad=0.02)
    cb2.set_label(r"$|\mathbf{B}_{\rm p}|$  (G)", fontsize=7)
    cb2.outline.set_linewidth(0.5)
    cb2.ax.tick_params(labelsize=6, width=0.5, length=2, pad=1.5)

    for ax in (axe, axm):
        ax.set_aspect("equal")
        # Cropped to just outside the star (R = 2.46e8 cm): the vacuum
        # carries nothing and only shrinks the part that does.
        ax.set_xlim(-3.0, 3.0)
        ax.set_ylim(-3.0, 3.0)
        ax.set_xticks([-2, 0, 2])
        ax.set_yticks([-2, 0, 2])
        ax.tick_params(labelsize=6.5, pad=1.5)
        for s_ in ax.spines.values():
            s_.set_color(C_MUTED)

    fig.savefig(HERE / "components.pdf")
    plt.close(fig)


if __name__ == "__main__":
    fig_window()
    fig_relaxation()
    fig_paired()
    fig_slices()
    fig_fieldlines()
    fig_components()
    print("wrote", *(p.name for p in sorted(HERE.glob("*.pdf"))))
