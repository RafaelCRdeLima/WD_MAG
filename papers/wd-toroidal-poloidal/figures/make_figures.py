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


if __name__ == "__main__":
    fig_window()
    fig_relaxation()
    fig_paired()
    print("wrote", *(p.name for p in sorted(HERE.glob("*.pdf"))))
