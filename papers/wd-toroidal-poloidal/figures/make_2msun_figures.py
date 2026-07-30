"""Figures for the 2 Msun toroidal-dominated configuration.

Run:  scf/.venv/bin/python3 papers/wd-toroidal-poloidal/figures/make_2msun_figures.py

Numbers are transcribed from the two investigation scripts that produced
them, investigations/mixed_2msun.py and the mesh study in
investigations/vector_potential_export.py's commit message, rather than
recomputed here: each SCF solve at the finest mesh takes minutes, and a
figure script that silently reruns them would drift from the values the
text quotes. If those investigations are rerun, update these tables.
"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = Path(__file__).resolve().parent
COL_IN = 88.0 / 25.4

C_FIELD = "#2a78d6"
C_SEED = "#52514e"
C_ALT = "#eb6834"
C_MUTED = "#898781"
C_GRID = "#e1e0d9"
C_BAND = "#f2f1ec"

# --- mesh convergence, domain fixed at 20 R_guess, K = 3.245e-3 --------
DR = np.array([8.219e6, 4.109e6, 2.056e6, 1.028e6])
NR = np.array([493, 986, 1971, 3942])
VE = np.array([7.162e-4, 3.702e-5, 1.305e-4, 1.724e-4])
MASS = np.array([2.0034, 2.0064, 2.0072, 2.0074])

# --- the poloidal trade-off, from investigations/mixed_2msun.csv -------
BT_BP = np.array([884.20, 98.24, 24.56, 8.84])       # energy ratio
VE_TOT = np.array([9.86e-5, 1.93e-3, 8.12e-3, 2.28e-2])
B_POLE = np.array([5.913e10, 1.774e11, 3.548e11, 5.913e11])
VE_GATE = 1.0e-3

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


def fig_convergence():
    fig, (ax, axv) = plt.subplots(
        2, 1, figsize=(COL_IN, COL_IN * 0.95), sharex=True,
        gridspec_kw=dict(height_ratios=[1.0, 1.0], hspace=0.10))
    for a in (ax, axv):
        style_axes(a)
        a.set_xscale("log")

    ax.plot(DR, MASS, color=C_FIELD, linewidth=1.4, marker="o", markersize=4.0,
            markeredgecolor="white", markeredgewidth=0.5, zorder=3)
    ax.axhline(MASS[-1], color=C_MUTED, linewidth=0.7, linestyle=(0, (4, 3)),
               zorder=1)
    ax.set_ylabel(r"$M$  ($M_\odot$)")
    ax.set_ylim(2.0020, 2.0086)
    for x, y, n in zip(DR, MASS, NR):
        ax.annotate(f"{n}", xy=(x, y), xytext=(0, -9),
                    textcoords="offset points", fontsize=6, color=C_MUTED,
                    ha="center")
    ax.text(1.1e6, 2.0079, "successive changes\n0.15%, 0.04%, 0.01%",
            fontsize=6.5, color=C_MUTED, ha="left", va="top")

    axv.axhspan(0, VE_GATE, color="#dce9fc", alpha=0.45, zorder=0)
    axv.plot(DR, VE, color=C_FIELD, linewidth=1.4, marker="o", markersize=4.0,
             markeredgecolor="white", markeredgewidth=0.5, zorder=3)
    axv.axhline(VE_GATE, color=C_MUTED, linewidth=0.7, linestyle=(0, (4, 3)),
                zorder=1)
    axv.text(9.0e6, VE_GATE * 1.25, "acceptance threshold", fontsize=6.5,
             color=C_MUTED, va="bottom")
    axv.set_yscale("log")
    axv.set_ylim(2e-5, 3e-3)
    axv.set_ylabel("VE")
    axv.set_xlabel(r"$\Delta r$  (cm)")
    axv.invert_xaxis()

    fig.savefig(HERE / "convergence_2msun.pdf")
    plt.close(fig)


def fig_tradeoff():
    fig, ax = plt.subplots(figsize=(COL_IN, COL_IN * 0.72))
    style_axes(ax)

    ax.axhspan(1e-5, VE_GATE, color="#dce9fc", alpha=0.45, zorder=0)
    ax.plot(BT_BP, VE_TOT, color=C_FIELD, linewidth=1.4, marker="o",
            markersize=4.0, markeredgecolor="white", markeredgewidth=0.5,
            zorder=3)
    ax.axhline(VE_GATE, color=C_MUTED, linewidth=0.7, linestyle=(0, (4, 3)),
               zorder=1)

    ax.plot([BT_BP[0]], [VE_TOT[0]], marker="o", markersize=7,
            markerfacecolor="none", markeredgecolor=C_ALT,
            markeredgewidth=1.1, zorder=4)
    ax.annotate("adopted", xy=(BT_BP[0], VE_TOT[0]), xytext=(-4, 14),
                textcoords="offset points", fontsize=7, color=C_ALT,
                ha="center")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.set_xlabel(r"$B_{\rm t}/B_{\rm p}$  (energy)")
    ax.set_ylabel("VE of the pair")
    ax.text(700, VE_GATE * 0.45, "certified", fontsize=6.5, color=C_MUTED,
            va="top")
    ax.text(700, VE_GATE * 1.6, "not certified", fontsize=6.5, color=C_MUTED,
            va="bottom")

    axb = ax.twiny()
    axb.set_xscale("log")
    axb.set_xlim(ax.get_xlim())
    axb.set_xticks(BT_BP)
    axb.set_xticklabels([f"{b:.0f}" for b in B_POLE / 1e10])
    axb.set_xlabel(r"exterior dipole $B_{\rm pole}$  ($10^{10}$ G)",
                   fontsize=7.5, labelpad=2)
    axb.tick_params(labelsize=6.5)
    for s in axb.spines.values():
        s.set_color(C_MUTED)

    fig.savefig(HERE / "tradeoff_2msun.pdf")
    plt.close(fig)


if __name__ == "__main__":
    fig_convergence()
    fig_tradeoff()
    print("wrote convergence_2msun.pdf, tradeoff_2msun.pdf")
