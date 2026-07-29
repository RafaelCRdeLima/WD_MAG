"""Figures for the purely-toroidal paper.

Run:  scf/.venv/bin/python3 papers/wd-toroidal/figures/make_figures.py

Every number plotted here is a measured result already recorded in the
repository -- nothing is recomputed, estimated or smoothed:

  Fig. 1  investigations/rho_c_1e9_M2_configurations.csv (read at run time,
          mode == "1_toroidal_only")
  Fig. 2  docs/teoria.md Sec 6.2d, the certified M(rho_c) table
  Fig. 3  docs/teoria.md Sec 6.1, the Omega_c continuation table

Figs. 2 and 3 are transcribed below rather than read from a file because
the sweeps behind them were run interactively and only their results were
kept; the tables in teoria.md are the source of record. If those tables
are ever corrected, correct them here too.
"""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
CSV = REPO / "investigations" / "rho_c_1e9_M2_configurations.csv"

# A&A single-column text width is 88 mm.
COL_IN = 88.0 / 25.4

# Two series maximum per panel, so the first two categorical hues suffice.
# Grayscale-safe by construction: every series also differs in dash pattern
# and marker shape, which is what survives a black-and-white print run.
C_FIELD = "#2a78d6"   # magnetized sequence
C_FREE = "#52514e"    # field-free reference
C_MUTED = "#898781"   # thresholds, annotations
C_GRID = "#e1e0d9"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "legend.fontsize": 7,
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


# ----------------------------------------------------------------------
# Fig. 1 -- equilibrium mass along the toroidal sequence at rho_c = 1e9
# ----------------------------------------------------------------------

def fig_mass_vs_field():
    rows = [r for r in csv.DictReader(CSV.open())
            if r["mode"] == "1_toroidal_only"]
    seq = sorted(
        (r for r in rows if r["is_interpolated_M2_crossing"] == "False"),
        key=lambda r: float(r["K_toroidal"]),
    )
    crossing = next(r for r in rows
                    if r["is_interpolated_M2_crossing"] == "True")

    b = [float(r["B_tor_max_G"]) / 1e13 for r in seq]
    m = [float(r["M_Msun"]) for r in seq]
    cert = [r["VE_certified"] == "True" for r in seq]

    fig, ax = plt.subplots(figsize=(COL_IN, COL_IN * 0.78))
    style_axes(ax)

    ax.axhline(2.0, color=C_MUTED, linewidth=0.7, linestyle=(0, (4, 3)),
               zorder=1)

    ax.plot(b, m, color=C_FIELD, linewidth=1.4, zorder=2)
    ax.plot([x for x, c in zip(b, cert) if c],
            [y for y, c in zip(m, cert) if c],
            linestyle="none", marker="o", markersize=4.0,
            color=C_FIELD, markeredgecolor="white", markeredgewidth=0.5,
            zorder=3)
    ax.plot([x for x, c in zip(b, cert) if not c],
            [y for y, c in zip(m, cert) if not c],
            linestyle="none", marker="o", markersize=4.0,
            markerfacecolor="white", markeredgecolor=C_FIELD,
            markeredgewidth=0.9, zorder=3)

    bx = float(crossing["B_tor_max_G"]) / 1e13
    ax.plot([bx], [2.0], marker="*", markersize=10, color=C_FIELD,
            markeredgecolor="white", markeredgewidth=0.7, zorder=5)

    # Both leaders run nearly horizontally through the empty upper-left
    # block, so they cross neither the sequence nor each other.
    ax.annotate(
        r"$M = 2\,M_\odot$ at $B_{\varphi,\max} = 4.3\times10^{13}$ G",
        xy=(bx, 2.0), xytext=(0.15, 2.13), fontsize=7, color="#0b0b0b",
        ha="left", va="center",
        arrowprops=dict(arrowstyle="-", linewidth=0.5, color=C_MUTED,
                        shrinkA=3, shrinkB=4),
    )
    ax.annotate(
        r"open symbol: not certified (VE $\geq 10^{-3}$)",
        xy=(b[-1], m[-1]), xytext=(0.15, 2.45), fontsize=7, color=C_MUTED,
        ha="left", va="center",
        arrowprops=dict(arrowstyle="-", linewidth=0.5, color=C_MUTED,
                        shrinkA=3, shrinkB=4),
    )

    ax.set_xlabel(r"$B_{\varphi,\max}$  ($10^{13}$ G)")
    ax.set_ylabel(r"$M$  ($M_\odot$)")
    ax.set_xlim(-0.15, 4.9)
    ax.set_ylim(1.28, 2.58)
    fig.savefig(HERE / "mass_vs_field.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------
# Fig. 2 -- the mass gain is independent of central density
# docs/teoria.md Sec 6.2d.  m = 1, K = 3e-3, no rotation.
# ----------------------------------------------------------------------

RHO_C = [1.0e10, 1.5e10, 1.0e12]
M_FREE = [1.4108, 1.4159, 1.4323]
M_TOR = [2.0722, 2.0874, 2.1426]
E_TOR_W = [0.1809, 0.1805, 0.1790]
RHO_NEUTRON = 1.94e10   # 16O, mu_e = 2; scf/eos.py, Boshkayev et al. 2013


def fig_mass_vs_rhoc():
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(COL_IN, COL_IN * 1.15), sharex=True,
        gridspec_kw=dict(height_ratios=[2.1, 1.0], hspace=0.08),
    )
    for a in (ax, axr):
        style_axes(a)
        a.axvspan(RHO_NEUTRON, 3e12, color="#f2f1ec", zorder=0)
        a.axvline(RHO_NEUTRON, color=C_MUTED, linewidth=0.6,
                  linestyle=(0, (1, 2)), zorder=1)

    ax.axhline(2.0, color=C_MUTED, linewidth=0.7, linestyle=(0, (4, 3)),
               zorder=1)
    ax.plot(RHO_C, M_TOR, color=C_FIELD, linewidth=1.4, marker="o",
            markersize=4.0, markeredgecolor="white", markeredgewidth=0.5,
            label=r"$E_{\rm tor}/|W| \simeq 0.18$", zorder=3)
    ax.plot(RHO_C, M_FREE, color=C_FREE, linewidth=1.4, marker="s",
            markersize=3.6, linestyle=(0, (5, 2)), markeredgecolor="white",
            markeredgewidth=0.5, label="field-free", zorder=3)

    ax.text(2.4e10, 2.42, "electron capture\n(not a white dwarf)",
            fontsize=6.5, color=C_MUTED, ha="left", va="top")
    ax.annotate("", xy=(1.0e10, M_TOR[0]), xytext=(1.0e10, M_FREE[0]),
                arrowprops=dict(arrowstyle="<->", linewidth=0.6,
                                color="#0b0b0b", shrinkA=1.5, shrinkB=1.5))
    ax.text(1.06e10, 1.73, r"$+47\%$", fontsize=7, color="#0b0b0b")

    ax.set_ylabel(r"$M$  ($M_\odot$)")
    ax.set_ylim(1.30, 2.50)
    ax.legend(loc="center", bbox_to_anchor=(0.62, 0.44), handlelength=2.4)

    # A +-1% band around the mean makes "constant to 0.9%" something the
    # reader can check against the axis instead of taking on faith.
    mean_e = sum(E_TOR_W) / len(E_TOR_W)
    axr.axhspan(mean_e * 0.99, mean_e * 1.01, color="#e8eefb", zorder=1)
    axr.plot(RHO_C, E_TOR_W, color=C_FIELD, linewidth=1.4, marker="o",
             markersize=4.0, markeredgecolor="white", markeredgewidth=0.5,
             zorder=3)
    axr.set_ylabel(r"$E_{\rm tor}/|W|$")
    axr.set_xlabel(r"$\rho_c$  (g cm$^{-3}$)")
    axr.set_xscale("log")
    axr.set_xlim(7e9, 2.2e12)
    axr.set_ylim(0.1705, 0.1895)
    axr.set_yticks([0.175, 0.180, 0.185])
    axr.text(2.0e12, mean_e * 1.01 + 0.0004, r"mean $\pm 1\%$",
             fontsize=6.5, color=C_MUTED, ha="right", va="bottom")

    fig.savefig(HERE / "mass_vs_rhoc.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------
# Fig. 3 -- rotation terminates the sequence earlier, not later
# docs/teoria.md Sec 6.1.  rho_c = 1e12, nr = 129, ntheta = 65.
# The Omega_c = 26.5, K = 0 point is the spurious branch discussed in the
# text and is deliberately not plotted.
# ----------------------------------------------------------------------

OMEGA_FREE = [0, 10, 15, 20, 24, 26]
LOSS_FREE = [0.000, 0.017, 0.040, 0.074, 0.112, 0.135]
VE_OMEGA_FREE = [0, 10, 20, 24, 26]
VE_FREE = [1.78e-3, 1.77e-3, 1.77e-3, 1.77e-3, 1.76e-3]

OMEGA_TOR = [0, 10, 15, 20, 24, 26, 26.5, 32]
LOSS_TOR = [0.000, 0.139, 0.460, 1.871, 2.525, 2.844, 2.922, 3.722]
VE_TOR = [5.18e-4, 5.27e-4, 5.37e-4, 2.41e-3, 1.90e-2, 3.82e-2, 4.44e-2,
          1.61e-1]


def fig_rotation():
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(COL_IN, COL_IN * 1.20), sharex=True,
        gridspec_kw=dict(hspace=0.10),
    )
    for a in (ax1, ax2):
        style_axes(a)
        a.axvspan(17, 18, color="#f2f1ec", zorder=0)

    ax1.axhline(1.0, color=C_MUTED, linewidth=0.7, linestyle=(0, (4, 3)),
                zorder=1)
    ax1.plot(OMEGA_TOR, LOSS_TOR, color=C_FIELD, linewidth=1.4, marker="o",
             markersize=4.0, markeredgecolor="white", markeredgewidth=0.5,
             label=r"$K = 3\times10^{-3}$", zorder=3)
    ax1.plot(OMEGA_FREE, LOSS_FREE, color=C_FREE, linewidth=1.4, marker="s",
             markersize=3.6, linestyle=(0, (5, 2)), markeredgecolor="white",
             markeredgewidth=0.5, label=r"$K = 0$", zorder=3)
    ax1.set_ylabel("mass-loss ratio")
    ax1.set_ylim(-0.25, 4.1)
    ax1.text(0.6, 1.12, "breakup", color=C_MUTED, fontsize=7)
    ax1.legend(loc="upper left", handlelength=2.4, borderaxespad=0.4,
               bbox_to_anchor=(0.0, 0.92))

    ax2.axhline(1e-3, color=C_MUTED, linewidth=0.7, linestyle=(0, (4, 3)),
                zorder=1)
    ax2.plot(OMEGA_TOR, VE_TOR, color=C_FIELD, linewidth=1.4, marker="o",
             markersize=4.0, markeredgecolor="white", markeredgewidth=0.5,
             zorder=3)
    ax2.plot(VE_OMEGA_FREE, VE_FREE, color=C_FREE, linewidth=1.4,
             marker="s", markersize=3.6, linestyle=(0, (5, 2)),
             markeredgecolor="white", markeredgewidth=0.5, zorder=3)
    ax2.set_yscale("log")
    ax2.set_ylabel("VE")
    ax2.set_xlabel(r"$\Omega_c$  (rad s$^{-1}$)")
    ax2.set_xlim(-1.5, 34)
    ax2.set_ylim(2e-4, 6e-1)
    ax2.text(33.4, 9.2e-4, "acceptance threshold", color=C_MUTED,
             fontsize=6.5, ha="right", va="top")
    ax2.annotate("practical limit", xy=(17.4, 5.5e-3), xytext=(9.2, 1.7e-2),
                 fontsize=6.5, color=C_MUTED, ha="left", va="center",
                 arrowprops=dict(arrowstyle="-", linewidth=0.5,
                                 color=C_MUTED, shrinkA=3, shrinkB=2))

    fig.savefig(HERE / "rotation_termination.pdf")
    plt.close(fig)


if __name__ == "__main__":
    fig_mass_vs_field()
    fig_mass_vs_rhoc()
    fig_rotation()
    print("wrote", *(p.name for p in sorted(HERE.glob("*.pdf"))))
