"""Central density of the non-rotating, field-supported 2.007 Msun star.

Phase 1 revisited (DIARIO 34): three magnetized runs and the field-free
control, all on the same sci-com binary (jobs 46067-46072). Data are the
TIME, MASS, KIN. ENERGY and MAXIMUM DENSITY columns of each run's
grid_diag.out, copied to data/rhomax_<run>.csv as t,mass,E_kin,rho_max.

The control stops at t = 4.92 s because it was killed from outside (the
sci-com epilog, DIARIO 34), not because it failed.

Run:  scf/.venv/bin/python3 investigations/phase1_revisit/plot_collapse.py
"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "plots"
COL_IN = 88.0 / 25.4

C_INK = "#0b0b0b"
C_MUTED = "#898781"
C_GRID = "#e1e0d9"
RHO_NEUTRON = 1.94e10     # mu_e = 2 threshold, as in the paper's Sec. 3

RUNS = [  # file tag, label, colour, linestyle, width
    ("p1ref192", r"$2.007\,M_\odot$, $192^3$", "#1b4f9c", "-", 1.1),
    ("p1ref96", r"$2.007\,M_\odot$, $96^3$", "#2a78d6", "-", 0.8),
    ("p1tight72", r"$2.007\,M_\odot$, smaller box", "#8fb8ea", (0, (3, 1.5)), 0.8),
    ("p1ctl208", r"field-free $1.346\,M_\odot$", C_INK, (0, (1, 1.2)), 0.9),
]

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8,
    "legend.fontsize": 6.2, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "legend.frameon": False, "figure.dpi": 200,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def main():
    OUT.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(COL_IN, COL_IN * 0.75))
    ax.grid(color=C_GRID, linewidth=0.4)
    for tag, lab, col, ls, lw in RUNS:
        t, _, _, rho = np.loadtxt(DATA / f"rhomax_{tag}.csv", delimiter=",", unpack=True)
        ax.plot(t, rho, color=col, linestyle=ls, linewidth=lw, label=lab)
        if tag != "p1ctl208":
            ax.plot(t[-1], rho[-1], marker="x", color=col, markersize=4, markeredgewidth=0.9)
    ax.axhline(RHO_NEUTRON, color=C_MUTED, linewidth=0.5, linestyle=(0, (2, 2)))
    ax.text(5.0, RHO_NEUTRON * 1.15, "neutronization", fontsize=6, color=C_MUTED, ha="right")
    ax.set_yscale("log")
    ax.set_xlim(0, 5.1)
    ax.set_ylim(2e8, 6e10)
    ax.set_xlabel("t (s)")
    ax.set_ylabel(r"$\rho_{\rm max}$  (g cm$^{-3}$)")
    ax.legend(loc="upper left")
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_collapse.{ext}")
    print("wrote", OUT / "fig_collapse.pdf")


if __name__ == "__main__":
    main()
