"""The equation of state itself: Eqs. 1 and 2 of the paper.

Run:  scf/.venv/bin/python3 papers/wd-toroidal-poloidal/figures/make_eos_figure.py

Kept apart from make_figures.py because it shares nothing with it: those
figures read cached output from real runs, this one is the analytic
equation of state and needs no data at all.

P(x) and rho(x) come from scf/eos.py -- the same module the SCF solver and
the Castro model writer call -- so what is plotted is the equation of
state the results were computed with, not a second implementation of the
same formulae that could drift from it.
"""

import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "scf"))

from eos import (A_CONST, B_of_mu_e, density,          # noqa: E402
                 neutronization_threshold_rho_c, pressure)

COL_IN = 88.0 / 25.4
C_FIELD = "#2a78d6"
C_SEED = "#52514e"
C_ALT = "#eb6834"
C_MUTED = "#898781"
C_GRID = "#e1e0d9"
C_BAND = "#f2f1ec"

MU_E = 2.0
RHO_C_STAR = 9.883938495e8      # the background star of Sect. 6

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


def main():
    B = B_of_mu_e(MU_E)
    x = np.logspace(-2.0, 2.0, 2000)
    rho = density(x, MU_E)
    P = pressure(x)

    # Gamma = dlnP/dlnrho, analytic rather than differenced: eos.py gives
    # dP/dx = 8A x^4 / sqrt(1+x^2), and drho/dx = 3B x^2, so
    # Gamma = (8A/3) x^5 / (P sqrt(1+x^2)).
    gamma = (8.0 * A_CONST / 3.0) * x**5 / (P * np.sqrt(1.0 + x**2))

    rho_n = neutronization_threshold_rho_c(MU_E)
    x_c = (RHO_C_STAR / B) ** (1.0 / 3.0)
    g_c = float(np.interp(RHO_C_STAR, rho, gamma))

    fig, (ax, axg) = plt.subplots(
        2, 1, figsize=(COL_IN, COL_IN * 1.05), sharex=True,
        gridspec_kw=dict(height_ratios=[1.35, 1.0], hspace=0.10))

    for a in (ax, axg):
        style_axes(a)
        a.axvspan(rho_n, 1e13, color=C_BAND, zorder=0)
        a.axvline(RHO_C_STAR, color=C_FIELD, linewidth=0.7,
                  linestyle=(0, (4, 2)), zorder=2)

    # Each limit normalised on the curve in the regime where it holds, so
    # the dashed lines lie on it there instead of being fitted by eye.
    nr = x < 0.05
    rel = x > 20.0
    ax.plot(rho, P, color="#0b0b0b", linewidth=1.4, zorder=4,
            label="Eqs. 1 and 2")
    ax.plot(rho, P[nr][-1] * (rho / rho[nr][-1]) ** (5.0 / 3.0),
            color=C_ALT, linewidth=0.9, linestyle=(0, (4, 2)), zorder=3,
            label=r"$P \propto \rho^{5/3}$ ($n = 3/2$)")
    ax.plot(rho, P[rel][0] * (rho / rho[rel][0]) ** (4.0 / 3.0),
            color=C_SEED, linewidth=0.9, linestyle=(0, (1, 1.5)), zorder=3,
            label=r"$P \propto \rho^{4/3}$ ($n = 3$)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(P.min() * 0.4, P.max() * 6)
    ax.set_ylabel(r"$P$  (dyn cm$^{-2}$)")
    ax.legend(loc="upper left", handlelength=2.4, borderaxespad=0.4)
    ax.text(rho_n * 1.6, P.min() * 1.5, "electron\ncapture", fontsize=6.5,
            color=C_MUTED, va="bottom")

    axg.plot(rho, gamma, color="#0b0b0b", linewidth=1.4, zorder=4)
    for v, lab in ((5.0 / 3.0, r"$5/3$"), (4.0 / 3.0, r"$4/3$")):
        axg.axhline(v, color=C_MUTED, linewidth=0.7, linestyle=(0, (4, 3)),
                    zorder=1)
        axg.text(4.0, v + 0.012, lab, fontsize=6.5, color=C_MUTED)
    axg.plot([RHO_C_STAR], [g_c], marker="o", markersize=4.0, color=C_FIELD,
             markeredgecolor="white", markeredgewidth=0.5, zorder=5)
    axg.annotate(r"$\rho_c$ of the star of Sect. 6" "\n"
                 rf"$\Gamma = {g_c:.3f}$",
                 xy=(RHO_C_STAR, g_c), xytext=(3.0e3, 1.47), fontsize=6.5,
                 color="#0b0b0b", ha="left",
                 arrowprops=dict(arrowstyle="-", linewidth=0.5,
                                 color=C_MUTED, shrinkA=1, shrinkB=3))
    axg.set_xlabel(r"$\rho$  (g cm$^{-3}$)")
    axg.set_ylabel(r"$\Gamma = \mathrm{d}\ln P/\mathrm{d}\ln\rho$")
    axg.set_ylim(1.29, 1.73)
    axg.set_xlim(rho.min(), 1e13)

    print(f"mu_e = {MU_E}, B = {B:.4e} g/cm^3")
    print(f"x(rho_c) = {x_c:.3f}  ->  Gamma(rho_c) = {g_c:.4f}")
    print(f"neutronization threshold = {rho_n:.4e} g/cm^3")
    fig.savefig(HERE / "eos.pdf")
    plt.close(fig)
    print("wrote eos.pdf")


if __name__ == "__main__":
    main()
