"""Figure for the braking sequence: where a 2 Msun WD crosses the window.

Run:  scf/.venv/bin/python3 investigations/plot_braking_sequence.py

Reads braking_sequence.csv. Two panels sharing the density axis: how much
magnetic support the star needs, and what peak field that support implies
against the Landau critical field. The second panel is the point of the
figure -- the window where the interesting nuclear physics happens sits
entirely above B_c for a space-filling toroidal field.
"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = Path(__file__).resolve().parent
COL_IN = 88.0 / 25.4

C_FIELD = "#2a78d6"
C_ALT = "#eb6834"
C_MUTED = "#898781"
C_GRID = "#e1e0d9"

B_C = 4.414e13
RHO_CC, RHO_NE, RHO_O = 3.0e9, 9.6e9, 1.94e10
CONFINE_GAIN = 2.71     # measured in confinement_cost.py at rho_c = 1e9:
                        # 4.275e13 space-filling vs 1.579e13 confined, same
                        # E_tor/|W|. Applied here as a constant, which it is
                        # only approximately -- the closed-region geometry
                        # drifts along the sequence.

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8,
    "legend.fontsize": 6.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "legend.frameon": False, "figure.dpi": 200,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def main():
    rows = []
    with (HERE / "braking_sequence.csv").open() as f:
        for line in f:
            if line.startswith("#") or line.startswith("rho_c"):
                continue
            p = line.strip().split(",")
            rows.append((float(p[0]), float(p[3]), float(p[4])))
    rho = np.array([r[0] for r in rows])
    emag = np.array([r[1] for r in rows])
    bt = np.array([r[2] for r in rows])

    fig, (ax, axb) = plt.subplots(
        2, 1, figsize=(COL_IN, COL_IN * 1.15), sharex=True,
        gridspec_kw=dict(height_ratios=[1.0, 1.15], hspace=0.10))

    for a in (ax, axb):
        a.set_xscale("log")
        a.grid(True, axis="y", color=C_GRID, linewidth=0.4, zorder=0)
        a.set_axisbelow(True)
        for side in ("top", "right"):
            a.spines[side].set_color(C_MUTED)
        # the three lines the track has to cross, in density order
        for x in (RHO_CC, RHO_NE, RHO_O):
            a.axvline(x, color=C_MUTED, linewidth=0.7, linestyle=(0, (4, 3)),
                      zorder=1)
        a.axvspan(RHO_CC, RHO_O, color="#f4f2ea", zorder=0)

    ax.plot(rho, emag, color=C_FIELD, linewidth=1.4, marker="o",
            markersize=3.4, markeredgecolor="white", markeredgewidth=0.5,
            zorder=3)
    ax.set_ylabel(r"$E_{\rm mag}/|W|$")
    ax.set_ylim(0.15, 0.225)
    ax.text(RHO_CC, 0.219, " C+C", fontsize=6.2, color=C_MUTED, ha="left")
    ax.text(RHO_NE, 0.219, r" $^{20}$Ne", fontsize=6.2, color=C_MUTED,
            ha="left")
    ax.text(RHO_O, 0.219, r" $^{16}$O", fontsize=6.2, color=C_MUTED,
            ha="right")
    ax.text(1.0e9, 0.163, "support varies by 21%\nacross a factor 25 in "
            r"$\rho_c$", fontsize=6.3, color=C_MUTED, va="bottom")

    axb.plot(rho, bt / B_C, color=C_FIELD, linewidth=1.4, marker="o",
             markersize=3.4, markeredgecolor="white", markeredgewidth=0.5,
             zorder=3, label="space-filling")
    axb.plot(rho, bt / B_C / CONFINE_GAIN, color=C_ALT, linewidth=1.2,
             linestyle=(0, (5, 2)), zorder=3, label="confined (estimate)")
    axb.axhline(1.0, color="#0b0b0b", linewidth=0.9, zorder=2)
    axb.text(8.5e8, 1.13, r"$B_c$", fontsize=7, color="#0b0b0b")
    axb.set_yscale("log")
    axb.set_ylim(0.2, 9.0)
    axb.set_ylabel(r"$\max|B_\varphi| \,/\, B_c$")
    axb.set_xlabel(r"$\rho_c$  (g cm$^{-3}$)")
    axb.legend(loc="upper left", fontsize=6.3)

    fig.savefig(HERE / "braking_sequence.pdf")
    plt.close(fig)
    print("wrote braking_sequence.pdf")
    for x, name in ((RHO_CC, "C+C"), (RHO_NE, "20Ne"), (RHO_O, "16O")):
        b = np.interp(np.log10(x), np.log10(rho), bt) / B_C
        e = np.interp(np.log10(x), np.log10(rho), emag)
        print(f"  at {name:5s} ({x:.2e}): E_mag/|W| = {e:.4f}, "
              f"B/B_c = {b:.2f} space-filling, {b / CONFINE_GAIN:.2f} confined")


if __name__ == "__main__":
    main()
