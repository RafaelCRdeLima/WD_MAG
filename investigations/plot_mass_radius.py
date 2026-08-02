"""Mass and shape of the 192^3 star over 12 s: does it keep what it has?

Run:  scf/.venv/bin/python3 investigations/plot_mass_radius.py

Reads bt_bp_192.csv, the columns tools/fbtbp.cpp measures over the cells above
rho_cut -- the star, not the box. The log's own MASS diagnostic answers a
different question: it integrates the whole domain, ambient included, so it
tests the scheme's conservation rather than the star's.

Three panels on the shared time axis of the energy figure, so the two can be
read side by side.

  (a) mass, against a +-1% band. The point of the panel is that the curve
      stays inside a band this narrow while the magnetic energy of the same
      star falls by 95%.
  (b) the radii. R_vol is drawn heavy and the two extents light, because
      R_eq and R_pol are set by the outermost single cell above the cut and
      a filament moves them; R_vol comes from the volume and cannot be moved
      by one cell. The grey band is one cell, 9.375e6 cm, around R_pol --
      that is 6% of it, and it is why no structure smaller than that in the
      R_pol curve means anything.
  (c) the oblateness, with the model's 0.383 marked. The measured value at
      t = 0 is 0.425 and the difference is discretisation, not evolution:
      R_pol lands on a cell boundary and R_eq is cut short by the density
      threshold biting into a shallow equatorial profile.
"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = Path(__file__).resolve().parent
COL_IN = 88.0 / 25.4

C_MASS = "#2a78d6"
C_VOL = "#eb6834"
C_EXT = "#898781"
C_MUTED = "#898781"
C_GRID = "#e1e0d9"
C_INK = "#0b0b0b"

DX_192 = 9.375e6 / 1.0e8      # one cell, in the 1e8 cm units of the file
OBLATE_MODEL = 0.3832
M_MODEL = 2.00515

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


def load():
    path = HERE / "bt_bp_192.csv"
    with open(path) as fh:
        n_comment = sum(1 for line in fh if line.startswith("#"))
    return np.genfromtxt(path, delimiter=",", names=True, skip_header=n_comment)


def main():
    d = load()
    t = d["t"]

    fig, (ax_m, ax_r, ax_o) = plt.subplots(
        3, 1, figsize=(COL_IN, COL_IN * 1.7), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.25, 0.85], "hspace": 0.12},
    )
    for ax in (ax_m, ax_r, ax_o):
        ax.grid(True, color=C_GRID, linewidth=0.4, alpha=0.9)
        ax.set_axisbelow(True)
        ax.set_xlim(0, 12)

    # (a) mass -----------------------------------------------------------
    m0 = d["M_Msun"][0]
    ax_m.axhspan(m0 * 0.99, m0 * 1.01, color=C_GRID, alpha=0.8, linewidth=0)
    ax_m.plot(t, d["M_Msun"], color=C_MASS, linewidth=1.1)
    ax_m.axhline(M_MODEL, color=C_MUTED, linewidth=0.6, linestyle=(0, (3, 2)))
    ax_m.text(1.15, M_MODEL - 0.0018, "modelo", fontsize=5.8, color=C_MUTED,
              ha="left", va="top")
    ax_m.text(11.8, m0 * 1.0104, r"$\pm 1\%$ do valor inicial", fontsize=5.8,
              color=C_MUTED, ha="right", va="bottom")
    ax_m.set_ylabel(r"$M_\star$  ($M_\odot$)")
    ax_m.set_ylim(m0 * 0.985, m0 * 1.015)
    ax_m.text(0.02, 0.92, "(a)", transform=ax_m.transAxes, fontsize=7, va="top")

    # (b) radii ----------------------------------------------------------
    ax_r.fill_between(t, d["R_pol_e8"] - DX_192 / 2, d["R_pol_e8"] + DX_192 / 2,
                      color=C_GRID, alpha=0.9, linewidth=0)
    ax_r.plot(t, d["R_eq_e8"], color=C_EXT, linewidth=0.8)
    ax_r.plot(t, d["R_pol_e8"], color=C_EXT, linewidth=0.8)
    ax_r.plot(t, d["R_vol_e8"], color=C_VOL, linewidth=1.3)

    ax_r.text(9.0, 4.45, r"$R_{\rm eq}$", color=C_INK, fontsize=7)
    ax_r.text(9.0, 3.55, r"$R_{\rm vol}$", color=C_VOL, fontsize=7.5)
    ax_r.text(9.0, 1.62, r"$R_{\rm pol}$", color=C_INK, fontsize=7)
    ax_r.annotate("uma celula", xy=(6.0, 2.06 + DX_192 / 2), xytext=(6.0, 2.62),
                  fontsize=5.6, color=C_MUTED, ha="center",
                  arrowprops=dict(arrowstyle="-", color=C_MUTED, linewidth=0.5))
    ax_r.set_ylabel(r"raio  ($10^8$ cm)")
    ax_r.set_ylim(1.2, 5.6)
    ax_r.text(0.02, 0.93, "(b)", transform=ax_r.transAxes, fontsize=7, va="top")

    # (c) oblateness ------------------------------------------------------
    ax_o.plot(t, d["Rpol_over_Req"], color=C_MASS, linewidth=1.1)
    ax_o.axhline(OBLATE_MODEL, color=C_MUTED, linewidth=0.6, linestyle=(0, (3, 2)))
    ax_o.text(11.8, OBLATE_MODEL - 0.012, "modelo  0.383", fontsize=5.8,
              color=C_MUTED, ha="right", va="top")
    ax_o.set_ylabel(r"$R_{\rm pol}/R_{\rm eq}$")
    ax_o.set_xlabel("t (s)")
    ax_o.set_ylim(0.33, 0.74)
    ax_o.text(0.02, 0.92, "(c)", transform=ax_o.transAxes, fontsize=7, va="top")

    fig.savefig(HERE / "mass_radius_192.pdf")
    fig.savefig(HERE / "mass_radius_192.png", dpi=200)

    M = d["M_Msun"]
    print(f"massa: {M.min():.5f} a {M.max():.5f} Msun, amplitude {100*(M.max()-M.min())/M[0]:.3f}%")
    print(f"R_vol: {d['R_vol_e8'].min():.3f} a {d['R_vol_e8'].max():.3f} e8 cm")
    print(f"wrote {HERE/'mass_radius_192.pdf'}")


if __name__ == "__main__":
    main()
