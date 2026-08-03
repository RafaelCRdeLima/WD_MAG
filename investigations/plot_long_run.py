"""The 192^3 run past 12 s: a steady pulsation, and a braking that is easing off.

Run:  scf/.venv/bin/python3 investigations/plot_long_run.py

Reads long_192.csv, the central density and axial angular momentum taken from
the TIME= diagnostics of run_rot192.log. The run is still going; everything
here is provisional to the last time in the file.

Two panels, and each carries a result that the first 12 s alone got wrong.

  (a) rho_max. Fitting the envelope over the first 12 s gave a damping
      e-folding of 12-20 s and the reading that the star was settling. Past
      the transient it is not: the oscillation locks to a period of 1.50 s at
      +-14.4%, and its envelope decays with an e-folding near 140 s. The
      first fit was measuring the initial-condition transient decaying on top
      of a mode that persists.

  (b) L_z, with the two slopes drawn. The braking rate halves once the field
      is gone -- 0.206%/s while the toroidal field lives, 0.099%/s after --
      which is the direct evidence that the transport was magnetic.

The dashed vertical marks t = 12 s, where the continuation restarted from
chk03436. It is a restart, not a discontinuity: the two segments are one
evolution written into one appended log.
"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = Path(__file__).resolve().parent
COL_IN = 88.0 / 25.4

C_RHO = "#2a78d6"
C_LZ = "#eb6834"
C_MUTED = "#898781"
C_GRID = "#e1e0d9"
C_INK = "#0b0b0b"

T_RESTART = 12.0
T_DYN = 0.4752

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
    path = HERE / "long_192.csv"
    with open(path) as fh:
        n = sum(1 for line in fh if line.startswith("#"))
    return np.genfromtxt(path, delimiter=",", names=True, skip_header=n)


def main():
    d = load()
    t, rho, lz = d["t"], d["rho_max"], d["Lz"]

    fig, (ax_r, ax_l) = plt.subplots(
        2, 1, figsize=(COL_IN, COL_IN * 1.25), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 0.85], "hspace": 0.12},
    )
    for ax in (ax_r, ax_l):
        ax.grid(True, color=C_GRID, linewidth=0.4, alpha=0.9)
        ax.set_axisbelow(True)
        ax.set_xlim(0, t.max())
        ax.axvline(T_RESTART, color=C_MUTED, linewidth=0.6, linestyle=(0, (2, 2)))

    # (a) central density -------------------------------------------------
    ax_r.plot(t, rho / 1e9, color=C_RHO, linewidth=1.0)
    ax_r.set_ylabel(r"$\rho_{\rm max}$  ($10^9$ g cm$^{-3}$)")
    ax_r.set_ylim(0.6, 3.2)
    ax_r.text(0.02, 0.93, "(a)", transform=ax_r.transAxes, fontsize=7, va="top")
    ax_r.text(T_RESTART - 0.3, 3.0, "reinicio", fontsize=5.8, color=C_MUTED,
              ha="right", va="top")

    # The band the oscillation settles into, from the peak/trough analysis.
    late = t >= 12.5
    hi, lo = 1.532, 1.145
    ax_r.axhspan(lo, hi, xmin=12.5 / t.max(), color=C_GRID, alpha=0.85, linewidth=0)
    ax_r.annotate(r"$P = 1.50$ s $= 3.15\,t_{\rm dyn}$, $\pm 14.4\%$",
                  xy=(17.0, hi), xytext=(0, 6), textcoords="offset points",
                  fontsize=6.2, color=C_INK, ha="center")

    # (b) angular momentum ------------------------------------------------
    ax_l.plot(t, lz / 1e50, color=C_LZ, linewidth=1.0)
    ax_l.set_ylabel(r"$L_z$  ($10^{50}$ g cm$^2$ s$^{-1}$)")
    ax_l.set_xlabel("t (s)")
    ax_l.set_ylim(1.96, 2.24)
    ax_l.text(0.02, 0.10, "(b)", transform=ax_l.transAxes, fontsize=7, va="bottom")

    # The two braking rates, fitted after the initial spin-up transient and
    # after the restart. Drawn as guide lines, offset so they sit clear of the
    # curve rather than on top of it.
    for lo_t, hi_t, lab, dy in ((1.5, 11.9, r"$-0.206\%$/s", -0.055),
                                (12.1, t.max(), r"$-0.099\%$/s", 0.035)):
        m = (t >= lo_t) & (t <= hi_t)
        s, c = np.polyfit(t[m], lz[m] / 1e50, 1)
        tt = np.linspace(lo_t, hi_t, 20)
        ax_l.plot(tt, s * tt + c, color=C_INK, linewidth=0.7, linestyle=(0, (1, 1.6)))
        ax_l.text(0.5 * (lo_t + hi_t), s * 0.5 * (lo_t + hi_t) + c + dy, lab,
                  fontsize=6.2, color=C_INK, ha="center")

    fig.savefig(HERE / "long_run_192.pdf")
    fig.savefig(HERE / "long_run_192.png", dpi=200)
    print(f"t ate {t.max():.2f} s; rho_max media em t>12.5: {rho[late].mean():.3e}")
    print(f"wrote {HERE/'long_run_192.pdf'}")


if __name__ == "__main__":
    main()
