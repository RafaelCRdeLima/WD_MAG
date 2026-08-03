"""192^3 against 256^3: what converges, what does not, and what was numerical.

Run:  scf/.venv/bin/python3 investigations/plot_convergence.py

Reads bt_bp_192.csv and bt_bp_256.csv. The 256^3 run only reaches t = 5.14 in
these files -- it died there on the first attempt and was restarted with a
smaller step -- but that window contains the whole growth phase, which is what
the comparison needs.

Three panels, one per test.

  (a) E_pol. The growth is exponential at both resolutions and the fits are
      clean, but the rate is NOT converged: 1.369/s at 192^3 against 2.007/s
      at 256^3, and the finer grid starts from a SMALLER seed. A rate measured
      at one resolution is not a measurement, and this is why.

  (b) V/V_0. The opposite verdict. The first peak falls from 2.51 to 1.79 and
      the second from an amplitude of 0.49 to 0.11, so the breathing that
      dominates the first seconds is largely the initial condition ringing on
      the grid.

  (c) B_tor/B_pol in peak amplitude, with equipartition marked. The finer grid
      carries the instability further: it crosses 1 at 256^3, meaning the
      poloidal field overtakes the toroidal in peak strength, which never
      happens at 192^3.

Mass is not plotted because there is nothing to see, which is itself the
result: 1.99971-2.00822 Msun at 192^3 and 2.00208-2.00730 at 256^3, both
inside half a percent with no trend.
"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = Path(__file__).resolve().parent
COL_IN = 88.0 / 25.4

C_192 = "#2a78d6"
C_256 = "#eb6834"
C_MUTED = "#898781"
C_GRID = "#e1e0d9"
C_INK = "#0b0b0b"

FIT_LO, FIT_HI = 1.3, 4.0
T_MAX = 5.6          # the 256^3 data stops at 5.14; no point drawing further

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


def load(name):
    path = HERE / name
    with open(path) as fh:
        n = sum(1 for line in fh if line.startswith("#"))
    return np.genfromtxt(path, delimiter=",", names=True, skip_header=n)


def main():
    runs = [("192$^3$", load("bt_bp_192.csv"), C_192),
            ("256$^3$", load("bt_bp_256.csv"), C_256)]

    fig, (ax_e, ax_v, ax_b) = plt.subplots(
        3, 1, figsize=(COL_IN, COL_IN * 1.75), sharex=True,
        gridspec_kw={"height_ratios": [1.15, 0.95, 0.9], "hspace": 0.12},
    )
    for ax in (ax_e, ax_v, ax_b):
        ax.grid(True, color=C_GRID, linewidth=0.4, alpha=0.9)
        ax.set_axisbelow(True)
        ax.set_xlim(0, T_MAX)

    # (a) poloidal energy, with the fitted growth ---------------------------
    ax_e.set_yscale("log")
    for lab, d, c in runs:
        m = d["t"] <= T_MAX
        ax_e.plot(d["t"][m], d["E_pol"][m], color=c, linewidth=1.1)
        f = (d["t"] >= FIT_LO) & (d["t"] <= FIT_HI)
        s, b = np.polyfit(d["t"][f], np.log(d["E_pol"][f]), 1)
        tt = np.linspace(FIT_LO, FIT_HI + 0.4, 20)
        ax_e.plot(tt, np.exp(b + s * tt), color=C_INK, linewidth=0.6,
                  linestyle=(0, (1, 1.6)))
        # the two guide lines nearly meet at the right edge, so the labels are
        # pushed apart vertically rather than sitting on top of each other
        dy = 7 if c == C_256 else -9
        ax_e.annotate(rf"$\sigma = {s:.2f}$/s", xy=(FIT_HI + 0.15, np.exp(b + s*(FIT_HI+0.15))),
                      xytext=(4, dy), textcoords="offset points",
                      fontsize=6.2, color=c, ha="left", va="center")
    ax_e.text(0.55, 6e47, "192$^3$", color=C_192, fontsize=7.5)
    ax_e.text(0.55, 1.4e47, "256$^3$", color=C_256, fontsize=7.5)
    ax_e.set_ylabel(r"$E_{\rm pol}$  (erg)")
    ax_e.set_ylim(4e43, 3e49)
    ax_e.text(0.02, 0.93, "(a)", transform=ax_e.transAxes, fontsize=7, va="top")

    # (b) volume ------------------------------------------------------------
    for lab, d, c in runs:
        m = d["t"] <= T_MAX
        ax_v.plot(d["t"][m], (d["R_vol_e8"][m] / d["R_vol_e8"][0]) ** 3,
                  color=c, linewidth=1.1)
    ax_v.axhline(1.0, color=C_MUTED, linewidth=0.6, linestyle=(0, (3, 2)))
    ax_v.annotate(r"$\times 2.51$", xy=(1.81, 2.51), xytext=(4, 0),
                  textcoords="offset points", fontsize=6.2, color=C_192, va="center")
    ax_v.annotate(r"$\times 1.79$", xy=(1.74, 1.79), xytext=(-6, 9),
                  textcoords="offset points", fontsize=6.2, color=C_256,
                  ha="right", va="center")
    ax_v.set_ylabel(r"$V/V_0$")
    ax_v.set_ylim(0.6, 2.9)
    ax_v.text(0.02, 0.93, "(b)", transform=ax_v.transAxes, fontsize=7, va="top")

    # (c) amplitude ratio ---------------------------------------------------
    ax_b.set_yscale("log")
    for lab, d, c in runs:
        m = d["t"] <= T_MAX
        ax_b.plot(d["t"][m], d["Bt_over_Bp_amp"][m], color=c, linewidth=1.1)
    ax_b.axhline(1.0, color=C_MUTED, linewidth=0.6, linestyle=(0, (3, 2)))
    ax_b.text(0.15, 1.12, "equiparti\u00e7\u00e3o", fontsize=5.8, color=C_MUTED,
              ha="left", va="bottom")
    ax_b.set_ylabel(r"$B_{\rm tor}^{\max}/B_{\rm pol}^{\max}$")
    ax_b.set_xlabel("t (s)")
    ax_b.set_ylim(0.6, 1.2e3)
    ax_b.text(0.02, 0.93, "(c)", transform=ax_b.transAxes, fontsize=7, va="top")

    fig.savefig(HERE / "convergence_192_256.pdf")
    fig.savefig(HERE / "convergence_192_256.png", dpi=200)

    for lab, d, c in runs:
        f = (d["t"] >= FIT_LO) & (d["t"] <= FIT_HI)
        s, _ = np.polyfit(d["t"][f], np.log(d["E_pol"][f]), 1)
        v = (d["R_vol_e8"] / d["R_vol_e8"][0]) ** 3
        print(f"{lab}: sigma={s:.4f}/s  E_pol(0)={d['E_pol'][0]:.3e}  "
              f"V_max={v.max():.2f}  Bt/Bp_min={d['Bt_over_Bp_amp'].min():.3f}")
    print(f"wrote {HERE/'convergence_192_256.pdf'}")


if __name__ == "__main__":
    main()
