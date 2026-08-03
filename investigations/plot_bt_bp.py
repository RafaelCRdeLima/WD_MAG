"""Figure for the 192^3 run: the toroidal field is destroyed by its own instability.

Run:  scf/.venv/bin/python3 investigations/plot_bt_bp.py

Reads bt_bp_192.csv. Three panels sharing the time axis, in the order the
argument has to be made:

  (a) the energies, so the growth is seen against the t = 0 interpolation
      floor -- without that floor drawn, an exponential in E_pol is not
      evidence of anything, because the poloidal field on the grid starts as
      interpolation error;
  (b) the two ratios, toroidal over poloidal, one in energy and one in peak
      amplitude, with equipartition marked -- they saturate at different
      places and the difference is the point;
  (c) the peak field against the Landau critical field, which is where the
      run leaves the range its own zero-temperature EOS is valid in.

Separate panels rather than twin axes: three quantities with three scales on
one frame would be unreadable and the comparison it invites is false.
"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = Path(__file__).resolve().parent
COL_IN = 88.0 / 25.4

C_TOR = "#2a78d6"
C_POL = "#eb6834"
C_FIELD = "#4a3aa7"
C_MUTED = "#898781"
C_GRID = "#e1e0d9"
C_INK = "#0b0b0b"

B_C = 4.414e13
T_A = 2.3829           # Alfven time, from inputs.rot192
FIT_LO, FIT_HI = 1.3, 4.0

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
    # skip_header rather than comments="#": with names=True genfromtxt takes the
    # FIRST line as the header even when it is commented, and the provenance
    # block at the top of the csv has commas in it.
    path = HERE / "bt_bp_192.csv"
    with open(path) as fh:
        n_comment = sum(1 for line in fh if line.startswith("#"))
    return np.genfromtxt(path, delimiter=",", names=True, skip_header=n_comment)


def growth_fit(t, e_pol):
    """Exponential fit over the growth phase; returns (sigma, e-folding of B)."""
    m = (t >= FIT_LO) & (t <= FIT_HI)
    sigma, intercept = np.polyfit(t[m], np.log(e_pol[m]), 1)
    return sigma, intercept, m


def main():
    d = load()
    t = d["t"]
    sigma, intercept, fitmask = growth_fit(t, d["E_pol"])

    fig, (ax_e, ax_r, ax_b) = plt.subplots(
        3, 1, figsize=(COL_IN, COL_IN * 1.75), sharex=True,
        gridspec_kw={"height_ratios": [1.15, 1.0, 0.85], "hspace": 0.12},
    )

    for ax in (ax_e, ax_r, ax_b):
        ax.grid(True, color=C_GRID, linewidth=0.4, alpha=0.9)
        ax.set_axisbelow(True)
        ax.set_xlim(0, 12)

    # (a) energies -------------------------------------------------------
    ax_e.set_yscale("log")
    ax_e.plot(t, d["E_tor"], color=C_TOR, linewidth=1.1)
    ax_e.plot(t, d["E_pol"], color=C_POL, linewidth=1.1)

    # The floor is the whole reason the growth is credible: E_pol at t = 0 is
    # interpolation error, so it is the noise level this measurement sits on.
    floor = d["E_pol"][0]
    ax_e.axhline(floor, color=C_MUTED, linewidth=0.6, linestyle=(0, (3, 2)))
    ax_e.text(11.8, floor * 1.5, "numerical floor (t = 0)", fontsize=5.8,
              color=C_MUTED, ha="right", va="bottom")

    tf = np.linspace(FIT_LO, FIT_HI + 0.6, 40)
    ax_e.plot(tf, np.exp(intercept + sigma * tf), color=C_INK,
              linewidth=0.7, linestyle=(0, (1, 1.6)))
    ax_e.annotate(rf"$\tau_B = {2/sigma:.2f}$ s $= {2/sigma/T_A:.2f}\,t_A$",
                  xy=(FIT_HI + 0.5, np.exp(intercept + sigma * (FIT_HI + 0.5))),
                  xytext=(3, 2), textcoords="offset points",
                  fontsize=6.2, color=C_INK, ha="left", va="bottom")

    ax_e.text(9.2, 8e48, r"$E_{\rm tor}$", color=C_TOR, fontsize=7.5)
    ax_e.text(9.2, 5.5e46, r"$E_{\rm pol}$", color=C_POL, fontsize=7.5)
    ax_e.set_ylabel("magnetic energy (erg)")
    ax_e.set_ylim(5e43, 3e50)
    ax_e.text(0.02, 0.93, "(a)", transform=ax_e.transAxes, fontsize=7,
              color=C_INK, va="top")

    # (b) ratios ---------------------------------------------------------
    ax_r.set_yscale("log")
    ax_r.plot(t, d["Et_over_Ep"], color=C_TOR, linewidth=1.1)
    ax_r.plot(t, d["Bt_over_Bp_amp"], color=C_POL, linewidth=1.1)
    ax_r.axhline(1.0, color=C_MUTED, linewidth=0.6, linestyle=(0, (3, 2)))
    ax_r.text(0.3, 1.15, "equipartition", fontsize=5.8, color=C_MUTED,
              ha="left", va="bottom")

    ax_r.text(8.0, 40.0, r"$E_{\rm tor}/E_{\rm pol}$", color=C_TOR, fontsize=7)
    ax_r.text(7.6, 2.9, r"$B_{\rm tor}^{\max}/B_{\rm pol}^{\max}$",
              color=C_POL, fontsize=7)
    ax_r.set_ylabel("toroidal / poloidal ratio")
    ax_r.set_ylim(0.7, 1e6)
    ax_r.text(0.02, 0.93, "(b)", transform=ax_r.transAxes, fontsize=7,
              color=C_INK, va="top")

    # (c) peak field vs the critical field --------------------------------
    bbc = d["B_over_Bc"]
    ax_b.axhspan(1.0, 1.7, color=C_GRID, alpha=0.7, linewidth=0)
    ax_b.plot(t, bbc, color=C_FIELD, linewidth=1.1)
    ax_b.axhline(1.0, color=C_MUTED, linewidth=0.6, linestyle=(0, (3, 2)))
    ax_b.text(11.8, 1.44, "unquantised EOS out of range", fontsize=5.8,
              color=C_MUTED, ha="right", va="center")
    ax_b.text(0.35, 0.62, r"$B_c$", fontsize=7, color=C_INK)
    ax_b.set_ylabel(r"$B^{\max}/B_c$")
    ax_b.set_xlabel("t (s)")
    ax_b.set_ylim(0.55, 1.7)
    ax_b.text(0.02, 0.93, "(c)", transform=ax_b.transAxes, fontsize=7,
              color=C_INK, va="top")

    fig.savefig(HERE / "bt_bp_192.pdf")
    fig.savefig(HERE / "bt_bp_192.png", dpi=200)

    frac = np.count_nonzero(bbc > 1.0) / bbc.size
    print(f"sigma_E   = {sigma:.4f} /s   (fit t = {FIT_LO}-{FIT_HI} s)")
    print(f"tau_B     = {2/sigma:.3f} s = {2/sigma/T_A:.2f} t_A")
    print(f"E_pol max = {d['E_pol'].max():.3e} erg at t = {t[np.argmax(d['E_pol'])]:.2f} s"
          f"  ({d['E_pol'].max()/d['E_pol'][0]:.0f}x the floor)")
    print(f"E_tor     = {d['E_tor'][0]:.3e} -> {d['E_tor'][-1]:.3e} erg"
          f"  ({d['E_tor'][-1]/d['E_tor'][0]:.2%} left)")
    print(f"above B_c : {frac:.0%} of the snapshots, peak {bbc.max():.2f} B_c")
    print(f"wrote {HERE/'bt_bp_192.pdf'}")


if __name__ == "__main__":
    main()
