"""192^3 against 256^3 over the whole 12 s: what converges and what does not.

Run:  scf/.venv/bin/python3 investigations/plot_convergence.py

Reads bt_bp_192.csv, bt_bp_256_long.csv and rotation_192.csv. An earlier
version of this figure stopped at t = 5.14 s, because that was as far as the
256^3 run had got; it now reaches stop_time = 12 s on both grids, so the
comparison covers growth, saturation AND decay rather than growth alone.

The verdict splits cleanly in two, and the split is the point of the figure.

  CONVERGED -- the star. Peak field, volume, mass and above all the rotation
  agree to a few percent between the grids. Averaged over 9-12 s (two
  pulsation periods, which removes the +-14% breathing that makes any single
  snapshot meaningless): Omega_core 6.65 against 6.73 rad/s, Omega_out 2.29
  against 2.36, and their ratio 0.345 against 0.350 -- 1.4% apart. The result
  that the star keeps its differential rotation is therefore a result and not
  a resolution.

  NOT CONVERGED -- the field, in energy and in rate. The finer grid grows
  faster (sigma 1.37 -> 2.01/s) and decays slower (gamma 0.282 -> 0.141/s),
  so it holds about twice the magnetic energy at t = 12. The two rates move
  in OPPOSITE directions, which is what identifies both as numerical: there
  is no explicit resistivity in these runs, so the decay IS the numerical
  resistivity, and halving it while dx falls by 1.33 puts it at dx^2.4 --
  consistent with the second-order scheme and extrapolating to zero.

Four panels.

  (a) E_tor and E_pol. Solid is toroidal, dashed poloidal. The fitted rates
      are drawn over their fit windows.
  (b) V/V_0. The first peak falls from 2.51 to 1.79; the EXCESS over unity
      falls as dx^2.3, so most of the initial breathing is the initial
      condition ringing on the grid rather than a physical pulsation. What
      survives is the late oscillation, which both grids settle into at the
      same amplitude.
  (c) Omega_out/Omega_core. The convergence result. Neither grid climbs
      towards 1; both sit near 0.35 and if anything become slightly MORE
      differential.
  (d) B_max/B_c. The field is above the Landau critical field for about half
      the run at both resolutions, which is a caveat on the microphysics, not
      a numerical artefact -- it converges.

Mass is not plotted because there is nothing to see, which is itself the
result: 2.0028 against 2.0042 Msun averaged over 9-12 s.
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

FIT_LO, FIT_HI = 1.3, 4.0        # the exponential growth of E_pol
DEC_LO, DEC_HI = 6.0, 12.0       # the decay of E_tor, after saturation
T_MAX = 12.0
P_PULSE = 1.498

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


def smooth(t, y, window):
    """Running mean over a fixed time window, for unevenly spaced samples."""
    return np.array([y[np.abs(t - ti) <= window / 2].mean() for ti in t])


def main():
    d192, d256 = load("bt_bp_192.csv"), load("bt_bp_256_long.csv")
    r192 = load("rotation_192.csv")
    # rotation at 192^3 lives in its own file (it runs to 60 s); at 256^3 it
    # is in the same table. Restrict the former to the common window.
    runs = [("192$^3$", d192, r192, C_192), ("256$^3$", d256, d256, C_256)]

    fig, (ax_e, ax_v, ax_o, ax_b) = plt.subplots(
        4, 1, figsize=(COL_IN, COL_IN * 2.25), sharex=True,
        gridspec_kw={"height_ratios": [1.4, 0.9, 0.75, 0.8], "hspace": 0.12},
    )
    for ax in (ax_e, ax_v, ax_o, ax_b):
        ax.grid(True, color=C_GRID, linewidth=0.4, alpha=0.9)
        ax.set_axisbelow(True)
        ax.set_xlim(0, T_MAX)

    # (a) the two energies, with the growth and decay rates ------------------
    ax_e.set_yscale("log")
    for lab, d, _, c in runs:
        m = d["t"] <= T_MAX + 0.01
        ax_e.plot(d["t"][m], d["E_tor"][m], color=c, linewidth=1.1)
        ax_e.plot(d["t"][m], d["E_pol"][m], color=c, linewidth=1.0,
                  linestyle=(0, (2.5, 1.5)))

        f = (d["t"] >= FIT_LO) & (d["t"] <= FIT_HI)
        s, b = np.polyfit(d["t"][f], np.log(d["E_pol"][f]), 1)
        tt = np.linspace(FIT_LO, FIT_HI + 0.3, 20)
        ax_e.plot(tt, np.exp(b + s * tt), color=C_INK, linewidth=0.6,
                  linestyle=(0, (1, 1.6)))
        # at the MIDDLE of the fit window, not its end: the two guide lines
        # end where E_tor runs, and a label there lands on top of it
        tm = 0.5 * (FIT_LO + FIT_HI)
        dy = 15 if c == C_256 else -15
        ax_e.annotate(rf"$\sigma = {s:.2f}$/s",
                      xy=(tm, np.exp(b + s * tm)),
                      xytext=(0, dy), textcoords="offset points",
                      fontsize=6.0, color=c, ha="center", va="center")

        g = (d["t"] >= DEC_LO) & (d["t"] <= DEC_HI)
        sd, bd = np.polyfit(d["t"][g], np.log(d["E_tor"][g]), 1)
        tt = np.linspace(DEC_LO - 0.3, DEC_HI, 20)
        ax_e.plot(tt, np.exp(bd + sd * tt), color=C_INK, linewidth=0.6,
                  linestyle=(0, (1, 1.6)))
        dy = 9 if c == C_256 else -9
        ax_e.annotate(rf"$\gamma = {sd:.2f}$/s",
                      xy=(DEC_LO + 1.0, np.exp(bd + sd * (DEC_LO + 1.0))),
                      xytext=(0, dy), textcoords="offset points",
                      fontsize=6.0, color=c, ha="center", va="center")

    ax_e.text(6.2, 1.5e50, r"$E_{\rm tor}$", fontsize=6.5, color=C_MUTED, ha="center")
    ax_e.text(9.0, 3.2e46, r"$E_{\rm pol}$", fontsize=6.5, color=C_MUTED, ha="center")
    ax_e.text(0.35, 2.0e45, "192$^3$", color=C_192, fontsize=7.5)
    ax_e.text(0.35, 2.6e44, "256$^3$", color=C_256, fontsize=7.5)
    ax_e.set_ylabel(r"$E$  (erg)")
    ax_e.set_ylim(3e43, 6e50)
    ax_e.text(0.02, 0.96, "(a)", transform=ax_e.transAxes, fontsize=7, va="top")

    # (b) volume -------------------------------------------------------------
    for lab, d, _, c in runs:
        m = d["t"] <= T_MAX + 0.01
        ax_v.plot(d["t"][m], (d["R_vol_e8"][m] / d["R_vol_e8"][0]) ** 3,
                  color=c, linewidth=1.0)
    ax_v.axhline(1.0, color=C_MUTED, linewidth=0.6, linestyle=(0, (3, 2)))
    ax_v.annotate(r"$\times 2.51$", xy=(1.81, 2.51), xytext=(4, 0),
                  textcoords="offset points", fontsize=6.2, color=C_192, va="center")
    ax_v.annotate(r"$\times 1.79$", xy=(1.74, 1.79), xytext=(-5, 8),
                  textcoords="offset points", fontsize=6.2, color=C_256,
                  ha="right", va="center")
    ax_v.set_ylabel(r"$V/V_0$")
    ax_v.set_ylim(0.5, 2.9)
    ax_v.text(0.02, 0.94, "(b)", transform=ax_v.transAxes, fontsize=7, va="top")

    # (c) the differential rotation -- the panel that converges ---------------
    for lab, _, r, c in runs:
        m = r["t"] <= T_MAX + 0.05
        q = (r["Om_out"] / r["Om_core"])[m]
        tq = r["t"][m]
        ok = np.isfinite(q) & (q > 0)
        ax_o.plot(tq[ok], q[ok], color=c, linewidth=0.5, alpha=0.55)
        ax_o.plot(tq[ok], smooth(tq[ok], q[ok], P_PULSE), color=c, linewidth=1.2)
    ax_o.axhline(1.0, color=C_MUTED, linewidth=0.6, linestyle=(0, (3, 2)))
    ax_o.text(T_MAX * 0.98, 1.02, "uniform rotation", fontsize=5.8,
              color=C_MUTED, ha="right", va="bottom")
    # No arrow towards "more differential" here: over these 12 s the ratio
    # does not fall, it drifts up by a few percent on both grids. The 192^3
    # run does turn downwards, but only after t = 20 s -- see rotation_192.
    ax_o.text(6.2, 0.44, r"both grids: $0.345$ against $0.350$ over $9$--$12$ s",
              fontsize=5.8, color=C_INK, ha="center", va="bottom")
    ax_o.set_ylabel(r"$\Omega_{\rm out}/\Omega_{\rm core}$")
    ax_o.set_ylim(0, 1.18)
    ax_o.text(0.02, 0.94, "(c)", transform=ax_o.transAxes, fontsize=7, va="top")

    # (d) the peak field against the Landau critical field --------------------
    for lab, d, _, c in runs:
        m = d["t"] <= T_MAX + 0.01
        ax_b.plot(d["t"][m], d["B_over_Bc"][m], color=c, linewidth=1.0)
    ax_b.axhline(1.0, color=C_MUTED, linewidth=0.6, linestyle=(0, (3, 2)))
    ax_b.text(11.7, 1.04, r"$B_c$", fontsize=6.2, color=C_MUTED,
              ha="right", va="bottom")
    ax_b.set_ylabel(r"$B^{\max}/B_c$")
    ax_b.set_xlabel("t (s)")
    ax_b.set_ylim(0.55, 2.05)
    ax_b.text(0.02, 0.94, "(d)", transform=ax_b.transAxes, fontsize=7, va="top")

    fig.savefig(HERE / "convergence_192_256.pdf")
    fig.savefig(HERE / "convergence_192_256.png", dpi=200)

    # the numbers quoted in the docstring and in the report
    print(f"{'':8s} {'sigma':>8s} {'gamma':>8s} {'Vmax':>7s} "
          f"{'<E_tor>':>10s} {'<E_pol>':>10s} {'<Om_c>':>7s} {'<Om_o>':>7s} {'<ratio>':>8s}")
    for lab, d, r, c in runs:
        f = (d["t"] >= FIT_LO) & (d["t"] <= FIT_HI)
        s, _ = np.polyfit(d["t"][f], np.log(d["E_pol"][f]), 1)
        g = (d["t"] >= DEC_LO) & (d["t"] <= DEC_HI)
        sd, _ = np.polyfit(d["t"][g], np.log(d["E_tor"][g]), 1)
        v = (d["R_vol_e8"] / d["R_vol_e8"][0]) ** 3
        late = (d["t"] >= 9.0) & (d["t"] <= 12.01)
        lr = (r["t"] >= 9.0) & (r["t"] <= 12.05)
        print(f"{lab:8s} {s:8.3f} {sd:8.3f} {v.max():7.2f} "
              f"{d['E_tor'][late].mean():10.3e} {d['E_pol'][late].mean():10.3e} "
              f"{r['Om_core'][lr].mean():7.3f} {r['Om_out'][lr].mean():7.3f} "
              f"{(r['Om_out'][lr] / r['Om_core'][lr]).mean():8.4f}")
    print(f"wrote {HERE/'convergence_192_256.pdf'}")


if __name__ == "__main__":
    main()
