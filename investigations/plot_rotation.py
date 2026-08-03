"""Angular momentum and the survival of differential rotation, 192^3 to 60 s.

Run:  scf/.venv/bin/python3 investigations/plot_rotation.py

Reads rotation_192.csv: t, Lz_star, Om_mean, Om_core, Om_mid, Om_out, as
measured by tools/fbtbp.cpp over the cells above the density cut. Omega is
mass-weighted inside three cylindrical shells -- below 0.15 R_eq, 0.45 to 0.55,
and 0.65 to 0.75 -- all chosen to stay inside the star at every phase of the
1.5 s pulsation. An earlier version used 0.85 to 1.15 R_eq and was empty in a
third of the snapshots, because the star breathes across that radius.

Three panels.

  (a) L_z of the star. Two slopes: the braking is twice as fast while the
      toroidal field lives as it is afterwards, which is what identifies the
      transport as magnetic rather than anything else in the scheme -- there
      is no explicit viscosity in these runs, and no explicit resistivity.

  (b) Omega in the three shells. The initial profile spans a factor of three
      from core to outer shell.

  (c) Omega_out/Omega_core, the differential rotation itself. If the field had
      done its work this would climb towards 1. It does not, and that is why
      the star survives: it never loses the differential rotation that holds
      2 Msun above the Chandrasekhar mass.

The curves carry the 1.498 s pulsation, which aliases badly against any coarse
sampling; the running mean over one period is drawn on top of (c) so the trend
can be read without reading the ripple as signal.
"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = Path(__file__).resolve().parent
COL_IN = 88.0 / 25.4

C_LZ = "#eb6834"
C_CORE = "#2a78d6"
C_MID = "#1baf7a"
C_OUT = "#4a3aa7"
C_MUTED = "#898781"
C_GRID = "#e1e0d9"
C_INK = "#0b0b0b"

P_PULSE = 1.498          # s, the measured pulsation period

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
    path = HERE / "rotation_192.csv"
    with open(path) as fh:
        n = sum(1 for line in fh if line.startswith("#"))
    return np.genfromtxt(path, delimiter=",", names=True, skip_header=n)


def smooth(t, y, window):
    """Running mean over a fixed time window, for unevenly spaced samples."""
    out = np.empty_like(y)
    for i, ti in enumerate(t):
        m = np.abs(t - ti) <= window / 2
        out[i] = y[m].mean()
    return out


def main():
    d = load()
    t = d["t"]

    fig, (ax_l, ax_o, ax_r) = plt.subplots(
        3, 1, figsize=(COL_IN, COL_IN * 1.7), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0, 0.95], "hspace": 0.12},
    )
    for ax in (ax_l, ax_o, ax_r):
        ax.grid(True, color=C_GRID, linewidth=0.4, alpha=0.9)
        ax.set_axisbelow(True)
        ax.set_xlim(0, t.max())

    # (a) angular momentum ------------------------------------------------
    lz = d["Lz_star"] / 1e50
    ax_l.plot(t, lz, color=C_LZ, linewidth=1.0)
    for lo, hi, lab, dy in ((1.5, 11.9, r"$-0.21\%$/s", -0.045),
                            (12.1, t.max(), r"$-0.10\%$/s", 0.035)):
        m = (t >= lo) & (t <= hi)
        if m.sum() < 3:
            continue
        s, c = np.polyfit(t[m], lz[m], 1)
        tt = np.linspace(lo, hi, 20)
        ax_l.plot(tt, s * tt + c, color=C_INK, linewidth=0.7, linestyle=(0, (1, 1.6)))
        ax_l.text(0.5 * (lo + hi), s * 0.5 * (lo + hi) + c + dy, lab,
                  fontsize=6.2, color=C_INK, ha="center")
    ax_l.set_ylabel(r"$L_z$  ($10^{50}$ g cm$^2$ s$^{-1}$)")
    ax_l.text(0.06, 0.10, "(a)", transform=ax_l.transAxes, fontsize=7, va="bottom")

    # (b) the three shells -------------------------------------------------
    # Labels in the gaps between the curves, not on them: the shells are well
    # separated in Omega and that separation is the point of the panel.
    for key, col, lab, ytext in (
            ("Om_core", C_CORE, r"core, $\varpi < 0.15\,R_{\rm eq}$", 8.35),
            ("Om_mid", C_MID, r"mid, $0.45$--$0.55$", 4.35),
            ("Om_out", C_OUT, r"outer, $0.65$--$0.75$", 1.66)):
        ax_o.plot(t, d[key], color=col, linewidth=0.9)
        ax_o.text(28.0, ytext, lab, color=col, fontsize=6.0, ha="center")
    ax_o.set_ylim(1.35, 9.0)
    ax_o.set_ylabel(r"$\Omega$  (rad s$^{-1}$)")
    ax_o.text(0.02, 0.92, "(b)", transform=ax_o.transAxes, fontsize=7, va="top")

    # (c) the differential rotation ----------------------------------------
    ratio = d["Om_out"] / d["Om_core"]
    ok = np.isfinite(ratio) & (d["Om_core"] > 0) & (d["Om_out"] > 0)
    ax_r.plot(t[ok], ratio[ok], color=C_MUTED, linewidth=0.6)
    ax_r.plot(t[ok], smooth(t[ok], ratio[ok], P_PULSE), color=C_OUT, linewidth=1.3)
    ax_r.axhline(1.0, color=C_MUTED, linewidth=0.6, linestyle=(0, (3, 2)))
    ax_r.text(t.max() * 0.98, 1.02, "uniform rotation", fontsize=5.8,
              color=C_MUTED, ha="right", va="bottom")
    ax_r.annotate("", xy=(58, 0.266), xytext=(58, 0.40),
                  arrowprops=dict(arrowstyle="->", color=C_INK, linewidth=0.6))
    ax_r.text(56.5, 0.44, "more\ndifferential", fontsize=5.8, color=C_INK,
              ha="right", va="bottom")
    ax_r.set_ylabel(r"$\Omega_{\rm out}/\Omega_{\rm core}$")
    ax_r.set_xlabel("t (s)")
    ax_r.set_ylim(0, 1.15)
    ax_r.text(0.02, 0.92, "(c)", transform=ax_r.transAxes, fontsize=7, va="top")

    fig.savefig(HERE / "rotation_192.pdf")
    fig.savefig(HERE / "rotation_192.png", dpi=200)

    r0, r1 = ratio[ok][:5].mean(), ratio[ok][-5:].mean()
    print(f"Lz_star: {d['Lz_star'][0]:.4e} -> {d['Lz_star'][-1]:.4e} "
          f"({100*(d['Lz_star'][-1]/d['Lz_star'][0]-1):+.2f}%)")
    print(f"Om_out/Om_core: {r0:.3f} -> {r1:.3f}  ({100*(r1/r0-1):+.1f}%)")
    print(f"wrote {HERE/'rotation_192.pdf'}")


if __name__ == "__main__":
    main()
