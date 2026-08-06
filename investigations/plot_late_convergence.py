"""The 46 s two-grid comparison: what converges late, and what does not.

Run:  scf/.venv/bin/python3 investigations/plot_late_convergence.py

Reads rotation_192.csv (192^3, to 60 s), bt_bp_192.csv (192^3 field, only to
12 s) and bt_bp_256_long.csv (256^3, everything to 58.2 s).

The first report compared the grids at t = 12 s. This figure is the comparison
over the whole evolution, and it separates cleanly into two panels that
converge and one that does not.

  (a) L_z, normalised at t = 12 s. THREE braking regimes, each about half the
      previous, and the two grids agree on all three: -0.200/-0.214 %/s while
      the toroidal field lives, -0.087/-0.094 after it starts to go, and
      -0.045/-0.052 in the residual phase. With no explicit viscosity or
      resistivity in the scheme, nothing else brakes anything, so the halving
      is what identifies the transport as magnetic -- and it now does so at two
      resolutions rather than one.

  (b) Omega_out/Omega_core. The 256^3 onset is about 6 s late; after that it
      steepens faster, and over the full baseline the two land within 10%.
      Measured only to t = 26.5 s the same two curves appear to disagree by a
      factor of eight, which is why the shaded region is drawn: it is the
      window that produced two wrong conclusions in earlier drafts.

  (c) E_mag. The panel that does not converge, and the only one where the
      longer baseline made the disagreement worse rather than better. The
      192^3 field data stops at t = 12 s -- the single point at t = 38.5 s is
      the one late instant that was processed with the field diagnostics, and
      it is drawn as a marker rather than a curve because that is all it is.
      The 256^3 plateau sits a factor 25 above it and is still rising.

Normalising (a) at t = 12 s rather than at t = 0 for the same reason as
elsewhere in this work: before that the star is in the initial transient,
breathing by a factor of 2.5 in volume, and any reference value taken there
depends on which phase of the pulsation it lands in.
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
C_SHADE = "#f0efec"

P_PULSE = 1.498
T_MAX = 58.5
REGIMES = [(1.5, 11.9), (12.1, 30.0), (30.0, 58.2)]

# The one late 192^3 field instant that was processed, from Table 4 of the
# first report. A single snapshot, not a curve, and drawn as such.
T_192_LATE, EMAG_192_LATE = 38.5, 2.782e46

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
    a = np.genfromtxt(path, delimiter=",", names=True, skip_header=n)
    o = np.argsort(a["t"])
    return {k: a[k][o] for k in a.dtype.names}


def smooth(t, y, window=P_PULSE):
    return np.array([y[np.abs(t - ti) <= window / 2].mean() for ti in t])


def main():
    r192, f192 = load("rotation_192.csv"), load("bt_bp_192.csv")
    d256 = load("bt_bp_256_long.csv")

    fig, (ax_l, ax_o, ax_e) = plt.subplots(
        3, 1, figsize=(COL_IN, COL_IN * 1.95), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 0.95, 1.0], "hspace": 0.12},
    )
    for ax in (ax_l, ax_o, ax_e):
        ax.grid(True, color=C_GRID, linewidth=0.4, alpha=0.9)
        ax.set_axisbelow(True)
        ax.set_xlim(0, T_MAX)
        for tb in (12.0, 30.0):
            ax.axvline(tb, color=C_MUTED, linewidth=0.5, linestyle=(0, (2, 2)))

    # (a) angular momentum, with the three braking regimes -------------------
    for lab, d, c in (("192", r192, C_192), ("256", d256, C_256)):
        t, lz = d["t"], d["Lz_star"]
        ref = lz[np.abs(t - 12.0) <= 1.5].mean()
        ax_l.plot(t, lz / ref, color=c, linewidth=1.0)
        for lo, hi in REGIMES:
            m = (t >= lo) & (t <= hi)
            if m.sum() < 4:
                continue
            s, b = np.polyfit(t[m], lz[m] / ref, 1)
            tt = np.linspace(lo, hi, 12)
            ax_l.plot(tt, s * tt + b, color=C_INK, linewidth=0.6,
                      linestyle=(0, (1, 1.6)))
    for (lo, hi), lab, yy in zip(REGIMES,
                                 (r"$-0.200\,|\,-0.214$", r"$-0.087\,|\,-0.094$",
                                  r"$-0.045\,|\,-0.052$"),
                                 (0.9565, 0.9565, 0.9565)):
        # along the bottom, not the top: the top-left corner belongs to the
        # panel letter and the first label landed on top of it
        ax_l.text(0.5 * (lo + hi), yy, lab, fontsize=5.8, color=C_INK,
                  ha="center", va="bottom")
    ax_l.set_ylabel(r"$L_z\,/\,L_z(12\,\mathrm{s})$")
    ax_l.set_ylim(0.950, 1.040)
    ax_l.text(0.02, 0.94, "(a)", transform=ax_l.transAxes, fontsize=7, va="top")

    # (b) the differential rotation ------------------------------------------
    ax_o.axvspan(12.0, 26.5, color=C_SHADE, zorder=0)
    ax_o.text(19.2, 0.266, "the $14$ s window\nthat misled", fontsize=5.4,
              color=C_MUTED, ha="center", va="bottom")
    for lab, d, c in (("192$^3$", r192, C_192), ("256$^3$", d256, C_256)):
        t = d["t"]
        q = d["Om_out"] / d["Om_core"]
        ok = np.isfinite(q) & (q > 0)
        ax_o.plot(t[ok], q[ok], color=c, linewidth=0.45, alpha=0.5)
        ax_o.plot(t[ok], smooth(t[ok], q[ok]), color=c, linewidth=1.25)
    ax_o.text(52.0, 0.2655, "192$^3$", color=C_192, fontsize=6.5, ha="center")
    ax_o.text(52.0, 0.3105, "256$^3$", color=C_256, fontsize=6.5, ha="center")
    ax_o.set_ylabel(r"$\Omega_{\rm out}/\Omega_{\rm core}$")
    ax_o.set_ylim(0.26, 0.395)
    ax_o.text(0.02, 0.94, "(b)", transform=ax_o.transAxes, fontsize=7, va="top")

    # (c) magnetic energy, the panel that does not converge -------------------
    ax_e.set_yscale("log")
    ax_e.plot(f192["t"], f192["E_tor"] + f192["E_pol"], color=C_192, linewidth=1.0)
    ax_e.plot(d256["t"], d256["E_tor"] + d256["E_pol"], color=C_256, linewidth=1.0)
    ax_e.plot([T_192_LATE], [EMAG_192_LATE], marker="o", markersize=3.2,
              color=C_192, markeredgewidth=0)
    ax_e.annotate("192$^3$: one processed\ninstant, $t=38.5$ s",
                  xy=(T_192_LATE, EMAG_192_LATE), xytext=(-4, 10),
                  textcoords="offset points", fontsize=5.4, color=C_192,
                  ha="right", va="bottom")
    ax_e.annotate("", xy=(46, 3.0e46), xytext=(46, 5.6e47),
                  arrowprops=dict(arrowstyle="<->", color=C_INK, linewidth=0.6))
    ax_e.text(47.2, 1.3e47, r"$\times 25$", fontsize=6.2, color=C_INK, va="center")
    ax_e.text(30.0, 2.6e48, "256$^3$ plateau, still rising",
              fontsize=5.6, color=C_256, ha="center")
    ax_e.set_ylabel(r"$E_{\rm mag}$  (erg)")
    ax_e.set_xlabel("t (s)")
    ax_e.set_ylim(8e45, 2e50)
    ax_e.text(0.02, 0.94, "(c)", transform=ax_e.transAxes, fontsize=7, va="top")

    fig.savefig(HERE / "late_convergence.pdf")
    fig.savefig(HERE / "late_convergence.png", dpi=200)

    print(f"{'regime':>16s} {'192^3':>10s} {'256^3':>10s}")
    for lo, hi in REGIMES:
        out = f"{lo:5.1f}-{hi:5.1f} s |"
        for d in (r192, d256):
            m = (d["t"] >= lo) & (d["t"] <= hi)
            s, _ = np.polyfit(d["t"][m], d["Lz_star"][m], 1)
            out += f" {100 * s / d['Lz_star'][m].mean():+9.4f}"
        print(out)
    print(f"wrote {HERE/'late_convergence.pdf'}")


if __name__ == "__main__":
    main()
