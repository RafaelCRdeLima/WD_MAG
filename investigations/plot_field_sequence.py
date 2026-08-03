"""Time sequence of the toroidal and poloidal field, in the meridional plane.

Run:  scf/.venv/bin/python3 investigations/plot_field_sequence.py [192|256]

Reads slices/*_eq.txt and slices/*_mer.txt (192^3), or the slices256/
equivalents, written by tools/fslice.cpp.

EACH COMPONENT IN THE PLANE WHERE IT IS IN-PLANE. This is the whole design and
an earlier version got it wrong: both rows were meridional, so the in-plane
field on both was the POLOIDAL one, and the B_phi row ended up carrying
dipolar-looking lines over a toroidal colour map. A toroidal field drawn that
way does not look toroidal.

  top     equatorial cut, z = 0. Here B_phi is IN the plane, so the streamlines
          of (B_x, B_y) are the toroidal field lines themselves -- circles
          while the star is axisymmetric, and their departure from circles is
          the instability. Colour is B_phi.
  bottom  meridional cut, y = 0. Here the in-plane (B_x, B_z) IS the poloidal
          field, and its streamlines are the poloidal loops. Colour is |B_pol|.

Two rows over the same times.

B_phi is signed and gets a diverging scale symmetric about zero, so grey means
no field rather than mid-range. |B_pol| is a magnitude and gets one hue, light
to dark. Both start at the same 1e10 G, because two maps meant to be compared
have to share a threshold in gauss.

The streamlines are integrated from B, not drawn as contours of a flux
function -- there is no flux function once axisymmetry goes, which is the
result.

Both rows share their colour scale across all times and across the two
resolutions, so panels can be compared by eye. That is the point of fixing
the limits rather than letting each panel autoscale.
"""

from pathlib import Path
import re
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm, SymLogNorm  # noqa: E402

HERE = Path(__file__).resolve().parent
FULL_IN = 180.0 / 25.4

C_NEG, C_MID, C_POS = "#2a78d6", "#f0efec", "#e34948"    # diverging pair + grey
POL_RAMP = ["#efeaf7", "#8c7fd0", "#4a3aa7", "#241a54"]  # one hue, light to dark
C_MUTED, C_INK, C_LINE = "#898781", "#0b0b0b", "#52514e"
C_STREAM_TOR = "#6b6a66"   # lighter, so it does not fight the diverging map

RHO_SURFACE = 1.0e6
B_SCALE = 7.5e13          # covers the peak reached at 256^3 (6.2e13 G)
B_LINTHRESH = 1.0e10

# The two rows share a threshold, and it has to be the same NUMBER of gauss or
# they cannot be compared. An earlier version had the poloidal map start at
# 1e9 while the toroidal one went grey below 1e11: the poloidal panel then gave
# colour to field a hundred times weaker than anything the toroidal panel
# showed, both rows painted the same stellar silhouette, and the figure read as
# if the two components were one map in two colour schemes. They are not --
# inside the star the median |B_tor|/|B_pol| is 1878 at t = 0, 25 at t = 5.4
# and 9 at t = 21.
#
# 1e10 rather than 1e11, because matching them at the higher value emptied the
# poloidal row almost completely and hid the growth that is half the result.
# At 1e10 the t = 0 poloidal panel is still blank -- the field there really is
# 5.5e8, below the floor -- while the later panels show it filling in.
BPOL_LO, BPOL_HI = B_LINTHRESH, B_SCALE

RUNS = {
    # 192^3 spans the whole run, so it gets the two late instants the 256^3
    # has not reached: the field is essentially gone by then and the star is
    # in the steady 1.5 s pulsation.
    # 192^3 spans the whole run and gets the late instants the 256^3 has not
    # reached: by t = 18-21 the field is mostly gone and the star is in the
    # steady 1.5 s pulsation.
    "192": ("slices",    ["plt00000", "plt00400", "plt00700", "plt01150",
                          "plt01550", "plt03436", "plt05012", "plt05805"]),
    "256": ("slices256", ["plt00000", "plt00550", "plt00950",
                          "plt01700", "plt02150"]),
}

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

CMAP_TOR = LinearSegmentedColormap.from_list("bphi", [C_NEG, C_MID, C_POS])
CMAP_POL = LinearSegmentedColormap.from_list("bpol", POL_RAMP)


def read(path):
    with open(path) as fh:
        hdr = [fh.readline() for _ in range(6)]
    t = float(re.search(r"= ([0-9.eE+-]+) s", hdr[1]).group(1))
    d = np.loadtxt(path)
    n = int(round(np.sqrt(d.shape[0])))
    r = [d[:, i].reshape(n, n) for i in range(6)]
    return t, r[0], r[1], r[2], r[3], r[4], r[5]     # x, z, rho, Bx, By, Bz


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "192"
    subdir, plots = RUNS[which]
    paths = [HERE / subdir / f"{p}_mer.txt" for p in plots]
    missing = [p for p in paths if not p.exists()]
    if missing:
        raise SystemExit("faltam: " + ", ".join(str(m.name) for m in missing))

    ncol = len(paths)
    fig, axes = plt.subplots(2, ncol, figsize=(FULL_IN, FULL_IN * 2.0 / ncol))
    fig.subplots_adjust(wspace=0.06, hspace=0.06, right=0.885)

    im_t = im_p = None
    for col, path in enumerate(paths):
        # --- top: equatorial, where B_phi is in the plane -------------------
        te, xe, ye, rhoe, bxe, bye, bze = read(Path(str(path).replace("_mer", "_eq")))
        varpi = np.hypot(xe, ye)
        bphi = np.where(varpi > 0, (-ye * bxe + xe * bye) / np.where(varpi > 0, varpi, 1), 0.0)
        xs_e, ys_e = xe / 1e8, ye / 1e8
        im_t = axes[0, col].pcolormesh(
            xs_e, ys_e, bphi, cmap=CMAP_TOR, shading="auto",
            norm=SymLogNorm(linthresh=B_LINTHRESH, vmin=-B_SCALE, vmax=B_SCALE, base=10))
        axes[0, col].streamplot(xs_e[0], ys_e[:, 0], bxe, bye, color=C_STREAM_TOR,
                                linewidth=0.3, density=0.6, arrowsize=0.3,
                                broken_streamlines=False)

        # --- bottom: meridional, where (B_x, B_z) is the poloidal field -----
        t, x, z, rho, bx, by, bz = read(path)
        xs, zs = x / 1e8, z / 1e8
        bpol = np.hypot(bx, bz)
        im_p = axes[1, col].pcolormesh(
            xs, zs, np.maximum(bpol, BPOL_LO), cmap=CMAP_POL, shading="auto",
            norm=LogNorm(vmin=BPOL_LO, vmax=BPOL_HI))
        axes[1, col].streamplot(xs[0], zs[:, 0], bx, bz, color=C_LINE,
                                linewidth=0.3, density=0.6, arrowsize=0.3,
                                broken_streamlines=False)

        for row, cx, cy, rr in ((0, xs_e, ys_e, rhoe), (1, xs, zs, rho)):
            ax = axes[row, col]
            ax.contour(cx, cy, rr, levels=[RHO_SURFACE], colors=[C_INK],
                       linewidths=0.45, linestyles=[(0, (2, 1.5))])
            ax.set_xlim(-5.5, 5.5)
            ax.set_ylim(-5.5, 5.5)
            ax.set_aspect("equal")
            ax.set_xticks([-4, 0, 4])
            ax.set_yticks([-4, 0, 4])
            if col > 0:
                ax.set_yticklabels([])
        axes[0, col].set_xticklabels([])
        axes[0, col].set_title(rf"$t = {t:.2f}$ s", fontsize=7.5, pad=3)
        axes[1, col].set_xlabel(r"$x$  ($10^8$ cm)")

    axes[0, 0].set_ylabel(r"$B_\phi$,  $z=0$" "\n" r"$y$  ($10^8$ cm)")
    axes[1, 0].set_ylabel(r"$|B_{\rm pol}|$,  $y=0$" "\n" r"$z$  ($10^8$ cm)")

    cax_t = fig.add_axes([0.895, 0.545, 0.010, 0.335])
    cb_t = fig.colorbar(im_t, cax=cax_t)
    cb_t.set_label(r"$B_\phi$  (G)", fontsize=7)
    cb_t.ax.tick_params(labelsize=5.5)
    cb_t.outline.set_linewidth(0.6)

    cax_p = fig.add_axes([0.895, 0.135, 0.010, 0.335])
    cb_p = fig.colorbar(im_p, cax=cax_p)
    cb_p.set_label(r"$|B_{\rm pol}|$  (G)", fontsize=7)
    cb_p.ax.tick_params(labelsize=5.5)
    cb_p.outline.set_linewidth(0.6)

    fig.text(0.5, 0.955, rf"${which}^3$", fontsize=9, ha="center", color=C_INK)

    out = HERE / f"field_sequence_{which}"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=200)
    print(f"wrote {out}.pdf  ({ncol} instantes)")


if __name__ == "__main__":
    main()
