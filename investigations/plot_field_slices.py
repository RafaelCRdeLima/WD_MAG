"""The field of the 192^3 run losing its axisymmetry.

Run:  scf/.venv/bin/python3 investigations/plot_field_slices.py

Reads slices/plt*_mer.txt and slices/plt*_eq.txt, written by tools/fslice.cpp.

Two rows, four times. The top row is the meridional cut, where at y = 0 and
x > 0 the out-of-plane component B_y IS B_phi and the in-plane (B_x, B_z) IS
the poloidal field, so nothing has to be rotated to be read. The bottom row is
the equatorial cut, which is where a non-axisymmetric mode shows itself: an
m = 1 kink displaces the toroidal column sideways and the colour map goes
lopsided.

Colour is B_phi, which is signed, so the map is diverging -- blue and red
poles with a neutral grey midpoint, and the scale is symmetric about zero so
that grey means zero field and not "middle of the range". It is a symmetric
log scale because B_phi spans four decades between the axis and the surface.

Streamlines rather than arrows for the in-plane field: a quiver over 192^2
cells is unreadable, and the in-plane field is what carries the geometry.

The streamlines mean different things in the two rows, and the labels say so.
In the meridional plane the in-plane field IS the poloidal field. In the
equatorial plane B_phi is in-plane, so the streamlines there are the toroidal
field lines seen face-on -- nearly circles while the configuration is
axisymmetric, and their departure from circles is the kink itself.
They are integrated from B itself, not drawn as contours of a flux function.
That distinction is the point of the figure -- the Grad-Shafranov figures in
papers/wd-toroidal-poloidal can use contours because that configuration is
axisymmetric and a poloidal line is exactly a contour of u. Here there is no
flux function after t ~ 1.3 s, which is the result.

The dotted contour is rho = 1e6 g/cm^3, a stand-in for the stellar surface.
"""

from pathlib import Path
import re

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm   # noqa: E402

HERE = Path(__file__).resolve().parent
SLICES = HERE / "slices"
FULL_IN = 180.0 / 25.4        # two-column width

C_NEG = "#2a78d6"             # diverging pair: blue <-> red, grey midpoint
C_POS = "#e34948"
C_MID = "#f0efec"
C_MUTED = "#898781"
C_INK = "#0b0b0b"
C_LINE = "#52514e"

PLOTFILES = ["plt00000", "plt00700", "plt01550", "plt03436"]
RHO_SURFACE = 1.0e6
B_SCALE = 3.3e13              # symmetric colour limit, just above max|B_phi| at t=0
B_LINTHRESH = 1.0e11

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8,
    "legend.fontsize": 6.5, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "legend.frameon": False, "figure.dpi": 200,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

CMAP = LinearSegmentedColormap.from_list("bphi", [C_NEG, C_MID, C_POS])


def read(path):
    """Return (time, X, Y, rho, B0, B1, Bout) on a 2-d grid.

    B0, B1 are the in-plane components in the order the file's axes are given,
    and Bout is the out-of-plane one.
    """
    with open(path) as fh:
        hdr = [fh.readline() for _ in range(6)]
    t = float(re.search(r"= ([0-9.eE+-]+) s", hdr[1]).group(1))
    axes = hdr[5].split()[1:3]                     # e.g. ['x', 'z']
    d = np.loadtxt(path)
    n = int(round(np.sqrt(d.shape[0])))
    c0 = d[:, 0].reshape(n, n)
    c1 = d[:, 1].reshape(n, n)
    rho = d[:, 2].reshape(n, n)
    b = {"x": d[:, 3].reshape(n, n), "y": d[:, 4].reshape(n, n),
         "z": d[:, 5].reshape(n, n)}
    out = ({"x", "y", "z"} - set(axes)).pop()
    return t, c0, c1, rho, b[axes[0]], b[axes[1]], b[out], axes


def panel(ax, path, equatorial):
    t, c0, c1, rho, b0, b1, bout, axes = read(path)

    if equatorial:
        # Out of the equatorial plane B_z is poloidal, so the colour has to be
        # B_phi computed from the in-plane components rather than taken as one
        # of them.
        varpi = np.hypot(c0, c1)
        safe = np.where(varpi > 0, varpi, 1.0)
        colour = np.where(varpi > 0, (-c1 * b0 + c0 * b1) / safe, 0.0)
    else:
        # At y = 0 with x > 0, B_y is exactly B_phi; for x < 0 it is -B_phi,
        # and leaving it signed is what makes the initial condition read as
        # antisymmetric rather than as two different fields.
        colour = bout

    im = ax.pcolormesh(c0 / 1e8, c1 / 1e8, colour, cmap=CMAP, shading="auto",
                       norm=SymLogNorm(linthresh=B_LINTHRESH, vmin=-B_SCALE,
                                       vmax=B_SCALE, base=10))

    speed = np.hypot(b0, b1)
    if speed.max() > 0:
        ax.streamplot(c0[0] / 1e8, c1[:, 0] / 1e8, b0, b1,
                      color=C_LINE, linewidth=0.35, density=0.7,
                      arrowsize=0.35, broken_streamlines=False)

    ax.contour(c0 / 1e8, c1 / 1e8, rho, levels=[RHO_SURFACE],
               colors=[C_INK], linewidths=0.5, linestyles=[(0, (2, 1.5))])

    lim = 5.5
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xticks([-4, 0, 4])
    ax.set_yticks([-4, 0, 4])
    return t, im


def main():
    fig, axes = plt.subplots(2, len(PLOTFILES),
                             figsize=(FULL_IN, FULL_IN * 0.56))
    fig.subplots_adjust(wspace=0.08, hspace=0.06, right=0.90)

    im = None
    for col, p in enumerate(PLOTFILES):
        t, im = panel(axes[0, col], SLICES / f"{p}_mer.txt", equatorial=False)
        panel(axes[1, col], SLICES / f"{p}_eq.txt", equatorial=True)
        axes[0, col].set_title(rf"$t = {t:.2f}$ s", fontsize=7.5, pad=3)

    for row, label in ((0, r"$z$  ($10^8$ cm)"), (1, r"$y$  ($10^8$ cm)")):
        axes[row, 0].set_ylabel(label)
    for col in range(len(PLOTFILES)):
        axes[1, col].set_xlabel(r"$x$  ($10^8$ cm)")
        axes[0, col].set_xticklabels([])
    for row in range(2):
        for col in range(1, len(PLOTFILES)):
            axes[row, col].set_yticklabels([])

    axes[0, 0].text(0.04, 0.95, "meridional  ($y = 0$)", transform=axes[0, 0].transAxes,
                    fontsize=6.3, color=C_INK, va="top")
    axes[1, 0].text(0.04, 0.95, "equatorial  ($z = 0$)", transform=axes[1, 0].transAxes,
                    fontsize=6.3, color=C_INK, va="top")

    cax = fig.add_axes([0.915, 0.13, 0.012, 0.74])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r"$B_\phi$  (G)", fontsize=7.5)
    cb.ax.tick_params(labelsize=6)
    cb.outline.set_linewidth(0.6)

    # The streamlines are not the same field in the two rows, so the note has
    # to say which is which; per-panel labels collided with the row labels.
    fig.text(0.915, 0.10,
             "linhas de campo no plano:\n"
             r"  cima $\rightarrow$ poloidal" "\n"
             r"  baixo $\rightarrow$ toroidal" "\n\n"
             r"tracejado: $\rho = 10^6$ g cm$^{-3}$",
             fontsize=5.8, color=C_MUTED, ha="left", va="top")

    fig.savefig(HERE / "field_slices_192.pdf")
    fig.savefig(HERE / "field_slices_192.png", dpi=200)
    print(f"wrote {HERE/'field_slices_192.pdf'}")


if __name__ == "__main__":
    main()
