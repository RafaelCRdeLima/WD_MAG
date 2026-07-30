"""The field of the adopted configuration, inside the star.

Run:  scf/.venv/bin/python3 papers/wd-toroidal-poloidal/figures/make_field_figures.py

Two figures, and the axisymmetry of this configuration decides the form of
both.

Meridional slice. For an axisymmetric field the poloidal field lines ARE
the contours of the flux function -- exactly, not approximately, since
B_pol = (-u_z, u_varpi)/varpi is perpendicular to grad u. So they are drawn
as contours rather than by integrating streamlines, which removes the
tracing error entirely. B_phi is out of this plane and appears as colour,
with arrows for the in-plane poloidal field.

3D field lines. Integrated in cylindrical coordinates,

    dvarpi/ds = B_varpi/|B|,  dz/ds = B_z/|B|,  dphi/ds = B_phi/(varpi |B|)

rather than by interpolating onto a Cartesian grid: it keeps the field
lines on their flux surfaces by construction, and u along a traced line is
a free check on the integration, reported below.

The field is toroidal-dominated with B_t/B_p = 26 in amplitude, so a line
winds about 26 times around the axis per poloidal circuit -- a tightly
wound helix on a nested torus, which is what the figure has to show.
"""

import sys
import warnings
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt                              # noqa: E402
from mpl_toolkits.mplot3d.art3d import Line3DCollection      # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
for p in (REPO / "scf", REPO / "dashboard"):
    sys.path.insert(0, str(p))
warnings.filterwarnings("ignore")

import diagnostics as diag                                   # noqa: E402
from gradshafranov import solve_gradshafranov                # noqa: E402
from scipy.interpolate import RegularGridInterpolator as RGI  # noqa: E402
from seed import r_guess                                     # noqa: E402
from sweep_worker import _solve_toroidal_certified           # noqa: E402
from terms.toroidal_sc import ToroidalSC                     # noqa: E402

RHO_C, MU_E, K_TOR, K0_POL = 1.0e9, 2.0, 3.245e-3, 1.0e-13
COL_IN = 88.0 / 25.4
WIDE_IN = 180.0 / 25.4
C_MUTED, C_GRID = "#898781", "#e1e0d9"

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8,
    "legend.fontsize": 6.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "legend.frameon": False, "figure.dpi": 200,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

BLUES = ["#f4f8fe", "#cde2fb", "#9ec5f4", "#5598e7", "#256abf", "#184f95",
         "#0d366b"]


def solve():
    res, r, th, _ = _solve_toroidal_certified(
        rho_c=RHO_C, R_guess=r_guess(RHO_C), K_tor=K_TOR, m_tor_sc=1.0,
        rotation=None, mu_e=MU_E, Nr_base=129, Ntheta=129, lmax=16,
        tol=1e-8, max_iter=200)
    rho, H = res["rho"], res["H"]
    omega2 = (r[:, None] * np.sin(th)[None, :]) ** 2
    u = solve_gradshafranov(-4.0 * np.pi * omega2 * rho * K0_POL, r, th, lmax=16)
    Bphi = ToroidalSC(K=K_TOR, m=1.0).B_phi(rho, np.sqrt(omega2))
    Br, Bth = diag.poloidal_field(u, r, th)
    R_eq, R_pol = diag.equatorial_polar_radii(H, r, th)
    return r, th, rho, H, u, Bphi, Br, Bth, max(R_eq, R_pol)


def to_meridional(r, th, fields, rmax, n=401):
    """Resample (r, theta) fields onto a (varpi, z) grid."""
    vp = np.linspace(0.0, rmax, n)
    zz = np.linspace(-rmax, rmax, 2 * n - 1)
    VP, ZZ = np.meshgrid(vp, zz, indexing="ij")
    RR = np.sqrt(VP**2 + ZZ**2)
    TT = np.arccos(np.clip(ZZ / np.maximum(RR, 1e-30), -1.0, 1.0))
    out = [RGI((r, th), f, bounds_error=False, fill_value=0.0)((RR, TT))
           for f in fields]
    return vp, zz, VP, ZZ, out


def main():
    from matplotlib.colors import LinearSegmentedColormap, LogNorm

    r, th, rho, H, u, Bphi, Br, Bth, Rstar = solve()
    sin_t, cos_t = np.sin(th)[None, :], np.cos(th)[None, :]
    Bvp = Br * sin_t + Bth * cos_t          # cylindrical components from
    Bz = Br * cos_t - Bth * sin_t           # the spherical ones
    Req, Rpol = diag.equatorial_polar_radii(H, r, th)
    print(f"R_eq = {Req:.4e}  R_pol = {Rpol:.4e}  R_pol/R_eq = {Rpol/Req:.4f}")

    rmax = 1.05 * Rstar
    vp, zz, VP, ZZ, (rho_m, H_m, u_m, bphi_m, bvp_m, bz_m) = to_meridional(
        r, th, (rho, H, u, Bphi, Bvp, Bz), rmax)
    # the surface is where the enthalpy crosses zero; a density cut at
    # 1e-3 of the peak sits well inside it, because this star carries a
    # very extended low-density envelope
    inside = H_m > 0.0
    dense = rho_m > rho.max() * 1e-3
    S8 = 1e8

    # ---------------- meridional slice ----------------
    cmap = LinearSegmentedColormap.from_list("b", BLUES)
    cmap.set_bad("#ffffff")
    fig, ax = plt.subplots(figsize=(COL_IN, COL_IN * 1.02))
    bp = np.ma.masked_where(~inside | (np.abs(bphi_m) <= 0), np.abs(bphi_m))
    vmax = float(np.percentile(np.abs(bphi_m)[inside], 99.5))
    mesh = ax.pcolormesh(VP / S8, ZZ / S8, bp, cmap=cmap,
                         norm=LogNorm(vmin=vmax / 1e3, vmax=vmax),
                         shading="gouraud", rasterized=True, zorder=1)
    # poloidal field lines = contours of u, exactly
    lev = np.linspace(u_m[inside].min(), u_m[inside].max(), 14)[1:-1]
    ax.contour(VP / S8, ZZ / S8, np.ma.masked_where(~inside, u_m),
               levels=lev, colors="#0b0b0b", linewidths=0.5, zorder=3)
    # poloidal direction, unit arrows: the amplitude spans orders of
    # magnitude and scaled arrows would show only the innermost region
    st = 34
    Vs, Zs = VP[::st, ::st], ZZ[::st, ::st]
    Us, Ws = bvp_m[::st, ::st].copy(), bz_m[::st, ::st].copy()
    nrm = np.hypot(Us, Ws)
    keep = inside[::st, ::st] & (nrm > 0)
    Us = np.where(keep, Us / np.where(nrm > 0, nrm, 1.0), np.nan)
    Ws = np.where(keep, Ws / np.where(nrm > 0, nrm, 1.0), np.nan)
    ax.quiver(Vs / S8, Zs / S8, Us, Ws, color="#eb6834", width=0.006,
              scale=26, pivot="mid", zorder=4)
    ax.contour(VP / S8, ZZ / S8, dense.astype(float), levels=[0.5],
               colors=C_MUTED, linewidths=0.6, linestyles="dashed",
               zorder=5)
    ax.contour(VP / S8, ZZ / S8, inside.astype(float), levels=[0.5],
               colors="#0b0b0b", linewidths=0.9, zorder=6)
    ax.set_aspect("equal")
    ax.set_xlim(0, rmax / S8)
    ax.set_ylim(-rmax / S8, rmax / S8)
    ax.set_xlabel(r"$\varpi$  ($10^8$ cm)")
    ax.set_ylabel(r"$z$  ($10^8$ cm)")
    cb = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(r"$|B_\varphi|$  (G)", fontsize=7.5)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=6.5, width=0.5, length=2, pad=1.5)
    for s in ax.spines.values():
        s.set_color(C_MUTED)
    fig.savefig(HERE / "field_meridional.pdf")
    plt.close(fig)
    print("wrote field_meridional.pdf")

    # ---------------- 3D field lines ----------------
    f_bvp = RGI((vp, zz), bvp_m, bounds_error=False, fill_value=0.0)
    f_bz = RGI((vp, zz), bz_m, bounds_error=False, fill_value=0.0)
    f_bphi = RGI((vp, zz), bphi_m, bounds_error=False, fill_value=0.0)
    f_u = RGI((vp, zz), u_m, bounds_error=False, fill_value=0.0)

    def trace(vp0, z0, nstep=26000, ds=None):
        ds = ds or 0.0025 * Rstar
        w, z, ph = vp0, z0, 0.0
        pts, mags = [], []
        for _ in range(nstep):
            p = np.array([[w, z]])
            bw, bz_, bf = f_bvp(p)[0], f_bz(p)[0], f_bphi(p)[0]
            b = np.sqrt(bw * bw + bz_ * bz_ + bf * bf)
            if b <= 0 or w <= 0:
                break
            pts.append((w * np.cos(ph), w * np.sin(ph), z))
            mags.append(b)
            w += ds * bw / b
            z += ds * bz_ / b
            ph += ds * bf / (w * b)
        return np.array(pts), np.array(mags)

    fig = plt.figure(figsize=(WIDE_IN * 0.62, WIDE_IN * 0.42))
    ax3 = fig.add_subplot(111, projection="3d")
    starts = [(0.28, 0.0), (0.50, 0.0), (0.70, 0.0)]
    norm = None
    for i, (fw, fz) in enumerate(starts):
        line, mag = trace(fw * Rstar, fz * Rstar)
        if len(line) < 10:
            continue
        u0 = f_u(np.array([[fw * Rstar, fz * Rstar]]))[0]
        uu = f_u(np.stack([np.sqrt(line[:, 0]**2 + line[:, 1]**2),
                           line[:, 2]], axis=-1))
        drift = np.max(np.abs(uu - u0)) / max(abs(u0), 1e-30)
        print(f"  line {i}: {len(line)} pts, |u| drift along it = {drift:.2e}")
        segs = np.stack([line[:-1] / S8, line[1:] / S8], axis=1)
        if norm is None:
            norm = plt.Normalize(mag.min(), mag.max())
        lc = Line3DCollection(segs, cmap=LinearSegmentedColormap.from_list(
            "b", BLUES[1:]), norm=norm, linewidths=0.55, zorder=3)
        lc.set_array(0.5 * (mag[:-1] + mag[1:]))
        ax3.add_collection3d(lc)

    uu_ = np.linspace(0, 2 * np.pi, 41)
    vv_ = np.linspace(0, np.pi, 21)
    Rs = Rstar / S8
    ax3.plot_wireframe(Rs * np.outer(np.cos(uu_), np.sin(vv_)),
                       Rs * np.outer(np.sin(uu_), np.sin(vv_)),
                       Rs * np.outer(np.ones_like(uu_), np.cos(vv_)),
                       rstride=6, cstride=6, color=C_GRID, linewidth=0.25,
                       zorder=1)
    lim = 0.62 * Rs
    ax3.set_xlim(-lim, lim); ax3.set_ylim(-lim, lim); ax3.set_zlim(-lim, lim)
    ax3.set_box_aspect((1, 1, 1))
    ax3.view_init(elev=22, azim=-58)
    ax3.tick_params(labelsize=6, pad=-2)
    for a in (ax3.xaxis, ax3.yaxis, ax3.zaxis):
        a.pane.set_facecolor("white"); a.pane.set_edgecolor(C_GRID)
        a.line.set_color(C_MUTED)
        a._axinfo["grid"]["color"] = C_GRID
        a._axinfo["grid"]["linewidth"] = 0.3
    ax3.set_xlabel(r"$x$ ($10^8$ cm)", fontsize=7, labelpad=-6)
    ax3.set_ylabel(r"$y$", fontsize=7, labelpad=-6)
    ax3.set_zlabel(r"$z$", fontsize=7, labelpad=-6)
    fig.savefig(HERE / "field_lines_3d.pdf")
    plt.close(fig)
    print("wrote field_lines_3d.pdf")


if __name__ == "__main__":
    main()
