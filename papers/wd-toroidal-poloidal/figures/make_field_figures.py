"""The field of the adopted configuration, inside the star.

Run:  scf/.venv/bin/python3 papers/wd-toroidal-poloidal/figures/make_field_figures.py

Three figures, and the axisymmetry of this configuration decides the form of
all of them.

3D poloidal loops (Fig. 2 of the paper). A poloidal line has no azimuthal
component, so it stays in one meridional plane and is exactly a contour of
u there: a planar closed loop about the O-point, repeated at every azimuth.
Levels straddle the separatrix so the figure shows which loops close inside
the star and which close outside it.

Meridional slice. Not used by the paper any more -- the 3D poloidal figure
replaced it -- but kept as a diagnostic, since it is the quickest way to see
the B_phi distribution and the density envelope together. For an axisymmetric field the poloidal field lines ARE
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

The field is toroidal-dominated, so a line is a tightly wound helix on a
nested torus: the drawn ones wind 17 to 19 times around the axis per
meridional circuit. That winding is a property of each flux surface and is
NOT the amplitude ratio B_t/B_p of the table, which is a ratio of two
maxima taken at different places in the star.

Which surfaces exist is set by the flux at the stellar surface. Along the
equator u rises from 0 on the axis to u_axis at the O-point and then falls
to 0.775 u_axis at the surface -- it does not return to zero. Surfaces
above that value are closed tori inside the star; below it they intersect
the surface and their lines leave the star, where rho and hence B_phi
vanish and the line stops winding. The lines here are seeded by FLUX
inside the closed region, not by radius: equally spaced radii land on
unequally spaced surfaces, straddle the separatrix, and mix the two
families in one picture without saying so.
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


def surface_of_revolution(r, th, H, n_th=121, n_ph=121):
    """The true stellar surface, R(theta) from the enthalpy, revolved.

    Not a sphere: this star is prolate, R_pol/R_eq = 1.10, and a sphere of
    max(R_eq, R_pol) overstates the equatorial radius by that whole margin.
    """
    th_s = np.linspace(0.0, np.pi, n_th)
    R_s = np.array([r[H[:, np.argmin(np.abs(th - t))] > 0].max()
                    for t in th_s])
    ph_s = np.linspace(0.0, 2.0 * np.pi, n_ph)
    sin_s, cos_s = np.sin(th_s)[:, None], np.cos(th_s)[:, None]
    X = (R_s[:, None] * sin_s) * np.cos(ph_s)[None, :]
    Y = (R_s[:, None] * sin_s) * np.sin(ph_s)[None, :]
    Z = (R_s[:, None] * cos_s) * np.ones_like(ph_s)[None, :]
    return X, Y, Z, R_s


def draw_star(ax3, X, Y, Z, S8, elev, azim, stride=5):
    """Draw the surface split into far and near halves about the view.

    mplot3d has no occlusion: it either sorts whole artists by their mean
    depth (computed_zorder=True, the default, which silently overrides any
    zorder set on an artist) or obeys zorder literally. Neither can put half
    a surface behind the field lines and half in front, which is what a
    wireframe around an interior field has to do. So the surface is split
    here, by the sign of the dot product with the camera direction, and the
    caller draws field lines at a zorder between the two halves. The near
    half is drawn faint, as glass rather than as a wall.
    """
    e, a = np.deg2rad(elev), np.deg2rad(azim)
    cam = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    near = (X * cam[0] + Y * cam[1] + Z * cam[2]) > 0.0
    for mask, zo, alpha in ((~near, 0, 1.0), (near, 20, 0.45)):
        ax3.plot_wireframe(
            np.where(mask, X, np.nan) / S8, np.where(mask, Y, np.nan) / S8,
            np.where(mask, Z, np.nan) / S8, rstride=stride, cstride=stride,
            color="#c9c8c0", linewidth=0.3, alpha=alpha, zorder=zo)


def fig_poloidal_3d(r, th, H, u, Br, Bth, Rstar):
    """The poloidal field in three dimensions, to pair with the total field.

    A poloidal line has no azimuthal component, so it stays in one meridional
    plane and IS a contour of u there -- exactly. Each line is therefore a
    planar closed loop about the O-point, and the three-dimensional field is
    those loops repeated at every azimuth; a few azimuths are drawn.

    The levels are chosen to straddle the separatrix: the inner three are the
    same surfaces the total-field figure winds its helices on, and the outer
    two leave the star and close outside it, where rho and hence B_phi vanish.
    That contrast is the twisted-torus geometry, and it is why the toroidal
    field is confined while the poloidal one reaches the exterior.
    """
    from matplotlib.colors import LinearSegmentedColormap, LogNorm

    S8 = 1e8
    Bp = np.hypot(Br, Bth)
    rmax = 2.0 * Rstar
    vp, zz, VP, ZZ, (u_m, bp_m, H_m) = to_meridional(
        r, th, (u, Bp, H), rmax, n=601)
    inside = H_m > 0.0

    k_o = np.unravel_index(
        np.argmax(np.where(inside, np.abs(u_m), -np.inf)), u_m.shape)
    u_axis = u_m[k_o]
    j_eq = np.argmin(np.abs(zz))
    w_edge = vp[inside[:, j_eq]].max()
    u_sep = u_m[np.argmin(np.abs(vp - w_edge)), j_eq] / u_axis
    print(f"  separatrix at u/u_axis = {u_sep:.4f}")

    fracs = sorted([0.97, 0.90, 0.80, 0.68, 0.55])   # contour wants increasing
    tmp = plt.figure()
    cs = tmp.gca().contour(VP, ZZ, u_m, levels=[f * u_axis for f in fracs])
    segs = [list(s) for s in cs.allsegs]
    plt.close(tmp)

    f_bp = RGI((vp, zz), bp_m, bounds_error=False, fill_value=0.0)
    cmap = LinearSegmentedColormap.from_list("b", BLUES[1:])
    finite = bp_m[inside & (bp_m > 0)]
    norm = LogNorm(vmin=np.percentile(finite, 1.0),
                   vmax=np.percentile(finite, 99.9))

    fig = plt.figure(figsize=(WIDE_IN * 0.62, WIDE_IN * 0.42))
    ax3 = fig.add_subplot(111, projection="3d", computed_zorder=False)
    # four azimuths, not more: every loop passes close to the axis at high
    # |z|, where the poloidal field is strongest, and more planes pile up
    # into an unreadable dark column there
    phis = np.deg2rad([0.0, 90.0, 180.0, 270.0])
    reach = 0.0

    for frac, level_segs in zip(fracs, segs):
        for seg in level_segs:
            if len(seg) < 8:
                continue
            w, z = seg[:, 0], seg[:, 1]
            reach = max(reach, np.hypot(w, z).max())
            mag = f_bp(np.stack([w, z], axis=-1))
            for ph in phis:
                line = np.stack([w * np.cos(ph), w * np.sin(ph), z], axis=-1)
                pieces = np.stack([line[:-1] / S8, line[1:] / S8], axis=1)
                lc = Line3DCollection(pieces, cmap=cmap, norm=norm,
                                      linewidths=0.65, zorder=10)
                lc.set_array(0.5 * (mag[:-1] + mag[1:]))
                ax3.add_collection3d(lc)
        out = "closes outside the star" if frac < u_sep else "closed inside"
        print(f"  u/u_axis = {frac:.2f}: {len(level_segs)} loop(s), {out}")

    XS, YS, ZS, R_s = surface_of_revolution(r, th, H)
    draw_star(ax3, XS, YS, ZS, S8, 22.0, -58.0, stride=6)

    lim = 1.02 * reach / S8
    ax3.set_xlim(-lim, lim); ax3.set_ylim(-lim, lim); ax3.set_zlim(-lim, lim)
    ax3.set_box_aspect((1, 1, 1), zoom=1.14)
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
    print(f"  R_eq = {R_s[len(R_s)//2]/S8:.3f}e8, R_pol = {R_s[0]/S8:.3f}e8, "
          f"loops reach {reach/S8:.3f}e8, frame {lim:.3f}e8")

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax3, fraction=0.030, pad=0.12)
    cb.set_label(r"$|\mathbf{B}_{\rm p}|$  (G)", fontsize=7.5)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=6.5, width=0.5, length=2, pad=1.5)

    fig.savefig(HERE / "field_poloidal_3d.pdf")
    plt.close(fig)
    print("wrote field_poloidal_3d.pdf")


def main():
    from matplotlib.colors import LinearSegmentedColormap, LogNorm

    r, th, rho, H, u, Bphi, Br, Bth, Rstar = solve()
    sin_t, cos_t = np.sin(th)[None, :], np.cos(th)[None, :]
    Bvp = Br * sin_t + Bth * cos_t          # cylindrical components from
    Bz = Br * cos_t - Bth * sin_t           # the spherical ones
    Req, Rpol = diag.equatorial_polar_radii(H, r, th)
    print(f"R_eq = {Req:.4e}  R_pol = {Rpol:.4e}  R_pol/R_eq = {Rpol/Req:.4f}")

    fig_poloidal_3d(r, th, H, u, Br, Bth, Rstar)

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

    # The magnetic axis: the O-point of the poloidal field, where u peaks.
    # It is what labels the flux surfaces, and what the lines wind around.
    k_o = np.unravel_index(
        np.argmax(np.where(inside, np.abs(u_m), -np.inf)), u_m.shape)
    w_axis, z_axis, u_axis = VP[k_o], ZZ[k_o], u_m[k_o]

    def seed_at_flux(frac):
        """Equatorial radius outside the magnetic axis carrying u = frac*u_axis.

        Seeding by flux rather than by radius is the point: u labels the
        surfaces, radius does not. Equally spaced radii land on unequally
        spaced -- and non-monotonic -- surfaces, because a radius inside the
        axis and one outside it can carry the same flux."""
        w_end = vp[inside[:, np.argmin(np.abs(zz))]].max()
        ws = np.linspace(w_axis, w_end, 4000)
        us = f_u(np.stack([ws, np.zeros_like(ws)], axis=-1)) / u_axis
        return float(np.interp(frac, us[::-1], ws[::-1]))

    def trace(vp0, z0, max_step=400000, ds=None):
        """One full meridional circuit, so the lines are comparable objects."""
        ds = ds or 0.0015 * Rstar
        w, z, ph = vp0, z0, 0.0
        pts, mags = [], []
        prev = np.arctan2(z - z_axis, w - w_axis)
        swept = 0.0
        for _ in range(max_step):
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
            ang = np.arctan2(z - z_axis, w - w_axis)
            swept += (ang - prev + np.pi) % (2 * np.pi) - np.pi
            prev = ang
            if abs(swept) >= 2 * np.pi:
                break
        return np.array(pts), np.array(mags), abs(ph) / (2 * np.pi)

    fig = plt.figure(figsize=(WIDE_IN * 0.62, WIDE_IN * 0.42))
    ax3 = fig.add_subplot(111, projection="3d", computed_zorder=False)
    # Only surfaces above the separatrix are closed tori. On the equator u
    # runs 0 at the axis, up to u_axis, then DOWN to 0.775 u_axis at the
    # stellar surface -- it does not return to zero, so surfaces below that
    # value intersect the surface and their lines leave the star. Seed inside
    # the closed region; drawing an open line beside closed ones without
    # saying so is what made the first version of this figure confusing.
    fluxes = [0.97, 0.90, 0.80]
    print(f"  magnetic axis at varpi = {w_axis / S8:.3f}e8 cm")
    drawn = []
    for frac in fluxes:
        w0 = seed_at_flux(frac)
        line, mag, turns = trace(w0, 0.0)
        if len(line) < 10:
            continue
        u0 = f_u(np.array([[w0, 0.0]]))[0]
        uu = f_u(np.stack([np.sqrt(line[:, 0]**2 + line[:, 1]**2),
                           line[:, 2]], axis=-1))
        drift = np.max(np.abs(uu - u0)) / max(abs(u0), 1e-30)
        print(f"  u/u_axis = {frac:.2f}: seed varpi = {w0 / S8:.3f}e8 cm, "
              f"{turns:.1f} turns/circuit, |u| drift = {drift:.2e}")
        drawn.append((frac, line, mag, turns))

    # one colour scale across all lines, or the strength comparison between
    # surfaces is lost
    norm = plt.Normalize(min(m.min() for _, _, m, _ in drawn),
                         max(m.max() for _, _, m, _ in drawn))
    cmap = LinearSegmentedColormap.from_list("b", BLUES[1:])
    for frac, line, mag, turns in drawn:
        segs = np.stack([line[:-1] / S8, line[1:] / S8], axis=1)
        lc = Line3DCollection(segs, cmap=cmap, norm=norm, linewidths=0.55,
                              zorder=10)
        lc.set_array(0.5 * (mag[:-1] + mag[1:]))
        ax3.add_collection3d(lc)

    # The surface used to be a SPHERE of max(R_eq, R_pol), drawn at radius
    # R_pol inside a frame of 0.62 R_pol -- both the wrong shape, inflating
    # the equator by the whole 10% of prolateness, and clipped by the axes,
    # which is why its south pole ran off the box. Now the true surface of
    # revolution, with a frame that contains all of it.
    XS, YS, ZS, R_s = surface_of_revolution(r, th, H)
    draw_star(ax3, XS, YS, ZS, S8, 22.0, -58.0, stride=6)

    lim = 1.02 * R_s.max() / S8
    print(f"  R_eq = {R_s[len(R_s)//2]/S8:.3f}e8, R_pol = {R_s[0]/S8:.3f}e8, "
          f"R_pol/R_eq = {R_s[0]/R_s[len(R_s)//2]:.4f}, frame {lim:.3f}e8")
    ax3.set_xlim(-lim, lim); ax3.set_ylim(-lim, lim); ax3.set_zlim(-lim, lim)
    ax3.set_box_aspect((1, 1, 1), zoom=1.14)
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

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax3, fraction=0.030, pad=0.12)
    cb.set_label(r"$|\mathbf{B}|$  (G)", fontsize=7.5)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=6.5, width=0.5, length=2, pad=1.5)

    fig.savefig(HERE / "field_lines_3d.pdf")
    plt.close(fig)
    print("wrote field_lines_3d.pdf")


if __name__ == "__main__":
    main()


def measure_winding(starts=(0.28, 0.50, 0.70)):
    """Turns around the axis per meridional circuit, and the local
    |B_phi|/|B_pol| along the same lines.

    Written because the paper first claimed the winding number was the
    amplitude ratio of the configuration table. It is not: that ratio is
    max|B_phi| / max|B_pol| with the maxima taken at different places, a
    global number, while the winding is a property of each flux surface.
    Measured, the two disagree and move independently -- the local ratio
    varies by nearly an order of magnitude across these surfaces while the
    winding changes by less than a fifth.
    """
    from scipy.interpolate import RegularGridInterpolator as RGI

    r, th, rho, H, u, Bphi, Br, Bth, Rstar = solve()
    sin_t, cos_t = np.sin(th)[None, :], np.cos(th)[None, :]
    Bvp, Bz = Br * sin_t + Bth * cos_t, Br * cos_t - Bth * sin_t
    vp, zz, VP, ZZ, (H_m, u_m, bphi_m, bvp_m, bz_m) = to_meridional(
        r, th, (H, u, Bphi, Bvp, Bz), 1.05 * Rstar)
    inside = H_m > 0.0

    # the O-point of the poloidal loops, the centre to measure angles from
    k = np.unravel_index(np.argmax(np.where(inside, np.abs(u_m), -np.inf)),
                         u_m.shape)
    wO, zO = VP[k], ZZ[k]

    def interp(A):
        return RGI((vp, zz), A, bounds_error=False, fill_value=0.0)
    f_w, f_z, f_f = interp(bvp_m), interp(bz_m), interp(bphi_m)

    print(f"O-point at varpi = {wO:.4e}, z = {zO:.4e} cm")
    print("  start/R    turns per circuit    mean local |B_phi|/|B_pol|")
    for frac in starts:
        w, z, phi = frac * Rstar, 0.0, 0.0
        ds = 0.0015 * Rstar
        prev = np.arctan2(z - zO, w - wO)
        swept, ratios = 0.0, []
        for _ in range(400000):
            p = np.array([[w, z]])
            bw, bz_, bf = f_w(p)[0], f_z(p)[0], f_f(p)[0]
            b = np.hypot(np.hypot(bw, bz_), bf)
            if b <= 0 or w <= 0:
                break
            bpol = np.hypot(bw, bz_)
            if bpol > 0:
                ratios.append(abs(bf) / bpol)
            w += ds * bw / b
            z += ds * bz_ / b
            phi += ds * bf / (w * b)
            ang = np.arctan2(z - zO, w - wO)
            step = (ang - prev + np.pi) % (2 * np.pi) - np.pi
            swept += step
            prev = ang
            if abs(swept) >= 2 * np.pi:
                print(f"    {frac:.2f}          {abs(phi) / (2 * np.pi):8.1f}"
                      f"              {np.mean(ratios):10.1f}")
                break
        else:
            print(f"    {frac:.2f}          circuit did not close")
