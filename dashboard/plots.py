"""Figure functions (Tab 1). Figures only — all physics comes from scf.* (rule R1).

Display units (R4): radial axes in km, field in gauss with colorbar in
scientific notation — conversion always via units.py, never a loose factor
here."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import units


def _meridional_grid(r, theta):
    """varpi, z in KM (not cm) — grid for the figure axes."""
    R, TH = np.meshgrid(r, theta, indexing="ij")
    varpi = units.cm_to_km(R * np.sin(TH))
    z = units.cm_to_km(R * np.cos(TH))
    return varpi, z


def _gauss_colorbar(fig, pc, ax, label):
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-2, 2))  # force scientific notation
    return fig.colorbar(pc, ax=ax, label=label, format=fmt)


def plot_density(rho, r, theta, H=None):
    """Color map of rho in the meridional plane, with the H=0 boundary highlighted."""
    varpi, z = _meridional_grid(r, theta)
    fig, ax = plt.subplots(figsize=(4.5, 5.5))
    pc = None
    for sign in (1, -1):
        pc = ax.pcolormesh(sign * varpi, z, rho, shading="auto", cmap="inferno")
    fig.colorbar(pc, ax=ax, label=r"$\rho$ (g/cm$^3$)")
    if H is not None:
        for sign in (1, -1):
            ax.contour(sign * varpi, z, H, levels=[0.0], colors="cyan", linewidths=1.2)
    ax.set_xlabel(r"$\varpi$ (km)")
    ax.set_ylabel(r"$z$ (km)")
    ax.set_aspect("equal")
    ax.set_title("Density (cyan: H=0)")
    fig.tight_layout()
    return fig


def plot_flux_contours(u, r, theta, rho=None, u_c=None, H=None, n_levels=7):
    """Contours of u=omega*A_phi — poloidal field lines.

    The region u>u_c is the object of interest (D6): it is entirely
    confined inside the star and is where the toroidal field is imposed —
    so it is filled with a translucent color, not just outlined. Field
    lines outside the stellar surface (rho<=0) are drawn light gray/dashed
    to de-emphasize them relative to the interior structure — nearly all
    of them are topologically closed loops too (u -> 0 on the axis and far
    from the star, so any level set closes around the single interior
    maximum of u, like ordinary dipole field lines), which is a different
    sense of "closed" from the torus boundary and was a source of
    confusion when both were drawn the same way.

    The H=0 curve (cyan) is kept for context, but tangency between it and
    the torus boundary (u=u_c) should be checked NUMERICALLY, not by eye
    — see toroidal.check_uc_tangency() / scf/tests/test_toroidal.py.
    Contour plotting uses 2D marching-squares interpolation on the curved
    (r,theta)->(varpi,z) grid, which can visually offset two independently
    drawn contour lines by up to about one grid cell even when they are
    mathematically exactly tangent (verified: the visual gap shrinks
    roughly linearly with grid spacing and vanishes at high resolution)."""
    varpi, z = _meridional_grid(r, theta)
    fig, ax = plt.subplots(figsize=(4.5, 6.1))
    umin, umax = np.min(u), np.max(u)
    levels = np.linspace(umin, umax, n_levels) if umax > umin else [umin]
    inside_star = rho > 0 if rho is not None else None

    for sign in (1, -1):
        x = sign * varpi
        if inside_star is not None:
            u_ext = np.where(inside_star, np.nan, u)
            ax.contour(x, z, u_ext, levels=levels, colors="lightgray",
                       linewidths=0.6, linestyles="dashed")
            u_int = np.where(inside_star, u, np.nan)
            ax.contour(x, z, u_int, levels=levels, colors="k", linewidths=0.7)
        else:
            ax.contour(x, z, u, levels=levels, colors="k", linewidths=0.7)

        if H is not None:
            ax.contour(x, z, H, levels=[0.0], colors="cyan", linewidths=0.9, alpha=0.6)

        if u_c is not None and umin < u_c < umax:
            ax.contourf(x, z, u, levels=[u_c, umax], colors=["tab:red"], alpha=0.18)
            ax.contour(x, z, u, levels=[u_c], colors="darkred", linewidths=1.4)

    ax.set_xlabel(r"$\varpi$ (km)")
    ax.set_ylabel(r"$z$ (km)")
    ax.set_aspect("equal")
    ax.set_title("Poloidal field lines")

    if u_c is not None:
        handles = [Patch(facecolor="tab:red", edgecolor="darkred", alpha=0.4,
                          label="torus boundary\n(last field line confined to the star)")]
        if H is not None:
            handles.append(Line2D([0], [0], color="cyan", alpha=0.6, lw=0.9,
                                   label="stellar surface (H=0)"))
        # figure-level legend in a reserved bottom strip (rect below), not
        # anchored to the axes — ax.legend(bbox_to_anchor=...) below the
        # xlabel collides with it because tight_layout() does not reliably
        # reserve room for legend artists placed outside the axes bbox.
        fig.tight_layout(rect=(0.0, 0.16, 1.0, 1.0))
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
                   fontsize=6.5, frameon=False, handlelength=1.4)
    else:
        fig.tight_layout()
    return fig


def plot_toroidal(Bphi, r, theta):
    """Color map of B_phi — shows the torus."""
    varpi, z = _meridional_grid(r, theta)
    fig, ax = plt.subplots(figsize=(4.5, 5.5))
    vmax = np.max(np.abs(Bphi))
    vmax = vmax if vmax > 0 else 1.0
    pc = None
    for sign in (1, -1):
        pc = ax.pcolormesh(sign * varpi, z, Bphi, shading="auto", cmap="RdBu_r",
                            vmin=-vmax, vmax=vmax)
    _gauss_colorbar(fig, pc, ax, r"$B_\phi$ (G)")
    ax.set_xlabel(r"$\varpi$ (km)")
    ax.set_ylabel(r"$z$ (km)")
    ax.set_aspect("equal")
    ax.set_title("Toroidal field")
    fig.tight_layout()
    return fig


def plot_density_profile(rho, r, theta):
    """rho(r) from the center to the surface, along the pole and the
    equator (closest grid line to theta=pi/2). Shows any density-peak
    migration away from the center in the strong-field regime (see
    docs/teoria.md §6)."""
    r_km = units.cm_to_km(r)
    i_pole = 0
    i_eq = int(np.argmin(np.abs(theta - np.pi / 2)))
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.plot(r_km, rho[:, i_eq], label="equator", color="tab:blue")
    ax.plot(r_km, rho[:, i_pole], label="pole", color="tab:orange")
    ax.set_xlabel("r (km)")
    ax.set_ylabel(r"$\rho$ (g/cm$^3$)")
    ax.set_title("Density profile (center to surface)")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend()
    fig.tight_layout()
    return fig


def plot_convergence(history):
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.semilogy(history, marker=".")
    ax.axhline(1e-6, color="gray", linestyle="--", linewidth=0.8, label="default tol")
    ax.set_xlabel("iteration")
    ax.set_ylabel(r"max|$\Delta\rho$|/$\rho_c$")
    ax.set_title("SCF convergence")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_virial_history(ve_history):
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.semilogy(ve_history, marker=".", color="darkred")
    ax.axhline(1e-3, color="gray", linestyle="--", linewidth=0.8, label="V3 limit (1e-3)")
    ax.set_xlabel("iteration")
    ax.set_ylabel("virial error (VE)")
    ax.set_title("Virial error per iteration")
    ax.legend()
    fig.tight_layout()
    return fig
