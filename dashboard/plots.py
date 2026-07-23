"""Figure functions (Tab 1). Figures only — all physics comes from scf.* (rule R1).

Display units (R4): radial axes in km, field in gauss with colorbar in
scientific notation — conversion always via units.py, never a loose factor
here."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

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


def plot_flux_contours(u, r, theta, u_c=None, n_levels=14):
    """Contours of u=omega*A_phi — poloidal field lines — with the last
    closed line (u_c) highlighted, if provided."""
    varpi, z = _meridional_grid(r, theta)
    fig, ax = plt.subplots(figsize=(4.5, 5.5))
    umin, umax = np.min(u), np.max(u)
    levels = np.linspace(umin, umax, n_levels) if umax > umin else [umin]
    for sign in (1, -1):
        ax.contour(sign * varpi, z, u, levels=levels, colors="k", linewidths=0.6)
        if u_c is not None and umin < u_c < umax:
            ax.contour(sign * varpi, z, u, levels=[u_c], colors="red", linewidths=1.8)
    ax.set_xlabel(r"$\varpi$ (km)")
    ax.set_ylabel(r"$z$ (km)")
    ax.set_aspect("equal")
    title = "Poloidal field lines"
    if u_c is not None:
        title += " (red: last closed line, $u_c$)"
    ax.set_title(title)
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
