"""Funcoes de figura (Aba 1). So' figuras — toda fisica vem de scf.* (regra R1).

Unidades de exibicao (R4): eixos radiais em km, campo em gauss com colorbar
em notacao cientifica — conversao sempre via units.py, nunca um fator solto
aqui."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

import units


def _meridional_grid(r, theta):
    """varpi, z em KM (nao cm) — malha para os eixos das figuras."""
    R, TH = np.meshgrid(r, theta, indexing="ij")
    varpi = units.cm_to_km(R * np.sin(TH))
    z = units.cm_to_km(R * np.cos(TH))
    return varpi, z


def _gauss_colorbar(fig, pc, ax, label):
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-2, 2))  # forca notacao cientifica
    return fig.colorbar(pc, ax=ax, label=label, format=fmt)


def plot_density(rho, r, theta, H=None):
    """Mapa de cor de rho no plano meridional, com a fronteira H=0 destacada."""
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
    ax.set_title("Densidade (ciano: H=0)")
    fig.tight_layout()
    return fig


def plot_flux_contours(u, r, theta, u_c=None, n_levels=14):
    """Contornos de u=omega*A_phi — linhas de campo poloidal — com a ultima
    linha fechada (u_c) destacada, se fornecida."""
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
    title = "Linhas de campo poloidal"
    if u_c is not None:
        title += " (vermelho: última linha fechada, $u_c$)"
    ax.set_title(title)
    fig.tight_layout()
    return fig


def plot_toroidal(Bphi, r, theta):
    """Mapa de cor de B_phi — mostra o toro."""
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
    ax.set_title("Campo toroidal")
    fig.tight_layout()
    return fig


def plot_convergence(history):
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.semilogy(history, marker=".")
    ax.axhline(1e-6, color="gray", linestyle="--", linewidth=0.8, label="tol padrão")
    ax.set_xlabel("iteração")
    ax.set_ylabel(r"max|$\Delta\rho$|/$\rho_c$")
    ax.set_title("Convergência do SCF")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_virial_history(ve_history):
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.semilogy(ve_history, marker=".", color="darkred")
    ax.axhline(1e-3, color="gray", linestyle="--", linewidth=0.8, label="limite V3 (1e-3)")
    ax.set_xlabel("iteração")
    ax.set_ylabel("erro virial (VE)")
    ax.set_title("Erro virial por iteração")
    ax.legend()
    fig.tight_layout()
    return fig
