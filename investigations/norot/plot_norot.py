"""The non-rotating 192^3 run to 60 s: every figure, from the files in data/.

Run:  scf/.venv/bin/python3 investigations/norot/plot_norot.py

Four figures, written to plots/ as .pdf and .png:

  norot_field          (a) E_tor/E_pol and (b) E_mag/E_mag(0), against the
                       rotating rot192. Without rotation the ratio crosses 1
                       near 3 s and stays between 0.19 and 0.48 to 60 s,
                       inside the mixed (Braithwaite) band. The energy does
                       NOT plateau past 26 s, as DIARIO 29 read it at 26 s: it
                       decays with an e-folding near 30 s while the ratio holds,
                       so the configuration loses energy without changing shape.

  norot_star           (a) rho_max, (b) R_pol/R_eq, (c) mass drift. The star
                       keeps its mass to 0.2%, but its radial pulsation GROWS:
                       (max-min)/(max+min) of rho_max goes from 0.15 at 5-10 s
                       to 0.50 at 55-60 s, e-folding ~35 s, period 1.69 s. The
                       rotating rot192 does the opposite (e-folding ~140 s decay,
                       plot_long_run.py). The oblateness is an oscillation, not a
                       settled state (DIARIO 29 quoted 0.823 at 26 s; the last
                       20 s sit at 0.95-0.99).

  norot_timestep       dt and the rejected advances. The first run (44987,
                       no ambient refill) stalls at 26 s as the ambient
                       evacuates; the second (45194, fill_ambient_bc = 1)
                       crosses the same crisis and recovers.

  norot_energy_budget  E_int, |E_grav|, E_kin and E_mag. E_int grows twelvefold,
                       but under ztwd the pressure depends on density alone, so
                       it is thermodynamically inert; E_grav barely moves.

The dashed line in the rotating curves is rot192, which differs in rotation AND
in the star (2.005 against 1.409 Msun): not a one-variable pair (DIARIO 29, 31).
"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "plots"
INV = HERE.parent
COL_IN = 88.0 / 25.4

C_NEW = "#2a78d6"     # norot192b, 45194, with ambient refill: the run
C_OLD = "#898781"     # norot192, 44987, died at 26 s
C_ROT = "#eb6834"     # rot192
C_MUTED = "#898781"
C_GRID = "#e1e0d9"
C_INK = "#0b0b0b"
C_BAND = "#e8e6df"

T_CRISIS = 26.1
FIT = (26.0, 60.0)

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


def read(path, delimiter=","):
    with open(path) as fh:
        n = sum(1 for line in fh if line.startswith("#"))
    return np.genfromtxt(path, delimiter=delimiter, names=True, skip_header=n)


def rot192():
    """rot192 field, 0-12 s from bt_bp_192.csv and 12-60 s from the late file."""
    a = read(INV / "bt_bp_192.csv")
    b = read(INV / "bt_bp_192_late.csv")
    late = b["t"] > a["t"][-1]
    t = np.r_[a["t"], b["t"][late]]
    e = np.r_[a["E_tor"] + a["E_pol"], (b["E_tor"] + b["E_pol"])[late]]
    r = np.r_[a["Et_over_Ep"], (b["E_tor"] / b["E_pol"])[late]]
    return t, e, r


def style(ax):
    ax.grid(True, color=C_GRID, linewidth=0.4, alpha=0.9)
    ax.set_xlim(0, 60)


def label(ax, s, y=0.93, va="top"):
    ax.text(0.02, y, s, transform=ax.transAxes, fontsize=7, va=va)


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=200)
    plt.close(fig)


def fig_field(new, old):
    rt, re, rr = rot192()
    fig, (ax_r, ax_e) = plt.subplots(2, 1, figsize=(COL_IN, COL_IN * 1.3), sharex=True,
                                     gridspec_kw={"hspace": 0.08})
    for ax in (ax_r, ax_e):
        style(ax)

    ax_r.axhspan(0.1, 10, color=C_BAND, linewidth=0, zorder=0)
    ax_r.text(58.5, 3.2, "mixed branch", fontsize=6.2, color=C_MUTED, ha="right")
    ax_r.plot(rt, rr, color=C_ROT, linewidth=1.0, linestyle=(0, (4, 2)), label="rotating (rot192)")
    ax_r.plot(old["t"], old["Et_over_Ep"], color=C_OLD, linewidth=0.9, label="non-rotating, first run (44987)")
    ax_r.plot(new["t"], new["Et_over_Ep"], color=C_NEW, linewidth=1.1, label="non-rotating (45194)")
    ax_r.set_yscale("log")
    ax_r.set_ylim(0.08, 1e6)
    ax_r.set_ylabel(r"$E_{\rm tor}/E_{\rm pol}$")
    ax_r.legend(loc="upper right", ncol=1)
    label(ax_r, "(a)", y=0.07, va="bottom")

    e_new = new["E_tor"] + new["E_pol"]
    e_old = old["E_tor"] + old["E_pol"]
    ax_e.plot(rt, re / re[0], color=C_ROT, linewidth=1.0, linestyle=(0, (4, 2)))
    ax_e.plot(old["t"], e_old / e_old[0], color=C_OLD, linewidth=0.9)
    ax_e.plot(new["t"], e_new / e_new[0], color=C_NEW, linewidth=1.1)
    s = (new["t"] >= FIT[0]) & (new["t"] <= FIT[1])
    k, c = np.polyfit(new["t"][s], np.log(e_new[s] / e_new[0]), 1)
    tt = np.linspace(*FIT, 50)
    ax_e.plot(tt, np.exp(k * tt + c), color=C_INK, linewidth=0.7, linestyle=(0, (1, 1.6)))
    ax_e.text(44, 0.16, rf"$e$-folding {abs(1 / k):.0f} s", fontsize=6.2, color=C_INK, ha="center")
    ax_e.set_yscale("log")
    ax_e.set_ylim(2e-4, 1.5)
    ax_e.set_ylabel(r"$E_{\rm mag}/E_{\rm mag}(0)$")
    ax_e.set_xlabel("t (s)")
    label(ax_e, "(b)", y=0.07, va="bottom")
    save(fig, "norot_field")
    return abs(1 / k)


def fig_star(new, gd, gd_old):
    fig, (ax_d, ax_q, ax_m) = plt.subplots(3, 1, figsize=(COL_IN, COL_IN * 1.55), sharex=True,
                                           gridspec_kw={"hspace": 0.08})
    for ax in (ax_d, ax_q, ax_m):
        style(ax)
    ax_d.plot(gd_old["t"], gd_old["rho_max"] / 1e9, color=C_OLD, linewidth=0.6)
    ax_d.plot(gd["t"], gd["rho_max"] / 1e9, color=C_NEW, linewidth=0.6)
    # envelope in 5 s windows: the pulsation grows instead of damping
    w = np.arange(0, 60, 5)
    hi = np.array([gd["rho_max"][(gd["t"] >= a) & (gd["t"] < a + 5)].max() for a in w]) / 1e9
    lo = np.array([gd["rho_max"][(gd["t"] >= a) & (gd["t"] < a + 5)].min() for a in w]) / 1e9
    amp = (hi - lo) / (hi + lo)
    s = w >= 5
    k = np.polyfit(w[s] + 2.5, np.log(amp[s]), 1)[0]
    ax_d.text(0.5, 0.88, rf"amplitude {amp[1]:.2f} $\to$ {amp[-1]:.2f}, growth $e$-folding {1 / k:.0f} s",
              transform=ax_d.transAxes, fontsize=6.2, color=C_INK, ha="center")
    ax_d.set_ylim(1.5, 6.4)
    ax_d.set_ylabel(r"$\rho_{\rm max}$  ($10^9$ g cm$^{-3}$)")
    label(ax_d, "(a)")

    ax_q.plot(new["t"], new["Rpol_over_Req"], color=C_NEW, linewidth=1.0, marker="o", markersize=1.6)
    ax_q.axhline(1.0, color=C_MUTED, linewidth=0.5, linestyle=(0, (2, 2)))
    ax_q.set_ylabel(r"$R_{\rm pol}/R_{\rm eq}$")
    label(ax_q, "(b)", y=0.07, va="bottom")

    ax_m.plot(gd_old["t"], 100 * (gd_old["mass"] / gd_old["mass"][0] - 1), color=C_OLD, linewidth=0.9,
              label="box, first run")
    ax_m.plot(gd["t"], 100 * (gd["mass"] / gd["mass"][0] - 1), color=C_NEW, linewidth=1.0, label="box")
    ax_m.plot(new["t"], 100 * (new["M_Msun"] / new["M_Msun"][0] - 1), color=C_NEW, linewidth=1.0,
              linestyle=(0, (1, 1.4)), label=r"star ($\rho > 10^5$)")
    ax_m.set_ylabel(r"$\Delta M/M_0$  (%)")
    ax_m.set_xlabel("t (s)")
    ax_m.legend(loc="upper left", bbox_to_anchor=(0.08, 1.0))
    label(ax_m, "(c)")
    save(fig, "norot_star")


def fig_timestep(dt, dt_old, rt):
    fig, ax = plt.subplots(figsize=(COL_IN, COL_IN * 0.62))
    style(ax)
    ax.plot(dt_old["t"], dt_old["dt"], color=C_OLD, linewidth=0.7, label="first run (44987), no ambient refill")
    ax.plot(dt["t"], dt["dt"], color=C_NEW, linewidth=0.7, label="45194, fill_ambient_bc = 1")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 2e-2)
    y = np.full(len(rt), 1.6e-5)
    ax.plot(rt["t_start"], y, "|", color=C_INK, markersize=4, markeredgewidth=0.4)
    ax.text(rt["t_start"].max() + 1.0, 1.35e-5, f"rejected advances ({int(rt['n_retry'].sum())})",
            fontsize=6.0, color=C_INK)
    ax.axvline(T_CRISIS, color=C_MUTED, linewidth=0.6, linestyle=(0, (2, 2)))
    ax.text(T_CRISIS + 0.6, 1.2e-2, "ambient crisis", fontsize=6.0, color=C_MUTED)
    ax.set_ylabel(r"$\Delta t$  (s)")
    ax.set_xlabel("t (s)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1)
    save(fig, "norot_timestep")


def fig_budget(new, gd):
    fig, ax = plt.subplots(figsize=(COL_IN, COL_IN * 0.7))
    style(ax)
    ax.plot(gd["t"], gd["E_int"], color=C_ROT, linewidth=1.0, label=r"$E_{\rm int}$ (inert under ztwd)")
    ax.plot(gd["t"], -gd["E_grav"], color=C_INK, linewidth=1.0, label=r"$|E_{\rm grav}|$")
    ax.plot(gd["t"], gd["E_kin"], color="#1baf7a", linewidth=0.8, label=r"$E_{\rm kin}$")
    ax.plot(new["t"], new["E_tor"] + new["E_pol"], color=C_NEW, linewidth=1.0, label=r"$E_{\rm mag}$ ($\rho > 10^5$)")
    ax.set_yscale("log")
    ax.set_ylabel("energy (erg)")
    ax.set_xlabel("t (s)")
    ax.legend(loc="lower right")
    save(fig, "norot_energy_budget")


def main():
    new = read(DATA / "bt_bp_norot192b.csv")
    old = read(DATA / "bt_bp_norot192.csv")
    gd = read(DATA / "grid_diag_norot192b.csv")
    gd_old = read(DATA / "grid_diag_norot192.csv")
    dt = read(DATA / "dt_norot192b.csv")
    dt_old = read(DATA / "dt_norot192.csv")
    rt = read(DATA / "retries_norot192b.csv")
    tau = fig_field(new, old)
    fig_star(new, gd, gd_old)
    fig_timestep(dt, dt_old, rt)
    fig_budget(new, gd)
    print(f"plots/ escrito; e-folding de E_mag em {FIT[0]:.0f}-{FIT[1]:.0f} s: {tau:.1f} s")


if __name__ == "__main__":
    main()
