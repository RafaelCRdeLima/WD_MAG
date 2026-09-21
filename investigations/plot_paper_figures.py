"""The three parameter-space figures for the paper.

Same house style as plot_bt_bp.py -- 88 mm column, serif, thin marks, inward
ticks -- and the same palette, which passes the colorblind checks (worst
adjacent pair dE 24.7 under protanopia against a target of 8).

Data are the 257^2 sweep (investigations/res257/) rather than the 129^2 one.
The two agree to 0.01% on every headline (DIARIO 25), but the converged set is
what the paper quotes.

  fig_ceiling.pdf   beta_min against E_tor/|W|, one curve per central density,
                    with the literature's requirement marked. The whole
                    argument of Section 3 in one panel.
  fig_frontier.pdf  the mass frontier converging to M_Ch from below.
  fig_tradeoff.pdf  mass against comparability, which replaces the single
                    quoted maximum (DIARIO 26).
  fig_decay.pdf     the m=1 destruction of the ordered field, BOTH meshes on
                    one axis. The residue differs by a factor of 25 between
                    them and the figure shows that rather than hiding it: the
                    destruction is the result, its late-time level is not.
  fig_probe.pdf     the time step of the probe pair, which is how the
                    comparable-energy configuration dies in the box we have.

Run:  scf/.venv/bin/python3 investigations/plot_paper_figures.py
"""

import csv
import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "paper" / "figures"
COL_IN = 88.0 / 25.4

C_A, C_B, C_C = "#2a78d6", "#eb6834", "#4a3aa7"
C_MUTED, C_GRID, C_INK = "#898781", "#e1e0d9", "#0b0b0b"
M_CH = 1.44
LIT = 0.2026          # E_tor/|W| the equilibrium literature requires
RHO_NEUTRON = 1.940e10

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


def load(pattern):
    rows = []
    for f in sorted(glob.glob(str(pattern))):
        head = open(f).readline()
        meta = dict(t.split("=", 1) for t in head.lstrip("# ").split()
                    if "=" in t)
        for r in csv.DictReader(l for l in open(f) if not l.startswith("#")):
            r["rho_c"] = float(meta["rho_c"])
            rows.append(r)
    return rows


def style(ax):
    ax.grid(True, color=C_GRID, linewidth=0.4, alpha=0.9)
    ax.set_axisbelow(True)


def fig_ceiling(rows):
    """beta_min against the field the star is asked to carry.

    One curve per central density, each the envelope over the poloidal shape
    parameters. The five curves collapse onto one another, and that collapse is
    the result rather than a plotting accident: the ceiling is a property of
    the confinement, not of the star. They are therefore drawn in one colour
    with the range annotated, instead of spending five hues on a single line.
    """
    fig, ax = plt.subplots(figsize=(COL_IN, COL_IN * 0.78))
    style(ax)
    by_rho = {}
    for r in rows:
        e, b = float(r["Etor_over_W"]), float(r["beta_min_tor_only"])
        d = by_rho.setdefault(r["rho_c"], {})
        d[e] = max(d.get(e, 0.0), b)

    for rho, d in sorted(by_rho.items()):
        e = np.array(sorted(d))
        b = np.array([d[x] for x in e])
        ax.plot(e, b, color=C_A, linewidth=0.9, alpha=0.85,
                marker="o", markersize=2.2, markeredgewidth=0)

    ax.text(1.05e-3, 0.45,
            "five central densities,\n"
            r"$5\times10^{8}$ to $5\times10^{9}$ g cm$^{-3}$" "\n"
            "(indistinguishable)",
            fontsize=6.0, color=C_A, ha="left", va="top", linespacing=1.3)
    # o bloco fica sob a linha beta=1, onde a curva ja' subiu para fora

    ax.axhline(1.0, color=C_INK, linewidth=0.7, linestyle=(0, (3, 2)))
    ax.text(0.165, 1.28, r"$\beta_{\rm min}=1$", fontsize=6.4,
            color=C_INK, ha="right", va="bottom")

    # o teto e o fator que separa da literatura
    ax.annotate("", xy=(LIT, 3.2e-2), xytext=(0.0101, 3.2e-2),
                arrowprops=dict(arrowstyle="<->", color=C_B, linewidth=0.7,
                                shrinkA=0, shrinkB=0))
    ax.text(np.sqrt(LIT * 0.0101), 4.2e-2, r"$\times 20$", fontsize=7.0,
            color=C_B, ha="center", va="bottom")
    ax.plot([0.0101], [1.0], "o", color=C_INK, markersize=4.0,
            markeredgewidth=0.8, markeredgecolor="white", zorder=4)
    ax.text(0.0115, 2.6, r"ceiling, $E_{\rm tor}/|W| = 0.010$",
            fontsize=6.2, color=C_INK, ha="left", va="bottom")

    ax.axvline(LIT, color=C_B, linewidth=0.8, linestyle=(0, (3, 2)))
    ax.text(LIT * 0.86, 12.0, "required by\nequilibrium models", fontsize=6.0,
            color=C_B, ha="right", va="top", linespacing=1.25)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$E_{\rm tor}/|W|$")
    ax.set_ylabel(r"$\beta_{\rm min}$")
    ax.set_xlim(9e-4, 0.62); ax.set_ylim(3e-3, 40)
    fig.savefig(OUT / "fig_ceiling.pdf")
    plt.close(fig)
    print("  fig_ceiling.pdf")


def fig_frontier(rows, hi):
    """The heaviest star that holds its own confined toroidal field, against
    central density, with the field-free sequence for scale. The gap between
    the two curves is what the field buys -- about 0.017 Msun -- and the upper
    curve approaches M_Ch without reaching it before neutronization ends the
    sequence."""
    fig, ax = plt.subplots(figsize=(COL_IN, COL_IN * 0.72))
    style(ax)

    def frontier(rs):
        out = {}
        for r in rs:
            if float(r["beta_min_tor_only"]) >= 1.0:
                m = float(r["M_msun"])
                out[r["rho_c"]] = max(out.get(r["rho_c"], 0.0), m)
        return out

    def floor_mass(rs):
        out = {}
        for r in rs:
            m = float(r["M_msun"])
            k = r["rho_c"]
            out[k] = min(out.get(k, 1e9), m)
        return out

    f = {**frontier(rows), **frontier(hi)}
    g = {**floor_mass(rows), **floor_mass(hi)}
    x = np.array(sorted(f)); y = np.array([f[k] for k in x])
    xg = np.array(sorted(g)); yg = np.array([g[k] for k in xg])

    ax.plot(xg, yg, color=C_MUTED, linewidth=1.0, linestyle=(0, (3, 2)),
            marker="o", markersize=2.2, markeredgewidth=0)
    ax.plot(x, y, color=C_A, linewidth=1.2, marker="o", markersize=3.0,
            markeredgewidth=0)
    ax.axhline(M_CH, color=C_INK, linewidth=0.7, linestyle=(0, (1, 2)))

    ax.text(6e8, M_CH + 0.004, r"$M_{\rm Ch}=1.44\,M_\odot$", fontsize=6.2,
            color=C_INK, ha="left", va="bottom")
    ax.text(2.2e9, 1.408, "holds its field", fontsize=6.2, color=C_A,
            ha="center", va="bottom")
    ax.text(xg[-1], yg[-1] - 0.012, "no field", fontsize=6.2, color=C_MUTED,
            ha="right", va="top")
    ax.axvspan(RHO_NEUTRON, 3e10, color=C_GRID, alpha=0.75, linewidth=0)
    ax.text(RHO_NEUTRON * 1.06, 1.33, "neutronization", fontsize=6.0,
            color=C_MUTED, rotation=90, ha="left", va="bottom")

    ax.set_xscale("log")
    ax.set_xlabel(r"$\rho_c$ (g cm$^{-3}$)")
    ax.set_ylabel(r"$M$ ($M_\odot$)")
    ax.set_xlim(4e8, 2.8e10); ax.set_ylim(1.29, 1.46)
    fig.savefig(OUT / "fig_frontier.pdf")
    plt.close(fig)
    print("  fig_frontier.pdf")


def fig_tradeoff():
    """Mass against comparability. Both rise with rotation, so criterion (iii)
    binds before the bar mode: the quoted maximum depends on where the band
    edge is drawn, which is why the curve replaces the point."""
    rows = list(csv.DictReader(open(HERE / "tradeoff_257.csv")))
    e = np.array([float(r["Etor_over_Epol"]) for r in rows])
    m = np.array([float(r["M_msun"]) for r in rows])
    ok = np.array([r["all_four_gates"] == "True" for r in rows])

    fig, ax = plt.subplots(figsize=(COL_IN, COL_IN * 0.72))
    style(ax)
    ax.axvspan(0.1, 10.0, color=C_A, alpha=0.07, linewidth=0)
    ax.plot(e, m, color=C_A, linewidth=1.1, zorder=2)
    ax.plot(e[ok], m[ok], "o", color=C_A, markersize=4.0, markeredgewidth=0.8,
            markeredgecolor="white", zorder=3)
    ax.plot(e[~ok], m[~ok], "o", color="white", markersize=4.0,
            markeredgewidth=0.9, markeredgecolor=C_B, zorder=3)

    ax.axhline(M_CH, color=C_INK, linewidth=0.7, linestyle=(0, (1, 2)))
    ax.text(3.05, M_CH + 0.018, r"$M_{\rm Ch}$", fontsize=6.4, color=C_INK,
            ha="left", va="bottom")
    ax.axvline(10.0, color=C_B, linewidth=0.7, linestyle=(0, (3, 2)))
    ax.text(10.6, 1.62, "comparable-energy\nband edge", fontsize=6.0,
            color=C_B, ha="left", va="center", linespacing=1.25)
    ax.text(4.35, 1.6484 - 0.03, r"$1.65\,M_\odot$", fontsize=6.4, color=C_A,
            ha="left", va="top")
    ax.text(8.6, 2.175, r"$2.13\,M_\odot$", fontsize=6.4, color=C_A,
            ha="right", va="bottom")

    ax.set_xscale("log")
    # ticks explicitos: os minor labels do log colidiam em "3x10^0 4x10^0"
    ax.set_xticks([3, 5, 7, 10, 15, 20])
    ax.get_xaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.minorticks_off()
    ax.set_xlabel(r"$E_{\rm tor}/E_{\rm pol}$")
    ax.set_ylabel(r"$M$ ($M_\odot$)")
    ax.set_xlim(2.8, 21); ax.set_ylim(1.38, 2.46)
    fig.savefig(OUT / "fig_tradeoff.pdf")
    plt.close(fig)
    print("  fig_tradeoff.pdf")


def fig_decay(_=None):
    """The ordered field is destroyed, on both meshes, and then stops decaying.

    Plotting 192^3 and 256^3 together is deliberate. The destruction is robust:
    both lose most of the magnetic energy over the same few Alfven times. What
    is not robust is where the field settles -- the minima differ by a factor
    of 25 between the meshes -- and the figure shows that rather than leaving
    it to the text, so that no reader takes the residue for a measurement.

    The 192^3 data live in two files, 0-12 s and 12-60 s, and are concatenated
    here; the 256^3 run reaches 78 s in one.
    """
    def read(name):
        n_hdr = sum(1 for l in open(HERE / name) if l.startswith("#"))
        return np.genfromtxt(HERE / name, delimiter=",", names=True,
                             skip_header=n_hdr)

    a1, a2 = read("bt_bp_192.csv"), read("bt_bp_192_late.csv")
    t_a = np.concatenate([a1["t"], a2["t"][1:]])
    e_a = np.concatenate([(a1["E_tor"] + a1["E_pol"]),
                          (a2["E_tor"] + a2["E_pol"])[1:]])
    b = read("bt_bp_256_long.csv")
    t_b, e_b = b["t"], b["E_tor"] + b["E_pol"]

    fig, ax = plt.subplots(figsize=(COL_IN, COL_IN * 0.72))
    style(ax)
    ax.plot(t_a, e_a, color=C_A, linewidth=1.0)
    ax.plot(t_b, e_b, color=C_B, linewidth=1.0)
    ax.text(58.5, e_a[-1] * 1.9, r"$192^3$", color=C_A, fontsize=6.8,
            ha="right", va="bottom")
    ax.text(76.5, e_b[-1] * 1.9, r"$256^3$", color=C_B, fontsize=6.8,
            ha="right", va="bottom")

    # os dois minimos, que e' onde as malhas discordam
    ia, ib = int(np.argmin(e_a)), int(np.argmin(e_b))
    for i, t, e, c in ((ia, t_a, e_a, C_A), (ib, t_b, e_b, C_B)):
        ax.plot([t[i]], [e[i]], "o", color=c, markersize=4.0,
                markeredgewidth=0.8, markeredgecolor="white", zorder=4)

    ax.annotate("", xy=(12.0, 1.25e50), xytext=(0.5, 1.25e50),
                arrowprops=dict(arrowstyle="<->", color=C_MUTED,
                                linewidth=0.7, shrinkA=0, shrinkB=0))
    ax.text(6.2, 1.45e50, r"$m=1$ disruption", fontsize=6.2, color=C_MUTED,
            ha="center", va="bottom")
    ax.text(78.0, 5.2e45,
            "minima differ by a factor of 25:\nthe destruction converges, "
            "the residue does not",
            fontsize=6.0, color=C_MUTED, ha="right", va="bottom",
            linespacing=1.3)

    ax.set_yscale("log")
    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$E_{\rm mag}$ (erg)")
    ax.set_xlim(0, 80); ax.set_ylim(3.5e45, 3.4e50)
    fig.savefig(OUT / "fig_decay.pdf")
    plt.close(fig)
    print("  fig_decay.pdf")


def fig_probe():
    """The time step of the probe pair.

    The control, carrying the literature's weak exterior dipole, holds a step
    of a few times 1e-3 s and reaches the target. The comparable-energy
    configuration peaks at t = 0.05 s and then loses the step in plateaus until
    the subcycle limit ends the run at t = 0.205 s. Nothing about stability is
    visible here and none is claimed: the pair is asymmetric by construction,
    and what the figure shows is the ambient constraint of Section 5.2.
    """
    rows = [l.split(",") for l in open(HERE / "probe_times.csv")
            if not l.startswith(("#", "run"))]
    fig, ax = plt.subplots(figsize=(COL_IN, COL_IN * 0.72))
    style(ax)
    for run, c, lab in (("tt192ctl", C_A, "literature geometry"),
                        ("tt192", C_B, "comparable energies")):
        t = np.array([float(r[2]) for r in rows if r[0] == run])
        t = np.sort(t)
        dt = np.diff(t) / 20.0
        ax.plot(t[1:], dt, color=c, linewidth=1.1, marker="o",
                markersize=2.4, markeredgewidth=0, label=lab)

    ax.plot([0.2046], [3.82e-4], "x", color=C_B, markersize=6,
            markeredgewidth=1.3, zorder=4)
    ax.annotate("abort: too many subcycles", xy=(0.2046, 3.82e-4),
                xytext=(0.30, 1.55e-4), fontsize=6.0, color=C_B,
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=C_B, linewidth=0.6,
                                shrinkA=1, shrinkB=3))
    ax.legend(loc="upper left", handlelength=1.4)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$\Delta t$ (s)")
    ax.set_xlim(1e-3, 4.0)
    fig.savefig(OUT / "fig_probe.pdf")
    plt.close(fig)
    print("  fig_probe.pdf")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = load(HERE / "res257" / "tt_sweep_*.csv")
    hi = load(HERE / "tt_sweep_hi_*.csv")
    print(f"{len(rows)} linhas em 257, {len(hi)} fora da grade")
    fig_ceiling(rows)
    fig_frontier(rows, hi)
    fig_tradeoff()
    fig_decay()
    fig_probe()


if __name__ == "__main__":
    main()
