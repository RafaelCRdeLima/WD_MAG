#!/usr/bin/env python3
"""Figures for the SNOOPY Pm scan.

Three panels:
  (a) the stress time series, so the saturation window can be judged rather
      than asserted;
  (b) alpha against Pm with the fitted power law and the extrapolation to our
      star, drawn as an extrapolation and not as data;
  (c) the Maxwell/Reynolds ratio, which is the setup's own sanity check --
      MRI turbulence sits near 3-5 and a value far off means the run is not
      doing what we think.

Pm is an ordered quantity, so the five runs get a sequential single-hue ramp
light-to-dark rather than five categorical colours.

Runs on partial data: anything without a saturated window is drawn in the time
series and omitted from the fit.
"""

import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyse_pm_scan import (RUNS, B0, T_SAT, ORBIT, load, mean,
                             block_error, drifting)

PM_STAR = 746.0

# sequential ramp, light -> dark, one hue
RAMP = ["#BBD3E8", "#8CB3D6", "#5B8FC0", "#31699F", "#173F6B"]
INK, MUTED, ACCENT = "#1A1A1A", "#6B6B6B", "#B0561A"


def collect():
    out = []
    for d in sorted(RUNS.glob("pm*"), key=lambda p: int(re.sub(r"\D", "", p.name))):
        tv = d / "timevar"
        if not tv.exists():
            continue
        rows = load(tv)
        if not rows:
            continue
        pm = int(re.sub(r"\D", "", d.name))
        sat = [r for r in rows if r[0] >= T_SAT]
        rec = dict(pm=pm, rows=rows, saturated=len(sat) >= 10)
        if rec["saturated"]:
            w = [r[3] - r[4] for r in sat]
            rec["W"] = mean(w) / B0**2
            dw, bm = block_error([r[0] for r in sat], w)
            rec["dW"] = dw / B0**2
            rec["drift"] = drifting(bm)
            rey, mx = mean([r[3] for r in sat]), mean([-r[4] for r in sat])
            rec["ratio"] = mx / rey if rey else float("nan")
        out.append(rec)
    return out


def fit(pts):
    xs = [math.log10(p["pm"]) for p in pts]
    ys = [math.log10(p["W"]) for p in pts]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return slope, my - slope * mx


def main():
    data = collect()
    if not data:
        print(f"no runs under {RUNS}")
        return
    done = [d for d in data if d.get("saturated")]

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7.0, 10.4))
    colors = {d["pm"]: RAMP[min(i, len(RAMP) - 1)] for i, d in enumerate(data)}

    # ---------------------------------------------------------- (a) time series
    for d in data:
        t = [r[0] / ORBIT for r in d["rows"]]
        w = [(r[3] - r[4]) / B0**2 for r in d["rows"]]
        ax1.plot(t, w, lw=1.3, color=colors[d["pm"]], label=f"Pm = {d['pm']}")
        if t:
            ax1.annotate(f"{d['pm']}", (t[-1], w[-1]), fontsize=8,
                         color=colors[d["pm"]], va="center",
                         xytext=(4, 0), textcoords="offset points")
    ax1.axvspan(T_SAT / ORBIT, 1e3, color="0.93", zorder=0)
    ax1.text(T_SAT / ORBIT + 0.4, 0.04, "averaging window", fontsize=8,
             color=MUTED, transform=ax1.get_xaxis_transform())
    ax1.set_yscale("log")
    ax1.set_xlim(0, max((r[0] for d in data for r in d["rows"]), default=1) / ORBIT)
    ax1.set_xlabel("orbits")
    ax1.set_ylabel(r"$W/B_0^2$")
    ax1.set_title("(a)  Angular-momentum flux, per orbit", fontsize=10, loc="left")
    ax1.grid(alpha=0.25, which="both", lw=0.5)
    ax1.legend(fontsize=8, ncol=2, framealpha=0.95, loc="lower right")

    # ---------------------------------------------------------- (b) alpha(Pm)
    if len(done) >= 2:
        pms = [d["pm"] for d in done]
        ws = [d["W"] for d in done]
        es = [d["dW"] for d in done]
        ax2.errorbar(pms, ws, yerr=es, fmt="o", ms=7, color="#31699F",
                     ecolor=MUTED, elinewidth=1.2, capsize=3, zorder=3,
                     label="measured")
        if len(done) >= 3:
            slope, inter = fit(done)
            gx = [10 ** (math.log10(min(pms)) + i * (math.log10(PM_STAR) -
                  math.log10(min(pms))) / 100) for i in range(101)]
            gy = [10**inter * x**slope for x in gx]
            infit = [(x, y) for x, y in zip(gx, gy) if x <= max(pms)]
            beyond = [(x, y) for x, y in zip(gx, gy) if x >= max(pms)]
            ax2.plot([p[0] for p in infit], [p[1] for p in infit], "-",
                     lw=1.6, color=ACCENT, zorder=2,
                     label=rf"fit, slope ${slope:.2f}$")
            ax2.plot([p[0] for p in beyond], [p[1] for p in beyond], ":",
                     lw=1.6, color=ACCENT, zorder=2, label="extrapolation")
            ax2.axvline(PM_STAR, color=INK, lw=1.1, ls="--")
            ax2.annotate(f"our star\nPm $\\approx$ {PM_STAR:.0f}",
                         (PM_STAR, min(ws)), fontsize=8, color=INK,
                         ha="right", va="bottom",
                         xytext=(-6, 0), textcoords="offset points")
            ax2.axvspan(max(pms), 2 * PM_STAR, color="0.95", zorder=0)
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_xlim(0.7, 2 * PM_STAR)
        ax2.set_xlabel(r"magnetic Prandtl number  $\mathrm{Pm}$")
        ax2.set_ylabel(r"$\langle W\rangle/B_0^2$")
        ax2.set_title("(b)  Transport against Pm — the shaded region is "
                      "extrapolated, not measured", fontsize=10, loc="left")
        ax2.grid(alpha=0.25, which="both", lw=0.5)
        ax2.legend(fontsize=8, framealpha=0.95, loc="upper left")
    else:
        ax2.text(0.5, 0.5, "not enough saturated runs yet",
                 ha="center", va="center", color=MUTED, transform=ax2.transAxes)
        ax2.set_xticks([]); ax2.set_yticks([])

    # ---------------------------------------------------------- (c) Maxwell/Reynolds
    if done:
        ax3.axhspan(3, 5, color="0.93", zorder=0)
        ax3.text(0.98, 4, "MRI turbulence, 3–5", fontsize=8, color=MUTED,
                 ha="right", va="center", transform=ax3.get_yaxis_transform())
        ax3.plot([d["pm"] for d in done], [d["ratio"] for d in done],
                 "o-", ms=7, lw=1.5, color="#31699F")
        ax3.set_xscale("log")
        ax3.set_xlabel(r"$\mathrm{Pm}$")
        ax3.set_ylabel("Maxwell / Reynolds")
        ax3.set_title("(c)  Sanity check on the setup, not a result",
                      fontsize=10, loc="left")
        ax3.grid(alpha=0.25, which="both", lw=0.5)
    else:
        ax3.text(0.5, 0.5, "no saturated runs yet", ha="center", va="center",
                 color=MUTED, transform=ax3.transAxes)
        ax3.set_xticks([]); ax3.set_yticks([])

    fig.tight_layout()
    out = Path(__file__).parent / "pm_scan"
    fig.savefig(f"{out}.pdf")
    fig.savefig(f"{out}.png", dpi=150)
    print(f"wrote {out}.pdf and .png  ({len(done)}/{len(data)} runs saturated)")


if __name__ == "__main__":
    main()
