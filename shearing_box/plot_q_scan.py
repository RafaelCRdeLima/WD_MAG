#!/usr/bin/env python3
"""Figures for the q scan.

  (a) transport against shear;
  (b) the first direct test of a MInIT coefficient against our own data.

Panel (b) is the point. MInIT takes its Maxwell coefficient from theory,
alpha^MRI = 1 - 4/q (Pessah & Chan 2008), which fixes the Maxwell-to-Reynolds
ratio as a function of shear alone. Plotting the measured ratio against
|1 - 4/q| tests that prediction: a slope of 1 means the shear dependence is
right, and any offset from the diagonal is a normalisation question.

The slope is convention-independent -- a constant factor between our stress
definition and theirs shifts the intercept and leaves the slope alone -- so it
is the slope that carries the physics here.
"""

import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyse_pm_scan import B0, T_SAT, ORBIT, load, mean, std

HERE = Path(__file__).parent
INK, MUTED, ACCENT, CAUTION = "#1A1A1A", "#6B6B6B", "#31699F", "#B0561A"


def collect():
    pts = []
    srcs = [(float(re.sub(r"^q", "", d.name).replace("p", ".")), d / "timevar")
            for d in sorted((HERE / "runs" / "q_scan").glob("q*"))]
    srcs.append((1.5, HERE / "runs" / "pm_scan" / "pm04" / "timevar"))
    for q, tv in srcs:
        if not tv.exists():
            continue
        rows = load(tv)
        sat = [r for r in rows if r[0] >= T_SAT]
        if len(sat) < 10:
            continue
        w = [r[3] - r[4] for r in sat]
        rey, mx = mean([r[3] for r in sat]), mean([-r[4] for r in sat])
        pts.append(dict(q=q, W=mean(w) / B0**2,
                        dW=std(w) / math.sqrt(len(w)) / B0**2,
                        ratio=mx / rey))
    pts.sort(key=lambda p: p["q"])
    return pts


def loglinfit(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sl = (sum((x - mx) * (y - my) for x, y in zip(xs, ys))
          / sum((x - mx) ** 2 for x in xs))
    return sl, my - sl * mx


def main():
    pts = collect()
    if len(pts) < 3:
        print("not enough saturated runs")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.5))

    # ---------------------------------------------------------- (a) W(q)
    q = [p["q"] for p in pts]
    w = [p["W"] for p in pts]
    e = [p["dW"] for p in pts]
    ax1.errorbar(q, w, yerr=e, fmt="o-", ms=7, lw=1.5, color=ACCENT,
                 ecolor=MUTED, capsize=3)
    s, i = loglinfit([math.log10(x) for x in q], [math.log10(y) for y in w])
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel(r"shear parameter  $q = -\,d\ln\Omega/d\ln\varpi$")
    ax1.set_ylabel(r"$\langle W\rangle/B_0^2$")
    ax1.set_title(rf"(a)  Transport rises with shear, $\propto q^{{{s:.2f}}}$",
                  fontsize=10, loc="left")
    ax1.grid(alpha=0.25, which="both", lw=0.5)
    ax1.annotate("our star spans\nthis whole range",
                 xy=(0.5, min(w)), fontsize=8, color=MUTED,
                 xytext=(4, 10), textcoords="offset points")

    # ---------------------------------------------------------- (b) MInIT test
    x = [abs(1.0 - 4.0 / p["q"]) for p in pts]
    y = [p["ratio"] for p in pts]
    s2, i2 = loglinfit([math.log10(v) for v in x], [math.log10(v) for v in y])

    gx = [10 ** (math.log10(min(x)) + k * (math.log10(max(x)) - math.log10(min(x))) / 60)
          for k in range(61)]
    ax2.plot(gx, [10**i2 * v**s2 for v in gx], "-", lw=1.6, color=CAUTION,
             label=rf"fit: slope ${s2:.2f}$, offset ${10**i2:.1f}\times$")
    ax2.plot(gx, [10**i2 * v for v in gx], "--", lw=1.2, color=MUTED,
             label="slope 1 (MInIT's shear law)")
    ax2.plot(x, y, "o", ms=8, color=ACCENT, zorder=3)
    for p, xv, yv in zip(pts, x, y):
        ax2.annotate(f"q={p['q']:g}", (xv, yv), fontsize=8, color=INK,
                     xytext=(6, -3), textcoords="offset points")
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel(r"$|\alpha^{\rm MRI}| = |1 - 4/q|$   (Pessah & Chan 2008)")
    ax2.set_ylabel("measured Maxwell / Reynolds")
    ax2.set_title("(b)  MInIT's shear law against our data", fontsize=10, loc="left")
    ax2.grid(alpha=0.25, which="both", lw=0.5)
    ax2.legend(fontsize=8, framealpha=0.95, loc="upper left")

    fig.tight_layout()
    out = HERE / "q_scan"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=150)
    print(f"wrote {out}.pdf and .png")
    print(f"  (a) W/B0^2 ~ q^{s:.2f}")
    print(f"  (b) Max/Rey = {10**i2:.2f} |1-4/q|^{s2:.2f}   "
          f"(MInIT predicts exponent 1)")


if __name__ == "__main__":
    main()
