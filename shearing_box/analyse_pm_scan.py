#!/usr/bin/env python3
"""alpha(Pm) from the SNOOPY Pm scan, with the literature check.

Reads runs/pm_scan/pm*/timevar and time-averages the turbulent stresses over
the saturated window.

Conventions
-----------
SNOOPY works in Alfven units, so b IS v_A and the magnetic energy is b^2/2.
The total angular-momentum flux is

    W = <v_x v_y> - <b_x b_y>            (Reynolds + Maxwell)

with the Maxwell term entering as MINUS <bxby>, which is why bxby comes out
negative in a run that is transporting outwards.

There is no sound speed in an incompressible box, so the usual
alpha = W / c_s^2 has no meaning here. The net-flux convention is used
instead: normalise by the imposed field, W / B0^2 with B0 = bz0.

What the answer is for
----------------------
Our star sits at Pm ~ 750 (investigations/magnetic_prandtl.py), which no DNS
reaches. The scan is therefore a BOUND, not a measurement: transport rises with
Pm in the disc literature, so the largest Pm we can afford gives a lower bound
on the transport at 750. The fitted slope is what gets extrapolated, and the
extrapolation has to be quoted with the result.
"""

import math
import re
from pathlib import Path

RUNS = Path(__file__).parent / "runs" / "pm_scan"
B0 = 0.1                    # bz0 in the stock MRI config
T_SAT = 60.0                # discard before this; saturation is by t ~ 100 (16 orbits)
ORBIT = 2.0 * math.pi

COLS = {"t": 0, "ev": 1, "em": 2, "vxvy": 9, "bxby": 16}


def load(path):
    rows = []
    with path.open() as fh:
        for line in fh:
            if line.startswith("t\t") or not line.strip():
                continue
            f = line.split()
            if len(f) < 17:
                continue
            try:
                rows.append([float(f[COLS[k]]) for k in ("t", "ev", "em", "vxvy", "bxby")])
            except ValueError:
                continue
    return rows


def mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")


def std(xs):
    if len(xs) < 2:
        return float("nan")
    m = mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def main():
    dirs = sorted(RUNS.glob("pm*"), key=lambda p: int(re.sub(r"\D", "", p.name)))
    if not dirs:
        print(f"no runs under {RUNS}")
        return

    print(f"Saturated average over t > {T_SAT:.0f} "
          f"({T_SAT/ORBIT:.1f} orbits), B0 = {B0}\n")
    print(f"{'Pm':>4} {'orbits':>7} {'n':>5} {'<em>':>10} {'Maxwell':>10} "
          f"{'Reynolds':>10} {'W':>10} {'W/B0^2':>9} {'Max/Rey':>8} {'+-':>8}")

    pts = []
    for d in dirs:
        tv = d / "timevar"
        if not tv.exists():
            print(f"{d.name:>4}  (no timevar yet)")
            continue
        rows = load(tv)
        if not rows:
            print(f"{d.name:>4}  (empty)")
            continue
        pm = int(re.sub(r"\D", "", d.name))
        t_end = rows[-1][0]
        sat = [r for r in rows if r[0] >= T_SAT]
        if len(sat) < 10:
            print(f"{pm:>4} {t_end/ORBIT:>7.1f} {len(sat):>5}   still spinning up")
            continue
        em = mean([r[2] for r in sat])
        rey = mean([r[3] for r in sat])
        max_ = mean([-r[4] for r in sat])
        w_series = [r[3] - r[4] for r in sat]
        w, dw = mean(w_series), std(w_series) / math.sqrt(len(w_series))
        print(f"{pm:>4} {t_end/ORBIT:>7.1f} {len(sat):>5} {em:>10.4g} "
              f"{max_:>10.4g} {rey:>10.4g} {w:>10.4g} {w/B0**2:>9.4g} "
              f"{max_/rey if rey else float('nan'):>8.2f} {dw/B0**2:>8.3g}")
        pts.append((pm, w / B0**2, dw / B0**2))

    if len(pts) >= 3:
        # least squares on log10(W/B0^2) vs log10(Pm)
        xs = [math.log10(p[0]) for p in pts]
        ys = [math.log10(p[1]) for p in pts if p[1] > 0]
        if len(ys) == len(xs):
            n = len(xs)
            mx, my = mean(xs), mean(ys)
            sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
            sxx = sum((x - mx) ** 2 for x in xs)
            slope = sxy / sxx
            inter = my - slope * mx
            print(f"\npower law:  W/B0^2 = {10**inter:.3g} * Pm^{slope:.3f}")
            print("literature for Keplerian net-flux boxes: transport rises with Pm,")
            print("Lesur & Longaretti (2007) and Fromang et al. (2007) find slopes")
            print("of order 0.5-1 over Pm = 1-16. A slope near zero or negative here")
            print("would mean our setup, not the physics, is wrong.")
            pm_star = 746.0
            print(f"\nEXTRAPOLATION to our star's Pm = {pm_star:.0f}: "
                  f"W/B0^2 ~ {10**inter * pm_star**slope:.4g}")
            print(f"  which is {pm_star**slope:.1f}x the Pm=1 value. This is an")
            print("  extrapolation of ~1.7 decades beyond the fitted range and must")
            print("  be quoted as such -- it is a lower bound only if the trend does")
            print("  not saturate, which these data cannot establish.")


if __name__ == "__main__":
    main()
