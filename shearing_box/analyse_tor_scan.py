#!/usr/bin/env python3
"""Transport against net toroidal flux, at our star's field ratio.

Reads runs/tor_scan/by*/timevar, with by0 = 0 taken from runs/pm_scan/pm04,
which is the same Re, Pm and q with a purely vertical field.

Read the scope note in scan_toroidal.sh before using these numbers. The
azimuthal MRI wavelength is ~6.5 v_Ay, which exceeds Lz = 1 at every ratio in
this scan, so the non-axisymmetric toroidal-field mode cannot grow here. What
is measured is how a net toroidal field modifies the VERTICAL-field MRI --
which MInIT has no term for, since its k_MRI is built from v_Az alone -- and
not which instability would dominate in the star.
"""

import math
import re
from pathlib import Path

from analyse_pm_scan import B0, T_SAT, ORBIT, load, mean, std

HERE = Path(__file__).parent
BZ = 0.1


def summarise(tv):
    rows = load(tv)
    if not rows:
        return None
    sat = [r for r in rows if r[0] >= T_SAT]
    if len(sat) < 10:
        return dict(partial=True, t_end=rows[-1][0])
    w = [r[3] - r[4] for r in sat]
    rey, mx = mean([r[3] for r in sat]), mean([-r[4] for r in sat])
    return dict(partial=False, t_end=rows[-1][0],
                W=mean(w) / B0**2, dW=std(w) / math.sqrt(len(w)) / B0**2,
                em=mean([r[2] for r in sat]),
                ratio=mx / rey if rey else float("nan"))


def main():
    srcs = [(0.0, HERE / "runs" / "pm_scan" / "pm04" / "timevar")]
    for d in sorted((HERE / "runs" / "tor_scan").glob("by*")):
        by = float(re.sub(r"^by", "", d.name).replace("p", "."))
        srcs.append((by, d / "timevar"))

    pts = []
    for by, tv in srcs:
        if not tv.exists():
            continue
        r = summarise(tv)
        if r:
            pts.append((by, r))
    pts.sort(key=lambda p: p[0])
    if not pts:
        print("no runs yet")
        return

    print(f"Re = 1000, Pm = 4, q = 1.5, bz0 = {BZ}. Average over t > {T_SAT:.0f}.")
    print("by0 = 0 is runs/pm_scan/pm04, the vertical-only reference.\n")
    print(f"{'by0':>6} {'B_tor/B_pol':>12} {'orbits':>7} {'W/B0^2':>16} "
          f"{'<em>':>9} {'Max/Rey':>8}")
    for by, r in pts:
        if r["partial"]:
            print(f"{by:>6.2f} {by/BZ:>12.1f} {r['t_end']/ORBIT:>7.1f}   still spinning up")
            continue
        print(f"{by:>6.2f} {by/BZ:>12.1f} {r['t_end']/ORBIT:>7.1f} "
              f"{r['W']:>9.4g} +-{r['dW']:>4.2f} {r['em']:>9.4g} {r['ratio']:>8.2f}")

    done = [(by, r) for by, r in pts if not r["partial"]]
    if len(done) >= 2:
        base = next((r["W"] for by, r in done if by == 0.0), None)
        if base:
            print(f"\nrelative to the vertical-only reference ({base:.1f}):")
            for by, r in done:
                if by:
                    print(f"  ratio {by/BZ:>4.1f}:  {r['W']/base:>5.2f}x")
            print("\nA ratio near 1 across the scan would mean the toroidal field")
            print("neither helps nor hinders the vertical-field MRI, and MInIT's")
            print("omission of it is harmless at our field geometry. A strong")
            print("departure is a term the closure is missing.")


if __name__ == "__main__":
    main()
