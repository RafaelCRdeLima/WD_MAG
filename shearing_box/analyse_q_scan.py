#!/usr/bin/env python3
"""alpha(q) from the q scan, with MInIT's own coefficient alongside.

q enters our problem twice and this script reports both:

  - the MEASURED transport, W/B0^2 from the box;
  - MInIT's THEORETICAL Maxwell coefficient, alpha^MRI = 1 - 4/q from
    Pessah & Chan (2008), which the subgrid model uses directly.

Having them side by side is the point. If the measured transport tracks the
theoretical coefficient across our star's range of q, MInIT is being fed a
shear dependence that its own closure already knows about; if it does not,
that is a discrepancy worth chasing before the closure is trusted.

Our star's rotation law is Komatsu j-constant, so q = 2w^2/(A^2+w^2) runs from
0 on the axis to 2 far out, with q = 1 at w = A = 0.468 R_eq.

q = 1.5 is not re-run here; it comes from runs/pm_scan/pm04, the same Re and Pm.
"""

import math
import re
from pathlib import Path

from analyse_pm_scan import (B0, T_SAT, ORBIT, load, mean, std,
                             block_error, drifting)

HERE = Path(__file__).parent
QRUNS = HERE / "runs" / "q_scan"
PM04 = HERE / "runs" / "pm_scan" / "pm04"       # this is q = 1.5, Pm = 4


def summarise(tv):
    rows = load(tv)
    if not rows:
        return None
    sat = [r for r in rows if r[0] >= T_SAT]
    if len(sat) < 10:
        return dict(partial=True, t_end=rows[-1][0], n=len(sat))
    w = [r[3] - r[4] for r in sat]
    dw, bm = block_error([r[0] for r in sat], w)
    rey, mx = mean([r[3] for r in sat]), mean([-r[4] for r in sat])
    return dict(partial=False, t_end=rows[-1][0], drift=drifting(bm), blocks=bm, n=len(sat),
                W=mean(w) / B0**2, dW=dw / B0**2,
                em=mean([r[2] for r in sat]),
                ratio=mx / rey if rey else float("nan"))


def main():
    pts = []
    for d in sorted(QRUNS.glob("q*")):
        q = float(re.sub(r"^q", "", d.name).replace("p", "."))
        tv = d / "timevar"
        if tv.exists():
            r = summarise(tv)
            if r:
                pts.append((q, r))
    if PM04.joinpath("timevar").exists():
        r = summarise(PM04 / "timevar")
        if r:
            pts.append((1.5, r))
    pts.sort(key=lambda p: p[0])

    if not pts:
        print(f"no runs under {QRUNS}")
        return

    print(f"Re = 1000, Pm = 4, net vertical flux. Average over t > {T_SAT:.0f}.\n")
    print(f"{'q':>5} {'orbits':>7} {'W/B0^2':>16} {'<em>':>9} {'Max/Rey':>8} "
          f"{'1-4/q':>9}  MInIT")
    for q, r in pts:
        a_mri = 1.0 - 4.0 / q if q else float("-inf")
        if r["partial"]:
            print(f"{q:>5.2f} {r['t_end']/ORBIT:>7.1f}   still spinning up")
            continue
        print(f"{q:>5.2f} {r['t_end']/ORBIT:>7.1f} "
              f"{r['W']:>9.4g} +-{r['dW']:>5.2f} {r['em']:>9.4g} "
              f"{r['ratio']:>8.2f} {a_mri:>9.3f}"
              + ("   DRIFTING" if r.get("drift") else ""))

    done = [(q, r) for q, r in pts if not r["partial"]]
    if len(done) >= 3:
        xs = [math.log10(q) for q, _ in done]
        ys = [math.log10(r["W"]) for _, r in done if r["W"] > 0]
        if len(ys) == len(xs):
            n = len(xs)
            mx, my = mean(xs), mean(ys)
            slope = (sum((x - mx) * (y - my) for x, y in zip(xs, ys))
                     / sum((x - mx) ** 2 for x in xs))
            print(f"\nW/B0^2 ~ q^{slope:.2f} over q = {done[0][0]} to {done[-1][0]}")
            print("Expect a positive slope: the MRI feeds on shear, and the linear")
            print("growth rate is gamma = q Omega/2, so transport should rise with q.")

        print("\nWhat this feeds:")
        print("  MInIT's alpha^MRI = 1 - 4/q is STRONGLY q-dependent and changes")
        print("  sign at q = 4, far outside our range, so it stays negative")
        print("  throughout and diverges as q -> 0. Near the rotation axis, where")
        print("  our star has q -> 0, the closure is therefore singular and the")
        print("  published implementation zeroes the MRI term for q <= 0. Whether")
        print("  the measured transport also collapses there is exactly what the")
        print("  low-q end of this scan tests.")


if __name__ == "__main__":
    main()
