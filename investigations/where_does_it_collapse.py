"""Where, spatially, does the runaway start?

Run:  scf/.venv/bin/python3 investigations/where_does_it_collapse.py <rundir>

Every scalar we have read from the logs is flat against every physical
parameter varied. Four configurations spanning a factor 2.8 in E_tor/|W| died
between 2.24 and 2.52 s; three resolutions from 64^3 to 128^3 died between
2.18 and 2.70 s; damped and undamped died 0.15 s apart. Mass is conserved to
0.02%, the centre of mass does not drift, and the CFL constraint sits within
about ten cells of the centre.

A death time that ignores every knob is either something common that has not
been identified or a property of the core, which is nearly identical in all
four stars: same rho_c, same equation of state, and B_phi = K rho varpi
vanishing on the axis. Scalars cannot separate those. Geometry can.

Three shapes, three different causes:

  spherical, centred          the core itself runs away -- structural, and
                              the field is irrelevant since it vanishes there
  a ring or column near the   the toroidal field's hoop tension squeezing
  axis                        toward the axis, which is what a mis-resolved
                              B_phi would do first
  off-centre or lopsided      an asymmetry the Cartesian mesh is imposing on
                              an axisymmetric configuration

So this measures, per plotfile: where the density maximum sits in (varpi, z),
how the density on the axis compares with the density on the equator at the
same spherical radius, and how much of the mass has moved inward.
"""

import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

import yt                                            # noqa: E402

yt.set_log_level(50)


def analyse(ds):
    lev0 = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                            dims=ds.domain_dimensions)
    rho = np.asarray(lev0["boxlib", "density"])

    lo = np.asarray(ds.domain_left_edge.to("cm"))
    hi = np.asarray(ds.domain_right_edge.to("cm"))
    n = np.asarray(ds.domain_dimensions)
    ax = [lo[d] + (hi[d] - lo[d]) / n[d] * (np.arange(n[d]) + 0.5)
          for d in range(3)]
    X = ax[0][:, None, None]
    Y = ax[1][None, :, None]
    Z = ax[2][None, None, :]
    VP = np.sqrt(X ** 2 + Y ** 2) + 0.0 * Z
    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)

    k = np.unravel_index(np.argmax(rho), rho.shape)
    rho_max = float(rho[k])
    vp_max = float(np.sqrt(ax[0][k[0]] ** 2 + ax[1][k[1]] ** 2))
    z_max = float(ax[2][k[2]])

    # shape: density on the axis against density on the equator, at matched
    # spherical radius. > 1 means prolate (matter piled along the axis, what
    # hoop tension does); < 1 means oblate.
    dx = (hi[0] - lo[0]) / n[0]
    shell = 4.0 * dx
    on_axis = VP < 1.5 * dx
    on_eq = np.abs(Z) < 1.5 * dx
    ratios = []
    for r0 in (2.0 * dx, 4.0 * dx, 6.0 * dx, 8.0 * dx):
        m = np.abs(R - r0) < shell
        a, e = rho[m & on_axis], rho[m & on_eq]
        if a.size > 2 and e.size > 2 and e.mean() > 0:
            ratios.append(a.mean() / e.mean())
    shape = float(np.mean(ratios)) if ratios else np.nan

    # how concentrated: fraction of the mass inside 1/4 of the domain radius
    cell = float(np.prod((hi - lo) / n))
    inner = R < 0.25 * float(hi[0])
    frac_in = float((rho[inner].sum() * cell) / (rho.sum() * cell))

    return rho_max, vp_max, z_max, shape, frac_in


def main():
    rundir = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    plts = sorted(p for p in rundir.glob("plt?????") if p.is_dir())
    if not plts:
        raise SystemExit(f"no plotfiles in {rundir}")

    print(f"{len(plts)} plotfiles from {rundir}\n")
    print("   t (s)   rho_max      at varpi     at z       axis/equator  "
          "M(r<R/4)")
    for p in plts:
        ds = yt.load(str(p))
        t = float(ds.current_time.to("s"))
        rho_max, vp, z, shape, frac = analyse(ds)
        print(f"  {t:6.3f}  {rho_max:.4e}  {vp:.3e}  {z:+.3e}  "
              f"{shape:11.3f}   {frac:.4f}")


if __name__ == "__main__":
    main()
