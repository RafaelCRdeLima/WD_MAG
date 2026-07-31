"""Phase 1 analysis: does the toroidal-dominated field survive?

Run:  scf/.venv/bin/python3 investigations/analyze_phase1_stability.py <rundir>

Three curves decide the Phase 1 question, and they answer different parts of
it, so all three are needed.

1. rho_c(t). The validity window. The initial condition is an SCF equilibrium
   sampled onto a Cartesian mesh with an IMPOSED poloidal component, so it is
   not exactly in equilibrium there and settles. Any statement about the field
   is only as good as the interval over which the star itself is still the
   star -- this is the same criterion as docs/teoria.md Sec 6.9.

2. E_tor/E_mag(t). Whether the field rearranges. A Tayler-unstable toroidal
   field converts toroidal energy into poloidal and kinetic energy; a stable
   one does not. Computed here from the plotfile's cell-centred B_x, B_y,
   B_z rather than from Castro's emag_density derive: that derive reads
   dat(i,j,k+1,2) with too few ghost cells and returns garbage (values of
   1e90) once the domain is decomposed in z across ranks. It is correct on a
   single rank, which is why the initial-condition check did not catch it.
   The bug is inherited from the wd_braithwaite problem and affects it too.
   etor_density is unaffected -- it never touches component 2 -- and is used
   here as an independent cross-check on the toroidal energy.

3. P_1/P_0(t), the azimuthal mode power. This is the actual signature. The
   initial condition is axisymmetric, so it has only m = 0. The Tayler
   instability is m = 1. Its growth against time, in Alfven units, is the
   measurement Phase 1 exists to make -- and if it grows, the e-folding time
   read off the exponential phase is the number that decides whether the
   configuration can sit still long enough to brake for Myr.

The mode power is computed by interpolating B onto a cylindrical grid and
transforming over phi, rather than by binning cells into annuli: binning
aliases the Cartesian mesh's own four-fold symmetry into m = 4 and, worse,
leaks it into m = 1 near the axis where annuli hold few cells.
"""

import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

import yt                                                  # noqa: E402
from scipy.interpolate import RegularGridInterpolator      # noqa: E402

yt.set_log_level(50)

N_VARPI, N_PHI, N_Z = 48, 64, 96      # cylindrical sampling for the transform
M_MAX = 4


def cylindrical_modes(ds):
    """Return (E_pol, E_tor, P_m[0..M_MAX]) from one plotfile."""
    lev0 = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                            dims=ds.domain_dimensions)
    bx = np.asarray(lev0["boxlib", "B_x"])
    by = np.asarray(lev0["boxlib", "B_y"])
    bz = np.asarray(lev0["boxlib", "B_z"])

    lo = np.asarray(ds.domain_left_edge.to("cm"))
    hi = np.asarray(ds.domain_right_edge.to("cm"))
    n = np.asarray(ds.domain_dimensions)
    ax = [lo[d] + (hi[d] - lo[d]) / n[d] * (np.arange(n[d]) + 0.5)
          for d in range(3)]

    interps = [RegularGridInterpolator(ax, f, bounds_error=False,
                                       fill_value=0.0)
               for f in (bx, by, bz)]

    # cylindrical grid inside the domain's inscribed cylinder
    vp_max = 0.85 * min(abs(lo[0]), abs(hi[0]), abs(lo[1]), abs(hi[1]))
    z_max = 0.85 * min(abs(lo[2]), abs(hi[2]))
    vp = np.linspace(vp_max / N_VARPI, vp_max, N_VARPI)
    ph = np.linspace(0.0, 2.0 * np.pi, N_PHI, endpoint=False)
    zz = np.linspace(-z_max, z_max, N_Z)
    VP, PH, ZZ = np.meshgrid(vp, ph, zz, indexing="ij")
    pts = np.stack([(VP * np.cos(PH)).ravel(), (VP * np.sin(PH)).ravel(),
                    ZZ.ravel()], axis=-1)

    BX, BY, BZ = [f(pts).reshape(VP.shape) for f in interps]
    cos_p, sin_p = np.cos(PH), np.sin(PH)
    B_vp = BX * cos_p + BY * sin_p
    B_ph = -BX * sin_p + BY * cos_p

    # azimuthal transform of the toroidal component -- the one the
    # instability feeds on
    coef = np.fft.rfft(B_ph, axis=1) / N_PHI
    w = VP[:, 0, :]                       # varpi weight for the volume element
    power = []
    for m in range(M_MAX + 1):
        power.append(float(np.sum(np.abs(coef[:, m, :]) ** 2 * w)))

    # energies from the plotfile fields, in Castro's Heaviside-Lorentz
    # convention (energy density = B^2/2, numerically equal to the Gaussian
    # B^2/8pi since B' = B/sqrt(4 pi))
    dv = float((ds.domain_width[0] / ds.domain_dimensions[0]) ** 3)
    X = ax[0][:, None, None]
    Y = ax[1][None, :, None]
    vpg = np.sqrt(X ** 2 + Y ** 2)
    b_tor = np.where(vpg > 0, (-Y * bx + X * by) / np.maximum(vpg, 1e-30), 0.0)
    E_mag = 0.5 * float(np.sum(bx ** 2 + by ** 2 + bz ** 2)) * dv
    E_tor = 0.5 * float(np.sum(b_tor ** 2)) * dv

    ad = ds.all_data()
    E_tor_castro = float(ad.quantities.total_quantity(
        ("boxlib", "etor_density"))) * dv
    rho_c = float(ad.quantities.extrema(("boxlib", "density"))[1])
    return E_mag, E_tor, E_tor_castro, rho_c, np.array(power)


def main():
    rundir = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    plts = sorted(rundir.glob("plt?????"))
    if not plts:
        raise SystemExit(f"no plotfiles in {rundir}")

    t_alfven = 1.2567          # s, reported by the problem at initialisation
    print(f"{len(plts)} plotfiles, t_alfven = {t_alfven} s\n")
    print("   t (s)    t/t_A   rho_c/rho_c0  E_tor/E_mag  E_tor/E_tor0  "
          "cross   P1/P0      P2/P0")

    rows = []
    rho0 = None
    for p in plts:
        ds = yt.load(str(p))
        t = float(ds.current_time.to("s"))
        E_mag, E_tor, E_tor_castro, rho_c, power = cylindrical_modes(ds)
        if rho0 is None:
            rho0, etor0 = rho_c, E_tor
        p0 = max(power[0], 1e-300)
        cross = E_tor_castro / max(E_tor, 1e-300)
        rows.append((t, rho_c, E_mag, E_tor, E_tor / max(E_mag, 1e-300),
                     power[1] / p0, power[2] / p0))
        print(f"  {t:7.3f}  {t/t_alfven:6.3f}   {rho_c/rho0:10.4f}  "
              f"{E_tor/max(E_mag,1e-300):10.4f}  {E_tor/etor0:11.4f}  "
              f"{cross:5.3f}  {power[1]/p0:.3e}  {power[2]/p0:.3e}")

    out = rundir / "phase1_stability.csv"
    with out.open("w") as f:
        f.write("# Phase 1 stability analysis, t_alfven = "
                f"{t_alfven} s\n")
        f.write("time_s,rho_c,E_mag,E_tor,E_tor_over_Emag,P1_over_P0,"
                "P2_over_P0\n")
        for r in rows:
            f.write(",".join(f"{v:.6e}" for v in r) + "\n")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
