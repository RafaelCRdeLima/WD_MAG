#!/usr/bin/env python3
"""Spin-down time from an actual torque, replacing the H ~ R_eq envelope.

What this replaces
------------------
Report III estimates the MRI braking time as R^2/(alpha c_s H) with H ~ R_eq,
gets ~130 s against a 78 s run, and concludes the differential rotation should
already have been erased. H ~ R_eq was the crudest input in that chain, and the
conclusion is sharp enough to be worth doing properly.

What this does instead
----------------------
The turbulent stress exerts a torque across the cylinder that separates the
inner half of the mass from the outer half. The inner region spins down on

    tau = L_z,inner / G(r_half),     G = 2 pi r_half^2 * Z(r_half) * W

with Z the vertical extent of the star at that radius and W the stress. Both
L_z,inner and the shell definitions come from the diagnostic itself
(tools/fbtbp.cpp): Om_core is mass-weighted inside 0.15 R_eq, Om_mid over
0.45-0.55, Om_out over 0.65-0.75, and r_half is the cylindrical radius
enclosing half the mass.

The stress comes from the box, converted out of Alfven units:

    W = (W/B0^2)_box * B_z^2 / (4 pi)

and (W/B0^2)_box is taken at the LOCAL shear rather than at Keplerian, which
is the second improvement: q at r_half is about 1 for the Komatsu law, where
the q scan measured 22.8 rather than the 25.2 of q = 1.5.

Assumptions that remain, and their direction
--------------------------------------------
1. B_z uniform, from the volume-typical value sqrt(8 pi E_pol/V)/sqrt(2). The
   field is certainly not uniform; if it is weaker at r_half than on average
   the torque falls and tau rises.
2. The star is an ellipsoid of semi-axes R_eq and R_pol, so
   Z(r) = 2 R_pol sqrt(1 - (r/R_eq)^2).
3. r_half is taken from the 192^3 diagnostic, which carries it; the 256^3 CSV
   does not. The two stars are the same object at two resolutions, so this is
   a small error, but it is an import.
4. The box value is at Pm = 4 against the star's ~750. The Pm scan says that is
   a factor 1.3-1.9, which is folded in as a range rather than a point.
"""

import csv
import math
from pathlib import Path

HERE = Path(__file__).parent
MSUN = 1.98892e33
C_S = 3.0e9


def load(fn):
    with (HERE / fn).open() as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    return [{k: float(v) for k, v in r.items()} for r in rows]


def main():
    d256 = load("bt_bp_256_long.csv")
    r192 = load("rotation_192.csv")
    s = d256[-1]                      # t = 78 s
    t = s["t"]

    R_eq = s["R_eq_e8"] * 1e8
    R_pol = s["Rpol_over_Req"] * R_eq
    R_vol = s["R_vol_e8"] * 1e8
    V = 4.0 / 3.0 * math.pi * R_vol**3
    rho = s["M_Msun"] * MSUN / V

    # r_half: the 256^3 CSV does not carry it; take the 192^3 value at the
    # nearest time and scale by that run's R_eq if needed.
    late192 = min(r192, key=lambda r: abs(r["t"] - 60.0))
    r_half = late192["r_half"] * 1e8

    B_pol = math.sqrt(8.0 * math.pi * s["E_pol"] / V)
    B_z = B_pol / math.sqrt(2.0)

    # Komatsu j-constant: q = 2 w^2/(A^2 + w^2), A = 0.468 R_eq
    A = 0.468 * R_eq
    q_half = 2.0 * r_half**2 / (A**2 + r_half**2)

    # box transport at the local shear, from the q scan (Re=1000, Pm=4)
    QTAB = [(0.5, 18.10), (1.0, 22.84), (1.5, 25.24), (1.9, 30.73)]
    lo = max((p for p in QTAB if p[0] <= q_half), default=QTAB[0])
    hi = min((p for p in QTAB if p[0] >= q_half), default=QTAB[-1])
    if hi[0] == lo[0]:
        WB0 = lo[1]
    else:
        f = (q_half - lo[0]) / (hi[0] - lo[0])
        WB0 = lo[1] + f * (hi[1] - lo[1])

    print(f"State at t = {t:.1f} s, 256^3\n")
    print(f"  R_eq        {R_eq:.3e} cm")
    print(f"  R_pol       {R_pol:.3e} cm   (R_pol/R_eq = {s['Rpol_over_Req']:.3f})")
    print(f"  r_half      {r_half:.3e} cm   = {r_half/R_eq:.3f} R_eq   [from 192^3]")
    print(f"  rho_mean    {rho:.3e} g/cm^3")
    print(f"  B_pol,rms   {B_pol:.3e} G  ->  B_z = {B_z:.3e} G")
    print(f"  q(r_half)   {q_half:.3f}   -> box gives W/B0^2 = {WB0:.1f}")

    Z = 2.0 * R_pol * math.sqrt(max(0.0, 1.0 - (r_half / R_eq) ** 2))
    W = WB0 * B_z**2 / (4.0 * math.pi)
    G = 2.0 * math.pi * r_half**2 * Z * W
    Lz_in = s["Lz_inner"]

    print(f"\n  Z(r_half)   {Z:.3e} cm")
    print(f"  stress W    {W:.3e} erg/cm^3")
    print(f"  torque G    {G:.3e} erg")
    print(f"  L_z,inner   {Lz_in:.3e} g cm^2/s")

    tau = Lz_in / G
    print(f"\n  tau = L_z,inner / G = {tau:.3e} s = {tau/78.0:.1f} x the run\n")

    print("  with the Pm correction (Pm 4 -> ~750 is a factor 1.3 to 1.9):")
    for f, lab in ((1.3, "saturating"), (1.9, "trend continues")):
        print(f"    {lab:>16}: tau = {tau/f:.0f} s = {tau/f/78.0:.1f} x the run")

    # the Report III envelope, recomputed the way it was actually written
    P = 7.812e24                       # degenerate pressure at rho_mean
    beta = P / (B_z**2 / (8.0 * math.pi))
    alpha = 2.0 / beta * 32.0
    nu_t = alpha * C_S * R_eq
    print(f"\n  Report III envelope, R_eq^2/(alpha c_s H) with H = R_eq:")
    print(f"    beta = {beta:.2e}, alpha = {alpha:.2e}  ->  tau = "
          f"{R_eq**2/nu_t:.0f} s")
    print(f"    the torque calculation is {tau/(R_eq**2/nu_t):.1f}x longer.")

    # THE COMPARISON THAT DECIDES IT
    #
    # Not against the star's total angular-momentum loss -- that is L_z leaving
    # the rho > 1e5 region altogether, a different thing from internal
    # redistribution. The like-for-like comparison is the MRI torque against
    # the observed rate of change of the INNER region's angular momentum, since
    # both are the rate at which L_z crosses r_half.
    def fit(key, lo, hi):
        pts = [(r["t"], r[key]) for r in d256 if lo <= r["t"] <= hi]
        n = len(pts)
        mx = sum(p[0] for p in pts) / n
        my = sum(p[1] for p in pts) / n
        return sum((x - mx) * (y - my) for x, y in pts) / \
               sum((x - mx) ** 2 for x, y in pts)

    for lo, hi in ((30.0, 78.0), (60.0, 78.0)):
        din, dout = fit("Lz_inner", lo, hi), fit("Lz_outer", lo, hi)
        print(f"\n  Observed over {lo:.0f}-{hi:.0f} s:")
        print(f"    dL_z,inner/dt = {din:+.3e} erg"
              f"   ({'INWARD transport, steepening' if din > 0 else 'outward'})")
        print(f"    dL_z,outer/dt = {dout:+.3e} erg")
        print(f"    predicted MRI torque {G:.2e} is {abs(G/din):.1f}x that, "
              f"and {'OPPOSITE' if din > 0 else 'the same'} in sign")


if __name__ == "__main__":
    main()
