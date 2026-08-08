#!/usr/bin/env python3
"""Is the MRI resolved? lambda_MRI/dx from the VOLUME-TYPICAL vertical field.

The number this replaces
------------------------
The claim on record is that lambda_MRI/dx crosses 6 at t ~ 2-3 s and reaches 36
at 256^3. That was computed from the PEAK B_z. A peak is not what sets the
wavelength over the bulk of the star: the MRI has to fit where the field
typically is, not where it happens to be strongest.

Why no plotfiles are needed
---------------------------
E_pol is a volume integral of B_pol^2/8pi over the star, and it is already in
bt_bp_256_long.csv for all 210 samples. So

    B_pol,rms = sqrt(8 pi E_pol / V)

is exactly a volume-typical field, obtained without touching a slice. V comes
from R_vol in the same row, so the two are consistent by construction.

The assumptions, each of which is a real limitation
---------------------------------------------------
1. B_z from B_pol. The MRI wavelength is built from the field along the
   rotation axis, and the CSV carries only the poloidal energy, which mixes
   B_varpi and B_z. Three cases are reported rather than one:
       upper   B_z = B_pol         (all poloidal flux vertical)
       central B_z = B_pol/sqrt(2) (the two components comparable)
       lower   B_z = B_pol/sqrt(3) (isotropic poloidal field)
   The spread between them is a factor 1.73, which is smaller than the
   peak-vs-typical effect this script exists to measure.

2. rho_mean, not local rho. v_A = B/sqrt(4 pi rho) is evaluated at the mean
   density M/V. In the core rho is higher and v_A lower, so the central estimate
   is optimistic there and pessimistic in the envelope.

3. <B^2>^(1/2) / sqrt(<rho>) is not <B/sqrt(rho)>. This is a typical-value
   estimate, not an average of the local Alfven speed.

None of these is worth more precision until the answer is known to within a
factor of a few, which is all that is needed to decide whether Q >= 6 holds.
"""

import csv
import math
from pathlib import Path

MSUN = 1.98892e33
CSV = Path(__file__).parent / "bt_bp_256_long.csv"

# 256^3 over [-9.046875e8, 8.953125e8]
DX = 1.8e9 / 256

# Q thresholds from the MRI literature
Q_LINEAR = 6.0    # below this the linear mode is not represented at all
Q_TURB = 15.0     # converged turbulence is normally held to need 15-20

# B_z = B_pol / SPLIT
CASES = {"upper": 1.0, "central": math.sqrt(2.0), "lower": math.sqrt(3.0)}


def load():
    with CSV.open() as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    return [{k: float(v) for k, v in r.items()} for r in rows]


def lambda_over_dx(row, split, omega_key="Om_mean"):
    R = row["R_vol_e8"] * 1e8
    V = 4.0 / 3.0 * math.pi * R**3
    rho = row["M_Msun"] * MSUN / V
    B_pol = math.sqrt(8.0 * math.pi * row["E_pol"] / V)
    B_z = B_pol / split
    v_Az = B_z / math.sqrt(4.0 * math.pi * rho)
    lam = 2.0 * math.pi * v_Az / row[omega_key]
    return lam / DX, B_pol, rho, lam


def main():
    rows = load()
    print(f"256^3, dx = {DX:.3e} cm, {len(rows)} samples, "
          f"t = {rows[0]['t']:.2f} to {rows[-1]['t']:.2f} s\n")

    print(f"{'t':>7} {'B_pol,rms':>11} {'rho_mean':>10} {'lam(cm)':>10} "
          f"{'Q_low':>7} {'Q_ctr':>7} {'Q_up':>7}")
    marks = [0.0, 1.0, 2.0, 3.0, 5.0, 9.0, 13.5, 20.0, 32.5, 45.0, 60.0, 78.0]
    for target in marks:
        row = min(rows, key=lambda r: abs(r["t"] - target))
        q = {k: lambda_over_dx(row, s)[0] for k, s in CASES.items()}
        _, B_pol, rho, lam = lambda_over_dx(row, CASES["central"])
        print(f"{row['t']:>7.2f} {B_pol:>11.3e} {rho:>10.3e} {lam:>10.3e} "
              f"{q['lower']:>7.2f} {q['central']:>7.2f} {q['upper']:>7.2f}")

    print()
    for name, split in CASES.items():
        qs = [(lambda_over_dx(r, split)[0], r["t"]) for r in rows]
        above6 = [t for q, t in qs if q >= Q_LINEAR]
        above15 = [t for q, t in qs if q >= Q_TURB]
        qmax, tmax = max(qs)
        frac6 = 100.0 * len(above6) / len(qs)
        print(f"{name:>8} (B_z = B_pol/{split:.3f}): "
              f"max Q = {qmax:.1f} at t = {tmax:.1f} s; "
              f"Q>={Q_LINEAR:.0f} in {frac6:.0f}% of samples"
              + (f", first at t = {min(above6):.2f} s" if above6 else ", NEVER")
              + (f"; Q>={Q_TURB:.0f} in {100.0*len(above15)/len(qs):.0f}%"
                 if above15 else f"; Q>={Q_TURB:.0f} NEVER"))

    # sensitivity to which Omega is used
    print("\nsensitivity to Omega (central case, B_z = B_pol/sqrt2):")
    for key in ("Om_core", "Om_mean", "Om_mid", "Om_out"):
        qs = [lambda_over_dx(r, CASES["central"], key)[0] for r in rows]
        n6 = sum(1 for q in qs if q >= Q_LINEAR)
        print(f"  {key:>8}: max Q = {max(qs):.1f}, "
              f"Q>=6 in {100.0*n6/len(qs):.0f}% of samples")


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = load()
    t = [r["t"] for r in rows]
    q_lo = [lambda_over_dx(r, CASES["lower"])[0] for r in rows]
    q_ct = [lambda_over_dx(r, CASES["central"])[0] for r in rows]
    q_up = [lambda_over_dx(r, CASES["upper"])[0] for r in rows]
    b_pol = [lambda_over_dx(r, CASES["central"])[1] for r in rows]
    lam = [lambda_over_dx(r, CASES["central"])[3] for r in rows]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 7.0), sharex=True)

    ax1.axhspan(Q_TURB, 1e3, color="0.90", zorder=0)
    ax1.axhspan(Q_LINEAR, Q_TURB, color="0.96", zorder=0)
    ax1.fill_between(t, q_lo, q_up, color="#4878A8", alpha=0.30, lw=0,
                     label=r"$B_z$ between $B_{\rm pol}/\sqrt{3}$ and $B_{\rm pol}$")
    ax1.plot(t, q_ct, color="#1F4E79", lw=2.0,
             label=r"$B_z = B_{\rm pol}/\sqrt{2}$")
    ax1.axhline(Q_LINEAR, color="#B0561A", lw=1.4, ls="--")
    ax1.axhline(Q_TURB, color="#7A2E1E", lw=1.4, ls=":")
    ax1.text(79, Q_LINEAR, r"  $Q=6$: linear mode", va="center",
             fontsize=8, color="#B0561A")
    ax1.text(79, Q_TURB, r"  $Q=15$: converged turbulence", va="center",
             fontsize=8, color="#7A2E1E")
    ax1.set_yscale("log")
    ax1.set_ylim(0.02, 60)
    ax1.set_ylabel(r"$Q = \lambda_{\rm MRI}/\Delta x$")
    ax1.set_title(r"MRI quality factor from the volume-typical $B_z$, $256^3$",
                  fontsize=11)
    ax1.legend(loc="lower left", fontsize=8, framealpha=0.95)
    ax1.grid(alpha=0.25, which="both", lw=0.5)

    ax2.plot(t, b_pol, color="#1F4E79", lw=1.8,
             label=r"$B_{\rm pol,rms} = \sqrt{8\pi E_{\rm pol}/V}$")
    ax2.set_yscale("log")
    ax2.set_ylabel(r"$B_{\rm pol,rms}$  (G)", color="#1F4E79")
    ax2.tick_params(axis="y", labelcolor="#1F4E79")
    ax2.set_xlabel(r"$t$  (s)")
    ax2.grid(alpha=0.25, which="both", lw=0.5)

    ax2b = ax2.twinx()
    ax2b.plot(t, lam, color="#B0561A", lw=1.4, ls="--")
    ax2b.axhline(DX, color="0.35", lw=1.2, ls="-.")
    ax2b.text(79, DX, r"  $\Delta x$", va="center", fontsize=8, color="0.35")
    ax2b.set_yscale("log")
    ax2b.set_ylabel(r"$\lambda_{\rm MRI}$  (cm), dashed", color="#B0561A")
    ax2b.tick_params(axis="y", labelcolor="#B0561A")

    ax1.set_xlim(0, 78)
    fig.tight_layout()
    out = Path(__file__).parent / "mri_wavelength"
    fig.savefig(f"{out}.pdf")
    fig.savefig(f"{out}.png", dpi=150)
    print(f"\nwrote {out}.pdf and .png")


if __name__ == "__main__":
    main()
    plot()
