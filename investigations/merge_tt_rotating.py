"""Join the tt_rot_NNN.csv and say what the four gates left standing.

The field-only sweep (DIARIO 19) found the twisted torus constructible but
never super-Chandrasekhar: E_tor/|W| capped near 0.0101, the field worth
+0.017 Msun, the frontier saturating at M_Ch. This asks the same question of a
rotating star, where the mass can come from somewhere else.

Two things to read out, and they are different questions:

  1. Does any point pass all four gates at once -- super-Chandrasekhar, holds
     its field, comparable energies, below the secular bar mode and not
     shedding? If yes, the exclusion of DIARIO 19 was about the SOURCE of the
     mass, not about twisted tori.
  2. Does rotation move the E_tor/|W| ceiling? Centrifugal support lowers the
     gas pressure that holds the field, so the ceiling could fall as the mass
     rises. The omega_frac = 0 rows are the control: they must reproduce the
     field-only sweep.

Run:  scf/.venv/bin/python3 investigations/merge_tt_rotating.py
Writes tt_rotating.csv.
"""

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "tt_rotating.csv"
M_CH = 1.44


def main():
    parts = sorted(HERE.glob("tt_rot_[0-9][0-9][0-9].csv"))
    if not parts:
        raise SystemExit("no tt_rot_NNN.csv here -- did the sweep run?")

    rows, missing = [], []
    for p in parts:
        head = p.read_text().splitlines()
        if any("DID NOT CONVERGE" in ln for ln in head):
            missing.append(p.name)
            continue
        meta = {}
        for tok in head[0].lstrip("# ").split():
            if "=" in tok:
                k, v = tok.split("=", 1)
                meta[k] = v
        with p.open() as f:
            for row in csv.DictReader(ln for ln in f if not ln.startswith("#")):
                row["rho_c"] = meta.get("rho_c", "")
                row["omega_frac"] = meta.get("omega_frac", "")
                row["k_frac"] = meta.get("k_frac", "")
                rows.append(row)
    if missing:
        print(f"{len(missing)} point(s) did not converge: {', '.join(missing)}")
    if not rows:
        raise SystemExit("every point failed to converge")

    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    held = [r for r in rows if float(r["beta_min_tor_only"]) >= 1.0]
    passed = [r for r in rows if r["verdict"] == "ALL FOUR GATES"]
    print(f"\n{len(rows)} rows from {len(parts)} points; "
          f"{len(held)} hold the field; {len(passed)} pass all four gates")

    if passed:
        best = max(passed, key=lambda r: float(r["M_msun"]))
        print(f"\nheaviest configuration passing every gate:")
        print(f"  M          = {float(best['M_msun']):.4f} Msun "
              f"({float(best['M_msun']) / M_CH:.2f} x Chandrasekhar)")
        print(f"  rho_c      = {float(best['rho_c']):.2e} g/cm^3")
        print(f"  omega_frac = {float(best['omega_frac']):.2f}, "
              f"T/|W| = {float(best['T_over_W']):.4f}, "
              f"shedding = {float(best['shedding']):.3f}")
        print(f"  K_frac     = {float(best['k_frac']):.2f}, "
              f"E_tor/|W| = {float(best['Etor_over_W']):.5f}")
        print(f"  m = {float(best['m']):.2f}, zeta = {float(best['zeta']):.2f}")
        print(f"  E_tor/E_pol = {float(best['Etor_over_Epol']):.2f}, "
              f"B_t/B_p = {float(best['Bt_over_Bp']):.2f}")
        print(f"  beta in the Castro ambient = "
              f"{float(best['beta_ambient']):.2e}")

    # does rotation move the ceiling? group the held rows by omega_frac
    print("\nE_tor/|W| ceiling against rotation "
          "(omega_frac = 0 is the field-only control):")
    ceil = {}
    mmax = {}
    for r in held:
        k = float(r["omega_frac"])
        e, m = float(r["Etor_over_W"]), float(r["M_msun"])
        ceil[k] = max(ceil.get(k, 0.0), e)
        mmax[k] = max(mmax.get(k, 0.0), m)
    print("  omega_frac   max E_tor/|W| held   max M held")
    for k in sorted(ceil):
        print(f"    {k:5.2f}         {ceil[k]:.5f}           "
              f"{mmax[k]:.4f} Msun")

    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
