"""Join the tt_sweep_NNN.csv the array writes, and say what the grid found.

The question the sweep asks is whether ANY point holds beta_min >= 1 with the
toroidal confined to closed field lines, and if so what B_t/B_p it buys. The
answer is a region in (M, beta_min), not a number, so this prints the frontier
rather than a verdict: for each central density, the largest mass whose
confined toroidal the star can still hold.

Run:  scf/.venv/bin/python3 investigations/merge_tt_sweep.py
Writes tt_sweep.csv.
"""

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "tt_sweep.csv"


def collect(parts):
    """Read a set of shards into one list of rows, carrying the header meta."""
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
                row["k_frac"] = meta.get("k_frac", "")
                rows.append(row)
    return rows, missing


def write(rows, out):
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main():
    # The off-grid high-density points get their own artefact: they answer a
    # different question (does the frontier cross 1.44 above the grid?) and
    # mixing them into tt_sweep.csv would hide that they are not a grid.
    hi = sorted(HERE.glob("tt_sweep_hi_*.csv"))
    if hi:
        hrows, _ = collect(hi)
        if hrows:
            write(hrows, HERE / "tt_sweep_hidens.csv")
            print(f"wrote tt_sweep_hidens.csv ({len(hrows)} rows from "
                  f"{len(hi)} off-grid points)")

    parts = sorted(HERE.glob("tt_sweep_[0-9][0-9][0-9].csv"))
    if not parts:
        raise SystemExit("no tt_sweep_NNN.csv here -- did the sweep run?")

    rows, missing = collect(parts)
    if missing:
        print(f"{len(missing)} point(s) did not converge: {', '.join(missing)}")
    if not rows:
        raise SystemExit("every point failed to converge")

    write(rows, OUT)

    held = [r for r in rows if float(r["beta_min_tor_only"]) >= 1.0]
    print(f"\n{len(rows)} rows from {len(parts)} points; "
          f"{len(held)} hold the confined toroidal at beta_min >= 1")

    if not held:
        print("\nNo point in the grid holds it. The obstruction DIARIO 15")
        print("measured at M = 2 Msun is not special to that mass: confining")
        print("the toroidal to the closed-line region breaks beta wherever")
        print("the field carries enough of the virial to matter.")
        print("\nthe three closest, by beta_min:")
        for r in sorted(rows, key=lambda r: -float(r["beta_min_tor_only"]))[:3]:
            print(f"  rho_c={float(r['rho_c']):.2e} M={float(r['M_msun']):.4f} "
                  f"E_tor/|W|={float(r['Etor_over_W']):.4f} "
                  f"m={r['m']} zeta={r['zeta']} "
                  f"beta_min={float(r['beta_min_tor_only']):.3e}")
    else:
        # Braithwaite & Spruit ask for comparable ENERGIES, and both pure
        # branches are unstable. A row with E_t/E_p = 0.16 is poloidal
        # dominated, which is the other unstable branch, not a success. Take
        # the band 0.1 < E_t/E_p < 10 as comparable and report the rest apart.
        def band(r):
            e = float(r["Etor_over_Epol"])
            return 0.1 < e < 10.0
        tt = [r for r in held if band(r)]
        pol = [r for r in held if float(r["Etor_over_Epol"]) <= 0.1]
        tor = [r for r in held if float(r["Etor_over_Epol"]) >= 10.0]
        print(f"of those, {len(tt)} sit in the comparable band "
              f"(0.1 < E_t/E_p < 10); {len(pol)} are poloidal dominated and "
              f"{len(tor)} toroidal dominated -- both unstable branches:")
        for r in sorted(tt, key=lambda r: -float(r["M_msun"]))[:10]:
            print(f"  rho_c={float(r['rho_c']):.2e} M={float(r['M_msun']):.4f} "
                  f"m={r['m']} zeta={r['zeta']} "
                  f"E_t/E_p={float(r['Etor_over_Epol']):.2f} "
                  f"beta_amb={float(r['beta_ambient']):.2e}")
        if tt:
            mmax = max(float(r["M_msun"]) for r in tt)
            print(f"\n  heaviest star holding a COMPARABLE twisted torus: "
                  f"{mmax:.4f} Msun")
            print(f"  Chandrasekhar for mu_e = 2 is 1.44 Msun -- this is "
                  f"{'ABOVE' if mmax > 1.44 else 'BELOW'} it.")

    print("\nfrontier -- largest M holding the confined toroidal, per rho_c:")
    # held, not band-restricted: this frontier is about beta alone.
    by_rho = {}
    for r in held:
        k = float(r["rho_c"])
        if k not in by_rho or float(r["M_msun"]) > float(by_rho[k]["M_msun"]):
            by_rho[k] = r
    for k in sorted(by_rho):
        r = by_rho[k]
        print(f"  rho_c = {k:.2e}  ->  M = {float(r['M_msun']):.4f} Msun "
              f"(E_tor/|W| = {float(r['Etor_over_W']):.4f})")
    if not by_rho:
        print("  (empty: no point holds it)")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
