"""Compare rho_c(t) between Castro and FLASH on the same star.

Run:  scf/.venv/bin/python3 flash_crosscheck/compare_codes.py

Both sides are field-free, self-gravitating, 64^3 on the same half-shift
domain, with the same global velocity damping and exterior sponge, and the
same 1D white-dwarf structure. What differs is the code -- which is the
point.

Castro reference: the field-free calibration run behind the paper's Fig. 1,
cached in papers/wd-toroidal-poloidal/figures/series.npz by
extract_series.py. Its own damping window was [0, 20 t_dyn], reproduced on
the FLASH side.

Not identical, and the differences are listed here rather than buried:
  * EOS. Castro evolves ztwd; FLASH evolves Helmholtz on a ztwd structure
    at T = 1e7 K, where the terms Helmholtz adds are ~1e-3 of the pressure
    (measured in make_wd_model.py, not assumed).
  * The damping is operator-split in FLASH (no in-hydro source hook in its
    driver) and uses exp(-r dt) rather than Castro's explicit 1 - r dt.
  * The ambient is 1 g/cm^3 in FLASH against Castro's much lower floor,
    forced by Helmholtz refusing to invert a harder vacuum.
  * t_dyn differs slightly (0.28797 s against 0.27581 s) because the 1D
    integration's surface convention is not the Castro problem's; masses
    agree to 0.026%. Everything below is in units of each run's own t_dyn.
"""

import glob
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
# Castro reference, MATCHED IN GEOMETRY to the FLASH run being compared.
# run_interp3d_test is the symmetric / vertex-centered field-free run
# (t=0 rho_c = 952180729.1, i.e. -4.78% of the 1e9 target, Sec 6.6's
# baseline case), with the same damping window [0, 20 t_dyn] and the same
# star. The half-shift Castro run lives in series.npz and is NOT the right
# reference for a symmetric FLASH run.
CASTRO_LOG = (REPO / "castro" / "Exec" / "science" / "wd_braithwaite"
              / "run_interp3d_test.log")
T_DYN_CASTRO = 0.2758062098

T_DYN_FLASH = 0.2879670
WINDOW = (0.4, 1.128)     # the Castro validity window, docs/teoria.md Sec 6.9

# Only the symmetric FLASH run is compared quantitatively. The half-shift
# FLASH run is excluded on purpose: the half-shift geometry puts a cell
# centre at r=0, and FLASH develops a spurious flow in exactly that cell --
# the dt limiter sits at the origin from step 1 and dt settles 100x below
# the symmetric run's. Since rho_c IS the central cell, that run cannot
# measure the quantity being compared. Castro needed a core patch
# (g(r=0)=0, Sec 6.7) to survive the same geometry.
RUNS = [
    ("run_sponge", "FLASH (damping + sponge)", "#2a78d6", "-"),
]

COL_IN = 88.0 / 25.4


def flash_series(subdir):
    import yt
    yt.set_log_level(50)
    files = sorted(glob.glob(str(HERE / subdir / "wd_hdf5_plt_cnt_*")))
    if not files:
        return None, None
    t, rc = [], []
    for f in files:
        ds = yt.load(f)
        t.append(float(ds.current_time))
        rc.append(float(ds.all_data().max(("flash", "dens"))))
    return np.array(t) / T_DYN_FLASH, np.array(rc)


def castro_series():
    """rho_c(t) from Castro's own log, via the app's parser (R1)."""
    import sys
    sys.path.insert(0, str(REPO / "braithwaite_app"))
    from core.star_builder import parse_rho_c_log
    rows = parse_rho_c_log(CASTRO_LOG, T_DYN_CASTRO)
    t = np.array([r[0] for r in rows])
    rc = np.array([r[1] for r in rows])
    return t, 100.0 * (rc / rc[0] - 1.0)


def main():
    ct, cdev = castro_series()
    print(f"Castro reference: {CASTRO_LOG.name}, {len(ct)} samples "
          f"to t/t_dyn = {ct[-1]:.2f}\n")

    series = {}
    for sub, label, _, _ in RUNS:
        tt, rc = flash_series(sub)
        if tt is None:
            print(f"  ({sub}: no plotfiles, skipped)")
            continue
        series[sub] = (tt, 100.0 * (rc / rc[0] - 1.0), label)

    print(f"{'t/t_dyn':>8}  {'Castro':>9}", end="")
    for sub in series:
        print(f"  {series[sub][2][:22]:>22}", end="")
    print()
    grid = np.arange(0.1, 1.75, 0.1)
    for tg in grid:
        row = f"{tg:8.2f}  {np.interp(tg, ct, cdev):+8.2f}%"
        for sub in series:
            tt, dv, _ = series[sub]
            row += f"  {np.interp(tg, tt, dv, right=np.nan):+21.2f}%" \
                if tg <= tt[-1] else f"  {'--':>22}"
        print(row)

    hi = WINDOW[1]
    print(f"\nAt the window's upper bound, t/t_dyn = {hi}:")
    c_at = np.interp(hi, ct, cdev)
    print(f"  Castro                 {c_at:+.2f}%")
    for sub in series:
        tt, dv, label = series[sub]
        if tt[-1] >= hi:
            f_at = np.interp(hi, tt, dv)
            print(f"  {label:22s} {f_at:+.2f}%   "
                  f"(x{abs(f_at / c_at):.1f} the Castro drift)")
        else:
            print(f"  {label:22s} run ended at {tt[-1]:.2f} t_dyn")

    # ---- figure ----
    fig, ax = plt.subplots(figsize=(COL_IN, COL_IN * 0.72))
    ax.grid(True, axis="y", color="#e1e0d9", linewidth=0.4, zorder=0)
    ax.set_axisbelow(True)
    ax.axvspan(*WINDOW, color="#f2f1ec", zorder=0)
    ax.axhline(0.0, color="#898781", linewidth=0.5, zorder=1)
    for pct in (1.0, 2.0):
        ax.axhspan(-pct, pct, color="#dce9fc", alpha=0.30, zorder=0)

    ax.plot(ct, cdev, color="#0b0b0b", linewidth=1.5, label="Castro", zorder=4)
    for sub, _, color, ls in RUNS:
        if sub not in series:
            continue
        tt, dv, label = series[sub]
        ax.plot(tt, dv, color=color, linewidth=1.3, linestyle=ls,
                label=label, zorder=3)
        ax.plot([tt[-1]], [dv[-1]], marker="x", markersize=5,
                color=color, markeredgewidth=1.2, zorder=5)

    ax.text(np.mean(WINDOW), 1.6, "validity window", fontsize=6.5,
            color="#898781", ha="center", va="bottom")
    ax.set_xlim(0, 1.75)
    ax.set_xlabel(r"$t/t_{\rm dyn}$")
    ax.set_ylabel(r"$\rho_c$ deviation (%)")
    ax.legend(loc="lower left", frameon=False, fontsize=6.5, handlelength=2.4)
    for side in ("top", "right"):
        ax.spines[side].set_color("#898781")
    fig.savefig(HERE / "compare_codes.pdf", bbox_inches="tight",
                pad_inches=0.02)
    print("\nwrote compare_codes.pdf")

    with (HERE / "compare_codes.csv").open("w") as fh:
        fh.write("t_ttdyn,castro_dev_pct," +
                 ",".join(f"{s}_dev_pct" for s in series) + "\n")
        for tg in np.arange(0.05, 1.75, 0.05):
            vals = [f"{np.interp(tg, series[s][0], series[s][1]):.4f}"
                    if tg <= series[s][0][-1] else "" for s in series]
            fh.write(f"{tg:.2f},{np.interp(tg, ct, cdev):.4f}," +
                     ",".join(vals) + "\n")
    print("wrote compare_codes.csv")


if __name__ == "__main__":
    main()
