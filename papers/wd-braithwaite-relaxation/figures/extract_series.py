"""Cache the time series behind the figures of the mixed-field paper.

Run:
    scf/.venv/bin/python3 papers/wd-toroidal-poloidal/figures/extract_series.py

Reads, through the app's own modules (R1 -- no diagnostic is recomputed
here, and in particular Div_B and the energy volume integrals go through
the interior-masking path, never the raw derived value):

  1. The field-free background star's rho_c(t), from the log the
     measurement window was derived from in the first place
     (run_halfshift_interp3d_test.log, docs/teoria.md Sec 6.9).
  2. E_tor/E_mag(t) and E_mag/|W|(t) for the ten seeds of batch C, from
     the plotfiles those runs left behind.
  3. rho_c(t) for the three extended-baseline runs, which is what the
     rank-ordering check of Sec 6.10 was made against.

Writes series.npz. The plotfiles live under the (gitignored) castro
tree, so this cannot be re-run from a fresh clone -- which is why its
output is committed.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "braithwaite_app"))

from core.extraction import extract_full_series      # noqa: E402
from core.star_builder import parse_rho_c_log        # noqa: E402

RUN_DIR = REPO / "castro" / "Exec" / "science" / "wd_braithwaite"
ARCHIVE = REPO / "braithwaite_app" / "data" / "raw_logs_archive"

# From the batch-C run logs ("t_dyn = ... s" and "|W| = ... erg").
T_DYN_S = 0.2758062098
W_ABS_ERG = 2.25755333e51
RHO_C_IC = 988393849.5

# docs/teoria.md Sec 6.9: the window this star yields, and the field's
# own relaxation time used as the lower bound.
T_FIELD_RELAX = 0.4
X_1PCT = 0.573
X_2PCT = 1.128

SEEDS_C = list(range(11, 21))
EXT_RUNS = ["seed42ext", "seed55ext", "seed123ext"]

OUT = HERE / "series.npz"


def main():
    store = {}

    star_log = RUN_DIR / "run_halfshift_interp3d_test.log"
    t, rho_c = zip(*parse_rho_c_log(star_log, T_DYN_S))
    store["star_t"] = np.asarray(t)
    store["star_rho_c"] = np.asarray(rho_c)
    print(f"background star: {len(t)} samples, "
          f"t/t_dyn up to {max(t):.2f}")

    for seed in SEEDS_C:
        run_id = f"seed{seed}_caeade"
        series = extract_full_series(RUN_DIR, run_id, W_ABS_ERG, T_DYN_S)
        if not series:
            print(f"  seed {seed}: no plotfiles, skipped")
            continue
        store[f"seed{seed}_t"] = np.array([d["t_ttdyn"] for d in series])
        store[f"seed{seed}_ratio"] = np.array(
            [d["E_tor_over_Emag"] for d in series])
        store[f"seed{seed}_emag"] = np.array(
            [d["E_mag_over_W"] for d in series])
        print(f"  seed {seed}: {len(series)} plotfiles, "
              f"E_tor/E_mag {series[0]['E_tor_over_Emag']:.3f} -> "
              f"{series[-1]['E_tor_over_Emag']:.3f}, "
              f"E_mag/|W| {series[0]['E_mag_over_W']:.4f} -> "
              f"{series[-1]['E_mag_over_W']:.4f}")

    for run in EXT_RUNS:
        log = ARCHIVE / f"run_{run}.log"
        if not log.exists():
            print(f"  {run}: log missing, skipped")
            continue
        t, rho_c = zip(*parse_rho_c_log(log, T_DYN_S))
        store[f"{run}_t"] = np.asarray(t)
        store[f"{run}_rho_c"] = np.asarray(rho_c)
        print(f"  {run}: {len(t)} samples, t/t_dyn up to {max(t):.2f}")

    store["meta"] = json.dumps({
        "t_dyn_s": T_DYN_S, "W_abs_erg": W_ABS_ERG, "rho_c_ic": RHO_C_IC,
        "t_field_relax": T_FIELD_RELAX, "X_1pct": X_1PCT, "X_2pct": X_2PCT,
        "seeds_C": SEEDS_C, "ext_runs": EXT_RUNS,
    })
    np.savez_compressed(OUT, **store)
    print("wrote", OUT.name, f"({OUT.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
