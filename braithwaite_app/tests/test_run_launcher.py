import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.run_launcher import RunSpec, render_inputs, read_progress


def test_render_inputs_uses_half_shift_geometry_matching_star_builder():
    spec = RunSpec(run_id="t", n_cell=64, rng_seed=42, e_mag_over_w=0.15, stop_time_s=0.4)
    text = render_inputs(spec)
    assert "-4.9765625" in text  # prob_lo for n_cell=64, matches core.star_builder's tested value
    assert "4.8234375" in text   # prob_hi


def test_read_progress_parses_w_abs_from_real_log():
    log_path = Path(
        "/home/rafael/wd-magnetizada/castro/Exec/science/wd_braithwaite/run_seed42_repro.log"
    )
    if not log_path.exists():
        pytest.skip("reproduction run log not present")
    p = read_progress(log_path, "seed42_repro")
    assert p.w_abs_erg == pytest.approx(2.25755333e51, rel=1e-6)
    assert p.t_dyn_s == pytest.approx(0.2758062098, rel=1e-6)
