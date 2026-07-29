"""Tests for the three Passo 1.5 functions -- each checked against real
numbers/behavior from this session, not just internal self-consistency.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.star_builder import (
    GravityPatchMissing,
    find_measurement_window,
    half_shift_domain,
    parse_rho_c_log,
    verify_gravity_patch,
)

CASTRO_WD_BRAITHWAITE = Path(
    "/home/rafael/wd-magnetizada/castro/Exec/science/wd_braithwaite"
)
CASTRO_SOURCE_DIR = Path("/home/rafael/wd-magnetizada/castro")
FINTERIOR = Path(
    "/home/rafael/wd-magnetizada/castro/external/amrex/Tools/Plotfile/finterior.gnu.ex"
)
FLINE = Path(
    "/home/rafael/wd-magnetizada/castro/external/amrex/Tools/Plotfile/fline.gnu.ex"
)


# ---------------------------------------------------------------------
# (a) half_shift_domain
# ---------------------------------------------------------------------
def test_half_shift_domain_64_matches_session_values():
    prob_lo, prob_hi = half_shift_domain(64, 4.9e8)
    assert prob_lo == pytest.approx(-4.9765625e8, rel=1e-12)
    assert prob_hi == pytest.approx(4.8234375e8, rel=1e-12)


def test_half_shift_domain_128_matches_session_values():
    prob_lo, prob_hi = half_shift_domain(128, 4.9e8)
    assert prob_lo == pytest.approx(-4.93828125e8, rel=1e-12)
    assert prob_hi == pytest.approx(4.86171875e8, rel=1e-12)


def test_half_shift_domain_centers_a_cell_on_origin():
    for n_cell in (64, 128):
        prob_lo, prob_hi = half_shift_domain(n_cell, 4.9e8)
        dx = (prob_hi - prob_lo) / n_cell
        center_index = n_cell // 2
        cell_center = prob_lo + dx * (center_index + 0.5)
        assert cell_center == pytest.approx(0.0, abs=1e-6)


def test_half_shift_domain_rejects_odd_n_cell():
    with pytest.raises(ValueError, match="even"):
        half_shift_domain(65, 4.9e8)


# ---------------------------------------------------------------------
# (b) find_measurement_window
# ---------------------------------------------------------------------
def test_measurement_window_reproduces_seed42_crossings_from_real_log():
    log_path = CASTRO_WD_BRAITHWAITE / "run_halfshift_interp3d_test.log"
    if not log_path.exists():
        pytest.skip(f"{log_path} not present -- real regression data unavailable")

    t_dyn_s = 0.2758062098
    rho_c_ic = 988393849.5
    series = parse_rho_c_log(log_path, t_dyn_s)

    result = find_measurement_window(series, rho_c_ic, t_field_relax_ttdyn=0.4)

    # exact values computed this session via the ad hoc crossing script
    assert result["X_1pct"] == pytest.approx(0.573, abs=0.01)
    assert result["X_2pct"] == pytest.approx(1.128, abs=0.01)
    assert result["valid"] is True
    assert result["window"][0] == pytest.approx(0.4)
    assert result["window"][1] == pytest.approx(1.128, abs=0.01)


def test_measurement_window_invalid_when_star_drifts_before_field_relaxes():
    # synthetic: rho_c leaves the 2% band at t/t_dyn=0.2, well before the
    # default 0.4 field-relaxation timescale -- must report no valid window,
    # not silently report "valid until 0.2".
    rho_c_ic = 1.0e9
    series = [
        (0.0, rho_c_ic),
        (0.1, rho_c_ic * 0.995),
        (0.2, rho_c_ic * 0.975),  # already past 2% here
        (0.5, rho_c_ic * 0.90),
    ]
    result = find_measurement_window(series, rho_c_ic, t_field_relax_ttdyn=0.4)
    assert result["valid"] is False
    assert "reason" in result
    assert "sem janela" in result["reason"] or "deriva antes" in result["reason"]


def test_measurement_window_valid_when_star_stays_flat_past_field_relax():
    rho_c_ic = 1.0e9
    series = [(t, rho_c_ic * (1.0 - 0.001 * t)) for t in [0, 0.2, 0.4, 0.6, 0.8, 1.0]]
    result = find_measurement_window(series, rho_c_ic, t_field_relax_ttdyn=0.4)
    assert result["valid"] is True
    assert result["window"][0] == 0.4


# ---------------------------------------------------------------------
# (c) verify_gravity_patch -- real functional check, not mocked
# ---------------------------------------------------------------------
def test_verify_gravity_patch_passes_on_current_build(tmp_path):
    if not CASTRO_WD_BRAITHWAITE.joinpath("Castro3d.gnu.MPI.ex").exists():
        pytest.skip("wd_braithwaite executable not built in this checkout")

    # should not raise
    verify_gravity_patch(
        castro_source_dir=CASTRO_SOURCE_DIR,
        executable_path=CASTRO_WD_BRAITHWAITE / "Castro3d.gnu.MPI.ex",
        finterior_path=FINTERIOR,
        fline_path=FLINE,
        static_ic_inputs=CASTRO_WD_BRAITHWAITE / "inputs.halfshift_interp3d_test",
        scratch_dir=tmp_path,
    )


def test_verify_gravity_patch_static_check_catches_missing_patch(tmp_path):
    # simulate a reset checkout: source dir without the patch signature
    fake_source = tmp_path / "fake_castro"
    (fake_source / "Source" / "gravity").mkdir(parents=True)
    (fake_source / "Source" / "gravity" / "Gravity.cpp").write_text(
        "// no patch here\nvoid interpolate_monopole_grav() {}\n"
    )
    with pytest.raises(GravityPatchMissing, match="not found"):
        verify_gravity_patch(
            castro_source_dir=fake_source,
            executable_path=CASTRO_WD_BRAITHWAITE / "Castro3d.gnu.MPI.ex",
            finterior_path=FINTERIOR,
            fline_path=FLINE,
            static_ic_inputs=CASTRO_WD_BRAITHWAITE / "inputs.halfshift_interp3d_test",
            scratch_dir=tmp_path,
        )
