"""Verifies the app builds the SAME star as the science runs, not a
similar one -- the real SCF solver (scf.hachisu_scf), not a reused
stale file. This is the fix for the gap the review caught: rho_c was
cosmetic until this module called the actual solver.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core import scf_store

EXISTING_MODEL_DAT = Path(
    "/home/rafael/wd-magnetizada/castro/Exec/science/wd_braithwaite/model.dat"
)


def test_scf_solve_is_fast_enough_to_run_synchronously():
    import time
    t0 = time.time()
    scf_store.converge_field_free_star(1.0e9, 2.0)
    elapsed = time.time() - t0
    # measured 16-25ms this session; generous margin, still far below any
    # threshold that would require backgrounding
    assert elapsed < 2.0, f"SCF solve took {elapsed:.2f}s -- reconsider running inline"


def test_build_model_dat_reproduces_the_real_science_star(tmp_path):
    if not EXISTING_MODEL_DAT.exists():
        pytest.skip(f"{EXISTING_MODEL_DAT} not present -- nothing to compare against")

    out_path = tmp_path / "model.dat"
    manifest = scf_store.build_model_dat(rho_c=1.0e9, mu_e=2.0, out_path=out_path)

    assert manifest["R_star_cm"] == pytest.approx(244853796.37, rel=1e-6)

    existing_lines = EXISTING_MODEL_DAT.read_text().splitlines()
    built_lines = out_path.read_text().splitlines()
    assert len(existing_lines) == len(built_lines)

    # density column (index 1) must match on every line -- the only
    # physically meaningful column (temperature/pressure are placeholders,
    # see castro_model_writer.py's docstring)
    mismatches = 0
    for existing, built in zip(existing_lines[1:], built_lines[1:]):
        existing_rho = float(existing.split()[1])
        built_rho = float(built.split()[1])
        if existing_rho != pytest.approx(built_rho, rel=1e-9):
            mismatches += 1
    assert mismatches == 0, f"{mismatches} density values differ beyond float noise"


def test_build_model_dat_raises_on_non_convergence(tmp_path, monkeypatch):
    def fake_hachisu_scf(*args, **kwargs):
        return {"converged": False, "iterations": 5, "history": [1e-2], "rho": None}

    monkeypatch.setattr(scf_store.scf_mod, "hachisu_scf", fake_hachisu_scf)
    with pytest.raises(RuntimeError, match="did not converge"):
        scf_store.build_model_dat(1.0e9, 2.0, tmp_path / "model.dat")


def test_neutronization_gate_still_works_after_scf_wiring():
    # regression guard: adding the SCF bridge must not disturb the
    # existing gate
    assert scf_store.neutronization_check(1.0e9, 2.0)["ok"] is True
    assert scf_store.neutronization_check(1.0e12, 2.0)["ok"] is False
