"""castro_model_writer -- validates field-free/non-rotating + spherical
symmetry, then writes a Castro model_parser-format 1D profile the
wd_braithwaite problem reads via read_model_file()."""

import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import castro_model_writer as cmw


def _spherical_star(n_r=50, n_theta=8, rho_c=1e9, R_star=1e8):
    r = np.linspace(0, 2 * R_star, n_r)
    theta = np.linspace(0, np.pi, n_theta)
    rho_1d = np.where(r < R_star, rho_c * (1.0 - (r / R_star) ** 2), 0.0)
    rho = np.tile(rho_1d[:, None], (1, n_theta))
    return r, theta, rho, rho_c


def test_rejects_field_present():
    r, theta, rho, rho_c = _spherical_star()
    with pytest.raises(ValueError, match="field-free"):
        cmw.check_field_free_non_rotating({"rho_c": rho_c, "k0": 1.5})


def test_rejects_rotation_present():
    r, theta, rho, rho_c = _spherical_star()
    with pytest.raises(ValueError, match="field-free"):
        cmw.check_field_free_non_rotating({"rho_c": rho_c, "Omega_c": 0.3})


def test_accepts_field_free_non_rotating():
    cmw.check_field_free_non_rotating({"rho_c": 1e9, "k0": 0.0, "K_tor": 0.0, "Omega_c": 0.0})


def test_rejects_non_spherical():
    r, theta, rho, rho_c = _spherical_star()
    rho_broken = rho.copy()
    rho_broken[10, 3] *= 1.5  # inject a theta-dependent bump
    with pytest.raises(ValueError, match="spherically symmetric"):
        cmw.check_spherical_symmetry(rho_broken, rho_c)


def test_write_model_file_roundtrip(tmp_path):
    r, theta, rho, rho_c = _spherical_star()
    out_path = tmp_path / "model.dat"
    manifest = cmw.write_model_file(
        r, theta, rho, {"rho_c": rho_c, "k0": 0.0, "Omega_c": 0.0}, out_path,
        run_hash="testhash", git_commit="deadbeef",
    )

    assert out_path.exists()
    lines = out_path.read_text().splitlines()
    assert lines[0] == "# r density temperature pressure X"
    assert len(lines) == 1 + len(r)

    # every row is parseable and density is non-negative and floored
    for line in lines[1:]:
        r_i, rho_i, T_i, p_i, X_i = (float(x) for x in line.split())
        assert rho_i >= 0
        assert T_i == cmw.PLACEHOLDER_TEMPERATURE_K
        assert X_i == 1.0

    # density column matches the theta-mean of the input profile (floored)
    written_rho = np.array([float(l.split()[1]) for l in lines[1:]])
    expected = np.maximum(rho.mean(axis=1), rho_c * 1e-10)
    np.testing.assert_allclose(written_rho, expected, rtol=1e-9)

    manifest_path = tmp_path / "model.dat.manifest.json"
    assert manifest_path.exists()
    with open(manifest_path) as f:
        saved_manifest = json.load(f)
    assert saved_manifest == manifest
    assert manifest["source_run_hash"] == "testhash"
    assert manifest["git_commit_scf"] == "deadbeef"
    assert manifest["rho_c_gcm3"] == rho_c
    # R_star should land near the analytic star radius used to build the fixture
    assert abs(manifest["R_star_cm"] - 1e8) / 1e8 < 0.1


if __name__ == "__main__":
    import subprocess
    raise SystemExit(subprocess.call([sys.executable, "-m", "pytest", __file__, "-v"]))
