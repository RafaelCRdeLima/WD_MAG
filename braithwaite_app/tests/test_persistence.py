import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import core.persistence as persistence


@pytest.fixture
def isolated_store(tmp_path, monkeypatch):
    monkeypatch.setattr(persistence, "RESULTS_CSV", tmp_path / "results.csv")
    return tmp_path / "results.csv"


def test_star_cache_key_depends_on_full_definition():
    k1 = persistence.star_cache_key(1e9, 2.0, 64, "scfhash1", "gitscf1", "gitwd1")
    k2 = persistence.star_cache_key(1e9, 2.0, 64, "scfhash1", "gitscf1", "gitwd1")
    assert k1 == k2  # deterministic

    # different mu_e -> different key, even with same rho_c/resolution --
    # this is exactly the gap the review caught (cache key must include
    # composition, not just (rho_c, resolution))
    k3 = persistence.star_cache_key(1e9, 2.2, 64, "scfhash1", "gitscf1", "gitwd1")
    assert k1 != k3

    # different SCF params hash -> different key
    k4 = persistence.star_cache_key(1e9, 2.0, 64, "scfhash2", "gitscf1", "gitwd1")
    assert k1 != k4


def test_save_and_load_star_result_roundtrip(isolated_store):
    key = persistence.star_cache_key(1e9, 2.0, 64)
    window_result = {"valid": True, "window": (0.4, 1.128)}
    persistence.save_star_result(key, rho_c=988393849.5, mu_e=2.0, resolution=64,
                                  window_result=window_result, VE=1e-4)

    loaded = persistence.load_star_result(key)
    assert loaded is not None
    assert loaded["window_lo"] == pytest.approx(0.4)
    assert loaded["window_hi"] == pytest.approx(1.128)
    assert bool(loaded["window_valid"]) is True


def test_load_star_result_cache_miss_for_unknown_key(isolated_store):
    assert persistence.load_star_result("nonexistent") is None


def test_schema_version_mismatch_is_a_cache_miss(isolated_store, monkeypatch):
    key = persistence.star_cache_key(1e9, 2.0, 64)
    window_result = {"valid": True, "window": (0.4, 1.128)}
    persistence.save_star_result(key, rho_c=988393849.5, mu_e=2.0, resolution=64,
                                  window_result=window_result)
    assert persistence.load_star_result(key) is not None

    # simulate a schema bump -- old rows must stop being cache hits, same
    # convention as dashboard/store.py's SCHEMA_VERSION
    monkeypatch.setattr(persistence, "STAR_CACHE_SCHEMA_VERSION", 2)
    assert persistence.load_star_result(key) is None


def test_star_and_seed_rows_share_one_file_linked_by_cache_key(isolated_store):
    key = persistence.star_cache_key(1e9, 2.0, 64)
    window_result = {"valid": True, "window": (0.4, 1.128)}
    persistence.save_star_result(key, rho_c=988393849.5, mu_e=2.0, resolution=64,
                                  window_result=window_result)

    measurement = {"t_ttdyn": 0.6, "E_mag_over_W": 0.005, "E_tor_over_Emag": 0.324,
                   "divB_interior_max": 1e-10}
    persistence.save_seed_result(key, seed=42, resolution=64, e_mag_over_w_target=0.15,
                                  measurement=measurement, plotfile_path="/fake/path")

    all_results = persistence.load_all_results()
    assert set(all_results["row_type"]) == {"star", "seed"}
    assert isolated_store.exists()  # one file, both row types

    seeds = persistence.load_seeds_for_star(key)
    assert len(seeds) == 1
    assert seeds.iloc[0]["E_tor_over_Emag"] == pytest.approx(0.324)
