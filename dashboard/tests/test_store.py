"""Smoke test for store.py: save, find in cache, reload, index."""

import numpy as np
import sys
import os
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import store


def test_save_load_roundtrip():
    with tempfile.TemporaryDirectory() as tmp:
        params = {"rho_c": 1e9, "k0": 1e-13, "mu_e": 2.0}
        scalars = {"M_Msun": 1.35, "VE": 1.4e-4}
        fields = {"rho": np.ones((4, 4)), "r": np.linspace(0, 1, 4)}

        assert store.run_exists(params, runs_dir=tmp) is None
        h = store.save_run(params, scalars, fields, runs_dir=tmp)
        assert store.run_exists(params, runs_dir=tmp) == h

        loaded = store.load_run(h, runs_dir=tmp)
        assert loaded["params"] == params
        assert abs(loaded["scalars"]["M_Msun"] - 1.35) < 1e-12
        assert np.allclose(loaded["fields"]["rho"], fields["rho"])
        assert loaded["manifest"]["hash"] == h

        idx = store.load_index(runs_dir=tmp)
        assert len(idx) == 1
        assert idx.iloc[0]["hash"] == h

        store.mark_reference(h, True, runs_dir=tmp)
        idx2 = store.load_index(runs_dir=tmp)
        assert bool(idx2.iloc[0]["reference"]) is True


if __name__ == "__main__":
    test_save_load_roundtrip()
    print("OK")
