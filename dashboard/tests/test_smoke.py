"""Smoke test (dashboard acceptance criterion): with default parameters
and k0=0, Tab 1 (pages/1_equilibrium.py) must reproduce the Chandrasekhar
limit within 1%. Runs the real page via streamlit.testing — this is not a
visual check."""

from pathlib import Path

from streamlit.testing.v1 import AppTest

PAGE = Path(__file__).resolve().parent.parent / "pages" / "1_equilibrium.py"


def test_default_params_reproduce_chandrasekhar():
    at = AppTest.from_file(str(PAGE), default_timeout=180)
    at.run()

    assert not at.exception, f"page raised an exception: {at.exception}"

    tables = at.table
    assert len(tables) >= 1, "scalars table not found"
    scalars = dict(zip(tables[0].value["quantity"], tables[0].value["value"]))

    M_msun = float(scalars["M/M_sun"])
    rel_err = abs(M_msun - 1.44) / 1.44
    print(f"M = {M_msun:.4f} Msun, relative error to the Chandrasekhar limit: {rel_err:.3%}")
    assert rel_err < 0.01, f"M={M_msun:.4f} Msun, error {rel_err:.3%} above 1%"

    VE = float(scalars["VE"])
    assert VE < 1e-3, f"VE={VE:.3e} above the plan's V3"


if __name__ == "__main__":
    test_default_params_reproduce_chandrasekhar()
    print("OK")
