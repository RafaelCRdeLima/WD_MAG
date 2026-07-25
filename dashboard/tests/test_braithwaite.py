"""Tab 5 (Braithwaite) acceptance test: honest skeleton -- renders with
no runs saved, renders with a background star available, seed config
save/list round-trips, and every Castro-dependent control is disabled."""

from pathlib import Path

from streamlit.testing.v1 import AppTest

PAGE = Path(__file__).resolve().parent.parent / "pages" / "5_braithwaite.py"
EQUILIBRIUM_PAGE = Path(__file__).resolve().parent.parent / "pages" / "1_equilibrium.py"

DISABLED_LABELS = {
    "generate seed (requires Castro)",
    "launch Castro evolution run",
    "refresh run status",
}


def test_renders_with_no_runs():
    at = AppTest.from_file(str(PAGE), default_timeout=60)
    at.run()
    assert not at.exception, f"page raised an exception: {at.exception}"


def test_background_star_selector_and_seed_config_roundtrip():
    # save a field-free equilibrium first (Tab 1, default params: no field, no rotation)
    at1 = AppTest.from_file(str(EQUILIBRIUM_PAGE), default_timeout=180)
    at1.run()
    assert not at1.exception
    for b in at1.button:
        if b.label == "save this run":
            b.click().run()
    assert not at1.exception

    at = AppTest.from_file(str(PAGE), default_timeout=60)
    at.run()
    assert not at.exception, f"page raised an exception: {at.exception}"

    metrics = {m.label: m.value for m in at.metric}
    assert "M/M_sun" in metrics and metrics["M/M_sun"] != "—", \
        "background star selector did not find the saved field-free run"

    for b in at.button:
        if b.label == "save this seed configuration":
            b.click().run()
    assert not at.exception
    assert any("Configuration saved" in s.value for s in at.success)
    assert len(at.dataframe) >= 1, "saved seed config table did not render"


def test_castro_dependent_controls_are_disabled():
    at = AppTest.from_file(str(PAGE), default_timeout=60)
    at.run()
    assert not at.exception
    found = {b.label: b.disabled for b in at.button if b.label in DISABLED_LABELS}
    assert found == {label: True for label in DISABLED_LABELS}, (
        f"expected every Castro-dependent control disabled, got {found}"
    )


if __name__ == "__main__":
    test_renders_with_no_runs()
    test_background_star_selector_and_seed_config_roundtrip()
    test_castro_dependent_controls_are_disabled()
    print("OK")
