"""Tab 2 (sweep) end-to-end acceptance test for the item-4 grid extension
(rotation modes, self-consistent toroidal field, out-of-validity marking).
Runs the real page via streamlit.testing, small grid for speed -- not a
visual check."""

from pathlib import Path

from streamlit.testing.v1 import AppTest

PAGE = Path(__file__).resolve().parent.parent / "pages" / "2_sweep.py"


def _shrink_mesh(at):
    for sl in at.sidebar.select_slider:
        if sl.label == "Nr":
            sl.set_value(65)
        elif sl.label == "Ntheta":
            sl.set_value(33)
    for sl in at.sidebar.slider:
        if sl.label == "parallel processes":
            sl.set_value(2)


def test_toroidal_sc_sweep_marks_out_of_validity_points():
    """rho_c axis straddling the (mu_e=2) neutronization threshold
    (~1.94e10 g/cm^3): must produce BOTH rho_c_valid True and False rows,
    and the priority plot ('M vs rho_c, one curve per K') must render
    without exception."""
    at = AppTest.from_file(str(PAGE), default_timeout=300)
    at.run()
    assert not at.exception

    for r in at.sidebar.radio:
        if r.label == "field":
            r.set_value("toroidal (self-consistent)")
    at.run()
    assert not at.exception

    _shrink_mesh(at)
    for ni in at.sidebar.number_input:
        if ni.label == "rho_c min":
            ni.set_value(1e8)
        elif ni.label == "rho_c max":
            ni.set_value(1e12)
    for sl in at.sidebar.slider:
        if sl.label == "N (rho_c)":
            sl.set_value(3)
        elif sl.label == "N (K>0)":
            sl.set_value(1)
    at.run()
    assert not at.exception

    for b in at.button:
        if b.label == "run sweep":
            b.click().run()
    assert not at.exception, f"sweep run raised: {at.exception}"

    rows = at.session_state["sweep_rows"] if "sweep_rows" in at.session_state else []
    assert len(rows) >= 4, f"expected at least 4 converged points, got {len(rows)}"

    valid_flags = {bool(r["rho_c_valid"]) for r in rows}
    assert valid_flags == {True, False}, (
        f"expected both valid and out-of-validity rho_c points in this grid, got {valid_flags}"
    )

    # new item-4 scalars must be present on every row
    for r in rows:
        for key in ("B_tor,max (G)", "T/|W|", "Bt/Bp (energy)", "Bt/Bp (amplitude)",
                    "equatorial mass-loss ratio"):
            assert key in r, f"missing scalar {key!r} in sweep row"

    headers = [h.value for h in at.subheader]
    assert any("priority plot" in h for h in headers), "priority plot section not rendered"


def test_rigid_rotation_sweep_records_non_convergence_as_failure():
    """High Omega_c at fixed rho_c is expected to stop converging (known,
    documented failure mode -- see scf/terms/rotation.py and Tab 1) --
    those points must land in sweep_fail_rows, not be silently dropped or
    reported as converged."""
    at = AppTest.from_file(str(PAGE), default_timeout=300)
    at.run()

    for r in at.sidebar.radio:
        if r.label == "rotation":
            r.set_value("rigid")
    at.run()
    assert not at.exception

    _shrink_mesh(at)
    for ni in at.sidebar.number_input:
        if ni.label == "Omega_c max":
            ni.set_value(40.0)
    for sl in at.sidebar.slider:
        if sl.label == "N (rho_c)":
            sl.set_value(3)
        elif sl.label == "N (Omega_c)":
            sl.set_value(5)
    at.run()
    assert not at.exception

    for b in at.button:
        if b.label == "run sweep":
            b.click().run()
    assert not at.exception, f"sweep run raised: {at.exception}"

    fails = at.session_state["sweep_fail_rows"] if "sweep_fail_rows" in at.session_state else []
    rows = at.session_state["sweep_rows"] if "sweep_rows" in at.session_state else []
    assert len(fails) > 0, "expected at least one non-converging point at high Omega_c"
    assert len(rows) > 0, "expected at least one converging point (e.g. Omega_c=0)"


if __name__ == "__main__":
    test_toroidal_sc_sweep_marks_out_of_validity_points()
    test_rigid_rotation_sweep_records_non_convergence_as_failure()
    print("OK")
