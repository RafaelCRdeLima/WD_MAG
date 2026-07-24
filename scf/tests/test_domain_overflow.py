"""domain_overflow_check() -- item 11/12 gap (docs/teoria.md Sec 8):
flags R_eq/R_pol landing suspiciously close to the domain edge, checked
independently in each direction (oblate stars overflow at the equator,
prolate stars overflow at the pole)."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import diagnostics as diag


def test_no_overflow_when_comfortably_inside():
    r = diag.domain_overflow_check(R_eq=0.5, R_pol=0.4, r_max=1.0, tol=0.1)
    assert not r["overflow"]
    assert not r["overflow_eq"]
    assert not r["overflow_pol"]


def test_flags_equatorial_overflow_oblate_star():
    # oblate: R_eq close to r_max, R_pol comfortably inside
    r = diag.domain_overflow_check(R_eq=0.95, R_pol=0.5, r_max=1.0, tol=0.1)
    assert r["overflow"]
    assert r["overflow_eq"]
    assert not r["overflow_pol"]


def test_flags_polar_overflow_prolate_star():
    # prolate: R_pol close to r_max, R_eq comfortably inside -- the
    # direction a pure toroidal-field-dominated star deforms toward
    r = diag.domain_overflow_check(R_eq=0.5, R_pol=0.95, r_max=1.0, tol=0.1)
    assert r["overflow"]
    assert r["overflow_pol"]
    assert not r["overflow_eq"]


def test_tol_boundary():
    # exactly at the tol boundary should not flag (strict >, not >=)
    r = diag.domain_overflow_check(R_eq=0.90, R_pol=0.5, r_max=1.0, tol=0.1)
    assert not r["overflow_eq"]


if __name__ == "__main__":
    test_no_overflow_when_comfortably_inside()
    test_flags_equatorial_overflow_oblate_star()
    test_flags_polar_overflow_prolate_star()
    test_tol_boundary()
    print("OK")
