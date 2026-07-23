"""Valida eos.sound_speed contra diferenca finita de P(rho) (dP/drho numerico)."""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from eos import pressure, density, sound_speed


def test_sound_speed_matches_finite_difference():
    for x0 in [0.1, 1.0, 5.0, 20.0]:
        dx = x0 * 1e-6
        P1, P2 = pressure(x0 - dx), pressure(x0 + dx)
        rho1, rho2 = density(x0 - dx), density(x0 + dx)
        cs2_numeric = (P2 - P1) / (rho2 - rho1)
        cs_analytic = sound_speed(x0)
        rel_err = abs(cs_analytic**2 - cs2_numeric) / cs2_numeric
        assert rel_err < 1e-6, f"x={x0}: erro relativo {rel_err}"


if __name__ == "__main__":
    test_sound_speed_matches_finite_difference()
    print("OK")
