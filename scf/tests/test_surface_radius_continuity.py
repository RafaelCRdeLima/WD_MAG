"""Teste de regressao GENERICO para a classe de bug encontrada em
diagnostics.surface_radius() (ver seu docstring e docs/teoria.md
Sec 1.11): qualquer criterio de superficie que interpole em rho (em vez
de H) degenera, porque rho e' clipado a exatamente 0.0 alem da superficie
(eos.density_of_enthalpy) -- a "interpolacao" vira sempre o ponto de
grade cru, e o raio reportado fica preso numa escada de largura ~dr
conforme um parametro varia continuamente.

Este teste nao verifica UMA funcao especifica -- verifica a ASSINATURA do
bug (patamares) num raio de superficie reportado, varrendo Omega_c
continuamente (o mesmo cenario que revelou o bug originalmente: uma
sequencia de rotacao cujo R_pol/R_eq caia em razoes de inteiros pequenos).
Qualquer reintroducao futura deste bug (nesta funcao ou em outra) deveria
fazer este teste falhar.

Confirmado que o teste discrimina: rodado contra uma reimplementacao
standalone da formula ANTIGA (rho-based, nao importada de lugar nenhum --
so' para provar que o teste pegaria o bug) da' fracao de valores unicos de
0.23 (degraus); a versao atual (H-based) da' 1.0 (continuo)."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from scf import hachisu_scf, initial_guess
from terms.rotation import Rotation
import diagnostics as diag


def test_r_eq_varies_continuously_with_omega_c_no_grid_plateaus():
    rho_c = 1e10
    R_guess = 1.6e8
    nr, ntheta = 97, 49
    r = np.linspace(0, 1.8 * R_guess, nr)
    theta = np.linspace(0, np.pi, ntheta)
    rho_seed = initial_guess(r, theta, rho_c, R_guess)

    omegas = np.linspace(0.0, 4.4, 30)
    R_eqs = []
    for Omega_c in omegas:
        rotation = Rotation(Omega_c=float(Omega_c), A=float("inf")) if Omega_c > 0 else None
        result = hachisu_scf(rho_seed, r, theta, rho_c, rotation=rotation,
                              lmax=16, tol=1e-8, max_iter=300)
        assert result["converged"]
        rho_seed = result["rho"]
        R_eq, _ = diag.equatorial_polar_radii(result["H"], r, theta)
        R_eqs.append(R_eq)

    R_eqs = np.array(R_eqs)
    unique_frac = len(np.unique(np.round(R_eqs, 3))) / len(R_eqs)
    print(f"unique fraction: {unique_frac:.3f}  "
          f"(buggy rho-based reference measured at 0.233 on this same scan)")
    # generous margin above the measured buggy value (0.233) and well
    # below the fixed value (1.0) -- catches a plateau regression without
    # being fragile to ordinary floating-point noise
    assert unique_frac > 0.7, (
        f"R_eq shows only {unique_frac:.1%} unique values across a continuous "
        "Omega_c scan -- looks grid-quantized (the surface_radius bug class, "
        "see module docstring), not continuous"
    )

    # also check monotonic growth (physically expected here: rotation
    # inflates R_eq) with no exact-repeat run of length > 2 (a plateau)
    diffs = np.diff(R_eqs)
    assert np.all(diffs >= 0), "R_eq should grow monotonically with Omega_c in this range"
    run_len = 1
    max_run = 1
    for d in diffs:
        run_len = run_len + 1 if d == 0.0 else 1
        max_run = max(max_run, run_len)
    assert max_run <= 2, f"found a plateau of {max_run} consecutive identical R_eq values"


if __name__ == "__main__":
    test_r_eq_varies_continuously_with_omega_c_no_grid_plateaus()
    print("OK")
