"""V-R0 (D0, prompt de rotacao): com rotation=poloidal=toroidal=None, o
hachisu_scf() da arquitetura de termos plugaveis tem que reproduzir os
resultados de ANTES da refatoracao BIT A BIT -- nao "proximo", nao
"dentro de 1e-10", identico. Compara diretamente contra uma copia
congelada do scf.py pre-refatoracao (_frozen_scf_pre_refactor.py) para ter
uma referencia externa de verdade, nao apenas "o codigo concorda consigo
mesmo".

Este e' o teste que tem que passar ANTES de qualquer fisica nova (rotacao,
toroidal autoconsistente) ser escrita -- ver D0 no prompt de fisica.
"""

import importlib.util
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from scf import hachisu_scf, initial_guess
from terms.poloidal import Poloidal


def _load_frozen_old_scf():
    path = os.path.join(os.path.dirname(__file__), "_frozen_scf_pre_refactor.py")
    spec = importlib.util.spec_from_file_location("_frozen_scf_pre_refactor", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _compare(rho_c, R_guess, k0, nr=81, ntheta=33, lmax=16, tol=1e-8, max_iter=200):
    scf_old = _load_frozen_old_scf()
    r = np.linspace(0, 1.3 * R_guess, nr)
    theta = np.linspace(0, np.pi, ntheta)
    rho0 = initial_guess(r, theta, rho_c, R_guess)

    old = scf_old.hachisu_scf(rho0.copy(), r, theta, rho_c, k0=k0, lmax=lmax, tol=tol, max_iter=max_iter)
    poloidal = Poloidal(k0=k0, lmax=lmax) if k0 != 0.0 else None
    new = hachisu_scf(rho0.copy(), r, theta, rho_c, poloidal=poloidal, lmax=lmax, tol=tol, max_iter=max_iter)

    assert old["converged"] and new["converged"]
    assert old["iterations"] == new["iterations"], (
        f"k0={k0}: iteration count changed ({old['iterations']} -> {new['iterations']}) "
        "-- the refactor must not alter the numerical path, only repackage it"
    )
    for key in ("rho", "Phi", "H", "u"):
        assert np.array_equal(old[key], new[key]), (
            f"k0={k0}: field '{key}' is not bit-for-bit identical to the pre-refactor result "
            f"(max abs diff = {np.max(np.abs(old[key] - new[key])):.3e})"
        )
    assert old["C"] == new["C"], f"k0={k0}: Bernoulli constant C differs"
    assert old["history"] == new["history"], f"k0={k0}: convergence history differs"


def test_v_r0_no_field_bit_for_bit():
    _compare(rho_c=1e9, R_guess=3.0e8, k0=0.0)


def test_v_r0_poloidal_bit_for_bit():
    _compare(rho_c=1e9, R_guess=3.0e8, k0=1.5e-13)


if __name__ == "__main__":
    test_v_r0_no_field_bit_for_bit()
    test_v_r0_poloidal_bit_for_bit()
    print("OK")
