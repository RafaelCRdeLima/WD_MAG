"""Valida eos.neutronization_threshold_rho_c(): os pontos tabelados batem
com Boshkayev, Rueda, Ruffini & Siutsou 2013 (Tabela 2), a interpolacao
entre eles e' monotonica, e fora da faixa tabelada a funcao extrapola
plano (nao inventa uma tendencia)."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from eos import neutronization_threshold_rho_c


def test_reference_points_match_the_paper():
    # 16O, mu_e=2.0 -- o conservador usado para toda a faixa mu_e=2
    # (4He, 12C e 16O compartilham mu_e=2, valores diferentes de limiar;
    # ver docstring)
    assert abs(neutronization_threshold_rho_c(2.0) - 1.94e10) / 1.94e10 < 1e-6
    # 56Fe, mu_e=56/26
    assert abs(neutronization_threshold_rho_c(56 / 26) - 1.18e9) / 1.18e9 < 1e-6


def test_interpolation_is_monotonic_between_references():
    mus = [2.0 + 0.02 * i for i in range(8)]  # 2.0 .. 2.14, dentro da faixa tabelada
    vals = [neutronization_threshold_rho_c(mu) for mu in mus]
    assert all(v2 <= v1 for v1, v2 in zip(vals, vals[1:])), \
        "threshold should decrease monotonically from 16O (mu_e=2) to 56Fe (mu_e~2.154)"


def test_flat_extrapolation_outside_table_range():
    # abaixo do menor mu_e tabelado -> mesmo valor do extremo (16O)
    assert neutronization_threshold_rho_c(1.0) == neutronization_threshold_rho_c(2.0)
    # acima do maior mu_e tabelado -> mesmo valor do extremo (56Fe)
    assert neutronization_threshold_rho_c(3.0) == neutronization_threshold_rho_c(56 / 26)


def test_default_is_the_mu_e_2_conservative_value():
    # mu_e=2.0 e' o padrao do dashboard em toda parte deste projeto
    assert neutronization_threshold_rho_c() == neutronization_threshold_rho_c(2.0)


if __name__ == "__main__":
    test_reference_points_match_the_paper()
    test_interpolation_is_monotonic_between_references()
    test_flat_extrapolation_outside_table_range()
    test_default_is_the_mu_e_2_conservative_value()
    print("OK")
