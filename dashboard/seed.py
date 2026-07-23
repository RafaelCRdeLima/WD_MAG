"""Chute de raio para semear o perfil inicial do SCF — nao e' fisica (o SCF
(rho_c,k0) deixa o raio de verdade emergir da convergencia), so' precisa ter
a ordem de grandeza certa para o chute inicial nao ser absurdo. Compartilhado
entre as paginas para nao duplicar o ajuste."""


def r_guess(rho_c):
    """Ajuste grosseiro a partir da relacao M-R de anas brancas degeneradas
    (mu_e=2), calibrado contra tests/test_chandrasekhar_shooting.py."""
    return 1.09e9 * (rho_c / 1.0e6) ** (-0.2436)
