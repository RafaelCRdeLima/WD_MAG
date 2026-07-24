"""Termo de Bernoulli poloidal, M_pol(u) = k0*u, f(u)=k0 constante (Lander
& Jones 2009). Envolve gradshafranov.solve_gradshafranov() -- fisica
inalterada em relacao ao scf.hachisu_scf() anterior a' refatoracao de
termos plugaveis (D0); este modulo so' reempacota isso como um termo. Ver
docs/teoria.md secao 1.4 para a derivacao de Grad-Shafranov e o historico
do bug de expoentes na funcao de Green.

Sinal no Bernoulli combinado (ver terms/__init__.py): M_pol SOMA a H
(H = C - Phi + M_pol + ...), igual ao H = C - Phi + M_u de antes.
"""

import numpy as np


class Poloidal:
    def __init__(self, k0, lmax=16):
        self.k0 = float(k0)
        self.lmax = lmax
        self.u = None

    def update(self, rho, r, theta, **kwargs):
        from gradshafranov import solve_gradshafranov
        lmax = kwargs.get("lmax", self.lmax)
        if self.k0 == 0.0:
            self.u = np.zeros_like(rho)
            return
        omega2 = (r[:, None] * np.sin(theta)[None, :]) ** 2
        source = -4 * np.pi * omega2 * rho * self.k0
        self.u = solve_gradshafranov(source, r, theta, lmax=lmax)

    def potential(self, r, theta):
        if self.u is None:
            return np.zeros((len(r), len(theta)))
        return self.k0 * self.u

    def energy(self, rho, r, theta):
        """E_pol = int (Br^2+Bth^2)/8pi dV, do campo poloidal derivado de
        u (diagnostics.poloidal_field)."""
        import diagnostics as diag
        Br, Bth = diag.poloidal_field(self.u, r, theta)
        E_pol = diag.volume_integral((Br**2 + Bth**2) / (8 * np.pi), r, theta)
        return {"E_pol": E_pol, "Br": Br, "Bth": Bth, "u": self.u}
