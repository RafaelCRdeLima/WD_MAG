"""Termo centrifugo de Bernoulli para rotacao diferencial j-constante.

    Omega(varpi) = Omega_c * A^2 / (A^2 + varpi^2)

Rotacao rigida e' o limite A -> infinito. Implementado ANALITICAMENTE
(A=inf aciona a formula rigida direto), nunca por um A grande-mas-finito
numerico -- isso arriscaria inf/inf ou overflow em A^2.

Derivacao (balanco hidrostatico no referencial rotante): a forca centrifuga
por unidade de massa e' Omega^2(varpi) varpi, radialmente para fora em
varpi. Seu potencial de Bernoulli C_rot obedece
d(C_rot)/dvarpi = Omega^2(varpi) varpi, ou seja

    C_rot(varpi) = int_0^varpi Omega^2(varpi') varpi' dvarpi'

Para a lei j-constante, a substituicao w = varpi'^2 reduz a integral a uma
funcao racional elementar, dando a forma fechada:

    C_rot(varpi) = (Omega_c^2 A^2 / 2) * varpi^2 / (A^2 + varpi^2)

Quando A -> infinito: varpi^2/(A^2+varpi^2) -> varpi^2/A^2 (dividindo
numerador e denominador por A^2 e descartando o termo varpi^2/A^2, que vai
a zero, no denominador), entao C_rot -> (1/2) Omega_c^2 varpi^2 -- o
potencial de rotacao rigida, recuperado como limite analitico genuino, nao
como um substituto numerico de A grande.

Sinal no Bernoulli combinado (ver terms/__init__.py): C_rot SOMA a H
(H = C - Phi + C_rot + ...), igual a' integral de Bernoulli classica de
Hachisu (1986) para estrelas rotantes, H + Phi - C_rot = C.

LIMITACAO CONHECIDA (rotacao rigida, perto do fim da sequencia): impor
Omega_c diretamente (em vez de uma razao axial alvo, como Hachisu faz) leva
a uma terminacao NUMERICA prematura da sequencia bem antes do breakup
Kepleriano real -- ver scf/tests/test_rotation.py::test_v_r1_sequence_termination
para os numeros exatos e os suspeitos (parametrizacao por Omega_c,
resolucao radial). Rotacao DIFERENCIAL nao sofre disso na faixa testada
(scf/tests/test_differential_rotation.py, validada a 0.40% contra Yoon &
Langer 2005, longe do breakup) -- a diferenciacao mantem o envelope longe
da ruptura equatorial mesmo com o nucleo girando rapido.
"""

import numpy as np


class Rotation:
    def __init__(self, Omega_c, A=float("inf")):
        if Omega_c < 0:
            raise ValueError("Omega_c deve ser >= 0 (o perfil e' controlado por A, nao pelo sinal de Omega_c)")
        if A <= 0:
            raise ValueError("A deve ser > 0 (ou inf, para rotacao rigida)")
        self.Omega_c = float(Omega_c)
        self.A = float(A)
        self.rigid = np.isinf(self.A)

    def Omega(self, varpi):
        """Perfil de velocidade angular, lei j-constante (rigida se A=inf)."""
        if self.rigid:
            return np.full_like(varpi, self.Omega_c, dtype=float)
        return self.Omega_c * self.A**2 / (self.A**2 + varpi**2)

    def C_rot(self, varpi):
        """Potencial centrifugo em forma fechada (ver docstring do modulo)."""
        if self.rigid:
            return 0.5 * self.Omega_c**2 * varpi**2
        return 0.5 * self.Omega_c**2 * self.A**2 * varpi**2 / (self.A**2 + varpi**2)

    def update(self, rho, r, theta, **kwargs):
        """Nada a atualizar -- C_rot depende so' da geometria fixa da
        malha, nao de rho. Mantido por uniformidade de interface com os
        outros termos."""
        pass

    def potential(self, r, theta):
        varpi = r[:, None] * np.sin(theta)[None, :]
        return self.C_rot(varpi)

    def energy(self, rho, r, theta):
        """T = (1/2) int rho Omega^2(varpi) varpi^2 dV -- a energia
        cinetica rotacional que entra no virial como 2T (ver scf.py e
        diagnostics.virial_error)."""
        import diagnostics as diag
        varpi = r[:, None] * np.sin(theta)[None, :]
        Om = self.Omega(varpi)
        T = 0.5 * diag.volume_integral(rho * Om**2 * varpi**2, r, theta)
        return {"T": T}
