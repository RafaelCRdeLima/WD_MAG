"""Campo PURAMENTE TOROIDAL autoconsistente, B_phi = K rho^m varpi^(2m-1),
m >= 1 (forma em lei de potencia; ver p.ex. Kiuchi & Yoshida 2008). NAO e'
a mesma coisa que toroidal.py::impose_toroidal() (D6) -- ver a tabela no
final deste docstring. Este ramo nao tem componente poloidal e e' resolvido
DENTRO do loop do SCF, nao imposto depois da convergencia.

---------------------------------------------------------------------------
DERIVACAO (a partir da forca de Lorentz, verificada simbolicamente com
SymPy para m=1, 2, 3/2 -- residuo exatamente 0 -- e por checagem fisica de
sinal; ver scf/tests/test_toroidal_sc.py)
---------------------------------------------------------------------------

Campo axissimetrico puramente toroidal, B = B_phi(r,theta) e_phi. Em
unidades gaussianas, a forca de Lorentz por volume e'
f_L = (1/4pi) (rot B) x B. Com

    (rot B)_r     = (1/(r sin theta)) d(B_phi sin theta)/dtheta
    (rot B)_theta = -(1/r) d(r B_phi)/dr
    (rot B)_phi   = 0

o produto vetorial da' f_L = (1/4pi) ( (rot B)_theta B_phi , -(rot B)_r
B_phi , 0 ). Definindo chi = varpi * B_phi (varpi = r sin theta), as duas
componentes se reescrevem, apos substituir as derivadas de (r B_phi) e
(B_phi sin theta) por derivadas de chi:

    f_L,r     = -(1/(4 pi varpi)) B_phi d(chi)/dr
    f_L,theta = -(1/(4 pi r^2 sin(theta) varpi)) B_phi d(chi)/dtheta

que e', com B_phi = chi/varpi,

    f_L = -(1/(8 pi varpi^2)) grad(chi^2)     [nas componentes r, theta]

Isso e' um gradiente PERFEITO f_L = -rho grad(Psi) se e somente se
chi^2 = varpi^2 B_phi^2 for funcao apenas de s = rho*varpi^2 (a condicao
de barotropia citada no prompt de fisica) -- porque entao
d(chi^2) = g'(s) ds com s=s(r,theta), e

    dPsi = (1/(8 pi rho varpi^2)) d(chi^2) = [g'(s)/(8 pi s)] ds

e' integravel em s sozinho. Substituindo B_phi = K rho^m varpi^(2m-1):

    chi = varpi B_phi = K rho^m varpi^(2m) = K s^m   (usando s=rho*varpi^2,
                                                        varpi^(2m)=s^m/rho^m
                                                        so' funciona limpo
                                                        porque rho^m se
                                                        cancela -- conferir:
                                                        chi = K rho^m
                                                        varpi^(2m), e
                                                        s^m = rho^m
                                                        varpi^(2m), logo
                                                        chi = K s^m)
    chi^2 = K^2 s^(2m)  =>  g(s) = K^2 s^(2m) ...

(nota: g(s) aqui e' chi^2(s), e o m na formula final vem de integrar
g'(s)/(8 pi s) = 2m K^2 s^(2m-1)/(8 pi s) = (m K^2/(4 pi)) s^(2m-2))

Integrando em s:

    M_tor(s) = Psi(s) = m K^2 / (4 pi (2m-1)) * s^(2m-1) ,   m >= 1

(m>=1 evita a divergencia logaritmica em m=1/2; a normalizacao Psi(0)=0
casa com a convencao M_pol(0)=0 ja usada no termo poloidal).

Balanco de momento grad(P)/rho = -grad(Phi) + f_L/rho = -grad(Phi) -
grad(Psi) da' grad(H + Phi + M_tor(s)) = 0, ou seja

    H + Phi + M_tor(s) = C

-- batendo exatamente com a equacao mestra do prompt de fisica quando
C_rot=M_pol=0. M_tor entra portanto com sinal MENOS relativo a H
(H = C - Phi - M_tor(s) + ...), o OPOSTO do sinal de M_pol -- isso nao e'
erro de digitacao, e foi checado de duas formas independentes: (1)
verificacao simbolica SymPy de f_L = -rho grad(Psi) para m=1, 2, 3/2
(residuo exatamente 0); (2) a checagem fisica de sinal exigida pelo
prompt -- ligar o campo a rho_c fixo tem que AUMENTAR a massa (ver
test_toroidal_sc.py :: test_mass_increases_with_field).

NAO decorre de M_tor depender de rho (correcao de uma nota anterior deste
docstring, que atribuia a causa errada). rho entrar em M_tor e' o motivo
de o passo 9 do loop virar busca de raiz em vez de inversao direta da EOS
(uma questao de ESTRUTURA do algoritmo) -- nao tem nada a ver com o SINAL
do termo no Bernoulli ou na identidade do virial. O sinal e' convencao
pura, fixada no momento em que se escreve "+M_tor" (em vez de "-M_tor")
do lado esquerdo do Bernoulli mestre; se a convencao tivesse sido escrita
com o sinal oposto, a identidade do virial abaixo sairia com sinal
oposto tambem, mesmo M_tor continuando a depender de rho do mesmo jeito.

---------------------------------------------------------------------------
Identidade do virial magnetico (V-R4) -- sinal e por que
---------------------------------------------------------------------------
int rho grad(M_tor).r dV = - int B_phi^2/(8 pi) dV     (sinal NEGATIVO)

Derivacao: o tensor de tensoes de Maxwell da' a identidade GERAL (nao
depende de nenhuma convencao de Bernoulli) int r.f_L dV = int B^2/(8pi) dV
para um campo puramente toroidal, via divergencia do tensor e
trace(T)=-B^2/(8pi). Como f_L = -rho*grad(M_tor) (ver derivacao acima),
substituindo: int r.(-rho*grad(M_tor)) dV = int B_phi^2/(8pi) dV, ou seja
int rho*grad(M_tor).r dV = -int B_phi^2/(8pi) dV. O sinal negativo vem
inteiramente de f_L = -rho*grad(M_tor) ter um MENOS -- que por sua vez
vem da convencao "+M_tor" no Bernoulli mestre, nao de rho aparecer dentro
de M_tor. Verificado numericamente com convergencia limpa em 4 resolucoes
(2.27% -> 1.00% -> 0.36% -> 0.14%, ver test_toroidal_sc.py::test_v_r4).

---------------------------------------------------------------------------
Este ramo NAO usa Grad-Shafranov
---------------------------------------------------------------------------
Sem funcao de fluxo, sem operador Delta*, sem funcao de Green para u.
B_phi e' uma funcao ALGEBRICA do rho local, nao a solucao de uma EDP --
gradshafranov.py fica inteiramente fora deste caminho.

---------------------------------------------------------------------------
impose_toroidal() (D6, toroidal.py)      vs    ToroidalSC (aqui)
---------------------------------------------------------------------------
Quando:    depois da convergencia do SCF        dentro do loop do SCF
Para que:  dado inicial do Castro (twisted torus) sequencia de equilibrio
Poloidal:  obrigatorio                           ausente
Bt/Bp:     livre, imposto a um alvo              fixado pela barotropia
                                                  (nao e' um botao livre)
---------------------------------------------------------------------------
Nao funda os dois. Ver D6 em plano_wd_magnetizada.md.
"""

import numpy as np


class ToroidalSC:
    def __init__(self, K, m=1.0):
        if m < 1:
            raise ValueError("m >= 1 exigido (ver derivacao no docstring do modulo; "
                              "tambem evita singularidade 1/(2m-1) em m=1/2)")
        self.K = float(K)
        self.m = float(m)

    def _coef(self):
        return self.m * self.K**2 / (4 * np.pi * (2 * self.m - 1))

    def M_tor(self, s):
        """M_tor(s) = m K^2/(4 pi (2m-1)) * s^(2m-1), s = rho*varpi^2 >= 0."""
        s = np.clip(s, 0.0, None)
        return self._coef() * s ** (2 * self.m - 1)

    def dM_tor_ds(self, s):
        """d(M_tor)/ds = (2m-1) * coef * s^(2m-2) -- usado so' em testes/
        verificacao (V-R4); o loop principal nao precisa da derivada."""
        s = np.clip(s, 0.0, None)
        return (2 * self.m - 1) * self._coef() * s ** (2 * self.m - 2)

    def B_phi(self, rho, varpi):
        # REVISADO (varredura da classe de bug de surface_radius, ver
        # docs/teoria.md Sec 1.11): rho_safe aqui e' so' um clamp de
        # seguranca numerica (evita base negativa em ruido de ponto
        # flutuante elevada a uma potencia fracionaria m) -- NAO e'
        # deteccao de superficie, nenhuma interpolacao envolvida.
        rho_safe = np.where(rho > 0, rho, 0.0)
        with np.errstate(invalid="ignore"):
            return self.K * rho_safe ** self.m * varpi ** (2 * self.m - 1)

    def update(self, rho, r, theta, **kwargs):
        """Nada persistente a atualizar -- B_phi/M_tor sao funcoes
        algebricas de rho, recalculadas sob demanda (potential_of_rho,
        energy)."""
        pass

    def potential_of_rho(self, rho_trial, r, theta):
        """-M_tor(rho_trial*varpi^2) -- a contribuicao ao Bernoulli (ja
        com o sinal trocado) para um rho TENTATIVO, usada dentro da busca
        de raiz por ponto em scf.py (o potencial deste termo depende do
        proprio rho sendo resolvido, ao contrario de rotation/poloidal)."""
        varpi = r[:, None] * np.sin(theta)[None, :]
        s = rho_trial * varpi**2
        return -self.M_tor(s)

    def energy(self, rho, r, theta):
        import diagnostics as diag
        varpi = r[:, None] * np.sin(theta)[None, :]
        Bphi = self.B_phi(rho, varpi)
        E_tor = diag.volume_integral(Bphi**2 / (8 * np.pi), r, theta)
        return {"E_tor": E_tor, "Bphi": Bphi}
