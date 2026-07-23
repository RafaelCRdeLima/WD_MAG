# Teoria do dashboard `wd-magnetizada`

> Esta é a versão em texto simples (fonte de verdade, boa para `git diff`
> e leitura no terminal). Para PDF com tipografia de livro (equações reais,
> não texto verbatim), ver `teoria.tex` — compilar com
> `xelatex teoria.tex` (duas vezes, para o sumário) — ou abrir `teoria.pdf`
> já compilado. O conteúdo é o mesmo nos três; `teoria.tex` foi escrito à
> mão a partir deste arquivo, não gerado automaticamente, então mudanças de
> conteúdo devem ser replicadas manualmente nos dois.

Este documento explica o que o dashboard calcula: as equações por trás de
cada número e cada figura, e o código que as implementa. Não é uma
introdução a física estelar nem ao MHD — pressupõe que o leitor já conhece
isso. O que este documento faz é a ponte entre a teoria (seção 1) e o código
real (`scf/`, `dashboard/`).

Ver `plano_wd_magnetizada.md` para o plano de projeto completo (decisões
D1–D6, arquitetura, fases). Este documento é sobre a *física implementada*,
não sobre o plano de trabalho.

---

## Tabela de símbolos e unidades

| Símbolo | Significado | Unidade (CGS) | Nome no código |
|---|---|---|---|
| `r`, `θ` | coordenadas esféricas | cm, rad | `r`, `theta` |
| `ϖ` | raio cilíndrico, `ϖ = r sinθ` | cm | `omega` (ver nota abaixo) |
| `z` | altura, `z = r cosθ` | cm | `z` (só em `plots.py`) |
| `ρ` | densidade de massa | g cm⁻³ | `rho` |
| `x` | momento de Fermi normalizado, `x = p_F/(m_e c)` | adimensional | `x` |
| `P` | pressão | dyn cm⁻² = erg cm⁻³ | `pressure(x)` |
| `H` | entalpia específica | erg g⁻¹ = cm² s⁻² | `H` |
| `Φ` | potencial gravitacional | erg g⁻¹ = cm² s⁻² | `Phi` |
| `A_φ` | componente φ do potencial vetor | G cm | `A_phi` |
| `u` | função de fluxo, `u = ϖ A_φ` | G cm³ | `u` |
| `B_r, B_θ, B_φ` | componentes do campo magnético | G (gauss) | `Br`, `Bth`/`Btheta`, `Bphi` |
| `f(u)` | função de corrente poloidal, `f(u) = k₀` | g^(-1/2) cm^(1/2) s^-1 | `k0` |
| `M(u)` | potencial de Bernoulli do poloidal, `M(u) = k₀ u` | erg g⁻¹ | `M_u` (variável local em `scf.py`) |
| `β(u)` | função de corrente toroidal, `β = ϖ B_φ` | G cm | não calculada à parte — `B_φ` já sai pronta de `impose_toroidal()` |
| `ζ` | amplitude do toroidal imposto | (unidades de `u^{-m}` × G) | `zeta` |
| `m` | expoente do toroidal imposto | adimensional, inteiro ≥ 1 | `m_tor` |
| `u_c` | valor de `u` na última linha fechada | G cm³ | `u_c` |
| `C` | constante de Bernoulli | erg g⁻¹ | `C` |
| `ρc` | densidade central (parâmetro de entrada) | g cm⁻³ | `rho_c` |
| `k₀` | amplitude do campo poloidal (parâmetro de entrada) | g^(-1/2) cm^(1/2) s^-1 | `k0` |
| `μₑ` | peso molecular médio por elétron, `Y_e = 1/μₑ` | adimensional | `mu_e` |
| `M` | massa total da estrela | g (exibido em M☉) | `M` (⚠ colide com `M(u)` na notação da espinha — ver nota) |
| `R_eq`, `R_pol` | raios equatorial e polar (onde `H=0`) | cm (exibido em km) | `R_eq`, `R_pol` |
| `W` | energia gravitacional | erg | `W` |
| `Π` | energia interna, `∫P dV` | erg | `Pi` |
| `E_pol`, `E_tor`, `ℳ` | energias magnéticas poloidal, toroidal, total | erg | `E_pol`, `E_tor`, `E_mag` |
| `VE` | erro virial | adimensional | `VE` |
| `G` | constante gravitacional | 6,674×10⁻⁸ cm³ g⁻¹ s⁻² | `G_CONST` |
| `A`, `B` | constantes da EOS | ver §1.1 | `A_CONST`, `B_of_mu_e(mu_e)` |
| `l` | grau da expansão de Legendre | inteiro ≥ 0 (Poisson) ou ≥ 1 (GS) | `l` |
| `t_din`, `v_A`, `t_Alfvén` | escalas de tempo/velocidade derivadas | s, cm/s, s | `t_din`, `v_A`, `t_alf` |

**Notas de notação:**
- O código usa a variável `omega` para `ϖ` (raio cilíndrico), **não** para
  velocidade angular — não há rotação neste projeto (D3, estrela estática).
  É uma escolha de nome específica deste código; se ler `omega` em
  `scf.py`/`gradshafranov.py`, é sempre `ϖ = r sinθ`.
- `M` é usado na espinha teórica para dois objetos diferentes: a massa
  estelar total (§3.7 do prompt original) e o potencial `M(u)` do Bernoulli
  (§3.4–3.5). O código resolve a colisão nomeando a massa `M` (retornada por
  `scf.total_mass()`) e o potencial `M_u` (variável local dentro de
  `hachisu_scf()`, nunca exposta fora da função). Neste documento, sempre que
  houver risco de ambiguidade, o texto diz "massa" ou "potencial de
  Bernoulli" por extenso.
- `Bth` aparece como `Btheta` em `diagnostics.py` (parâmetro da função) e
  `Bth` nas páginas do dashboard (variável local) — mesmo objeto, dois nomes,
  por causa de `flake8`/legibilidade nos dois contextos.

---

## 1. Núcleo teórico comum

### 1.1 Equação de estado

Gás de elétrons completamente degenerado a T = 0:

```
P = A [ x(2x² − 3)(x² + 1)^{1/2} + 3 sinh⁻¹ x ]
ρ = B x³
x ≡ p_F / (m_e c)

A = 6,01 × 10²²  dyn cm⁻²
B = 9,82 × 10⁵ / Y_e  g cm⁻³     (Y_e = 0,5 → B ≈ 1,96 × 10⁶)
```

`A` é uma constante física fixa (combinação de `m_e`, `c`, `ℏ`); `B` depende
da composição só através de `Y_e` (número de elétrons por núcleon). O
dashboard parametriza por `μₑ = 1/Y_e` em vez de `Y_e` diretamente.

→ `scf/eos.py :: pressure()` (P(x)), `density()` (ρ(x)), `B_of_mu_e()` (B(μₑ))

Entalpia específica, analítica — é isso que torna o SCF possível:

```
H = ∫ dP/ρ = (8A/B) [ √(1 + x²) − 1 ]
```

→ `scf/eos.py :: enthalpy()`

Inversa, aplicada a cada iteração do SCF:

```
x = √[ (1 + HB/8A)² − 1 ]        (H ≤ 0 => ρ = 0)
```

→ `scf/eos.py :: x_of_enthalpy()`, `density_of_enthalpy()`

O termo `−1` normaliza H = 0 na superfície (fronteira ρ=0 da estrela).

**Índice politrópico efetivo.** Perto da superfície `x ≪ 1` dá `ρ ∝ H^{3/2}`,
isto é `n = 3/2`. No núcleo relativístico `x ≫ 1` dá `ρ ∝ H³`, isto é `n = 3`.
A EOS desliza entre os dois, e `n = 3` é o caso marginal onde a massa fica
independente da densidade central — a origem do limite de Chandrasekhar. Esse
fato é a causa direta da decisão de parametrização em §1.4/§1.5.

**Velocidade do som.** Adicionada depois da espinha original (não estava no
plano), por ser necessária para a amplitude da perturbação exportada na
Aba 3:

```
c_s = √(dP/dρ) ,   dP/dx = 8A x⁴/√(1+x²)   (derivada analítica de P(x))
```

→ `scf/eos.py :: sound_speed()`, validada contra diferença finita em
`scf/tests/test_sound_speed.py`

---

### 1.2 Gravidade própria

```
∇²Φ = 4πGρ
```

Resolvida por expansão em polinômios de Legendre. Com
`ρ(r,θ) = Σ_l ρ_l(r) P_l(cosθ)`:

```
Φ_l(r) = − (4πG / (2l+1)) [ r^{−(l+1)} ∫₀^r ρ_l(r′) r′^{l+2} dr′
                          + r^{l}     ∫_r^∞ ρ_l(r′) r′^{1−l} dr′ ]
```

Os pesos `r′^{l+2}` e `r′^{1−l}` carregam o `r′²` do elemento de volume.

→ `scf/poisson.py :: solve_poisson()` (implementação exata desta fórmula,
D_l/E_l calculados por soma cumulativa trapezoidal); `legendre_matrix()`
calcula `P_l(cosθ)`.

Validado contra a solução fechada da esfera uniforme e o teorema das cascas
de Newton em `scf/tests/test_poisson.py`.

`G = 6,674 × 10⁻⁸ cm³ g⁻¹ s⁻¹` — constante física (valor de referência
CODATA arredondado), não um parâmetro de projeto. → `scf/poisson.py ::
G_CONST`, repetida em `dashboard/units.py :: G_CONST` (mesmo valor, dois
módulos, porque `poisson.py` não pode depender do dashboard — R1).

---

### 1.3 Campo magnético e a função de fluxo

Com `ϖ = r sinθ`, define-se a função de fluxo

```
u = ϖ A_φ
```

e o campo axissimétrico separa em poloidal mais toroidal:

```
B = (1/ϖ) [ ∇u × ê_φ  +  β(u) ê_φ ]
```

Em componentes esféricas:

```
B_r = (1 / (r² sinθ)) ∂u/∂θ
B_θ = − (1 / (r sinθ)) ∂u/∂r
B_φ = β(u) / ϖ
```

→ `scf/diagnostics.py :: poloidal_field()` implementa `B_r`, `B_θ` (por
diferenças finitas de `u`, com proteção contra a singularidade de coordenada
em `ϖ→0`). `B_φ` não é calculado a partir de um `β(u)` intermediário; sai
diretamente de `scf/toroidal.py :: impose_toroidal()`, que já implementa a
forma funcional escolhida para `β` (ver §1.9) substituída na fórmula acima.

A forma geral `B = (1/ϖ)[∇u×ê_φ + β ê_φ]` como identidade vetorial não
aparece escrita em nenhum lugar do código — o código vai direto às três
componentes escalares.

**Interpretação geométrica: os contornos de `u` são as linhas de campo
poloidal.** É por isso que a Aba 1 os plota (§2).

---

### 1.4 Grad-Shafranov

```
Δ* u = − 4π ϖ² ρ f(u) − β β′(u)

Δ* = ∂²/∂ϖ² − (1/ϖ) ∂/∂ϖ + ∂²/∂z²
```

`f(u) = dM/du` é a função de corrente que gera o poloidal; `β(u)` gera o
toroidal. A escolha mais simples, e a usada aqui, é `f(u) = k₀` constante
(Lander & Jones 2009). Como o toroidal é imposto depois de o poloidal
convergir (§1.9, D6) e não entra na própria equação GS resolvida pelo SCF, o
termo `−ββ′(u)` **não aparece na fonte que o código de fato monta** — a fonte
implementada é só `−4πϖ²ρk₀`.

→ `scf/gradshafranov.py :: solve_gradshafranov()` — fonte passada por
`scf/scf.py :: hachisu_scf()` como `source = -4*np.pi*omega2*rho*k0`.

A densidade de corrente que corresponde a essa fonte é

```
J_φ = c ρ ϖ f(u)
```

Esta fórmula **não corresponde a nenhuma função no código** — `J_φ` nunca é
calculada como quantidade nomeada. Ela foi usada apenas analiticamente,
durante a depuração que achou o bug descrito no box abaixo, para montar o
lado direito da lei de Ampère integral (§1.4, testes sem derivada).

**Resolução por expansão de Legendre associada.** O operador `Δ*` separa em
funções `P_l¹(cosθ)`, **não** em `P_l(cosθ)`. A estrutura radial da função de
Green é análoga à de Poisson, mas os pesos das integrais **não são os
mesmos**:

```
pesos em Poisson (para Φ):         r′^{l+2}   e   r′^{1−l}
pesos em Grad-Shafranov (para u):  r′^{l+1}   e   r′^{−l}
```

→ `scf/gradshafranov.py :: solve_gradshafranov()`, linhas onde `D_l`/`E_l`
são montados com `r ** (l + 1)` e `r ** (-l)`.

> **Por que assim — os expoentes da função de Green (G4).**
> A diferença de uma potência vem do peso em `ϖ` da fonte da GS comparado ao
> `r′²` do elemento de volume em Poisson, e depende de a equação ser
> resolvida para `u` ou para `A_φ`. **Este foi um bug real neste projeto:**
> uma versão anterior de `gradshafranov.py` usava `r′^{l+2}` e `r′^{1−l}`
> (copiados por engano da estrutura de `poisson.py`, que tem uma potência de
> `r` a mais por causa da substituição `χ = rΦ_l` usada na redução do
> Laplaciano escalar — `Δ*` não precisa dessa substituição). O bug inflava o
> campo poloidal por um fator ~7000× e a energia magnética por ~5×10⁷.
>
> O bug sobreviveu a toda a validação peça-a-peça — inclusive uma forma
> fechada derivada "à mão" — porque essa forma fechada foi construída
> reusando a mesma equação indicial (`r^{l+1}`, `r^{−l}` para as soluções
> homogêneas, que **estavam certas**) e por isso herdava a mesma normalização
> das integrais internas, que estava errada. Só dois testes que não usam a
> função de Green nem uma segunda derivada — consistência via lei de Ampère —
> revelaram o fator fixo. Ver a nota completa (com a cadeia de raciocínio) no
> topo de `scf/gradshafranov.py`.
>
> **Lição operacional:** "resolve com a fórmula, confere com a mesma
> fórmula" não é validação independente, mesmo em arquivos/funções
> diferentes, se ambos herdam a mesma normalização de uma derivação
> compartilhada.

**Testes de normalização sem derivada** — a defesa contra essa classe de
erro:

```
Fluxo:    u(ϖ,z) = ∫₀^ϖ B_z(ϖ′,z) ϖ′ dϖ′
Ampère:   oint_C B·dl = 4π k₀ ∫_S ρ ϖ dA
```

- **Ampère**: implementado e no repositório → `scf/tests/test_gradshafranov.py
  :: test_ampere_law()`. Usa um caso sintético com fonte de modo `l=1` puro
  (não a EOS/SCF completos), calcula `oint B·dl` num laço retangular em `(r,θ)` e
  compara com `−∫∫ [source/sinθ] dr dθ` (forma equivalente derivada da lei de
  Ampère para esta geometria). Fecha a ~2%.
- **Fluxo**: foi *proposto* durante a depuração (ver histórico da conversa
  que corrigiu o bug) mas **não foi implementado como teste separado** — o
  teste de Ampère sozinho já revelou e confirmou a correção do bug, então o
  teste de consistência de fluxo nunca chegou a ser escrito. Fica registrado
  como lacuna em §7.

---

### 1.5 Bernoulli e a parametrização

```
H + Φ − M(u) = C ,        M(u) = ∫ f(u) du
```

Com `f = k₀` constante, `M(u) = k₀ u`.

A constante de integração é fixada **no centro**:

```
C = H(ρc) + Φc − M(u_c)
```

onde `Φc`, `u_c` aqui denotam os valores no centro (`r=0`) — **atenção**:
este uso de `u_c` (valor de `u` no centro) é diferente do `u_c` usado em
§1.9 (valor de `u` na última linha fechada, na superfície). São o mesmo
símbolo para dois pontos de avaliação diferentes; o código nunca precisa dos
dois ao mesmo tempo, mas o leitor deve notar a colisão. Como `ϖ = r sinθ` se
anula no centro para qualquer `θ`, `u_c(centro) = 0` sempre e o termo
`M(u_c)` desaparece — ele fica na fórmula por generalidade, não porque
contribua.

→ `scf/scf.py :: hachisu_scf()` — linha `C = H_c + Phi[0, 0] - M_u[0, 0]`.

> **Por que assim — a parametrização por (ρc, k₀) (G4).**
> A receita clássica de Hachisu fixa `C` impondo `H = 0` em dois pontos da
> superfície (polo e equador). Isso é mal-posto para esta EOS: como `n → 3`
> no núcleo (§1.1), a massa fica quase independente da densidade central e o
> raio despenca, de modo que resolver para `C` a partir do raio inverte um
> mapa quase singular. Testado e confirmado neste projeto: a iteração de
> Picard sob essa receita é linearmente instável em toda a faixa
> `ρc = 10⁶`–`10¹⁰` g/cm³, com ou sem sub-relaxação.
>
> A parametrização adotada é por **(ρc, k₀)**, duas entradas independentes,
> sem nenhuma condição de superfície — `C` é fixada por uma condição local
> (`H` no centro), não global (integral sobre toda a massa). Essa mudança
> por si só resolveu a instabilidade: convergência geométrica até precisão
> de máquina em ~12 iterações, contra divergência exponencial da receita
> antiga.
>
> **Ressalva:** no regime de campo forte, quando o pico de densidade migra
> para fora do centro, esta âncora também perde o sentido físico que a torna
> estável — ver §6.

---

### 1.6 O loop SCF

```
1.  ρ ← palpite inicial (politropo esférico n = 3)
2.  A_φ ← 0
3.  repetir:
4.      Φ   ← Poisson(ρ)
5.      A_φ ← GradShafranov(ρ, f)
6.      u   ← ϖ A_φ ;  M ← ∫ f du
7.      C   ← H(ρc) + Φc − M(u_c)
8.      H   ← C − Φ + M(u)
9.      ρ_novo ← EOS⁻¹(H)          [H ≤ 0 => ρ = 0]
10.     ρ   ← (1−ω) ρ + ω ρ_novo        ω ≈ 0,3
11. até max|Δρ|/ρc < tol
```

→ `scf/scf.py :: hachisu_scf()` implementa os passos 1, 3–9 e 11
literalmente (passo 2 é implícito: `u` começa em zero antes da primeira
iteração).

**Discrepância com o código real, passo 10:** a implementação **não faz
sub-relaxação**. A linha correspondente em `hachisu_scf()` é `rho = rho_new`
— substituição direta, equivalente a `ω = 1`, não `ω ≈ 0,3`. Não há parâmetro
`ω` na assinatura da função. Isso não é um descuido: é consequência direta da
mudança descrita no box de §1.5. A receita antiga (duas condições de
superfície) exigia sub-relaxação para tentar estabilizar uma iteração que era
instável de qualquer forma; a parametrização `(ρc, k₀)` converge por
substituição direta porque a instabilidade que a sub-relaxação tentava
mascarar foi removida na raiz. Ver a nota de projeto no topo de
`scf/scf.py` para o histórico completo (inclusive testes que mostraram que
sub-relaxar a receita antiga não resolvia o problema).

O critério de convergência (`tol`) é um parâmetro numérico, não físico — o
dashboard expõe valores de `1e-4` a `1e-8` (Aba 1), com `1e-6` como padrão.

---

### 1.7 Virial e diagnósticos de energia

```
W = ½ ∫ ρ Φ dV                      (gravitacional, negativa)
Π = ∫ P dV                          (interna)
E_pol = ∫ (B_r² + B_θ²) / 8π  dV
E_tor = ∫ B_φ² / 8π  dV
ℳ = E_pol + E_tor
```

→ `scf/diagnostics.py :: gravitational_energy()` (W), `pressure_integral()`
(Π), `magnetic_energies()` (E_pol, E_tor, ℳ — chamada de `E_mag` no
código). Todas usam `volume_integral()` como base (`dV = r² sinθ dr dθ dφ`,
`φ` integrado analiticamente para `2π`).

Teorema do virial escalar para configuração estática:

```
W + 3Π + ℳ = 0
```

Erro virial, usado como portão de qualidade:

```
VE = | W + 3Π + ℳ | / |W|          critério de aceite: VE < 10⁻³
```

→ `scf/diagnostics.py :: virial_error()`. O limite `10⁻³` é convenção de
projeto (V3 do plano), não uma constante física — aparece hard-coded como
comparação em `dashboard/pages/1_equilibrio.py` e `3_exportacao.py`, não em
`units.py` nem em `diagnostics.py` (nenhum módulo de física define esse
limite como constante nomeada — é uma decisão de UI/gate repetida em duas
páginas).

**Identidade do virial magnético.** Exata, e é o teste que liga o setor
magnético ao gravitacional:

```
∫ ρ ∇M(u) · r  dV  =  ∫ B²/8π  dV
```

Esta identidade **não corresponde a nenhuma função no código**. Foi
verificada numericamente de forma ad hoc durante a depuração do bug de
§1.4 (comparando `k0 * ∫ρ r ∂u/∂r dV` com `∫B²/8π dV` para uma solução
convergida), mas não existe uma função `scf/diagnostics.py` que a calcule
nem um teste comitado que a exercite. Fica em §7 como lacuna.

> **Nota conceitual obrigatória.** `M(u)` **não é** a energia magnética. É o
> potencial específico da força de Lorentz — apenas a parcela que entra no
> Bernoulli. Não existe identidade *local* entre `M(u)` e `B²/8π`. A
> identidade *global* acima é o que força as duas quantidades a concordarem
> em ordem de grandeza. Medida neste projeto, no regime linear (campo
> perturbativo), a razão `(E_mag/|W|) / (M_u/H_c)` vale ≈ 0,5 e é constante em
> `k₀` (confirmado dobrando `k₀`: ambas as razões quadruplicam, a razão entre
> elas não muda); ela desliza para ~0,38 conforme o campo deixa de ser
> perturbação pequena (`k₀ ≳ 10⁻¹²` na configuração `ρc=10⁹`, `R≈3×10⁸` cm).
> Essa foi a observação que, junto com o teste de Ampère, confirmou que o bug
> do box de §1.4 era real e não um artefato de unidades.

**Limite físico.** `ℳ < |W|` é vínculo rígido: `E_mag/|W| ≥ 1` é configuração
impossível (violaria o virial com `Π ≥ 0`). Use como âncora de sanidade ao
ler qualquer resultado — se o dashboard mostrar `E_mag/|W|` maior que
alguns décimos, desconfie antes de acreditar.

---

### 1.8 As duas razões Bt/Bp

```
Bt/Bp (energia)    = E_tor / E_pol
Bt/Bp (amplitude)  = max|B_φ| / max|B_pol|
```

→ `scf/toroidal.py :: bt_bp_ratios()`

Diferem por ordens de grandeza porque o toroidal está confinado a um volume
pequeno (§1.9). A literatura frequentemente não diz qual usa. **O dashboard
sempre mostra as duas, rotuladas** (Aba 1 e Aba 3).

---

### 1.9 O campo toroidal e o twisted torus

Na formulação barotrópica, `β = β(u)`: a função toroidal depende só da
função de fluxo. É essa condição que anula a componente φ da força de
Lorentz e permite a existência da integral de Bernoulli (§1.5).

Consequência geométrica: fora da estrela não há corrente, logo `B_φ = 0` lá;
como `β = β(u)`, ela precisa se anular em toda linha de fluxo que escapa da
superfície. **O toroidal fica automaticamente confinado à região de linhas
poloidais fechadas** — o toro torcido não é imposto por decreto geométrico,
ele cai da consistência da equação.

Forma funcional adotada:

```
β(u) = ζ (u − u_c)^{m+1} Θ(u − u_c) ,     m ≥ 1
```

com `u_c` o valor de `u` na última linha fechada (aqui, `u_c` = valor de `u`
na superfície — ver a distinção de notação em §1.5). O expoente `≥ 1` mantém
`ββ′` contínua na borda do toro.

→ `scf/toroidal.py :: impose_toroidal()` implementa esta forma **já
substituída em `B_φ = β/ϖ`**: `B_φ = ζ(u−u_c)^{m+1}/ϖ` para `u > u_c`, `0`
fora. `find_uc()` implementa a busca de `u_c` como o **máximo de `u` ao longo
de toda a superfície estelar** (`H=0`, varrendo `θ` do polo ao equador) — é
essa escolha específica de "última linha fechada" que o código usa; contornos
com `u` maior que esse máximo, por definição, não tocam a superfície em
nenhum ponto e ficam inteiramente internos.

`ζ` não é um parâmetro de entrada direto do usuário — o dashboard pede a
razão `Bt/Bp` (energia) alvo e resolve `ζ` para atingi-la, aproveitando que
`B_φ` é linear em `ζ` (logo `E_tor` é quadrático):

→ `scf/toroidal.py :: solve_zeta_for_energy_ratio()`

> **Por que o toroidal é imposto e não resolvido (G4).**
> A condição `β = β(u)` confina o toroidal a um toro de volume pequeno, e a
> razão de energias resultante — se `β` fosse extraída de uma condição de
> fechamento barotrópica geral em vez de escolhida livremente — sai em
> poucos por cento. **Não é possível atingir Bt/Bp ~ 1/2 por essa via.** Por
> isso o projeto resolve o SCF barotrópico apenas para o poloidal (§1.4–1.6,
> com `f(u)=k₀`) e impõe o toroidal por cima, com a razão desejada
> (`ζ` resolvida para o alvo), carregando no Castro fora de equilíbrio exato
> e relaxando com amortecimento (sponge, Aba 3). Para estudo dinâmico não é
> preciso equilíbrio exato — basta estar perto o suficiente para o transiente
> não destruir a topologia.

**Fração de volume do toro** (quanto da estrela tem `u > u_c`):
→ `scf/toroidal.py :: closed_torus_volume_fraction()`

**Espessura radial do toro** (usada para checar resolução de malha, Aba 3):
→ `scf/toroidal.py :: torus_radial_extent()`, medida ao longo do equador por
padrão.

---

### 1.10 Escalas de tempo e unidades

```
t_din    = √( R_eq³ / GM )
v_A      = ⟨B⟩ / √(4π ρ̄)
t_Alfvén = R_eq / v_A
```

→ `dashboard/units.py :: dynamical_time()`, `alfven_speed()`,
`alfven_time()`. `⟨B⟩` e `ρ̄` (médias volumétricas) são calculadas em
`dashboard/pages/3_exportacao.py` diretamente com `diagnostics.volume_integral()`
— não há uma função `mean_field()` dedicada em `scf/`.

A razão `t_Alfvén / t_din` mede o custo da simulação: no regime de campo
fraco ela chega a 10³–10⁴, que é proibitivo; no regime de campo forte deste
projeto, com `E_mag/|W| ~ 0,1`, ela fica em 2–3 (D4 do plano).

**Unidade natural de campo.** De `B²/8π ~ Gρ²R²`:

```
B_unit = R ρ √(8πG)
```

Para `ρc ~ 10⁹` e `R ~ 10⁸` cm isso dá `~10¹⁴` G, que é também a ordem do
campo virial de uma anã branca. **Esta fórmula não está implementada em
nenhum lugar do código** — foi usada apenas como estimativa de ordem de
grandeza durante a depuração do bug de §1.4, para checar se o campo relatado
fazia sentido físico antes de procurar o erro numérico. **Âncora de
sanidade:** com `E_mag/|W| ~ 0,1`, o campo exibido pelo dashboard deve estar
na casa de `10¹³`–`10¹⁴` G; se aparecer muito diferente disso, desconfie
antes de acreditar (foi exatamente essa desconfiança que achou o bug).

**Convenção do Castro:** o campo é carregado como `B′ = B / √(4π)`. O
dashboard exibe **sempre em gauss** e converte apenas na exportação (Aba 3).

→ `dashboard/units.py :: gauss_to_castro()`, `castro_to_gauss()`. Conferido
neste ciclo de trabalho: o fator bate com o plano (`√(4π) ≈ 3,5449`). **Nota
de estado:** `gauss_to_castro()` existe e está correta, mas **não é chamada
em nenhum lugar do pipeline de exportação atual**
(`dashboard/pages/3_exportacao.py` escreve `B_phi` no HDF5 diretamente em
gauss, com um atributo `units` explícito documentando isso). A conversão
para a convenção `B′` do Castro é responsabilidade do `problem_initialize.H`
do lado do Castro (ainda não escrito — Fase 0 do plano pendente), que lerá o
HDF5 em gauss e aplicará `gauss_to_castro()` (ou o equivalente em C++) ao
montar o estado interno. Ver §7.

Todas as conversões de unidade de exibição (gauss, km) — tanto o número
quanto a formatação da string — vivem em `dashboard/units.py`, ponto único
de verdade (regra R4 do dashboard): `cm_to_km()`, `g_to_msun()`,
`format_gauss()`, `format_km()`, `format_km_value()`.

---

## 2. Aba 1 — Equilíbrio

→ `dashboard/pages/1_equilibrio.py`

Executa uma corrida única do SCF (§1.6) e mostra o resultado. Usa todo o
núcleo teórico (§1.1–1.10).

**Parâmetros de entrada** (barra lateral): `ρc`, `μₑ`, `k₀` (opcional,
liga/desliga campo poloidal), razão `Bt/Bp` alvo e `m` (toroidal, opcional),
mais os parâmetros numéricos da malha (`Nr`, `Ntheta`, `l_max`, `tol`,
`max_iter`). `θ` sempre cobre `[0, π]` inteiro — sem simetria equatorial,
decisão do plano (D3) para não mascarar modos assimétricos (`m=1`).

**Faixa de `k₀`.** Não é conhecida a priori (depende de `ρc`, `R`, `μₑ` de
forma não trivial — ver o bug de §1.4, que por muito tempo fez a faixa
aparente parecer errada por 3–4 ordens de grandeza). O dashboard sonda
empiricamente subindo `k₀` geometricamente até `VE > 10⁻³` ou a SCF parar de
convergir, numa malha grosseira, e grava o resultado num cache em
`dashboard/k0_range_cache.json`, indexado por `(ρc, μₑ)`.
→ `dashboard/pages/1_equilibrio.py :: _estimate_k0_max()`

**Escalares exibidos** (tabela "Escalares"):

| Escalar na tela | Definição | Função |
|---|---|---|
| `M/M☉` | massa total / `M_SUN` | `scf.total_mass()`, `units.M_SUN` |
| `R_eq`, `R_pol` (km) | raio onde `ρ→0`, no equador e no polo | `diagnostics.equatorial_polar_radii()` |
| `R_pol/R_eq` | achatamento | razão direta |
| `ρc confirmado` | `ρ[r=0]` pós-convergência | deve bater com o `ρc` de entrada |
| `ρ média` | `M / volume do elipsoide (R_eq, R_eq, R_pol)` | aproximação geométrica, não integral exata |
| `W`, `E_int=Π`, `E_mag`, `E_pol`, `E_tor` | §1.7 | `diagnostics.virial_error()`, `magnetic_energies()` |
| `E_mag/|W|` | intensidade de campo adimensional | razão direta (ver âncora de sanidade em §1.7/1.10) |
| `B_pol,max`, `B_central`, `B_tor,max` | valores de campo, gauss | máximos/pontuais de `Br`,`Bth`,`Bphi` |
| `fração de volume do toro` | §1.9 | `toroidal.closed_torus_volume_fraction()` |
| `VE` | §1.7 | `diagnostics.virial_error()` |

**Convergência**: dois gráficos (`max|Δρ|/ρc` e `VE`, ambos por
iteração, escala log) → `plots.plot_convergence()`, `plots.plot_virial_history()`.
O histórico de `VE` só existe se `hachisu_scf(..., track_virial=True)` —
custa integrais extras a cada iteração, por isso é opcional (ver §1.6).

**Bt/Bp**: as duas razões de §1.8, lado a lado, sempre rotuladas.

**Figuras do plano meridional** (§1.3):
- **densidade**: mapa de cor de `ρ`, com a fronteira `H=0` (a superfície da
  estrela) sobreposta em ciano → `plots.plot_density()`
- **linhas de campo poloidal**: contornos de `u` — fisicamente, cada
  contorno *é* uma linha de campo poloidal (§1.3), porque `B` é tangente às
  curvas de `u` constante por construção (`B·∇u = 0` decorre diretamente das
  fórmulas de `B_r`,`B_θ` em termos de derivadas de `u`). A última linha
  fechada (`u_c`, §1.9) é destacada em vermelho quando há toroidal imposto →
  `plots.plot_flux_contours()`
- **campo toroidal**: mapa de cor de `B_φ`, mostra o toro confinado →
  `plots.plot_toroidal()`

Todos os eixos radiais em km, campo sempre em gauss com colorbar em notação
científica (regra R4, ver §1.10).

---

## 3. Aba 2 — Varredura

→ `dashboard/pages/2_varredura.py`, `dashboard/sweep_worker.py`

Roda o SCF (§1.6) numa grade de `(ρc, k₀)` em paralelo
(`ProcessPoolExecutor`), com cache por hash dos parâmetros
(`store.run_exists()`/`store.save_run()`, §5). Cada ponto da grade é uma
chamada independente a `scf.hachisu_scf()` seguida dos mesmos diagnósticos
da Aba 1 — nenhuma física nova, só orquestração.

→ `dashboard/sweep_worker.py :: run_one()` — a função picklable que cada
processo da grade executa.

**O que é uma sequência de equilíbrio.** Fixando `ρc` e variando `k₀` (ou
vice-versa), obtém-se uma sequência de configurações de equilíbrio. **A
sequência termina** quando o SCF para de convergir ou quando `VE` ultrapassa
o limite de aceite (§1.7) e não melhora com resolução — ver os números
medidos em §6. Uma terminação de sequência é, em si, um resultado físico
(não um erro a ser corrigido): sinaliza o limite de validade daquela família
de equilíbrios.

**Diagrama M-R.** `R_eq` (km) no eixo x, `M/M☉` no eixo y, colorido por
`k₀`. A reta horizontal em `1,44 M☉` marca o limite de Chandrasekhar
**sem campo** (`μₑ=2`) — não é o limite físico da sequência magnetizada, é
só uma referência de leitura. Overlay opcional de pontos da literatura de
`dashboard/data/referencias/bera_bhattacharya_2014.csv`, se o arquivo
existir (não existe atualmente neste repositório — nenhum dado foi
digitalizado; o dashboard funciona sem ele e avisa que falta).

**Por que o máximo de massa relevante é sobre o plano inteiro, não uma
fatia.** Este é o ponto que mais gera comparação errada com a literatura.
`M_max(k₀=0)` já depende de `ρc`: só se aproxima do limite de Chandrasekhar
assintoticamente, para `ρc` alto (`~10¹¹`–`10¹²` g/cm³ — ver
`scf/tests/test_scf_v1.py`, validado a 0,78% nessa faixa). Uma fatia de `ρc`
baixo (por exemplo `ρc=10⁹`, usada durante a depuração deste projeto) dá
`M(k₀=0) = 1,39 M☉` — **já abaixo** do limite de Chandrasekhar sem campo, e
qualquer massa medida ao longo dessa fatia (com ou sem campo) não é
comparável ao número de referência da literatura (`M_max ~ 1,9 M☉` em Bera &
Bhattacharya 2014), que é o máximo tomado sobre **todo** o plano `(ρc, k₀)`,
não sobre uma reta `k₀` variável com `ρc` fixo. A varredura desta aba é
exatamente o que permite fazer essa comparação corretamente — mapeando o
plano, não uma fatia.

**Convergência da grade.** Pontos que não convergem são registrados (não
descartados silenciosamente) e mostrados num expansor separado; a taxa de
convergência da grade é, ela mesma, informação sobre onde a família de
soluções termina.

**Mapa de calor de VE.** Sobre a grade `(ρc, k₀)`, em `log₁₀(VE)` — revela
onde o método está no limite da validade sem precisar ler número por número.

---

## 4. Aba 3 — Exportação

→ `dashboard/pages/3_exportacao.py`

Pega um equilíbrio já salvo (Aba 1 ou 2), impõe o toroidal (§1.9) e gera os
três artefatos de saída: HDF5 de dado inicial, `inputs` do Castro,
`run_manifest.json`.

**Imposição do toroidal.** Mesma função de §1.9
(`toroidal.solve_zeta_for_energy_ratio()`), com controles próprios desta
aba (a razão alvo pode ser diferente da usada quando o equilíbrio foi
salvo). A continuidade de `ββ′` na borda do toro é **garantida
analiticamente** para `m ≥ 1` (não é recalculada numericamente aqui — é uma
propriedade da forma funcional de §1.9, `(u-u_c)^{m+1}` e sua derivada vão a
zero em `u=u_c` para qualquer `m≥1`); a página só confirma que `m_tor≥1` foi
respeitado.

**Por que o campo é inicializado pelo potencial vetor, não por B nos
centros.** Isto é uma decisão do lado Castro (D3 do plano, seção 5), não
algo que o SCF resolve — mas o dashboard já prepara os dados nesse formato:
o HDF5 exportado contém `A_φ` (calculado como `u/ϖ` a partir da função de
fluxo convergida), **não** `B_r`, `B_θ` diretamente. O algoritmo de
transporte restrito do Castro mantém `B` nas faces da malha; inicializar via
`∇×A` interpolado nas arestas garante `∇·B = 0` na precisão de máquina por
construção. Inicializar com valores de `B` nos centros das células não dá
essa garantia. **O que falta:** a interpolação de `A_φ` para as arestas da
malha cartesiana e o cálculo de `∇×A` ali é responsabilidade do
`problem_initialize_mhd_data.H` do Castro, que ainda não foi escrito (Fase 0
do plano pendente) — o dashboard entrega `A_φ` numa malha esférica `(r,θ)`;
a interpolação para a malha cartesiana do Castro acontece do outro lado.

**Convenção B′ = B/√(4π).** Ver §1.10. O HDF5 exportado por esta aba
guarda `B_φ` em gauss puro, com um atributo `units` no cabeçalho dizendo
isso explicitamente — a conversão para a convenção do Castro ainda não
acontece neste pipeline (ver nota de estado em §1.10 e a lacuna em §7).

**Escalas derivadas e o custo da simulação.** `t_din`, `v_A`, `t_Alfvén` e a
razão `t_Alfvén/t_din` de §1.10, calculados a partir do `⟨B⟩` e `ρ̄` do
equilíbrio carregado. Um badge de sucesso aparece quando a razão está entre
0,3 e 3 — o regime de campo forte que torna a simulação barata (D4 do
plano).

**Checagem de resolução do toro.** Dado o número de células por lado da
caixa do Castro (128/256/384) e o tamanho da caixa (múltiplo de `R_eq`),
calcula quantas células atravessam a espessura radial do toro
(`toroidal.torus_radial_extent()`, medida no equador). Menos de 10 células
gera aviso — o toro é pequeno comparado à estrela e é ele que precisa ser
resolvido (plano, seção 6).

**Parâmetros da caixa.** `castro.small_dens` é escolhido como
`ρc × 10⁻⁸` — esta é uma **convenção de projeto deste dashboard**, não um
valor derivado da física nem citado no plano; não há justificativa teórica
documentada para o expoente `10⁻⁸` além de ser um piso de densidade
suficientemente baixo para não perturbar a estrela e suficientemente alto
para não gerar `v_A` absurdo no vácuo numérico (o "problema do fluff",
plano seção 5). O raio de início do sponge (`1,5 R_eq`), a duração do
amortecimento (`5 t_din`) e a amplitude da perturbação (`10⁻⁴ c_s`, usando
`eos.sound_speed()` no centro) vêm diretamente dos valores citados no plano
(seção 5), não são recalculados/otimizados por este dashboard.

**R5 — exportação bloqueada se `VE ≥ 10⁻³`**, sem opção de forçar. Implementado
como um `if`/`else` simples em torno do botão de exportação — não há
mecanismo de override em nenhuma camada.

---

## 5. Aba 4 — Registro

→ `dashboard/pages/4_corridas.py`, `dashboard/store.py`

Sem física — este módulo só persiste o que as outras abas já calcularam
(regra R1/R3 do dashboard).

**O que é gravado**, por corrida, em `dashboard/runs/<hash>/`:

| Arquivo | Conteúdo | Função |
|---|---|---|
| `params.json` | parâmetros de entrada completos | `store.save_run()` |
| `scalars.json` | escalares derivados (mesmos da tabela da Aba 1) | idem |
| `fields.npz` | `rho, Phi, u, H, Bphi, r, theta` na malha | idem |
| `manifest.json` | hash, timestamp, git, dependências | idem |

`dashboard/runs/index.csv` agrega hash + parâmetros + escalares de todas as
corridas, para a Aba 4 carregar rápido sem abrir cada diretório
(`store.load_index()`).

**Por que hash de git importa.** `manifest.json` grava
`git_commit_hash(REPO_ROOT)` — o commit de `scf/` **e** `dashboard/` (o
mesmo repositório git cobre os dois, inicializado especificamente para dar
proveniência a este dashboard; `amrex/`, `castro/`, `microphysics/` são
excluídos via `.gitignore`, são repositórios próprios). Um escalar sem o
commit associado não pode ser reproduzido com confiança — se o código mudar
(por exemplo, o bug de §1.4 sendo corrigido), resultados salvos antes e
depois **não são comparáveis**, mesmo com os mesmos parâmetros de entrada.
`git_dirty()` também é gravado — sinaliza se havia mudanças não commitadas
no momento da corrida, o que torna a reprodução exata impossível mesmo
sabendo o hash.

→ `dashboard/store.py :: git_commit_hash()`, `git_dirty()`,
`dependency_versions()` (versões de `numpy`, `scipy`, `streamlit`, `plotly`,
`h5py`, `python`).

**Cache.** O hash é `sha256(json(params, sort_keys=True))[:12]`
(`store.params_hash()`) — determinístico nos parâmetros, usado tanto para
nomear o diretório da corrida quanto para a Aba 2 pular pontos já
calculados.

**Versão de esquema:** **não existe** um campo `schema_version` (ou
equivalente) em `manifest.json` ou `scalars.json` atualmente. Quando o
conjunto de escalares mudou durante este projeto (por exemplo, ao adicionar
`B_pol,max (G)` ao schema da Aba 2), corridas salvas antes da mudança
ficaram com colunas ausentes em `index.csv`, quebrando gráficos que
esperavam a coluna nova — isso já aconteceu de verdade durante o
desenvolvimento e foi contornado apagando o cache antigo manualmente, não
por um mecanismo de versionamento. Fica registrado como lacuna em §7.

**Funcionalidades:** tabela filtrável/ordenável (`st.data_editor`,
formatada por coluna conforme a regra R4 — gauss em notação científica, km
com 2 casas), comparação lado a lado de duas corridas, recarregar uma
corrida na Aba 1 (via `st.session_state["reload_run_params"]` +
`st.switch_page`), marcar como referência (`store.mark_reference()` — hoje
só grava a flag; nenhuma outra aba lê `reference` ainda, ver §7).

---

## 6. Limitações conhecidas

Números medidos, não descrições vagas — todos na configuração
`ρc = 10⁹` g/cm³, `R ≈ 3×10⁸` cm, `nr=161`, `ntheta=65`, `l_max=16` salvo
onde indicado:

- **O critério `VE < 10⁻³` falha acima de `k₀ ≈ 2×10⁻¹²`** nesta
  configuração. O último ponto válido é `k₀ ≈ 1,6×10⁻¹²`, `M ≈ 1,50 M☉`,
  `VE ≈ 6,6×10⁻⁴`. Em `k₀ ≈ 2,3×10⁻¹²` (`M ≈ 2,02 M☉`), `VE ≈ 1,57×10⁻³` —
  acima do limite.
- **Nesse mesmo ponto (`k₀ ≈ 2,3×10⁻¹²`) o pico de densidade migra para fora
  do centro** — de `r_idx=1` (essencialmente a origem) para `r_idx≈25` (raio
  real, não artefato de grade) — e `R_pol/R_eq` cai a `~0,61`: evacuação
  equatorial genuína. A âncora em `ρc(r=0)` (§1.5) perde sentido físico
  exatamente aí, porque o centro deixou de ser o ponto de densidade máxima.
- **Está em aberto se a falha do `VE` acima é resolução insuficiente ou
  terminação de sequência genuína — mas há evidência forte para a segunda
  hipótese**: um estudo de convergência feito neste projeto (`l_max` de 16
  a 48, malha de `129²` a `385²`, quase 9× mais pontos) manteve `VE` em
  `1,2`–`1,6×10⁻³`, **sem tendência de queda** — um platô acima do limite,
  não uma curva convergindo para zero. Isso é a assinatura esperada de uma
  terminação física, não de sub-resolução.
- **Todos os resultados de sequência (`k₀` variável) documentados acima são
  de uma única fatia `ρc = 10⁹`**, onde `M(k₀=0) = 1,39 M☉` — abaixo do
  próprio limite de Chandrasekhar sem campo (`1,44 M☉`). Comparações com
  máximos de massa da literatura (`M_max ~ 1,9 M☉`, Bera & Bhattacharya 2014)
  exigem varredura no plano `(ρc, k₀)` inteiro (Aba 2), não uma fatia — ver
  §3.
- **A conversão `B′ = B/√(4π)` do Castro não é aplicada em nenhum lugar do
  pipeline de exportação atual** (§1.10, §4) — o HDF5 exportado guarda `B`
  em gauss puro, documentado via atributo. Fica responsabilidade do
  `problem_initialize.H` do Castro (não escrito ainda).
- **A Fase 0 do plano (compilar o Castro com `USE_MHD=TRUE`) está pendente**
  — dependências de sistema (`gfortran`, `libhdf5-openmpi-dev`,
  `libopenmpi-dev`) ainda não instaladas neste ambiente. Nada na Aba 3 foi
  testado contra um Castro real.

---

## 7. Questões em aberto

Lacunas identificadas ao escrever este documento — não preenchidas por
conta própria (regra G1):

1. **Teste de consistência de fluxo** (`u(ϖ,z) = ∫₀^ϖ B_z ϖ′ dϖ′`, §1.4) foi
   proposto durante a depuração do bug da função de Green mas nunca
   implementado como teste separado. O teste de Ampère sozinho bastou para
   achar e confirmar o bug; o de fluxo ficaria como segunda linha de defesa
   independente.
2. **Identidade do virial magnético** (`∫ρ∇M(u)·r dV = ∫B²/8π dV`, §1.7) não
   tem função dedicada em `diagnostics.py` nem teste comitado — só foi
   verificada numericamente de forma ad hoc numa sessão de depuração.
3. **A fórmula `J_φ = cρϖf(u)`** (§1.4) não corresponde a nenhuma função —
   nunca é calculada como quantidade nomeada no código.
4. **`B_unit = Rρ√(8πG)`** (§1.10, unidade natural de campo) não está
   implementada — foi usada só como estimativa de ordem de grandeza durante
   a depuração.
5. **Sem `schema_version`** em `manifest.json`/`scalars.json` (§5) — mudanças
   no conjunto de escalares já quebraram `index.csv` de corridas antigas
   uma vez durante este projeto.
6. **`store.mark_reference()`** grava a flag `reference` no índice, mas
   nenhuma outra aba (em particular a Aba 2, que segundo o prompt original
   deveria mostrar corridas de referência nos gráficos) lê essa flag ainda.
7. **`castro.small_dens = ρc × 10⁻⁸`** (§4) é uma convenção sem
   justificativa teórica documentada — funciona como ponto de partida, não
   como valor calibrado.
8. **A questão original que motivou toda a investigação do bug de §1.4**
   — se a razão `(E_mag/|W|)/(M_u/H_c) ≈ 0,5` (regime linear) tem uma
   derivação fechada, ou é só o que se observa numericamente nesta EOS e
   nesta faixa de parâmetros — não foi resolvida analiticamente, só
   confirmada como autoconsistente (constante sob `k₀ → 2k₀`).
9. **`gauss_to_castro()`/`castro_to_gauss()`** existem e estão corretas
   (conferido neste ciclo), mas não são chamadas por nenhum código de
   exportação ainda — ver limitação em §6.

---

## 8. Referências

**Método SCF e equilíbrios magnetizados**
- Hachisu, I. 1986, ApJS 61, 479 — o método SCF
- Tomimura, Y. & Eriguchi, Y. 2005, MNRAS — twisted torus, referência canônica
- Lander, S. K. & Jones, D. I. 2009, MNRAS — campos mistos, funções livres
- Lander, S. K. & Jones, D. I. 2012 — estabilidade de campos mistos

**Anãs brancas magnetizadas**
- Das, U. & Mukhopadhyay, B. 2014, MNRAS 445, 3951 — SCF para AB magnetizada
- Bera, P. & Bhattacharya, D. 2014, MNRAS — M-R autoconsistente com Lorentz
- Bera, P. & Bhattacharya, D. 2016, MNRAS 456, 3375 — geometria de campo
- Bera, P. & Bhattacharya, D. 2017, MNRAS 465, 4026 — estudo de perturbação
- Nityananda, R. & Konar, S. 2014 — crítica aos modelos super-Chandrasekhar
- Coelho, J. G. et al. 2014 — limites do virial

**Estabilidade**
- Markey, P. & Tayler, R. J. 1973 — a instabilidade m = 1
- Braithwaite, J. & Spruit, H. C. 2004 — relaxação para configuração estável
- Braithwaite, J. & Nordlund, Å. 2006

**Códigos**
- Castro: https://github.com/AMReX-Astro/Castro
- Documentação MHD: https://amrex-astro.github.io/Castro/docs/mhd.html
- XNS: https://www.arcetri.inaf.it/science/ahead/XNS/
