# Diário da campanha — anã branca ultramassiva magnetizada em MHD 3D

Registro do que foi feito, do que foi encontrado, do que deu errado e por quê.
Ordenado por assunto, não por data: a cronologia está nas mensagens de commit.

Mantido em português por ser documento interno. Os relatórios em `reports/` são
em inglês.

---

## 1. O que está sendo simulado

Uma anã branca de **2.005 M⊙** — 38% acima do limite de Chandrasekhar — que só
existe por causa da rotação diferencial, com um campo interior toroidal de
3.2×10¹³ G ao lado de um dipolo exterior de 10⁹ G.

| | |
|---|---|
| ρ_c | 3.0×10⁹ g/cm³ |
| R_eq, R_pol | 3.92×10⁸, 1.50×10⁸ cm (achatamento 0.383) |
| lei de rotação | Ω(ϖ) = Ω_c A²/(A²+ϖ²), A = 1.834×10⁸ cm |
| P na superfície | 4.3 s (três vezes mais rápido que a WD mais rápida observada) |
| T/\|W\| | 0.0993, abaixo dos limiares de barra (0.14 secular, 0.27 dinâmico) |
| E_tor/E_pol | 3.1×10⁵ na malha (2.2×10⁷ no modelo analítico) |

Castro 26.07 / AMReX, MHD ideal com transporte restrito, HLLD, malha cartesiana
única. EOS `ztwd`: elétrons degenerados a temperatura zero, **barotrópica**.

A configuração vem da família de Subramanian & Mukhopadhyay (2015) — mesma lei
de rotação, mesmo regime de massa, campo e período. É uma configuração tipo
remanescente de fusão, não modelo de estrela observada.

---

## 2. Os runs

| run | malha | alcance | estado |
|---|---|---|---|
| `dir_rot192` | 192³ | t = 60 s | completo, campo processado até 60 s |
| `dir_rot256` | 256³ | rumo a 100 s | processado até 64.5 s; corrente ativa |
| `dir_mixed192` | 192³ | — | campanha TT, morreu em t = 0.06 s, removido |

O 256³ foi encadeado em janelas de 3 h, quatorze submissões, ~35–40 h de CPU em
256 núcleos e ~460 GB de saída depois da poda de checkpoints.

---

## 3. Resultados

### 3.1 Convergidos

**A estrela sobrevive.** Massa conservada em 0.02% ao longo de 60 s enquanto a
energia magnética cai três ordens de grandeza. O que destrói o campo não
desliga a estrela.

**A rotação diferencial não é apagada.** Ω_out/Ω_core fica perto de 1/3 nas duas
malhas e em nenhuma caminha para rotação uniforme. Em t = 12 s as duas
concordam em 1.4%.

**O perfil afina com o tempo**, −22.6% contra −20.3% sobre 46 s de base. A
malha fina começa ~6 s depois e depois afina mais rápido.

**A frenagem de L_z cai pela metade em três regimes sucessivos**, e as duas
malhas concordam em todos:

| | 192³ | 256³ |
|---|---|---|
| t = 1.5–11.9 s (campo vivo) | −0.1996 %/s | −0.2144 %/s |
| t = 12–30 s (campo decaindo) | −0.0873 %/s | −0.0935 %/s |
| t = 30–58 s (resíduo) | −0.0446 %/s | −0.0523 %/s |

Sem viscosidade nem resistividade explícitas, nada mais no esquema freia coisa
alguma — o fator dois identifica o transporte como magnético, e agora em duas
resoluções.

**A instabilidade é m=1** e cresce em poucos tempos de Alfvén. Não é a MRI
axissimétrica: construída com o campo vertical, λ_MRI = 2.9×10⁴ cm, ou 0.003 de
uma célula. Sobra kink de Tayler ou MRI não-axissimétrica de campo toroidal, e o
número de modo sozinho não separa as duas.

### 3.2 Não convergidos

**As taxas, e em sentidos opostos.** O crescimento sobe 47% de 192³ para 256³, o
decaimento cai pela metade. A medida grossa exagera a violência nas duas pontas.

**O resíduo do campo.** E_mag estaciona em ~3×10⁴⁶ erg no 192³ e 7×10⁴⁷ no
256³ — fator 25.

**O nível do resíduo.** E_mag no mínimo é 2.5×10⁴⁶ erg no 192³ e 6.3×10⁴⁷ no
256³ — fator 25.

**A taxa de crescimento depois do mínimo**, +0.0485 /s no 192³ contra
+0.0196 /s no 256³. Cai 2.5× com o refinamento, mesma direção e mesmo fator do
decaimento antes dela. Taxa que cai pela metade quando a malha refina é
assinatura de processo dominado pela malha.

### 3.2b O campo para de decair e volta a crescer — nas duas malhas

O primeiro resultado positivo da campanha: algo que a estrela **faz**.

| | 192³ | 256³ |
|---|---|---|
| mínimo de E_mag | t = 41.1 s | t = 44.7 s |
| E_mag ali | 2.54×10⁴⁶ | 6.35×10⁴⁷ |
| crescimento depois | +0.0485 /s (e-fold 21 s) | +0.0196 /s (e-fold 51 s) |
| do mínimo ao fim | +114% | +21% |

Mínimos localizados por parábola em ln E_mag sobre t = 25–55 s, para a pulsação
de 1.5 s não decidir a resposta.

**A existência da virada converge**; nada de quantitativo nela converge. O
mecanismo está disponível e não é exótico: a rotação diferencial sobrevive e
continua afinando, então há energia de cisalhamento livre, e enrolar poloidal
residual em toroidal é o que cisalhamento faz com campo.

**A previsão registrada falhou nas duas direções.** Estava escrito que o 192³
viraria *mais tarde* e *mais fraco*, porque mais dissipação numérica atrasaria
o ponto em que a regeneração supera o decaimento. Ele vira 3.6 s **antes** e
cresce 2.5× **mais rápido**.

O modo do erro é informativo. Se fosse dínamo correndo contra dissipação
numérica, menos dissipação daria mais cedo e mais forte; observa-se o
contrário. Uma leitura compatível é que as duas malhas enrolam rumo a uma
saturação e a que parte 25× mais baixa tem mais caminho — mas a distância só
encolhe de 25× para 16× na janela disponível, então também não estão
visivelmente convergindo para um nível comum.

**Com o 256³ levado a t = 78 s, a subida virou exponencial estável.** Ajustes
de ln E_mag em janelas independentes de 10 s: +0.0354 /s em t = 55–65 e
+0.0349 /s em t = 65–78 — concordância de 1.5%, e-folding de 29 s. Do mínimo
ao fim, E_mag mais que dobra (+124%, contração explica 3.5%) e o campo de pico
vai de 4.05 para 8.37×10¹² G.

**A oscilação ficou implausível.** Trinta e três segundos de subida monótona
exigiriam de uma senoide período acima de 133 s: 280 tempos dinâmicos, 89
pulsações da própria estrela. Sem apontar um modo com esse período, não se
sustenta.

**E as duas malhas estão convergindo.** A distância em E_mag encolhe de forma
monótona — 24.5× no mínimo, 16.1× em t = 51, 14.1× em t = 58, que é onde o
192³ acaba. Nível fixado pela malha manteria a razão; nível para o qual as
duas sobem fecha a razão, e é isso que se observa. Extrapolando as duas taxas,
elas se encontrariam perto de t ≈ 190 s, muito além dos dois runs.

**O que segue não convergido é a taxa.** Na janela que as duas cobrem, a
grossa cresce ~2.4× mais rápido — mesma direção e fator parecido com o
decaimento antes dela. O fenômeno é robusto; a escala de tempo dele não é.

### 3.3 O diagnóstico que enquadra tudo

Extraindo a difusividade magnética numérica do decaimento medido: η_num ≈
0.6–0.9 c_s·dx, dando **Rm ≈ 3 no 192³ e ≈ 6 no 256³**. A resistividade física
da matéria degenerada dá Rm ~ 10¹⁸ e tempo ôhmico de 3×10¹⁰ anos.

Nossa resistividade é dezessete ordens de grandeza maior que a real. Rm ~ 5 não
é MHD quase-ideal — é regime difusivo, sem faixa inercial.

**Consequência:** MHD não ideal é inviável aqui. Para o η explícito dominar o
numérico ele precisa ser algumas vezes maior, o que derruba Rm para ~1–2. Para
ter η controlado *e* Rm ~ 10 seria preciso dx < 10⁶ cm, ou ~2000³. Não existe
janela. Isso é limitação quantificada, não omissão, e substitui a frase
"não implementamos resistividade explícita" nos relatórios.

**A separação que salva metade dos resultados:** a instabilidade de Tayler é
ideal. O crescimento, o caráter m=1 e o transporte de momento angular enquanto o
campo existe são robustos ao Rm baixo. O decaimento, o resíduo e o crescimento
novo são dissipativos e devem ser reportados como limitados por resolução.

---

## 4. Erros cometidos e como foram pegos

Registrados porque o modo de falha costuma ser mais transferível que o
resultado.

**Diagnósticos nativos do Castro corrompidos.** `emag_density`, `etor_density` e
`Div_B` empacotam componentes de faces com index types diferentes num FAB
cell-centred e leem além do fim da caixa; `Div_B` chegava a 3×10¹⁴⁶. Tudo foi
reconstruído de B_x, B_y, B_z. Também produziam `inf` nos plotfiles, o que
tornou 70 arquivos ilegíveis até serem removidos do `derive_plot_vars`.

**Campo lido em Gauss em vez de Heaviside-Lorentz.** Levou a afirmar que a malha
retinha 27% do pico toroidal do modelo. O correto é 96.0% a 96³ e 98.9% a 192³.
A correção também revelou que B excede B_c em ~44% do run.

**MRI afirmada com λ construído da velocidade de Alfvén total.** A MRI
axissimétrica exige o campo *vertical*, e com ele λ é 0.003 de uma célula.
Retratado.

**Semente m1/m0 ≈ 10⁻³ que era erro de interpolação em anéis.** A projeção
direta por célula dá 1.2×10⁻¹⁶ — a condição inicial é axissimétrica à precisão
de máquina.

**A = R_eq assumido na lei de rotação.** Deu Ω equatorial 2.6× errado e uma
gravidade equatorial 8× a de ponto-massa, que deveria ter sido o sinal. Correto:
P_eq = 4.3 s.

**Comparação de convergência sobre janela curta demais, duas vezes.** Sobre 14 s
o afinamento parecia artefato de malha grossa; sobre 28 s parecia real em
direção mas não em amplitude; sobre 46 s converge em 10%. A causa é o atraso de
~6 s no início da malha fina. **Uma janela menor que o atraso entre as malhas
mede o atraso, não a física.** Um efeito atrasado medido assim parece ausente —
e, se a curva atrasada deriva para o outro lado antes, parece ter sinal oposto.

**Previsão do 192³ errada nas duas direções.** Registrei que ele viraria mais
tarde e mais fraco que o 256³; virou antes e mais forte. A previsão ter sido
escrita antes é o que fez o fracasso significar algo — sem ela seria fácil
racionalizar qualquer um dos dois desfechos.

**B_c verificado, β esquecido.** Ao escolher a amplitude poloidal da campanha
TT conferi `max|B|/B_c` e ignorei a razão entre pressão magnética e do gás. Em
B_pole = 3×10¹² o β vale 0.088 no envelope e 1.3×10⁻⁴ no ambiente: o campo
domina, empurra o fluido, a densidade cruza zero em t = 0.06 s. O limite certo
não era o de Landau.

**Calibração de B_pole pelo máximo de |B_r| em vez de `surface_dipolarity`.**
Errou a varredura inteira por fator 7.6 e me levou a recomendar um ponto que
está a 1.90 B_c. A versão corrigida valida a si mesma: em B_pole = 10⁹
reproduz E_tor/E_pol = 2.18×10⁷, o valor documentado do modelo analítico.

**Abort antigo lido como novo.** O log é cumulativo entre janelas e o script
usava `head -1`, então reportava a falha de t = 5.146 s em toda janela desde
então. Custou um ciclo de diagnóstico e quase uma redução de CFL desnecessária.

---

## 5. Infraestrutura: o que quebrou e o padrão que resolveu

Detalhado em `cluster/cenapad/ONBOARDING.md`. Resumo dos incidentes:

**PBS mata o `mpirun` na parede**, então todo o encadeamento colocado depois
dele nunca executava. Dois runs saudáveis pararam assim. Resolvido com parada
graciosa via `dump_and_stop` e uma folga de 15 min.

**Jobs concorrentes no mesmo diretório** se sobrescreveram, e um reiniciou do
checkpoint meio escrito do outro. Resolvido com diretório por run e lock
`RUNNING`. O lock só funciona no mesmo nó — `kill -0` não enxerga outros.

**Alvo hardcoded no script** dessincronizou do arquivo de entrada ao estender o
run, e teria parado a corrente na primeira janela parecendo término limpo. Agora
é lido do `inputs`.

**Reinício exatamente no `stop_time`** faz o código tentar um passo de ~10⁻¹⁵ s e
abortar com mensagem enganosa. Guard adicionado.

**Queda de nó** (`PRTE has lost communication`) parece nada: o log para
no meio de um passo perfeitamente saudável. Custou uma janela e quebrou a
corrente, porque o guard de "nenhum checkpoint novo" disparou. O guard está
certo; o espaçamento estava errado. `check_int` 2000 → 1000.

**Jobs retidos** (`H`, "too many failed attempts to run") três vezes. Provável
limite de submissões por usuário, com as duas campanhas competindo por vagas.

**Ferramentas de análise falhando com `GLIBCXX`** — módulo não carregado na
sessão nova. E varreduras longas morrendo com o terminal por falta de `nohup`.

**Guard de corrente com furo, sete janelas queimadas.** A detecção de reinício
exclui `chk00000`, escrito na inicialização, e o guard de "nenhum checkpoint
novo" não excluía. Na primeira janela ele comparava `""` com `"chk00000"`,
concluía que houve progresso e ressubmetia — indefinidamente, num run que
morria em t = 0.06 s. Corrigido nos dois scripts.

**Fila `paralela` exige dois nós.** Pedir um dá `Job violates queue and/or
server resource limits`. Trabalho de nó único vai para `par128`.

**A corrente precisa de vaga livre a cada três horas** para se ressubmeter. A
do `wdrot256` morreu em silêncio às 21:10 de 6 de agosto quando seu `qsub`
interno foi recusado com quatro jobs de outra campanha na fila.

**Disco:** 657 GB no `dir_rot256`, dos quais 54 checkpoints. Podados para dois,
liberando 198 GB.

---

## 6. Literatura

Quatro artigos em `references/`, todos com Mukhopadhyay como autor. São um
programa em quatro estágios: equilíbrios GRMHD com rotação (2015),
observabilidade (2020), anisotropia de matéria e orientação de campo (2022),
canal de formação por evolução estelar 1D (2024).

**O ponto que importa:** os quatro calculam *estrutura*. Nenhum evolui o campo
em 3D, e o de 1D não poderia — num modelo esférico unidimensional um campo
toroidal entra como termo de pressão e não tem como ficar instável. Estabilidade
nessa literatura significa razão de energias, critério radial ou ponto de
retorno ao longo de uma sequência, nunca teste dinâmico.

É a lacuna que esta campanha ataca. A pergunta não é se dá para construir uma
anã magnetizada de 2 M⊙ — dá, por vários métodos — mas se ela sobrevive quando
se deixa evoluir.

---

## 7. Produtos

- `reports/report_rot192_rot256.pdf` — relatório I, 13 páginas, física primeiro,
  numérica ao fim.
- `reports/report_late_convergence.pdf` — relatório II, 6 páginas, os três
  testes sobre 46 s de base.
- `investigations/bt_bp_256_long.csv` — 178 linhas, t = 0 a 64.5 s, 18 colunas.
- `investigations/bt_bp_192_late.csv` — campo do 192³ de t = 12 a 60 s.
- `investigations/plot_late_convergence.py` — a figura dos três testes.
- `tools/fbtbp.cpp`, `fslice.cpp`, `fmodes.cpp` — diagnósticos reconstruídos.
- `cluster/cenapad/ONBOARDING.md` — o que a infraestrutura custou aprender.

---

## 8. Em aberto

**Crescimento ou oscilação?** Precisa de tempo físico, não de resolução nem de
análise adicional. Se E_mag virar para baixo até t ≈ 90 s é oscilação; se seguir
acelerando, é regeneração. O 256³ está rodando rumo a 100 s.

**A taxa de crescimento é física?** Ela cai 2.5× de 192³ para 256³, igual ao
decaimento antes dela. Ou é dominada pela malha, ou depende do nível de
partida, que difere por 25×. Duas malhas não separam as duas leituras; uma
terceira seria necessária, e 384³ custa ~5× o 256³.

### A campanha TT falhou na primeira tentativa, e por quê

Primeira submissão em 7 de agosto: a estrela explode em **t = 0.06 s**, com
densidade negativa (−5.3×10⁵ g/cm³) em células a ϖ ≈ z ≈ 1.6×10⁸ cm — o raio
polar, onde a densidade despenca e o dipolo imposto é mais forte.

A causa não é numérica. Verifiquei `max|B|/B_c` ao escolher B_pole = 3×10¹² e
ignorei a razão β entre pressão magnética e do gás:

| B (G) | β no envelope (ρ=10⁶) | β no ambiente (ρ=2×10⁴) |
|---|---|---|
| 10⁹ (config. atual) | 7.9×10⁵ | 1.2×10³ |
| 10¹² | 0.79 | 1.2×10⁻³ |
| 3×10¹² | **0.088** | **1.3×10⁻⁴** |

O campo domina a pressão do gás por uma ordem de grandeza dentro da estrela e
por quatro fora. Ele empurra o fluido, a densidade cruza zero, o Castro aborta.

**A raiz é geométrica.** O campo toroidal é confinado por construção
(`B_φ ∝ ρ`, some onde a estrela acaba). O poloidal imposto é um dipolo de
vácuo que se estende para fora indefinidamente — multiplicar sua amplitude por
3000 põe 10¹² G num ambiente de 2×10⁴ g/cm³. Os toros torcidos de Braithwaite
não são assim: neles o poloidal fecha *dentro* da estrela e só a parcela que
atravessa a superfície vira dipolo exterior, fraco.

O limite é severo: β = 1 no ambiente exige B ≲ 3.4×10¹⁰ G, o que deixa
E_tor/E_pol na casa de 10⁴. **Por este caminho o toro torcido é inalcançável**
— não por B_c, não pelo virial, mas porque um dipolo exterior não tem o que o
segure.

Custo do erro: sete janelas desperdiçadas, porque o guard de "nenhum checkpoint
novo" tinha um furo. A detecção de reinício exclui `chk00000` (escrito na
inicialização) mas o guard não excluía, então na primeira janela ele comparava
`""` com `"chk00000"`, concluía que houve progresso e ressubmetia — para
sempre, num run que morre logo após inicializar. Corrigido nos dois scripts.

**O caminho que resta é o de Braithwaite:** não construir o toro torcido
analiticamente, e sim partir de um campo qualquer e deixar a estrela relaxar
para o dele. `castro_problems/wd_braithwaite` já existe com essa
infraestrutura, e o rascunho `papers/wd-braithwaite-relaxation` é sobre isso.
Tem a vantagem de não exigir equilíbrio inicial nenhum — a relaxação é o
método, não uma concessão.

**Campanha TT — o toro torcido.** O run que falhou. Configuração mista com energias
poloidal e toroidal comparáveis, com rotação, a 192³. É pergunta de estabilidade,
portanto ideal, portanto robusta ao Rm ~ 6 que não temos como consertar — a única
coisa que não dá para consertar. Gerador em
`investigations/export_tt_model.py`.

O que se abre mão: com o SCF atual o campo poloidal é imposto sobre um
equilíbrio já convergido, não resolvido junto, então o par sai do balanço virial
pela própria energia poloidal. A B_t/B_p = 8.8 o erro virial é 2.3×10⁻², vinte
vezes o limiar. O modelo **não é equilíbrio** e não pode ser apresentado como
tal — é condição inicial a relaxar, que é o que o rascunho
`papers/wd-toroidal-poloidal` prescreve para essa faixa. O transiente inicial
passa a ser parte do experimento.
