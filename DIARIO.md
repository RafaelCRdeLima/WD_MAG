# Diário da campanha — anã branca ultramassiva magnetizada em MHD 3D

Registro do que foi feito, do que foi encontrado, do que deu errado e por quê.
Ordenado por assunto, não por data: a cronologia está nas mensagens de commit.

Mantido em português por ser documento interno. Os relatórios em `reports/` são
em inglês.

**Convenção (8 de agosto de 2026):** achados novos, bibliografia e decisões
entram aqui por padrão, sem precisar ser pedido. O valor deste registro está
tanto nas reversões quanto nos números finais — o afinamento da rotação
diferencial foi chamado de artefato e depois confirmado; o toro torcido foi
dado como construível e depois não. Previsões registradas antes e os erros que
as derrubaram fazem parte do que se escreve.

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

**Alocação do projeto proj503:** 520.134 unidades concedidas, 23.386 usadas até
8 de agosto — 4.5%. Compute não é a restrição desta campanha; vagas de
submissão são. Não há motivo para economizar resolução por medo de gastar.

**O recurso escasso é vaga na fila, não núcleo-hora.** Com duas campanhas ativas
(WD e Brana) chega-se facilmente a seis ou sete jobs entre rodando e
enfileirados, e a corrente precisa de **uma vaga livre no instante em que cada
janela fecha**. Ela não avisa quando não consegue — só para. Foi assim que o
`wdrot256` morreu às 21:10 de 6 de agosto. Antes de submeter uma segunda
corrente, vale conferir quantos jobs já estão na fila.

**Disco:** 657 GB no `dir_rot256`, dos quais 54 checkpoints. Podados para dois,
liberando 198 GB.

---

## 6. Literatura

### 6.1 Os quatro artigos em `references/`

Todos com Mukhopadhyay como autor. São um programa em quatro estágios:
equilíbrios GRMHD com rotação (Subramanian & Mukhopadhyay 2015, MNRAS 454,
752), observabilidade (Gupta, Mukhopadhyay & Tout 2020, MNRAS 496, 894),
anisotropia de matéria e orientação de campo (Deb, Mukhopadhyay & Weber 2022,
ApJ 926, 1), e canal de formação por evolução estelar 1D (Zuraiq et al. 2024,
ASSP).

Os quatro calculam **estrutura**. Nenhum evolui o campo em 3D, e o de 1D não
poderia — num modelo esférico unidimensional um campo toroidal entra como termo
de pressão e não tem como ficar instável.

### 6.2 A busca de agosto de 2026: o que existe e o que não existe

**Não existe MHD 3D de uma anã branca isolada, magnetizada, em rotação
diferencial, testada dinamicamente.** A lacuna que esta campanha ataca é real.
O que existe se divide em duas literaturas que não se tocam: equilíbrios (6.1) e
remanescentes de fusão, que fazem MHD 3D mas de outro objeto — núcleo quente
mais disco espesso, não uma configuração super-Chandrasekhar em equilíbrio.

**Ji & Fisher 2013** (arXiv:1302.5700), primeiras simulações multidimensionais
de fusão de anãs com campo. O disco é fortemente instável à MRI, o campo cresce
rápido até ~2×10⁸ G, e a MRI **freia** o remanescente. Quando a MRI é
resolvida, ela apaga rotação diferencial.

**Pakmor & Pelisoli 2024** (arXiv:2407.02566), fusão de anãs de hélio em alta
resolução, ~50 rotações. Duas fases de amplificação: dínamo de pequena escala
primeiro, depois **dínamo de grande escala dirigido pela MRI** produzindo campo
azimutal ordenado ao longo de dezenas de períodos rotacionais.

**Becerra, Rueda, Lorén-Aguilar & García-Berro 2018** (ApJ 857), o mais próximo
do nosso objeto — anãs super-Chandrasekhar magnetizadas pós-fusão. Mas é modelo
1D com torques, campo fixo entre 10⁶ e 10⁹ G, sem MHD 3D.

**Braithwaite & Spruit** (astro-ph/0510316) e trabalhos relacionados:
equilíbrios estáveis exigem intensidades poloidal e toroidal **comparáveis**;
puramente toroidal e puramente poloidal são ambos instáveis. A configuração
estável é o toro torcido — dipolo polar estabilizado por toroidal semelhante,
como Prendergast previu em 1956. Isso valida o desenho da campanha TT.

### 6.3 O que a busca faz com os nossos resultados

**A sobrevivência da rotação diferencial não se sustenta como resultado.** A
literatura mostra a MRI apagando rotação diferencial quando resolvida. Nós não
a resolvemos — λ_MRI = 0.004 de uma célula — e vemos a rotação sobreviver. Isso
é indistinguível de "nossa simulação não tem como destruí-la", e precisa ser
dito na primeira página dos relatórios, não nos caveats.

**O crescimento do campo ficou ambíguo.** Antes da busca eu o via como imune à
objeção do Rm, por ser efeito Ω puro. Não é: existe um dínamo de grande escala
conhecido neste contexto, e ele é dirigido pela MRI. Ou vemos efeito Ω, que é
mecanismo diferente e mais fraco, ou um análogo numérico de algo que precisa de
MRI para existir.

**Estender o 192³ perdeu o valor que eu havia atribuído.** Convergência entre
malhas não distingue efeito Ω de MRI mal resolvida: as duas malhas têm o mesmo
problema. O teste que distingue é rodar **sem rotação diferencial**. Se o
crescimento sumir, é cisalhamento; se persistir, é numérico.

### 6.4 A MRI é alcançável — o cálculo

λ_MRI = 2π v_Az/Ω depende do campo **vertical**, então ela é ajustável por
construção, não só por refinamento.

| B_z (G) | λ_MRI (cm) | λ/dx no 256³ | λ/R |
|---|---|---|---|
| 3.6×10⁹ (atual) | 2.9×10⁴ | 0.004 | 4×10⁻⁵ |
| 5×10¹² | 4.0×10⁷ | 5.7 | 0.13 |
| **10¹³** | **8.0×10⁷** | **11.4** | **0.27** |
| 2×10¹³ | 1.6×10⁸ | 22.8 | 0.53 |

Existe uma **janela**: λ = 10 dx exige B_z > 8.8×10¹² G no 256³, e a MRI é
suprimida quando λ > R, isto é B_z > 2.4×10¹⁴ G. Fator 27 entre os dois
limites.

E em B_z = 10¹³ G o campo total ao lado do toroidal de 3.2×10¹³ dá 0.76 B_c —
dentro da faixa de validade da EOS. No núcleo β = 792, folgado.

**O obstáculo não é resolução, é geometria.** β cai abaixo de 1 já em
ρ = 10⁷ g/cm³ com esse campo, então ele tem que ser **confinado** ao interior,
não um dipolo de vácuo. É exatamente o mesmo obstáculo que matou a campanha TT,
e a mesma correção resolve os dois: uma fonte de Grad–Shafranov concentrada
num toro em vez de ∝ ρϖ², recalibrada pelo pico interior em vez do dipolo de
superfície.

## 6.5 Campanha CT — corrente confinada em toro

Nome escolhido para distinguir da TT: lá o poloidal era **imposto** sobre um
equilíbrio pronto; aqui a corrente entra na fonte de Grad–Shafranov e o campo
nasce confinado.

### O que a busca encontrou sobre o esquema

**A técnica não é nova, e é padrão em estrelas de nêutrons.**

**Ciolfi & Rezzolla 2013** (arXiv:1306.2803) — a referência central. Na equação
de Grad–Shafranov há duas funções arbitrárias do fluxo ψ: β(ψ), que fixa o
campo toroidal e a corrente azimutal, e F(ψ), que fixa a fonte de corrente
poloidal. A escolha comum é F constante, e é ela que produz dipolo de vácuo.

O primeiro elemento da prescrição deles é justamente o que a CT precisa:
tornar F(ψ) **não constante**, concentrando as correntes perto do eixo de
simetria. Isso **amplia a região de linhas fechadas**, e eles observam que no
limite dessa ideia o campo poloidal fica **inteiramente confinado à estrela**,
citando Fujisawa et al. 2012.

O segundo elemento é um termo adicional que cancela a reação do toroidal sobre
as linhas poloidais. Para β(ψ) eles usam uma forma com função degrau em ψ̄, o
fluxo na última linha fechada, que confina o toroidal à região fechada.

Note que o problema deles é o **espelho** do nosso: queriam mais toroidal
(construções anteriores travavam abaixo de 10% e eles chegaram a 90%), nós
queremos mais poloidal. A alavanca é a mesma.

**Pili & Bucciantini 2014** (arXiv:1401.4308) resolvem configurações mistas de
toro torcido em regime não perturbativo, **no código XNS** — o mesmo que
Subramanian & Mukhopadhyay usaram para gerar a família de onde nossa
configuração vem. Vale ler o corpo do artigo antes de fixar a forma funcional.

**O que não foi encontrado:** aplicação disso a anãs brancas. A literatura de
equilíbrio de WD super-Chandrasekhar usa sequências ou puramente toroidais ou
puramente poloidais; o toro torcido de Braithwaite não aparece ali.

### O que a busca custou à proposta

Os requisitos de resolução da MRI são mais duros do que eu supus. A prática
estabelecida quer **Q ≳ 15 vertical e ~20 azimutal** para turbulência MRI
convergida; seis células por comprimento de onda é o mínimo apenas para a
**fase linear**.

Nosso Q ≈ 11 em B_z = 10¹³ G resolve a fase linear com folga e fica abaixo do
padrão para a saturada. E λ/R = 0.27 dá só ~4 comprimentos de onda na estrela.

**Portanto a CT não promete "incluir a MRI".** Ela promete *ver a MRI crescer*
— medir a taxa linear contra a predição analítica, que é resultado legítimo e
verificável. Afirmar coisa alguma sobre saturação, transporte ou dínamo MRI
exigiria Q ≳ 20, isto é B_z ≈ 1.8×10¹³ G, e aí λ/R = 0.48: duas ondas na
estrela. **A janela útil é mais estreita que o fator 27 sugeria**, porque
subir B_z melhora Q e piora o número de ondas ao mesmo tempo.

### A mudança técnica

Nosso solver resolve Δ* u = −4πϖ²ρ f(u) − ββ'(u), e hoje é chamado com
f(u) = k₀ constante e sem termo β (o toroidal vem à parte, de `ToroidalSC`).
F constante é exatamente o que gera o dipolo de vácuo.

A CT troca f por uma função de u que concentra a corrente para dentro, por
exemplo anulando-a abaixo de um fluxo limiar. **Isso torna a equação não
linear** — u aparece dos dois lados — e exige iteração: resolver, avaliar
f(u), resolver de novo. A forma da equação no módulo já antecipa isso; o uso
atual é que não exercita.

### Previsão registrada antes de implementar

Espero E_tor/E_pol entre 1 e 10, pico total abaixo de B_c, β > 1 em todo o
interior, e Q entre 8 e 15. Se Q sair abaixo de 6 a MRI nem linear aparece, e
a CT vira só a campanha do toro torcido — que ainda vale, por outro motivo.

---

## 6.6 O que a CT entregou, e o que ela matou

### A prescrição certa, e a errada que tentei primeiro

**Fujisawa, Yoshida & Eriguchi 2012** (arXiv:1204.5830) é a referência que
Ciolfi cita para o limite confinado, e é a que resolve o problema. Eles obtêm
campo central **duas ordens de grandeza** acima do de superfície.

Na formulação deles a densidade de corrente tem duas funções arbitrárias do
fluxo Ψ. A primeira, κ(Ψ), é a parte força-livre e dá diretamente o campo
toroidal; *essa* leva corte em Ψ_max, senão o toroidal vazaria para o vácuo. A
segunda, µ(Ψ), é a corrente toroidal não-força-livre e é ela que controla a
localização. Todos os trabalhos anteriores usavam **µ constante** — que é
exatamente o nosso `f = k₀` e é a razão de o nosso campo ser dipolo de vácuo.
Fujisawa usa **lei de potência com expoente negativo**, µ ∝ (Ψ + ε)^m: com
m < −1 o potencial magnético cresce sem limite conforme Ψ cai rumo ao eixo e as
linhas se concentram; m = 0 recupera o caso constante.

**Minha primeira tentativa fez o oposto da ideia:** um limiar que *desligava* a
corrente onde o fluxo era baixo. Deu B_int/B_ext ≈ 3, contra os ~100
necessários, porque corrente localizada de um sinal só ainda carrega momento de
dipolo e ainda produz dipolo exterior. **Concentrar corrente não é o mesmo que
cortá-la**, e a literatura concentra. Trocada a prescrição, B_int/B_ext vai a
973 e β_min de 0.34 a 11.7 — o confinamento funciona.

### A fronteira, e por que ela mata o toro torcido

Fixando para cada forma o maior campo que a estrela **aguenta** (β_min = 1):

| m | B_pol pico | E_tor/E_pol | \|B\|/B_c | λ/dx₂₅₆ |
|---|---|---|---|---|
| 0.00 | 5.8×10¹² | 37.5 | 0.74 | 4.9 |
| −0.50 | 9.2×10¹² | 43 | 0.74 | 6.4 |
| **−1.00** | **2.8×10¹³** | **52.9** | **0.74** | **8.9** |
| −1.30 | 7.1×10¹³ | 63 | 1.60 | 11.8 |
| −1.80 | 2.3×10¹⁴ | 96 | 5.20 | 20.8 |

De m = −1.3 em diante o campo total cruza B_c, então só as três primeiras
linhas são admissíveis — e nelas **E_tor/E_pol nunca desce de 37**.

**A tensão é estrutural.** O confinamento confina *estreitando*: reduz o dipolo
exterior, que salva β na superfície, ao custo do **volume** que o poloidal
ocupa — e energia poloidal é campo ao quadrado vezes volume. As três exigências
puxam em direções incompatíveis:

- toro torcido quer poloidal com energia comparável, logo espalhado;
- β > 1 na superfície quer poloidal confinado, logo estreito;
- MRI resolvível quer B_z alto, o que empurra o total contra B_c.

Com o toroidal de 3.2×10¹³ G fixo, **β ≥ 1 força domínio toroidal de pelo menos
38:1. O toro torcido de Braithwaite não é construível nesta estrela.**

### O que sobra é outra pergunta, e talvez melhor

A linha m = −1 é admissível e entrega o que nenhuma configuração anterior
entregava: **λ/dx = 8.9 com β = 1 e B = 0.74 B_c**. MRI na fase linear
resolvida, num campo que a estrela suporta.

Mas isso não é "toro torcido estável". É a mesma configuração
toroidal-dominada de sempre, agora com poloidal forte o bastante para a MRI
existir na malha. A pergunta muda de *"uma geometria estável sobrevive?"* para
*"quando a MRI pode crescer, ela apaga a rotação diferencial como a literatura
de fusão diz?"*.

É a objeção que a busca bibliográfica levantou contra o nosso resultado
principal, e agora é decidível: se com λ_MRI resolvido a rotação diferencial
for apagada, a sobrevivência que medimos era artefato de resolução; se não for,
o resultado se fortalece muito.

**Campanha ML — a MRI na fase linear.** O nome diz o que ela pode e o que não
pode afirmar: Q ≈ 9 resolve o crescimento linear e fica abaixo do Q ≳ 15–20 que
a turbulência MRI convergida exige. Nada sobre saturação ou transporte
dirigido por MRI sai deste run.

**Previsão registrada antes de implementar:** a MRI cresce em ~1/Ω ≈ 0.12 s, e a
rotação diferencial é apagada em algumas dezenas de segundos. Se nada crescer,
Q = 8.9 não bastou e a conclusão é sobre a malha, não sobre a estrela.

---

## 6.7 Campanha ML no ar — primeiras 700 iterações

Submetida em 8 de agosto (job 994485, 256³). O run **sobreviveu ao arranque**,
que era a falha esperada, mas roda apertado.

**O confinamento melhorou o problema por um fator 21, não o eliminou.** Campo
exterior de 6.6×10¹¹ G contra os 3.0×10¹² da TT, dando β no ambiente de
2.7×10⁻³ contra 1.3×10⁻⁴. Continua abaixo de 1, e é isso que gera os
`Invalid density`. A diferença é que agora o mecanismo de retry do Castro dá
conta, em vez de o run abortar em t = 0.06 s.

| ρ | β com o campo exterior da ML |
|---|---|
| 10⁹ (núcleo) | 1.8×10⁵ |
| 10⁷ | 84 |
| 10⁶ | 1.8 |
| 10⁵ | 0.039 |
| 2×10⁴ (ambiente) | 0.0027 |

**A taxa de retries é estável**, que é o critério que separa "apertado" de
"deteriorando": 1.40 por passo nos primeiros 293 passos, 1.39 nos 419
seguintes. E ρ_max caiu só 1.3% em t = 0.22 s — a estrela assenta devagar, não
toca.

**Errei a estimativa de custo no `inputs`.** Escrevi "uma ou duas janelas" sem
prever os retries, que derrubam o ritmo para ~60%. A t = 3 s exigiria ~7
janelas. Mas o alvo de 3 s é margem, não requisito: a MRI tem e-folding de
0.16 s, então quatro e-foldings — suficiente para um ajuste exponencial
defensável — chegam em **t ≈ 0.64 s**, duas janelas a partir daqui.

**Critério de deterioração, declarado antes:** se os retries por passo passarem
de 5, o run vai abortar e o próximo modelo sai a 80% da fronteira β = 1,
aceitando Q = 6.4 em vez de 8.0.

### A ML morreu em t = 0.221 s, e o critério que declarei era o errado

`ABORT in this window: too many subcycles`, passo 712, com `DT = 5.18e-05`.

Os retries por passo ficaram em **1.39 do começo ao fim** — o critério que eu
mandei vigiar nunca disparou. O sinal real era o **dt colapsando**: 3.6×10⁻⁴ no
passo 293, 5.2×10⁻⁵ no 712, fator 7 em 400 passos. Os retries estavam
funcionando, só que a um passo cada vez menor, até o subciclo estourar.

**Lição transferível: num run com retry ativo, a taxa de retry pode ser
constante enquanto o run morre. Vigiar `dt`.**

Registro também o que funcionou: o detector de abort por janela, escrito depois
do falso alarme do log cumulativo, imprimiu `ABORT in this window` com o STEP
correto. Sem ele eu teria lido de novo o abort histórico de t = 5.146 s.

### E a ML não é consertável nesta família — a conta

O ambiente de 2×10⁴ g/cm³ aguenta no máximo **3.4×10¹⁰ G** (β = 1). A MRI na
fase linear exige λ/dx > 6. Escalando a amplitude por α a partir do pico de
10¹³ G, as duas exigências viram:

α < 3.4×10⁻³ · (B_int/B_ext)   e   α > 6 / (λ/dx)

o que só é satisfeito se **B_int/B_ext × λ/dx > 1757**.

| m | B_int/B_ext | λ/dx | produto |
|---|---|---|---|
| 0.0 | 7.6 | 8.46 | 64 |
| −1.0 | 38.7 | 3.14 | 122 |
| −1.8 | 336 | 0.90 | 303 |
| −2.5 | 973 | 0.68 | **661** |

O melhor da família é 661. **Falta um fator 2.7, e nenhuma amplitude resolve** —
as duas condições se fecham em direções opostas sobre α. Confinar mais sobe a
razão mas transforma o campo numa agulha, derrubando o λ que vem do B_z médio.

**Conclusão: não existe configuração nesta família em que a MRI seja resolvível
e o ambiente aguente o campo.** A ML morre pela mesma razão estrutural que a
TT, agora medida em vez de suposta.

O que restaria: eliminar o ambiente como limitante — fronteira de vácuo em vez
de atmosfera preenchendo a caixa, ou densidade ambiente muito maior numa caixa
menor. Ambos são mudanças de setup, não de modelo, e nenhum é barato.

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

## 8. Em aberto, revisto depois da busca de literatura

**O teste que decide o crescimento: rodar sem rotação diferencial.** Se E_mag
parar de crescer, o mecanismo é cisalhamento e a afirmação fica limpa. Se
persistir, é numérico. Mais barato que estender o 192³ e responde o que a
comparação entre malhas não responde — porque as duas malhas compartilham o
mesmo defeito.

**Incluir a MRI é possível e não exige refinamento.** λ_MRI ∝ B_z, então basta
um campo poloidal interior de ~10¹³ G para trazê-la a ~11 células no 256³,
dentro da janela 8.8×10¹² < B_z < 2.4×10¹⁴ G e com o campo total em 0.76 B_c.
O obstáculo é o mesmo da campanha TT: o campo precisa ser **confinado**, e o
SCF atual só sabe impor dipolo de vácuo.

**A correção única que destrava as duas coisas:** trocar a fonte de
Grad–Shafranov por uma concentrada num toro e recalibrar pelo pico interior.
Isso dá simultaneamente o toro torcido (energias comparáveis, geometria estável
de Braithwaite) e a MRI resolvida. É a próxima peça de trabalho e é de
modelagem, não de máquina.

**Continua em aberto sem caminho barato:** se a taxa de crescimento é física.
Ela cai 2.5× de 192³ para 256³, igual ao decaimento antes dela. Uma terceira
malha custaria ~5× o 256³ e ainda assim não separaria efeito Ω de MRI mal
resolvida.

**Descartado:** MHD não ideal, por aritmética (seção 3.3). Estender o 192³ além
de 60 s, porque o ganho que eu atribuía a isso não existe (seção 6.3).
