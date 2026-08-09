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

## 6.8 O que fica depois da ML — a limitação que vai para a primeira página

Duas campanhas mortas na mesma parede, e a parede agora está medida.

**A conta que fecha as duas.** O ambiente de 2×10⁴ g/cm³ suporta 3.4×10¹⁰ G
antes de o campo dominar a pressão do gás. A MRI na fase linear exige
λ/dx > 6. Escalando a amplitude, as duas exigências só coexistem se
B_int/B_ext × λ/dx > 1757, e o máximo da família é 661. Confinar mais sobe a
razão e transforma o campo numa agulha, derrubando o B_z de volume que dá o
comprimento de onda. Não é ajuste malfeito; é a geometria do problema.

**A consequência para os relatórios, e ela é desconfortável.** A busca de
literatura (seção 6.2) mostrou que a MRI apaga rotação diferencial em
remanescentes de fusão quando é resolvida — Ji & Fisher 2013, Pakmor &
Pelisoli 2024. Nós não a resolvemos, e agora sabemos que **não conseguimos
resolvê-la neste setup a nenhum custo acessível**.

Portanto: *"a rotação diferencial sobrevive"* continua indistinguível de
*"nossa simulação não tem como destruí-la"*, e essa ambiguidade deixou de ser
provisória. Ela pertence à primeira página dos dois relatórios como limitação
declarada, não ao fim como caveat. Hoje os relatórios apresentam a
sobrevivência da rotação como resultado; depois da busca isso não se sustenta
sem a ressalva.

**O que continua de pé, e é bastante.** Nada disso toca a frenagem de L_z, que
converge em três regimes nas duas malhas e identifica o transporte como
magnético. Nem o afinamento do perfil, convergido em 10% sobre 46 s. Nem a
virada do campo em t ≈ 41–45 s nas duas malhas. O que cai é a interpretação da
rotação diferencial como *resultado sobre a estrela* em vez de *comportamento
deste setup*.

**O que restaria tecnicamente**, e nenhum é barato: fronteira de vácuo em vez de
atmosfera preenchendo a caixa, ou ambiente muito mais denso numa caixa menor.
Ambos são mudanças de setup, não de modelo, e ambos invalidariam a
comparabilidade com os dois runs já feitos.

### Infraestrutura consertada no mesmo dia

**Retry no `qsub` da corrente**, nos três scripts. Três correntes morreram
porque o script chamava `qsub` uma vez, a fila estava cheia com a outra
campanha naquele segundo, e a submissão era recusada — com a linha
`chain: submission N of M` já impressa, então o `.out` parecia entrega saudável
com o ID faltando. Custava um run *e* parecia sucesso. Seis tentativas de dez
em dez minutos cobrem uma janela inteira da campanha concorrente, e cada
tentativa fica registrada.

**`ONBOARDING.md` ganhou as duas assinaturas de falha do dia:** o retry no
`qsub`, e a de que taxa de retry constante pode esconder um run morrendo —
vigiar `dt`.

---

## 6.9 Refinar não resolve, mas a limitação é menor e mais comum do que eu disse

Duas correções à seção 6.8, ambas para melhor.

### Refinar é impossível, não caro

λ_MRI construído com o campo **inicial** contra o tamanho de célula:

| malha | λ/dx | custo relativo |
|---|---|---|
| 256³ (atual) | 0.004 | 1× |
| 1024³ | 0.017 | 256× |
| 4096³ | 0.066 | 65.536× |

Para λ/dx = 6 seriam **372.000 células por lado**, 5×10¹⁶ células, 4.5×10¹²
vezes o custo — cerca de 2×10¹⁰ anos de máquina. O custo escala como N⁴
(N³ células vezes N passos pelo CFL). Não é questão de esperar; é impossível
por doze ordens de grandeza.

### Mas a MRI passa a ser resolvida durante o próprio run

Eu vinha afirmando que a MRI nunca é resolvida. **Está errado.** λ_MRI escala
com o campo poloidal, e a instabilidade de Tayler fabrica campo poloidal.
Usando o max|B_pol| medido:

| t (s) | λ/dx no 192³ | λ/dx no 256³ |
|---|---|---|
| 0 | 0.08 | 0.10 |
| 1.0 | 2.2 | 3.0 |
| 2.0 | 3.5 | 13.6 |
| 3.0 | 16.2 | 36.4 |
| 5.2 | 32.7 | 35.9 |
| 12.0 | 10.3 | 36.5 |

**A partir de t ≈ 2–3 s a MRI está resolvida e assim permanece.**

A afirmação correta é: *não a resolvemos durante os primeiros ~2 s, e a
resolvemos depois.* O que fica de fora é exatamente a janela em que o campo é
ordenado e forte — que é quando a literatura de fusão a vê agir.

Duas ressalvas: usei o **pico** de B_pol, limite superior; o valor típico de
volume é menor por talvez 3 a 10, o que empurraria o cruzamento para t ≈ 3–4 s
e deixaria o regime tardio marginal em vez de folgado. E λ/R cresce junto — em
λ/dx = 36 o comprimento de onda é meia estrela, e aí a MRI global perde
sentido. **Medir o B_z típico de volume com as fatias em disco é tarefa
pendente e não precisa de cluster.**

### E a limitação é a padrão do campo, não um defeito nosso

Busca de 8 de agosto. [Subgrid modelling of MRI-driven turbulence in
differentially rotating neutron stars](https://arxiv.org/abs/2509.07081) abre
reconhecendo que **a maioria das simulações de fusão de estrelas de nêutrons
não resolve a MRI**, pelo mesmo motivo: comprimento de onda pequeno demais, e
resolução proibitiva para cobrir todas as escalas.

As três saídas da comunidade, e nenhuma é refinar:

- **modelos de sub-grade**, representando a turbulência não resolvida em termos
  das quantidades de larga escala;
- **large-eddy simulations**, já usadas em fusões binárias;
- **caixas de cisalhamento locais**, que resolvem a MRI num pedaço e servem
  para **calibrar** os coeficientes de transporte dos modelos globais.

Nossa densidade de 10⁹ g/cm³ torna o problema pior que num disco de fusão, onde
ρ é ordens de grandeza menor e λ_MRI proporcionalmente maior.

E o Pakmor & Pelisoli 2024 descreve a mesma estrutura que medimos: dínamo de
pequena escala primeiro, **depois** dínamo de grande escala dirigido pela MRI.
Nem eles resolvem a MRI a partir do campo semente — ela entra depois que outro
mecanismo amplificou o campo.

### A formulação que vai para os relatórios

Não *"não resolvemos a MRI"*, que soa como incompetência. E sim:

> A MRI não é resolvida durante a fase de campo ordenado, como em praticamente
> toda simulação global desta classe; a literatura contorna isso com modelos de
> sub-grade calibrados em caixas de cisalhamento, que não implementamos.

Preciso, situa o trabalho no estado da arte, e aponta o caminho em vez de só
admitir o buraco. **Um modelo de sub-grade é a rota concreta** e não exige
refinar nada — exige um termo de transporte turbulento no esquema com
coeficiente vindo da literatura. Mudança no Castro, não trivial, mas de outra
ordem de dificuldade que 372.000³.

---

## 6.10 Sub-grade e LES: qual das duas famílias serve, e por quê

Estudo de 8 de agosto, sem implementar nada.

### MInIT — viável, com coeficientes emprestados

[Miravet-Tenés et al. 2025](https://arxiv.org/abs/2509.07081). Três peças: duas
densidades de energia turbulenta evoluídas (MRI e instabilidades parasíticas),
cada uma com `∂_t e + ∇·(v e) = S`; tensores de Maxwell e Reynolds ligados a
elas por coeficientes constantes; e o transporte vindo essencialmente da
componente ϖφ, ou seja **um torque**, não um tensor viscoso completo.

A favor:

- **O gancho existe e já usamos.** `problem_source.H` tem 79 linhas e implementa
  o amortecimento da acomodação. Um torque entra na mesma estrutura.
- **A lei de rotação do artigo é a nossa** — Ω_c/(1 + ϖ²/A²), Komatsu, a
  j-constante. Os coeficientes foram derivados para esse perfil.
- **O modelo já contempla malha grossa:** o comprimento de decaimento é
  λ = min(Δ, λ_MRI), com Δ o tamanho de célula. Foi construído para o caso mal
  resolvido.
- Escalares advectados são nativos no AMReX.

Contra:

- **O Castro removeu difusão de velocidade** (só sobrou térmica), então não há
  infraestrutura de divergente de tensor. Contornável porque precisamos de uma
  componente só, escrita como fonte.
- **Os coeficientes parasíticos (−1.4, −0.8) são calibrados para estrela de
  nêutrons**, politrópica γ = 2 em densidade nuclear. Os da MRI são teóricos e
  transferem; os parasíticos, transferir para `ztwd` é **suposição**. O
  resultado dependeria de coeficientes de outro regime — defensável se
  declarado, mas transforma "medimos transporte por MRI" em "aplicamos modelo
  calibrado alhures".

Custo: semanas. Duas variáveis de estado, fontes locais, um torque, e um teste
de verificação reproduzindo a evolução das energias turbulentas deles.

### LES — não serve, e o motivo é instrutivo

[Viganò, Aguilera-Miret et al. 2020](https://arxiv.org/abs/2004.00870) estendem
para GRMHD o modelo **gradiente**: expansão de Taylor dos termos não lineares
nos fluxos, "fisicamente agnóstico". Vantagem real sobre o MInIT — **sem
coeficientes calibrados noutro regime**.

Mas a revisão de [Schmidt-Brückner 2025](https://arxiv.org/abs/2509.06801) diz
que modelos de sub-grade têm impacto prático limitado na maioria dos códigos
astrofísicos **porque os solvers já têm difusão numérica significativa**. É a
nossa situação em grau extremo: η_num ≈ 0.6–0.9 c_s·dx, Rm ≈ 5.

**LES estima o que ocorre abaixo da escala de filtro, pressupondo cascata
turbulenta continuando abaixo da célula. A Rm ≈ 5 não há cascata** — a
dissipação acontece *na* célula. Não há nada abaixo para modelar, e o termo de
sub-grade somaria a uma dissipação já excessiva.

### A distinção que eu não tinha feito

- **LES modela uma cascata**: energia que as escalas resolvidas já têm,
  transferida para baixo. Pressupõe faixa inercial.
- **MInIT modela uma instabilidade**: energia que a malha não consegue criar,
  porque o modo não cabe. Não pressupõe cascata; injeta o que falta.

Nossa física ausente é a segunda. A MRI é um modo que não cabe na malha, não uma
cascata mal amostrada. **O MInIT é a classe certa; o LES não é.**

Confirma isso onde o LES funciona: [Aguilera-Miret et al.
2020](https://arxiv.org/abs/2009.06669) capturam amplificação nos primeiros
10 ms após fusão, onde há turbulência violenta de Kelvin–Helmholtz nas escalas
resolvidas. Nossa estrela, assentada e pulsando suavemente, não está nesse
regime.

---

## 6.11 Caixas de cisalhamento — RETRATADO: não é barata, e a conta que faltava

> **AVISO, escrito no mesmo dia.** A conclusão original desta seção — que a
> caixa era "a rota mais barata" — **está errada**. Calculei o requisito
> ESPACIAL (tamanho da caixa, células, fator de qualidade), vi 128³ sem
> autogravidade e concluí "barato" sem calcular o número de PASSOS DE TEMPO.
> Custo é células × passos, fórmula que eu havia usado duas seções antes para
> mostrar que refinar escala como N⁴, e esqueci aqui. A correção está em
> 6.11b — que **também saiu errada**, para o lado oposto. O veredito final
> está em 6.11c. Leia as três na ordem; o percurso é o registro.
>
> Segundo erro do mesmo tipo em dois dias: também escrevi "uma ou duas janelas"
> para a ML sem contar os retries. **Estimar custo por uma dimensão só.**

Estudo de 8 de agosto, sem implementar nada.

### Os números espaciais são favoráveis — e só eles

**A aproximação local é satisfeita com folga de quatro ordens de grandeza.**
λ_MRI = 2.5×10⁴ cm contra R_eq = 3.9×10⁸ cm, razão 6×10⁻⁵. Caixa local exige
λ ≪ R; nosso caso é muito mais confortável que num disco fino, onde λ pode ser
fração apreciável da altura de escala.

**O cisalhamento está na faixa padrão.** A lei j-constante dá
q = −dlnΩ/dlnϖ = 2ϖ²/(A²+ϖ²):

| ϖ/R_eq | q |
|---|---|
| 0.25 | 0.44 |
| 0.468 (= A) | 1.00 |
| 0.80 | 1.49 |
| 1.00 | 1.64 |

A MRI opera em 0 < q < 2, kepleriano é 1.5. **A metade externa da estrela está
no regime padrão da literatura de discos**, com q ≈ 1.5 perto de 0.8 R_eq.

**O custo é trivial.** Caixa de 4λ de lado:

| malha | lado | dx | Q |
|---|---|---|---|
| 64³ | 9.9×10⁴ cm | 1539 cm | 16 |
| **128³** | 9.9×10⁴ cm | **770 cm** | **32** |
| 256³ | 9.9×10⁴ cm | 385 cm | 64 |

Q = 32 no 128³ é o padrão para turbulência MRI convergida. Uma malha de 128³
sobre um domínio de **1 km**, sem autogravidade, sem Poisson, sem gradientes
globais. A coisa mais barata que se discutiu nesta campanha.

### O que resolve

Exatamente a objeção que sobrou do MInIT (seção 6.10): em vez de importar
α = −1.4 e β = −0.8 de estrela de nêutrons a densidade nuclear sob politrópica
γ = 2, mede-se nas nossas condições.

E um refinamento que reduz o trabalho pela metade: **os coeficientes da MRI,
α = 1 − 4/q, dependem só do cisalhamento** — são função de q e transferem
trivialmente. Os parasíticos vêm da saturação não linear, que depende da
microfísica. Só essa metade precisa de calibração, e é a que a caixa mede bem.

### O que não resolve, e o obstáculo

A caixa é um pedaço: dá coeficientes, não resultados sobre a estrela. O caminho
completo continua sendo **caixa → coeficientes → MInIT no run global**.

**O Castro não faz caixa de cisalhamento** — exige contorno periódico com
deslizamento, que não existe lá. Seria outro código; Athena++, PLUTO, Snoopy e
similares trazem o problema como exemplo, por ser um dos testes mais
padronizados da área. Custo de aprendizado, não de máquina.

Aproximação a declarar: caixas padrão usam EOS isotérmica ou gamma-law, não
`ztwd`. Num domínio de 1 km a densidade é praticamente constante, então uma
caixa isotérmica com c_s casado ao valor degenerado provavelmente basta —
**hipótese a verificar, não fato**.

## 6.11b A conta que faltava: o custo é temporal — TAMBÉM ERRADA

### No Castro, o inventário de implementação

**Já existe:** referencial rotativo com Coriolis e centrífuga, inclusive
atualização implícita do Coriolis, que é o termo chato.

**Trivial:** o termo de maré 2qΩ²x é fonte algébrica local, entra no
`problem_source.H`.

**Desnecessário:** advecção orbital tipo FARGO. A variação de velocidade
orbital na caixa é 1.2×10⁶ cm/s contra c_s = 3×10⁹ — o passo já é limitado
pelo som.

**O trabalho real:** contorno periódico com deslizamento, exigindo *remap*
conservativo das células fantasma, e o tratamento das EMF na fronteira para
preservar div B = 0 sob transporte restrito. [Stone & Gardiner
2010](https://arxiv.org/abs/1006.0139): se as EMF nas duas faces radiais não
coincidem, o fluxo vertical líquido não se conserva, e a MRI é sensível a isso.

### Mas o custo mata antes

| | |
|---|---|
| altura de escala efetiva H = c_s/Ω | 3.7×10⁸ cm ≈ R_eq |
| λ_MRI/H | 6.7×10⁻⁵ |
| dx/H na caixa | 2.1×10⁻⁶ |
| **passos por órbita** | **1.0×10⁷** |
| numa caixa de disco típica | 1.3×10³ |

**7400 vezes piores por órbita.** Cem órbitas dão 10⁹ passos numa malha 128³:
**2100 vezes o trabalho do run global 256³**, ou ~2×10⁷ core-hours contra 520
mil unidades de alocação total.

**Causa raiz: v_A/c_s ≈ 10⁻⁵**, β ~ 10⁹ no interior. Caixas de disco rodam
β = 100 a 1000. A MRI evolui numa escala absurdamente longa comparada ao tempo
de travessia sonoro que fixa o passo.

A caixa escapa da resolução **espacial** e não da **temporal**. E isso independe
de código: Athena++ teria o mesmo problema.

### O que restaria

Rodar caixas em β acessível (10³–10⁴) e extrapolar α_SS(β) até 10⁹. Prática
padrão, mas reintroduz o que a caixa deveria eliminar — extrapolação para fora
do regime medido, agora de cinco ordens de grandeza em β.

## 6.11c O custo escala como Q⁴, e Q era escolha minha

> Terceira versão, mesmo dia. Rafael perguntou o óbvio — "qual o problema de ser
> mais demorado? Quanto tempo?" — e a conta em horas de parede desmontou 6.11b.

**O erro comum às duas versões anteriores:** tratei Q = 32 como dado do
problema. Era escolha minha, feita sem nota, e é ela que domina o custo.
Células ∝ Q³, passos ∝ Q (dt ∝ dx), logo **custo ∝ Q⁴**. Nunca variei o
parâmetro que mandava no resultado. Ontem isso deu "barato" porque olhei só
células; hoje deu "inviável" porque olhei só passos — sempre no mesmo ponto.

Calibração a partir do nosso próprio 256³: 1.9×10¹² cell-steps por ~1.4×10⁴
core-h, ou 1.3×10⁸ cell-steps/core-h.

| Q | malha | core-h/órbita | 100 órbitas | % alocação | órbitas que os 497k restantes compram |
|---|---|---|---|---|---|
| 8 | 33³ | 6.5×10² | 6.5×10⁴ | 13% | 761 |
| **16** | **65³** | **1.0×10⁴** | **1.0×10⁶** | **201%** | **48** |
| 32 | 130³ | 1.7×10⁵ | 1.7×10⁷ | 3214% | 3 |

**Q = 16 é viável, e é a barra certa** — Sano+2004 pede Q ≥ 15–20 para
turbulência MRI convergida. A alocação restante compra 48 órbitas; a literatura
mede α com 30–100.

**Tempo de parede: ~40 h por órbita em 2 nós.** Trinta órbitas ≈ 50 dias de
execução contínua, 2–3 meses reais com fila e janelas de 3 h. Não comprime com
mais nós: 65³ em 256 cores já são 4300 células/core, e integração temporal é
sequencial. **O tempo de parede, não a alocação, é o custo que dói.**

**E a objeção de 6.11b cai.** A caixa roda no β real da estrela — é exatamente
isso que dá λ/H = 6.7×10⁻⁵ e o custo alto. Não há extrapolação em β a fazer:
α sai medido no regime que nos interessa.

**O que continua valendo de 6.11b:** o inventário de implementação no Castro. O
referencial rotativo existe, o termo de maré é trivial, a advecção orbital é
desnecessária, e o trabalho real é o contorno deslizante com remap conservativo
e EMF consistente com CT ([Stone & Gardiner 2010](https://arxiv.org/abs/1006.0139)).
Semanas de implementação antes da primeira órbita.

### Ordenação final das três saídas

1. **Caixa de cisalhamento a Q = 16** — a única rota que MEDE α no β da
   estrela, sem coeficiente emprestado nem extrapolação. Custo: ~10⁴ core-h por
   órbita, 30 órbitas dentro da alocação, 2–3 meses de parede, mais semanas de
   implementação do contorno deslizante no Castro. Cara em tempo, não em
   alocação.
2. **MInIT** — a classe certa de modelo, gancho já existe no
   `problem_source.H`, barato. Depende de coeficientes de outro regime; os da
   MRI (α = 1 − 4/q) transferem, os parasíticos não. É o que a caixa
   alimentaria.
3. **LES** — descartado. Modela cascata, e a Rm ≈ 5 não há cascata.

As duas primeiras compõem: a caixa calibra, o MInIT aplica no run global. Isso
é o programa, se houver 2–3 meses para gastar.

**Enquanto não houver, a limitação da MRI permanece e deve ser declarada como
tal nos relatórios** — na primeira página, não nas ressalvas.

## 6.12 A saída: o problema era o conjunto de equações, não o código

> 8 de agosto. Rafael: "estude uma solução viável, mesmo que migremos para o
> PLUTO". A busca desmontou tanto a premissa da pergunta quanto o custo de
> 6.11c.

### PLUTO não resolveria

Volume finito compressível, como o Castro: mesmo `dt = cfl·dx/c_s`. O que ele
traz pronto é caixa de cisalhamento com FARGO
([Mignone+2012](https://www.aanda.org/articles/aa/full_html/2012/09/aa19557-12/aa19557-12.html)),
poupando as semanas do contorno deslizante. Mas FARGO acelera quando a
velocidade orbital domina, e aqui Δv_orb/c_s = 4×10⁻⁴; o ganho deles é 3.75× e
precisamos de 10⁴. **Trocar de código dentro da mesma classe não muda nada.**

### O que resolve: incompressível

Os três cálculos de custo em 6.11/b/c pressupunham compressibilidade sem que eu
percebesse que era pressuposto. O passo é limitado por c_s, que aqui não carrega
física: Mach = 4×10⁻⁴, v_A/c_s = 10⁻⁵. A aproximação incompressível vale quando
fluxo E Alfvén são << c_s — no nosso caso não é marginal, é ideal.

| | compressível (Q=16, 64³) | incompressível espectral |
|---|---|---|
| limite do passo | dx/c_s = 1.4×10⁻⁷ s | dx/v_turb = 8×10⁻³ s |
| passos por órbita | 5.3×10⁶ | ~100 |
| 100 órbitas | 677 dias | ~10⁴ passos, **horas** |

O fator 5×10⁴ é simplesmente c_s/v_turb.

**Código: [SNOOPY](https://ipag.osug.fr/~lesurg/snoopy.html)** (Lesur),
pseudo-espectral, MHD Boussinesq/incompressível em caixa de cisalhamento, GPL.
Decomposição em ondas de cisalhamento com remap periódico: contorno deslizante
resolvido, e sem o problema de EMF do CT porque é espectral.

**Precedente na nossa classe de objeto:**
[Guilet & Müller 2015](https://arxiv.org/pdf/1501.07636) — MRI em protoestrelas
de nêutrons, Boussinesq escolhido justamente por dar contorno limpo mantendo
flutuabilidade; [Rembiasz+2016](https://arxiv.org/pdf/1603.00466);
Reboul-Salze+2021/2022 (dínamo αΩ em PNS). Objeto compacto, alto β, rotação
diferencial.

### E os coeficientes já podem existir

[Miravet-Tenés+2025](https://arxiv.org/abs/2509.07081) (MNRAS 545, aceito
set/2025) implementou **MInIT em simulações globais newtonianas de estrelas de
nêutrons magnetizadas e diferencialmente rotativas** — nosso problema, um objeto
ao lado. Código Aenus, HLL + PPM + RK3, esférico axissimétrico.

Coeficientes publicados:
- MRI, de teoria (Pessah & Chan 2008): α^MRI_ϖφ = 1 − 4/q, β^MRI_ϖφ = 1
- Parasíticos, calibrados em caixa
  ([Miravet-Tenés+2022](https://arxiv.org/abs/2210.02173)): α^PI = −1.4,
  β^PI = −0.8

Tratam explicitamente o nosso caso: **quando λ_MRI < Δ, substituem k por 2π/Δ**,
e zeram γ_MRI onde q ≤ 0 ou q ≥ 4. A célula deles é ~10× λ_MRI — deliberadamente
não resolvida, como a nossa.

Resultado deles: achatamento do perfil nas regiões internas, transporte de
momento angular para fora, Ω_max decaindo mais rápido com campo mais forte
(10¹⁴ G contra 3.5×10¹³ G). **Isto é exatamente a previsão contra a qual o nosso
"a rotação diferencial sobrevive" precisa ser testado.**

Nosso q pela lei de Komatsu: q = 2ϖ²/(A²+ϖ²), de 0 a 2 — dentro da faixa
0 < q < 4 em que o MInIT liga a MRI. A física está no domínio de validade.

### Ressalvas que eles próprios registram

- Subgrade **só no momento, não na indução**: transporta momento angular, não
  gera campo. Sem dínamo MRI de grande escala.
- Run deles axissimétrico; o nosso é 3D, o que joga a favor.
- Estimativa deles para 3D plenamente resolvido: ~10⁶ CPU-anos. Ninguém faz.

### Programa

1. **MInIT no Castro** — semanas, sem cluster novo. Tensores no
   `problem_source.H`. Validar reproduzindo os resultados publicados antes de
   aplicar à anã branca.
2. **Caixas SNOOPY no nosso β e q** — dias, estação de trabalho. Verificar se os
   coeficientes parasíticos transferem de PNS para anã branca. Era o item de 2–3
   meses de 6.11c; deixou de ser bloqueio e virou controle.
3. **256³ global com MInIT ligado** — mesmo custo dos runs atuais.

### Boris não serve, e a razão importa

[Matsumoto+2019](https://arxiv.org/abs/1902.02810), Boris-HLLD, aparece em toda
busca por "reduzir passo em MHD". Ele limita a velocidade de **Alfvén** e serve
a **baixo β**, onde v_A domina. O nosso é o oposto: alto β, c_s domina, v_A é
10⁻⁵ dela. Não ganha nada. Registro porque é a armadilha óbvia da busca.

A alternativa compressível seria RSST — reduzir c_s artificialmente
([Hotta+2014](https://iopscience.iop.org/article/10.1088/0004-637X/786/1/24),
[Iijima+2019](https://arxiv.org/abs/1812.04135)), padrão em convecção solar pelo
mesmo motivo. Daria ~10² em vez de 10⁴, e exigiria implementação. Incompressível
é estritamente melhor aqui, e já está escrito.

### A lição das quatro versões

6.11 a 6.12, mesmo dia. Erro em 6.11: custo por células só. Em 6.11b: custo por
passos só, com Q fixo sem perceber que era escolha. Em 6.11c: Q variado, mas
compressibilidade fixa, também sem perceber que era escolha. **Três vezes o
mesmo padrão — um parâmetro tratado como dado do problema quando era premissa
minha.** A correção veio de fora nas três, por Rafael perguntar "quanto tempo?"
e "e se migrarmos?".

### O teste científico que isto destrava

Nossos dois grids concordam que a rotação diferencial **acentua**: −22.6% (256³)
e −20.3% (192³). Miravet-Tenés+2025, com MInIT ligado, encontram **achatamento**
do perfil interno, transporte de momento angular para fora, e Ω_max decaindo
mais rápido quanto mais forte o campo (10¹⁴ G contra 3.5×10¹³ G).

**São previsões de sinal oposto.** Se ligarmos o MInIT e o perfil achatar, o
nosso acentuamento era a assinatura da MRI ausente. A Fase 1 deixa de ser
refinamento e vira o teste decisivo do resultado central dos dois relatórios.

## 6.13 MInIT: a formulação, e por que casa com o Castro

Equações extraídas de [Miravet-Tenés+2025](https://arxiv.org/html/2509.07081v1),
com o modelo original em [Miravet-Tenés+2022](https://arxiv.org/abs/2210.02173).
PDFs em `references/minit_2022_mri.pdf` e `references/minit_2025_ns.pdf`.

### As duas equações de evolução

    ∂_t e_MRI + ∇_j(v̄_j e_MRI) = 2 γ_MRI e_MRI − 2 γ_PI e_PI
    ∂_t e_PI  + ∇_j(v̄_j e_PI)  = 2 γ_PI  e_PI  − C e_PI^{3/2} / (√ρ̄ λ)

com C = 8.6 e λ = min[Δ, λ_MRI].

Taxas, calculadas do campo **resolvido** em cada célula:

    γ_MRI = (q/2) Ω                     q ≡ −d ln Ω / d ln ϖ
    γ_PI  = σ k_MRI √(2 e_MRI/ρ̄)        σ = 0.27
    k_MRI = √(1 − (2−q)²/4) · Ω / v̄_Az   v̄_Az = b̄_z/√ρ̄

### Os tensores, algébricos nas energias

    M̄_ij = α^MRI_ij e_MRI + α^PI_ij e_PI          (Maxwell)
    R̄_ij = (β^MRI_ij e_MRI + β^PI_ij e_PI)/ρ̄      (Reynolds)
    F̄_ij = γ^PI_ij e_PI / √ρ̄                      (Faraday)

Coeficientes: α^MRI_ϖφ = 1 − 4/q e β^MRI_ϖφ = 1 vêm de **teoria** (Pessah & Chan
2008); α^PI_ϖφ = −1.4 e β^PI_ϖφ = −0.8 vêm de **caixa** (Miravet-Tenés+2022).

Entram no momento:

    ∂_t p̄_i + ∇_j[ρ̄ v̄_i v̄_j + (P̄* + Tr M̄)δ_ij − b̄_i b̄_j + ρ̄ R̄_ij − M̄_ij] = f̄_i

### Por que não é LES, e por que a objeção da Rm não se aplica

Descartei LES em 6.10 porque a Rm ≈ 5 não há cascata inercial para modelar. O
MInIT **não modela cascata**. É um balanço de energia de instabilidade:

1. a MRI cresce a γ_MRI = qΩ/2, taxa **linear e analítica**, tirada do
   cisalhamento resolvido — não precisa de turbulência nenhuma;
2. as instabilidades parasitas (Kelvin–Helmholtz e tearing sobre os modos-canal)
   comem essa energia e saturam o crescimento — é o mecanismo de saturação de
   Goodman & Xu 1994 e Pessah 2010;
3. os tensores são função algébrica das duas energias.

Nunca é preciso resolver λ_MRI. Precisa-se de Ω(ϖ), q, ρ e B_z resolvidos — que
temos. **A objeção da cascata era contra o LES e continua válida; não alcança o
MInIT.**

### Encaixe no Castro

| Peça | Situação |
|---|---|
| e_MRI, e_PI como escalares advectados | nativo — `NumAdv` / estado auxiliar |
| divergência dos tensores no momento | `problem_source.H`, gancho que já usamos no damping |
| fontes rígidas (stiff) | padrão das reações; `do_react=0` hoje, mas a infra existe |
| B_z, ρ, v | direto do estado |
| EOS barotrópica | **simplifica** — sem equação de energia para acertar |

**O único trabalho real: q = −d ln Ω/d ln ϖ em coordenadas cartesianas.** O run
deles é esférico axissimétrico, onde ϖ é coordenada. Aqui Ω = (x v_y − y v_x)/ϖ²
e a derivada é direcional ao longo de ϖ̂ = (x,y,0)/ϖ, por diferenças finitas no
grid resolvido. Cuidado no eixo (ϖ→0) e no envelope de baixa densidade, onde Ω é
ruidoso.

Eles zeram γ_MRI onde q ≤ 0 ou q ≥ 4. Nossa lei de Komatsu dá
q = 2ϖ²/(A²+ϖ²), variando de 0 a 2 — **dentro da faixa em que o MInIT liga a
MRI.** E quando λ_MRI < Δ, trocam k por 2π/Δ; a célula deles é ~10× λ_MRI,
deliberadamente não resolvida, como a nossa.

## 6.14 SNOOPY: em incompressível, β desaparece — e sobra o Pm

### O que a caixa incompressível conhece

A pressão térmica é multiplicador de Lagrange impondo ∇·v = 0; c_s é infinita e
some. **β não é parâmetro da caixa incompressível.** Os grupos adimensionais que
restam:

    q          parâmetro de cisalhamento
    L/λ_MRI    quantos comprimentos de onda cabem
    Re = ΩL²/ν,  Rm = ΩL²/η,  Pm = ν/η

O β = 10⁹ que tornava a caixa compressível cara **não aparece na física**, só na
numérica que escolhemos. Isso responde à ressalva que eu tinha deixado aberta:
não é β que precisa casar entre PNS e anã branca.

Bônus de consistência: nossa EOS `ztwd` é barotrópica, logo N² = 0, logo sem
flutuabilidade. Guilet & Müller precisaram do modo Boussinesq para gradiente de
entropia; **para nós o modo incompressível puro é exatamente certo.**

### Mas o Pm é um problema, e talvez um resultado

Anã branca: ν ~ 3×10⁻² cm²/s, η ~ 6×10⁻² cm²/s, **Pm ~ 0.5–0.6**, Rm ~ 10¹⁴–10¹⁵
(zona convectiva de anã branca CO em cristalização —
[Fuentes+2024](https://iopscience.iop.org/article/10.3847/2041-8213/ad3100)).
Protoestrela de nêutrons: Pm enorme. **Os coeficientes parasitas foram
calibrados no regime errado para nós.**

E a literatura diz que isso importa muito. Fromang+2007 e Lesur & Longaretti 2007
acham lei de potência íngreme do transporte com Pm, e **Pm crítico ~ 2–4 abaixo
do qual a turbulência MRI morre**. Pm ~ 0.58 está abaixo disso.

Ressalva que impede conclusão apressada, em duas frentes:

- O Pm_crit ~ 2–4 é de caixa **sem fluxo líquido**, onde a MRI depende de dínamo.
  Temos **fluxo vertical líquido**: a MRI é instabilidade linear e cresce de
  qualquer jeito. [Simon & Hawley
  2011](https://iopscience.iop.org/article/10.1088/0004-637X/740/1/18) acham que
  o transporte sobrevive a Pm baixo desde que Rm supere um crítico — e o nosso
  Rm físico é 10¹⁴.
- O valor Pm ~ 0.58 é de anã branca CO fria em cristalização, não de remanescente
  de fusão a 2 M⊙ e quente. **Precisa ser recalculado para as nossas condições.**

De qualquer modo isto promove a caixa: deixa de ser tarefa de calibração e vira
**pergunta física própria** — a MRI transporta momento angular no interior de
uma anã branca a Pm de ordem unidade?

### Como usar o SNOOPY

[Página do Lesur](https://ipag.osug.fr/~lesurg/snoopy.html), GPL, C com FFTW3
e MPI/OpenMP. Adimensional: fixa-se Ω = 1 e o tamanho da caixa, e o que se
escolhe é q, L/λ_MRI, Re, Rm.

Plano mínimo, dias e não meses:

1. compilar e reproduzir um caso MRI padrão de disco (q = 1.5) contra a
   literatura — validação do nosso uso, não do código;
2. varrer q em 0.5, 1.0, 1.5, 2.0 cobrindo a faixa de Komatsu;
3. varrer Pm em torno de 1 (0.25 a 4) com Rm alto, medindo se o transporte
   sobrevive e como escala;
4. extrair α^PI e β^PI pelo procedimento de Miravet-Tenés+2022 e comparar com
   −1.4 e −0.8.

Se baterem, o MInIT roda com os coeficientes publicados. Se não, temos os
nossos, medidos no nosso regime — que é um resultado publicável por si.

### Riscos abertos

- SNOOPY V6.0 é de 2011; verificar se compila com FFTW3 e MPI atuais.
- Recalcular ν e η para 2 M⊙ quente e degenerado, não para anã branca fria.
- Procedimento de extração dos coeficientes está no 2022; ler antes de rodar.

## 6.15 Duas ressalvas do MInIT, e a caixa volta para primeiro

> Rafael: "o MInIT no Castro pode resolver esse problema da MRI? E por que não
> usamos antes?" Depois: "mas é possível que o SNOOPY seja muito melhor pra nós?"

### Por que não usamos antes: eu o coloquei em segundo

§6.10 é **de hoje de manhã**. Já cita o artigo de 2025, já conclui "viável", já
registra o gancho, os escalares nativos e o tratamento de malha grossa. Registra
inclusive algo mais forte do que eu afirmei à tarde: **a lei de rotação do artigo
deles é a nossa Komatsu j-constante, e os coeficientes foram derivados para esse
perfil** — a confirmar em `references/minit_2025_ns.pdf`, e isso decide se os
parasitas transferem.

Não foi descuido de meses; foram horas. Pus o MInIT em segundo pela objeção dos
coeficientes emprestados, disse que a caixa resolveria, e gastei o dia errando o
custo da caixa três vezes. **A caixa nunca foi pré-requisito; eu a transformei em
um.**

### Ressalva 1: o MInIT modela a MRI errada para o nosso campo

    k_MRI = √(1 − (2−q)²/4) · Ω / v̄_Az        v̄_Az = b̄_z/√ρ̄

Só a componente **vertical**.

> **CORRIGIDA no mesmo dia — ver §6.17.** Eu escrevi aqui que o campo é
> toroidal-dominado por 10⁷, número do **modelo analítico em t=0**. A estrela
> evoluída fica entre E_tor/E_pol = 4 e 90. A ressalva enfraquece muito.

O que sobra dela: o que destrói o campo nos nossos runs é o kink de Tayler m=1,
instabilidade diferente da que o modelo carrega, e o MInIT não a cobre. Ligar o
MInIT continua legítimo, mas a resposta significa "quanto a MRI axissimétrica
acrescentaria", não "o que a física ausente faz".

### Ressalva 2: com B_z fraco a saturação vira dependente de malha

λ_MRI cai muito abaixo da célula, o modelo entra permanentemente no ramo
k → 2π/Δ, e então γ_PI = σ(2π/Δ)√(2e/ρ). **A amplitude saturada passa a ser
fixada pela resolução.**

Isto nos favorece: temos 192³ e 256³ com diagnóstico idêntico. Ligar o MInIT nos
dois e testar convergência do transporte é teste que o artigo deles
(axissimétrico, resolução única) não faz.

### O que o MInIT não cobre de jeito nenhum

Os termos de subgrade entram **só no momento, não na indução**. Nosso outro
resultado principal — campo ordenado destruído em poucos tempos de Alfvén, três
ordens de energia magnética — fica intocado. O modelo não diz se um dínamo MRI
teria regenerado o campo. Metade da história segue sem cobertura.

### Reordenação: SNOOPY primeiro

Quatro razões, e nenhuma é preferência:

1. **Custo** — dias rodando código publicado, contra semanas escrevendo o nosso.
2. **Erro silencioso** — MInIT malimplementado não quebra o run; produz
   transporte plausível. Um bug de subgrade é indistinguível de física.
3. **Ordem lógica** — o MInIT pressupõe que a MRI axissimétrica opera. Se ela não
   sustenta turbulência no nosso Pm, ele dá resposta confiante e errada, e é a
   que iria para o artigo.
4. **Testa a ressalva 1** — uma caixa pode ter fluxo toroidal líquido mais
   vertical fraco na nossa proporção, e medir qual instabilidade domina. O MInIT
   não pode fazer esse teste sobre si mesmo.

**Mas SNOOPY é local e não substitui nada.** Não diz se a rotação diferencial da
nossa estrela sobrevive; isso exige o run global. O programa segue
caixa → coeficientes → run global. Muda só que a caixa é o primeiro passo.

### O preço honesto do SNOOPY: a extrapolação vai de β para Rm

DNS em caixa alcança Rm ~ 10³–10⁴; o nosso físico é ~10¹⁴. Fixa-se Pm = 0.6
exatamente, mas com Re e Rm muito abaixo dos reais — e se a dependência em Pm
medida a Rm = 10⁴ vale a Rm = 10¹⁴ é problema em aberto desde Fromang+2007.

Atenuante: com fluxo vertical líquido a MRI é instabilidade linear, bem menos
sensível a Rm que o dínamo sem fluxo.

Ainda assim é ganho: sai de "coeficientes emprestados de outro objeto em regime
desconhecido" para "coeficientes medidos no nosso q e no nosso Pm, extrapolados
em Rm" — **uma extrapolação em vez de três, e num parâmetro que sabemos
nomear.** Precisa constar do artigo.

## 6.16 Outros códigos: FLASH, Einstein Toolkit, e o princípio geral

> Rafael: "e outros programas, como Flash, Einstein Toolkit...".

### O princípio, antes da lista

**Nenhum código da classe volume-finito-compressível resolve a MRI.** FLASH,
PLUTO, Athena++, Einstein Toolkit, AREPO — mesma parede, porque o obstáculo é o
conjunto de equações mais o regime físico, não a implementação. Duas alavancas
só: mudar as equações (incompressível, §6.12) ou modelar o que falta (§6.13).
Trocar de código dentro da classe é movimento lateral.

### FLASH — já respondido, e negativamente

Investigação completa em `flash_crosscheck/README.md`: FLASH 4.8 com Helmholtz,
unidade `WDHydrostatic`, modelo ztwd validado a 0.026% em massa, damping e
sponge portados, correção de meia-célula da §6.6 reproduzida em dois códigos
independentes (−1.23% contra −1.16%).

**FLASH segura a estrela ~10× pior** — −9.03% contra −0.91% no limite da janela,
e com forma qualitativamente distinta: Castro sobe, vira e assenta; FLASH cai
monotonicamente e ainda caía ao morrer. Nenhuma configuração do FLASH produz
janela de medida válida; o bloqueio é o binning de zona interna do
`Multipole_new`, que não tolera célula de raio zero.

Resultado negativo caracterizado, não pergunta aberta. E não teria tocado a MRI.

### Einstein Toolkit — a compacidade decide

    GM/Rc² = 9.3×10⁻⁴        estrela de nêutrons ≈ 0.2, fator 215
    v_rot/c = 7.4×10⁻²

Correções de RG na casa de 0.1%, contra incerteza dominante de resistividade
numérica a Rm ≈ 5. Pagaríamos evolução da métrica (~25 variáveis extras por
célula) e MHD tipicamente com limpeza de divergência em vez de CT. Não.

### Onde trocar de ferramenta ajudaria

| Problema nosso | É de código? |
|---|---|
| EOS barotrópica, 6×10⁴⁹ erg sem destino | **Sim — e já temos.** `microphysics/EOS/helmholtz` está no repositório. Sem migração. |
| MHD do Castro é nível único, sem AMR | **Sim, limitação real.** FLASH tem AMR-MHD, mas segura a estrela 10× pior. |
| Resistividade numérica Rm ≈ 5 | Não. Só resolução. |
| MRI não resolvida | Não. Equações ou subgrade. |

O primeiro é o mais aproveitável: a ressalva do EOS barotrópico está nos dois
relatórios e a solução está instalada.

### Dois não considerados antes

- **Idefix** — GPU, do mesmo Lesur do SNOOPY, caixa de cisalhamento nativa. Se o
  CENAPAD tiver nós com GPU, é eixo de ganho que não avaliamos. **A perguntar.**
- **AREPO** — malha móvel, ferramenta da comunidade de fusão de anãs brancas;
  Pakmor & Pelisoli 2024 está citado no nosso `inputs.ml256`. MHD lá é
  Powell/oito ondas, não CT — recuo. Relevante para conectar com a literatura de
  remanescentes, não para a MRI.

### 6.16b GPU: alavanca errada, e a URL do SNOOPY estava morta

> Rafael: "SNOOPY usa GPU?"

**Correção de link.** `ipag-old.osug.fr` recusa conexão. O endereço vivo é
`https://ipag.osug.fr/~lesurg/snoopy.html`, que responde 200 mas está atrás de
proteção anti-bot (Anubis) — download tem de ser manual, por navegador.

**SNOOPY não usa GPU.** Nada na literatura nem na descrição menciona CUDA,
OpenCL ou Kokkos; é C com FFTW3, MPI e OpenMP, V6.0 de 2011. Não verificado na
fonte por causa do Anubis: tratar como muito provável, não confirmado. Há razão
estrutural além da data — pseudo-espectral faz FFT global a cada passo, e FFT
distribuída exige *all-to-all*; é o padrão que menos ganha com GPU, porque o
gargalo é rede, não aritmética.

**E GPU é a alavanca errada.** SNOOPY roda em horas; acelerar horas não muda
decisão. O ganho já foi capturado trocando as equações — incompressível deu
5×10⁴, GPU daria 10–100.

**Onde mudaria: Idefix.** [Mignone/Lesur+2023](https://arxiv.org/abs/2304.13746),
Kokkos, V100 e Mi250, com caixa de cisalhamento, CT, advecção orbital e MHD
não-ideal; o artigo cita 12 simulações em ~3000 GPU-horas. No nosso caso
compressível, 10⁴ core-h/órbita valendo uma GPU por ~100 cores dá ~100
GPU-h/órbita e **~3000 GPU-horas para 30 órbitas**. GPU *resgata* a rota
compressível da inviabilidade — mas a incompressível custa horas de CPU e
continua 10³–10⁴ à frente. **GPU conserta o método errado; o método certo
dispensa GPU.**

Consequência: a pergunta sobre GPU no CENAPAD (§6.16) perde urgência. Importaria
indo de Idefix; indo de SNOOPY, uma estação de trabalho basta.

**A ler antes de fixar parâmetros das caixas:**
[On the numerical convergence of MRI simulations](https://arxiv.org/pdf/2511.06022)
(nov/2025) — fala direto da nossa discussão de Q.

## 6.17 O campo NÃO é toroidal por 10⁷, e isso reordena o programa

> Rafael: "o SNOOPY vai servir para calibrarmos o impacto do campo toroidal no
> MRI?" Fui verificar a razão antes de responder e ela me desmentiu.

### O número que eu vinha repetindo é da condição inicial

| t (s) | E_tor/E_pol |
|---|---|
| 0.03 | 7.3×10⁵ |
| **13.5** | **4.0** ← mínimo |
| 32.5 | 20 |
| 78.0 | 89 |

O 2.18×10⁷ é do **modelo analítico**; na malha em t=0 já são 3.1×10⁵, e a tabela
do relatório sempre trouxe as duas colunas — o erro foi só meu, ao repetir o
número do modelo como se descrevesse a estrela evoluída. Em amplitude,
B_tor/B_pol fica entre 2 e 9.5. **O poloidal deixa de ser desprezível nos
primeiros segundos.**

Corrigido em `references/README.md` e na ressalva 1 da §6.15. A legenda do
`report_rot192_rot256.tex` está correta: é tabela de t=0 e traz modelo e malha
lado a lado.

### O que o SNOOPY pode e não pode sobre o toroidal

**Pode.** Uma caixa carrega fluxo azimutal líquido e vertical líquido ao mesmo
tempo, na razão que se quiser, e modos não-axissimétricos são k_y ≠ 0 — nativos
num código espectral. Dá para medir se o transporte é dominado pela MRI de campo
vertical (a que o MInIT modela) ou pela azimutal, e qual α sai de cada.

**Não pode: Tayler.** A instabilidade é movida pela curvatura das linhas
toroidais em torno do eixo, a tensão de aro. Caixa cartesiana não tem curvatura,
logo o m=1 não cresce ali por construção. E o kink de Tayler é o que destrói
nosso campo. **Nem SNOOPY nem MInIT tocam o segundo resultado principal** — não
por falta de ferramenta, mas por o problema ser global.

### A medida pendente virou a primeira coisa a fazer

Com B_tor/B_pol entre 2 e 9.5, v_Az é fração real de v_A e λ_MRI não é minúsculo.
O número que temos — λ/dx cruza 6 em t ≈ 2–3 s e chega a 36 no 256³ — foi
calculado do **pico** de B_z. Pico não é típico, e medir o típico em volume
ficou pendente há dias.

**Se o típico ficar em 6–36, já estamos resolvendo a MRI axissimétrica**, ao
menos marginalmente. Então "a rotação diferencial acentua" já inclui MRI, e nem
SNOOPY nem MInIT são o próximo passo.

Usa fatias que já estão em disco, não precisa de cluster, e **decide qual dos
três programas é sequer necessário.** Vai na frente de tudo.

### Ordem revisada

0. **λ_MRI/dx do B_z típico em volume**, das fatias em disco. Horas.
1. SNOOPY, se 0 disser que não resolvemos.
2. Coeficientes.
3. MInIT nos dois grids.

E segue sem plano para o campo, em qualquer ramo.

## 6.18 MEDIDO: não resolvemos a MRI. Q ≈ 0.4 onde o resultado é medido

`investigations/mri_wavelength.py`, figura em `mri_wavelength.pdf`.

### Não precisou de fatia nenhuma

Eu havia dito que a medida exigia plotfiles, que estão no lovelace — errado.
E_pol **é** a integral de volume de B_pol²/8π sobre a estrela e já está no
`bt_bp_256_long.csv` nos 210 instantes. Então

    B_pol,rms = √(8π E_pol / V)

é exatamente o campo típico em volume, sem tocar em fatia. V vem do R_vol da
mesma linha, logo os dois são consistentes por construção.

### O resultado

| t (s) | B_pol,rms (G) | λ_MRI (cm) | Q = λ/dx |
|---|---|---|---|
| 1.0 | 6.2×10¹⁰ | 3.1×10⁶ | 0.43 |
| 3.0 | 2.6×10¹¹ | 2.0×10⁷ | 2.80 |
| **4.5** | ~9×10¹¹ | ~4.6×10⁷ | **6.5 ← máximo** |
| 13.5 | 5.2×10¹¹ | 2.3×10⁷ | 3.30 |
| 32.5 | 1.1×10¹¹ | 3.8×10⁶ | 0.55 |
| 45.0 | 4.7×10¹⁰ | 1.9×10⁶ | 0.27 |
| 78.0 | 7.0×10¹⁰ | 2.6×10⁶ | 0.37 |

Três hipóteses para B_z a partir de B_pol, porque o CSV só traz a energia
poloidal:

| B_z | Q máximo | Q ≥ 6 | Q ≥ 15 |
|---|---|---|---|
| B_pol (tudo vertical) | 9.2 em t=4.5 | 10% das amostras, 1ª em t=3.50 | **nunca** |
| B_pol/√2 (central) | 6.5 em t=4.5 | 3%, 1ª em t=4.04 | **nunca** |
| B_pol/√3 (isotrópico) | 5.3 em t=4.5 | **nunca** | **nunca** |

Sensibilidade a qual Ω se usa: Q máximo vai de 3.8 (Om_core) a 10.7 (Om_out).

### O que isto decide

**O "36" estava errado por fator 5–7.** Vinha do pico de B_z; o típico em volume
é muito menor. A afirmação registrada de que λ/dx cruza 6 em t ≈ 2–3 s e chega a
36 não sobrevive.

**Não resolvemos a MRI**, exceto marginalmente numa janela estreita em
t ≈ 3.5–5 s, e mesmo aí só sob a hipótese mais favorável de B_z. É o momento em
que E_pol picou, durante a ruptura m=1. Três por cento do run.

**Q = 15 nunca é alcançado em hipótese nenhuma.** Turbulência MRI convergida
está fora de alcance por construção nesta malha.

**E o mais importante: onde o resultado principal é medido, Q ≈ 0.3–0.5.** O
acentuamento de −22.6% da rotação diferencial é medido em t = 40–78 s, faixa em
que estamos 12 a 20 vezes abaixo da barra linear. **A rotação diferencial
sobrevive num regime onde a MRI não tem como existir na malha.**

### Consequências

1. O programa dos três passos **está de pé**. SNOOPY e MInIT são necessários, e
   a medida não os dispensou como eu havia especulado.
2. Isto é exatamente o número quantitativo que a primeira página dos dois
   relatórios precisa: não "não resolvemos a MRI", mas **"Q ≈ 0.4 na janela em
   que a medida é feita, contra Q ≥ 6 para o modo linear e Q ≥ 15 para
   turbulência convergida"**.
3. A ressalva 1 da §6.15 volta a ter força, por caminho diferente do original.
   Não é que o campo seja toroidal por 10⁷ — é que o poloidal, embora
   comparável, é **fraco em termos absolutos**: B_pol,rms ~ 7×10¹⁰ G dá
   λ_MRI = 2.6×10⁶ cm contra dx = 7.0×10⁶ cm.

### Limitações da própria medida, declaradas no cabeçalho do script

- B_z inferido de B_pol; três casos reportados, espalhamento de 1.73 — menor que
  o efeito pico-vs-típico que a medida existe para capturar.
- ρ_mean, não ρ local. No núcleo ρ é maior e v_A menor.
- ⟨B²⟩^½/√⟨ρ⟩ não é ⟨B/√ρ⟩. Estimativa de valor típico, não média da velocidade
  de Alfvén local.

Nenhuma merece mais precisão: a resposta é fator 15 abaixo do limiar, e as três
juntas não movem isso.

## 6.19 Pm ≈ 750, não 0.58 — e isso é boa notícia

`investigations/magnetic_prandtl.py`.

### A reversão

Eu citei Pm ~ 0.58 a partir de um resultado de busca e construí em cima disso a
ressalva de que os coeficientes do MInIT não transferiam e de que estaríamos
**abaixo** do Pm crítico. Calculado para as nossas condições, **Pm ≈ 750**.

O 0.58 é da zona convectiva de uma anã branca CO **fria em cristalização** —
outro objeto, densidade e temperatura muito menores, e possivelmente com
viscosidade iônica dominando em vez da eletrônica.

### O cálculo

Transporte eletrônico em matéria degenerada, [Shternin
2008](https://arxiv.org/abs/0803.3893) Eq. (1)–(2):

    η_visc = n_e v_F p_F/(5 ν_e),   ν_e = ν_ee + ν_ei

A assimetria que fixa o Pm: **colisões elétron-elétron conservam o momento
total dos elétrons**, logo não degradam corrente e **não entram na
condutividade** — mas degradam tensão de cisalhamento e entram na viscosidade.

    σ = n_e e²/(m* ν_ei)          só ei
    η_mag = c²/(4πσ)
    Pm = (η_visc/ρ)/η_mag ∝ 1/ν_coll²

Como ν_ei ∝ n_i e tanto ν quanto σ carregam um n_e compensante, **a densidade
iônica cancela** e o Pm depende de T só via ν_ee e o logaritmo de Coulomb —
muito menos sensível a T do que eu supunha.

### Resultado, ρ = 4.8×10⁷ g/cm³, T = 10⁸ K

| | |
|---|---|
| x = p_F/m_e c | 2.91 (relativístico) |
| Γ (acoplamento iônico) | 13 — líquido, fusão em ~175 |
| ν_ee/ν_ei | 0.038 — **ei domina**, ee é correção |
| σ | 2.2×10²² s⁻¹ |
| ν | 2.4 cm²/s |
| η | 3.3×10⁻³ cm²/s |
| **Pm** | **746** |

Robusto: Pm > 190 em toda a faixa ρ = 5×10⁷–10⁹ e T = 10⁷–10⁹, e permanece 31
mesmo com logaritmo de Coulomb = 5. Para chegar a Pm = 1 seria preciso Λ ≈ 27,
impossível.

Escala estelar: Re = 2.9×10¹⁷, Rm = 2.2×10²⁰.

### Duas consequências, ambas a favor

1. **A preocupação com Pm_crit ~ 2–4 morre.** Estamos duas a três ordens
   **acima**, não abaixo. A turbulência MRI não corre risco de morrer por Pm
   baixo.
2. **Estamos no mesmo regime de alto Pm das caixas de protoestrela de nêutrons**
   onde os coeficientes parasitas do MInIT foram calibrados. A chance de
   transferirem é muito maior do que eu temia em §6.14.

### Mas a caixa continua sem alcançar

Pm = 750 exige Rm = 750·Re. Uma DNS chegando a Re ~ 10³ precisaria de
Rm ~ 10⁶ — fora de alcance, como o Rm = 10²⁰ físico.

**A caixa limita em vez de medir:** varrer Pm = 1, 2, 4, 8, 16 a Re fixo, medir
a tendência, e declarar a extrapolação. E aqui o sinal é conhecido — **o
transporte cresce com Pm na literatura de discos, logo uma caixa a Pm = 16 dá
um LIMITE INFERIOR do transporte a 750.** Limite inferior basta, se já apagar a
rotação diferencial.

A extrapolação sai de β e vai para (Rm, Pm) — mas num parâmetro cujo sinal
sabemos.

### Limitações declaradas no cabeçalho do script

- Λ_ei é o ponto fraco: Pm ∝ 1/Λ² e Λ é de ordem 1 mas não é conhecido sem
  tabelas de condutividade (condegin de Potekhin, Itoh et al.). Fica como
  parâmetro, com a sensibilidade impressa.
- Correlações íon-íon em líquido fortemente acoplado suprimem ν_ei abaixo da
  estimativa não correlacionada, o que **empurraria Pm para cima**.
- Viscosidade iônica e radiativa ignoradas. No interior degenerado os elétrons
  dominam; num envelope convectivo frio talvez não — a razão mais provável de o
  valor da literatura diferir.
- **Nossa EOS é barotrópica e não tem temperatura alguma.** T é hipótese sobre
  remanescente de fusão, não saída da simulação. Daí a varredura.

## 6.20 Varredura de Pm no SNOOPY — lançada

`shearing_box/scan_pm.sh`, análise em `analyse_pm_scan.py`.

### Por que q = 1.5 e não o nosso q, primeiro

Cisalhamento Kepleriano é onde Lesur & Longaretti 2007 e Fromang+2007
publicaram α(Pm). Rodar ali **mede a tendência que precisamos e ao mesmo tempo
confere o nosso uso do código** contra números de terceiros. Varrer o nosso
q = 0–2 vem depois que essa conferência passar; fazer na ordem inversa deixaria
qualquer discordância sem atribuição possível.

### A configuração padrão já é o nosso caso

Verificado antes de mexer: `bz0 = 0.1`, **fluxo vertical líquido** — a MRI é
instabilidade linear, não dínamo, que é a nossa situação. E λ_MRI/dz = 40,
folgadamente resolvido. `by0` existe para o toroidal líquido que vamos precisar
depois, na razão 2 a 9.5 da §6.17.

### Parâmetros

| | |
|---|---|
| malha | 64³, caixa (4,4,1) |
| q | 1.5 (Kepleriano, para validação) |
| Re | 1000, **fixo** — só a escala resistiva se move |
| Pm | 1, 2, 4, 8, 16 |
| t_final | 200 ≈ 32 órbitas; média sobre t > 60 |
| execução | 3 concorrentes × 4 threads = 12 cores |
| duração medida | ~87 min por run, ~2.9 h no total |

### Ressalva declarada antes de existirem números

**A 64³ o topo da varredura é marginal.** Rm = 16000 põe a escala resistiva
perto da malha, e as rodadas publicadas de alto Pm usam 128³. **Pm = 8 e 16 são
provisórios** até que um deles seja repetido a 128³ — o que custa 16× e fica
para depois.

### Convenções, para não errar sinal depois

SNOOPY trabalha em unidades de Alfvén: b **é** v_A, e a energia magnética é
b²/2. O fluxo total de momento angular é

    W = ⟨v_x v_y⟩ − ⟨b_x b_y⟩       (Reynolds + Maxwell)

com Maxwell entrando como **menos** ⟨bxby⟩ — daí `bxby` sair negativo numa
rodada que transporta para fora. Sem velocidade do som numa caixa
incompressível, o α = W/c_s² usual não existe; usa-se a normalização de fluxo
líquido, W/B0² com B0 = bz0.

### O que a varredura pode e não pode concluir

Nossa estrela está em Pm ≈ 750 (§6.19), que DNS nenhuma alcança. **A varredura é
limite, não medida.** O que se extrai é a *inclinação* de α(Pm), e ela é
extrapolada 1.7 décadas além da faixa ajustada. Isso tem de ser citado junto com
o resultado: é limite inferior **apenas se a tendência não saturar**, e estes
dados não estabelecem isso.

Critério de sanidade fixado agora: a literatura acha inclinação de ordem 0.5–1
sobre Pm = 1–16. **Inclinação perto de zero ou negativa significa que o nosso
arranjo está errado, não a física.**

## 6.21 Resultado da varredura: inclinação 0.19, e o topo achatou

Figura em `shearing_box/pm_scan.pdf`. Cinco rodadas completas a 31.8 órbitas,
1401 amostras saturadas cada, média sobre t > 60.

| Pm | ⟨E_mag⟩ | Maxwell | Reynolds | W/B0² | Max/Rey |
|---|---|---|---|---|---|
| 1 | 0.262 | 0.160 | 0.0335 | 19.32 ± 0.21 | 4.76 |
| 2 | 0.326 | 0.191 | 0.0362 | 22.67 ± 0.26 | 5.26 |
| 4 | 0.389 | 0.215 | 0.0374 | 25.24 ± 0.30 | 5.74 |
| 8 | 0.499 | 0.271 | 0.0479 | **31.84 ± 0.33** | 5.65 |
| 16 | 0.511 | 0.272 | 0.0473 | **31.88 ± 0.32** | 5.74 |

Ajuste: **W/B0² = 19.7 · Pm^0.194**

### O critério que registrei estava errado — mesmo erro de antes

Escrevi em §6.20 que a literatura dá inclinação 0.5–1. **Errado.** O resumo do
próprio [Lesur & Longaretti 2007](https://arxiv.org/abs/0704.2943) diz
**δ entre 0.25 e 0.5**, sobre 0.12 < Pm < 8 e 200 < Re < 6400.

Peguei "0.5–1" de um resumo de busca e registrei como se tivesse verificado —
**exatamente o que produziu o erro do Pm = 0.58 na §6.19**, dois dias seguidos.
A regra que falta é simples: *pré-registro só vale contra número que eu li na
fonte.* Um critério pré-registrado errado é pior que nenhum, porque dá falsa
autoridade ao veredito.

Com a faixa correta, nosso 0.19 fica **logo abaixo** de 0.25–0.5. Mesmo sinal,
mesma ordem. Aceitável para uma resolução não validada, mas **não é
concordância** e não deve ser relatado como tal.

### O achado real: o topo achatou

Pm = 8 e Pm = 16 concordam em 0.039, contra 2σ = 0.92. **Indistinguíveis.**

Duas leituras, e os dados não separam:

1. **α(Pm) satura de verdade** acima de Pm ~ 8. Note que LL07 parou em Pm = 8 —
   ninguém mediu além.
2. **A 64³ a resistividade numérica já domina**, e o código não distingue
   Rm = 8000 de Rm = 16000. Foi a ressalva que registrei em §6.20 antes de
   existirem números, e ela disparou.

O teste que separa: **Pm = 16 a 128³.** Se W subir, o ponto de 64³ era limitado
por resolução. Custo 16× o de 64³ — cerca de 17 h com 12 threads a t_final = 100.
**Não lancei; é a estação de trabalho do Rafael e são muitas horas de todos os
cores.** Decisão dele.

### O que dá para dizer sobre a estrela, hoje

Transporte no nosso Pm ≈ 746, relativo a Pm = 1:

- se a tendência continuar: **3.7×** (extrapolando 1.7 décadas)
- se já saturou em Pm ~ 8: **1.6×**

**Nos dois casos é fator de poucas unidades, não ordens de grandeza.** Isso é
resultado útil por si: a correção de alto Pm ao transporte MRI é modesta, então
o que decide o destino da rotação diferencial é o α em si, não o Pm.

Também favorável ao MInIT: se α varia tão pouco de Pm = 1 a 746, a diferença de
regime entre protoestrela de nêutrons e anã branca importa menos do que a §6.14
temia — os coeficientes parasitas têm boa chance de transferir.

### Aferições que passaram

- Maxwell/Reynolds entre 4.8 e 5.7, acima da faixa 3–5 de fluxo nulo mas
  esperado para fluxo líquido.
- Saturação em ~2 órbitas nas cinco, com 30 órbitas de turbulência sustentada.
- Nenhuma morreu — coerente com §6.19: estamos acima do Pm crítico.

## 6.22 α ≈ 8×10⁻⁴, e o tempo de frenagem bate com a duração do run

`reports/report_shearing_box.tex` (Relatório III, 8 páginas).

### α, o número que faltava nomear

Pm e Rm dizem em que regime o plasma está. **α é a resposta que se quer:**
eficiência de transporte, tensão turbulenta em unidades da pressão,

    α = W_ϖφ/P = [⟨ρ δv_ϖ δv_φ⟩ − ⟨B_ϖ B_φ⟩/4π]/P

entrando na dinâmica como viscosidade efetiva ν_t = α c_s H.

**Caixa incompressível não dá α direto** — não há pressão, o denominador não
existe. Ela dá W/B0², e a conversão precisa do campo e da pressão da estrela:

    α = (2/β)·(W/B0²),      β = P/(B_z²/8π)

Logo α é pequeno na nossa estrela **não porque a MRI seja ineficiente, mas
porque o campo vertical é fraco perante a pressão.**

### A conta que incomoda

Com B_pol,rms = 7×10¹⁰ G tardio (logo B_z ≈ 5×10¹⁰), P = 7.8×10²⁴ erg/cm³
degenerada a ρ = 4.8×10⁷, e W/B0² = 32 do topo da varredura:

| | |
|---|---|
| β no campo vertical | 8×10⁴ |
| **α implicado** | **8×10⁻⁴** |
| ν_t = α c_s H, H ~ R_eq | 7.6×10¹⁴ cm²/s |
| **tempo de frenagem R²/ν_t** | **~130 s** |
| **duração do nosso run mais longo** | **78 s** |

**O tempo de frenagem é da ordem da duração do run.** Se a MRI operasse no nível
que esta caixa mede, com o campo que a estrela de fato tem, ela teria apagado
fração substancial da rotação diferencial dentro da janela simulada. **Nossos
runs mostram o perfil acentuando 22.6%.**

É a tensão posta quantitativamente, e é por isso que a MRI não resolvida não é
ressalva acadêmica. Coerente em sinal com Miravet-Tenés+2025, que ligam subgrade
numa estrela de nêutrons diferencialmente rotativa e acham achatamento.

### O que enfraquece a estimativa

Envelope, e declarada como tal no relatório:

- **H ~ R_eq é grosseiro.** A divergência da tensão age na escala em que campo e
  rotação variam; H efetivo maior alonga a frenagem.
- **W/B0² = 32 é medido a Pm = 16 e q = 1.5**, não no nosso Pm nem na nossa
  faixa de q.
- **α ∝ B_z², e o campo usado é o tardio já decaído** — exatamente a quantidade
  que o Relatório II mostra não convergida.

Fator dez para qualquer lado não surpreenderia. **Fator 10⁴, que é o que "a
rotação diferencial está a salvo" exigiria, surpreenderia muito.**


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
