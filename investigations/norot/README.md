# Sem rotação — dados, saídas e figuras

Tudo do run sem rotação 192³ num lugar só. A leitura está no DIARIO §28, §29, §31 e **§32** (60 s).

| | primeiro run | run completo |
|---|---|---|
| job | 44987 | **45194** (+ 45832, retomada que saiu sem rodar) |
| diretório no sci-com | `~/wd-mag/runs/dir_norot192` | `~/wd-mag/runs/dir_norot192b` |
| ambiente | sem reposição | `castro.fill_ambient_bc = 1` |
| até | t = 26.1 s (morreu: ambiente evacuado, §28) | **t = 60.0 s**, passo 26735, 4 d 19 h no node1 |
| plotfiles | 184 (250 GB) | 135 (`plot_int = 200`, 210 GB), `chk26735` final |

Estrela: 1.409 M⊙, ρ_c = 3.0×10⁹ g/cm³, EOS `ztwd` (barotrópica), sem rotação; toroidal interior + dipolo exterior, como o `rot192`.

## Figuras — `plots/` (pdf e png)

Geradas por `plot_norot.py` a partir de `data/` (e de `../bt_bp_192*.csv` para o rotante):

```bash
scf/.venv/bin/python3 investigations/norot/plot_norot.py
```

| figura | o que mostra |
|---|---|
| `norot_field` | (a) E_tor/E_pol contra o `rot192`: cruza 1 perto de 3 s e fica em 0.19–0.48 até 60 s, no ramo misto. (b) E_mag/E_mag(0): sobe um pouco até 26 s e **decai com e-folding de 28 s** de 26 a 60 s; aos 60 s restam 3.4%, contra 0.09% com rotação |
| `norot_star` | (a) ρ_max: a pulsação radial **cresce**, amplitude 0.15 → 0.50, e-folding ~35–41 s, P = 1.69 s. (b) R_pol/R_eq oscila entre 0.74 e 1.06. (c) massa: +0.18% na caixa, +0.15% acima de 10⁵ |
| `norot_timestep` | dt: os dois runs caem juntos a partir de 16 s; o primeiro estagna aos 26 s, o 45194 atravessa e volta ao patamar em t ≈ 43 s. 454 avanços rejeitados, todos entre 16 e 44 s |
| `norot_energy_budget` | E_int cresce 12× (inerte sob `ztwd`), \|E_grav\| estável, E_kin ~10⁴⁹ com a pulsação, E_mag |

## Dados — `data/`

| arquivo | origem | conteúdo |
|---|---|---|
| `bt_bp_norot192b.csv` | `fbtbp` nos 135 plotfiles (job 46024) | t, E_tor, E_pol, razões, B de pico, massa e raios acima de 10⁵, L_z, Ω por região — 24 colunas nomeadas |
| `bt_bp_norot192.csv` | `fbtbp`, primeiro run (1 a cada 4 plotfiles) | idem, 0–26 s |
| `grid_diag_norot192b.csv`, `grid_diag_norot192.csv` | `grid_diag.out`, 1 a cada 10 passos | massa, momentos, L, E_kin, E_int, E_grav, E_tot, centro de massa, T_max, ρ_max |
| `dt_norot192b.csv`, `dt_norot192.csv` | `amr_diag.out`, 1 a cada 10 passos | dt, subciclos, tempo de parede por passo |
| `retries_norot192b.csv` | log completo do 45194 | avanços rejeitados por passo e motivo (densidade inválida / validade do passo) |

## Saídas brutas — `raw/`

- `fbtbp_raw_norot192b.txt`, `fbtbp_raw_norot192.txt`: saída direta do `fbtbp` (largura fixa; números negativos encostam no vizinho, leia com regex).
- `jobs_sacct.txt`: 44987, 45194, 45832, 46024 — nós, horários, estado.

Ficaram no sci-com, grandes demais para o repositório: os plotfiles e checkpoints (`dir_norot192b`, 210 GB), o log completo do passo a passo (`~/WD_MAG/cluster/scicom/logs/rafael-wdmag-pair-45194.out`, 36 MB) e os `.out` de diagnóstico completos.

## Configuração — `setup/`

`inputs.norot192b`, `job_pair.sh` (submissão), `job_pair_restart.sh` (retomada encadeada por `afterany`), `extract_norot192b.sh` (extração). Executável: `Castro/Exec/science/wd_scf_stability/Castro3d.gnu.MPI.ex` de 18/09/2026, `EOS_DIR = ztwd`, `USE_MHD = TRUE`.
