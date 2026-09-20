#!/bin/bash
#SBATCH --job-name=rafael-wdmag-resolution
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=ry-short
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Checagem de convergencia das duas manchetes, em 257x257 contra os 129x129 de
# tudo que foi rodado ate' agora. E' o que um referee pede primeiro, e este
# projeto ja' foi mordido por isso: a taxa de crescimento m=1 sobe 47% de 192^3
# para 256^3 e o residuo de E_mag muda por fator 25 (DIARIO 3.2).
#
# Medido antes de submeter: o custo sobe 1.8x de 129 para 257, nao 4x -- o
# solve nao e' dominado pela parte O(N^2). Um ponto isolado ja' mostrou beta_min
# convergindo a 0.03% em m = 0 e m = -1, que sao as formas das 40 configuracoes
# que passam os quatro portoes, e a 2.4% em m = -2, que nao carrega manchete
# nenhuma. Esta rodada confirma isso na fronteira, onde beta_min ~ 1 e a
# sensibilidade pode ser outra.
#
# Um job, sem array (ver job_tt_sweep.sh).

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source "$HOME/WD_MAG/cluster/scicom/env.sh"
mkdir -p logs
cd "$WDMAG"
mkdir -p sweep_logs res257
echo "no: $(hostname)  cpus: ${SLURM_CPUS_PER_TASK}"

echo "=== varredura de campo, 50 pontos em 257x257 ==="
seq 0 49 | xargs -P 8 -I{} sh -c \
  "$PY investigations/twisted_torus_sweep.py --index {} --nr 257 \
     --out res257/tt_sweep_{}.csv > sweep_logs/r257_{}.log 2>&1 \
   || echo 'FALHOU {}'"
echo "csv: $(ls res257/tt_sweep_*.csv 2>/dev/null | wc -l)/50"

echo "=== bisseccao do maximo de massa em 257x257 ==="
$PY investigations/max_mass_frontier.py --nr 257
