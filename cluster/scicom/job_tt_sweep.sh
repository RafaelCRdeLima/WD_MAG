#!/bin/bash
#SBATCH --job-name=rafael-wdmag-tt-sweep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=ry-short
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# UM job, nao um array. A versao em array (--array=0-49%12) perdia a onda
# inteira a cada rodizio: os fins do 44307 se agrupam em ondas de 12 -- o
# tamanho do throttle -- e dentro de cada onda TODAS terminam no mesmo segundo,
# as que completam e as que morrem com SIGKILL. Nao era memoria (pico de RSS
# medido: 0.10 GB), nao era no' doente, nao era preempcao (PreemptMode=OFF).
# Era o mecanismo de array. Um job unico com paralelismo interno nao o toca.
#
# 50 pontos (rho_c, K_tor) em 8 processos de uma thread cada.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source "$HOME/WD_MAG/cluster/scicom/env.sh"
mkdir -p logs
cd "$WDMAG"
echo "no: $(hostname)  cpus: ${SLURM_CPUS_PER_TASK}"

mkdir -p sweep_logs
seq 0 49 | xargs -P 8 -I{} sh -c \
  "$PY investigations/twisted_torus_sweep.py --index {} > sweep_logs/tt_{}.log 2>&1 \
   || echo 'FALHOU {}'"

echo "csv escritos: $(ls investigations/tt_sweep_[0-9][0-9][0-9].csv 2>/dev/null | wc -l)/50"
$PY investigations/merge_tt_sweep.py
