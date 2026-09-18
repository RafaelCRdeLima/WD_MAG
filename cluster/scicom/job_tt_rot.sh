#!/bin/bash
#SBATCH --job-name=rafael-wdmag-tt-rot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=ry-short
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Campo + rotacao: 75 pontos (rho_c, omega_frac, K_frac). Um job, sem array
# -- ver o cabecalho do job_tt_sweep.sh para por que.
#
# Aceita indices na linha de comando para rodar um subconjunto:
#   sbatch job_tt_rot.sh 42        (portao)
#   sbatch job_tt_rot.sh           (a grade toda)

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source "$HOME/WD_MAG/cluster/scicom/env.sh"
mkdir -p logs
cd "$WDMAG"
mkdir -p sweep_logs
echo "no: $(hostname)  cpus: ${SLURM_CPUS_PER_TASK}"

if [ "$#" -gt 0 ]; then IDX=$(printf "%s\n" "$@"); else IDX=$(seq 0 74); fi

echo "$IDX" | xargs -P 8 -I{} sh -c \
  "$PY investigations/twisted_torus_rotating.py --index {} \
     > sweep_logs/rot_{}.log 2>&1 || echo 'FALHOU {}'"

echo "csv: $(ls investigations/tt_rot_[0-9][0-9][0-9].csv 2>/dev/null | wc -l)"
echo "=== pontos que passam os quatro portoes ==="
# || true: com set -e, um grep sem match derruba o job DEPOIS do trabalho
# todo feito -- foi o que marcou o 44362 como FAILED com os dados intactos.
grep -h "ALL FOUR GATES" investigations/tt_rot_*.csv 2>/dev/null | wc -l || true
echo "=== os mais pesados que seguram o campo ==="
grep -h -v "^#" investigations/tt_rot_*.csv 2>/dev/null \
  | awk -F, '$7+0 >= 1.0 {print $3, $4, $9, $12}' | sort -rn | head -10 || true
echo "(colunas: M_msun  T/|W|  E_tor/E_pol  veredito)"
