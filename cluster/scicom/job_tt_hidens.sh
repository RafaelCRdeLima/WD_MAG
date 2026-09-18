#!/bin/bash
#SBATCH --job-name=rafael-wdmag-tt-hidens
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=ry-short
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Fecha a extrapolacao que a varredura principal deixou aberta. A fronteira
# beta_min >= 1 sobe com rho_c -- 1.3227 Msun em 5e8, 1.4186 em 5e9 -- e a
# grade parou em 5e9, enquanto a neutronizacao (mu_e = 2) so' entra em
# 1.94e10. Se a fronteira cruzar 1.44 acima da grade, a afirmacao "um toro
# torcido que a estrela segura nao e' super-Chandrasekhar" cai. Medir.
#
# Um job, sem array -- ver o cabecalho do job_tt_sweep.sh.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source "$HOME/WD_MAG/cluster/scicom/env.sh"
mkdir -p logs
cd "$WDMAG"
mkdir -p sweep_logs
echo "no: $(hostname)  cpus: ${SLURM_CPUS_PER_TASK}"

# 1.94e10 e' o limiar; 1.8e10 fica logo abaixo dele.
for RC in 8e9 1.2e10 1.6e10 1.8e10; do
  for KF in 0.05 0.08 0.11 0.15 0.20 0.30; do
    echo "$RC $KF"
  done
done | xargs -P 8 -n 2 sh -c \
  "$PY investigations/twisted_torus_sweep.py --rho-c \$0 --k-frac \$1 \
     --tag \$0_\$1 > sweep_logs/hi_\$0_\$1.log 2>&1 || echo 'FALHOU \$0 \$1'"

echo "=== fronteira acima da grade ==="
grep -h -v "^#" investigations/tt_sweep_hi_*.csv 2>/dev/null | \
  awk -F, '$6+0 >= 1.0 {print $4, $5, $6}' | sort -rn | head -5
echo "(colunas: M_msun  E_tor/|W|  beta_min)"
