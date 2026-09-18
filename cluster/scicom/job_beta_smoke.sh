#!/bin/bash
#SBATCH --job-name=rafael-wdmag-beta-smoke
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=ry-short
#SBATCH --time=01:00:00
#SBATCH --output=%x-slurm_job-%j.out
#SBATCH --error=%x-slurm_job-%j.err

# Smoke test: roda o portao de beta que ja' rodou na estacao e compara os
# numeros. Se beta_min nao bater com 6.05e-3, a maquina nao esta' reproduzindo
# a estacao e nada mais vale a pena. Portao antes de qualquer coisa cara.

set -euo pipefail
echo "$SLURM_SUBMIT_DIR"
cd "$SLURM_SUBMIT_DIR"
source "$HOME/WD_MAG/cluster/scicom/env.sh"

echo "no: $(hostname)   python: $($PY --version)"
cd "$WDMAG"
time $PY investigations/twisted_torus_beta.py
echo done
