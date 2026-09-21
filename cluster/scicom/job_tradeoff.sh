#!/bin/bash
#SBATCH --job-name=rafael-wdmag-tradeoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=ry-short
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source "$HOME/WD_MAG/cluster/scicom/env.sh"
mkdir -p logs
cd "$WDMAG"
echo "no: $(hostname)"
$PY investigations/tradeoff_curve.py --nr 257
