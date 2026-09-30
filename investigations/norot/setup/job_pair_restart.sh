#!/bin/bash
#SBATCH --job-name=rafael-wdmag-pair-rst
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --partition=z3-long
#SBATCH --time=10-00:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Retomada do par controlado a partir do ultimo checkpoint completo. Mesmo
# ambiente e mesmos cuidados do job_pair.sh (ler o cabecalho de la: --mem
# obrigatorio, NUNCA dois jobs seus no mesmo no').
#
# Uso tipico, encadeado no job que pode estourar a janela (27/09/2026, 45194
# em t = 34,7 de 60 com ETA incerta entre 02/10 e depois de 05/10):
#
#   sbatch --dependency=afterany:45194 job_pair_restart.sh norot192b
#
# afterany: roda se o 45194 acabar por TIMEOUT. Se ele chegar a stop_time, este
# job le o tempo do ultimo chk, ve que ja passou e sai sem rodar nada.
#
# Checkpoint: so diretorios chkNNNNN com Header. O AMReX grava em chkNNNNN.temp
# e renomeia no fim; um chk interrompido pelo SIGKILL do timeout fica .temp.

set -euo pipefail
TAG="${1:?passe norot192b ou rot192b}"
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

EU_UID=$(id -u); EU_NOME="${SLURM_JOB_USER:-$USER}"; NO=$(hostname -s)
MEUS=$(squeue -h -o "%i %u %N" 2>/dev/null \
       | awk -v a="$EU_UID" -v b="$EU_NOME" -v n="$NO" -v eu="$SLURM_JOB_ID" \
             '($2==a || $2==b) && index($3,n)>0 && $1!=eu {print $1}' | paste -sd,) || true
[ -n "$MEUS" ] && { echo "ABORTADO: job(s) $MEUS seu(s) ja rodam em $NO."; exit 1; }

source /opt/ohpc/pub/apps/spack/0.22.2/share/spack/setup-env.sh
H=$(spack find -l openmpi@5.0.3 %gcc@9.4.0 2>/dev/null | grep -oE "^[a-z0-9]{7}" | head -1)
spack load "/$H"
source /opt/rh/gcc-toolset-13/enable

MPI_LIB="$(readlink -f "$(dirname "$(command -v mpirun)")/../lib")"
export LD_LIBRARY_PATH="$MPI_LIB:${LD_LIBRARY_PATH:-}"

PROB="$HOME/wd-mag/Castro/Exec/science/wd_scf_stability"
RUN="$HOME/wd-mag/runs/dir_$TAG"
[ -d "$RUN" ] || { echo "FALTA $RUN: nada para retomar"; exit 1; }
cd "$RUN"

CHK=$(ls -d chk* 2>/dev/null | grep -E '^chk[0-9]+$' | sort -V \
      | while read -r c; do if [ -f "$c/Header" ]; then echo "$c"; fi; done | tail -1) || true
[ -n "$CHK" ] || { echo "Nenhum checkpoint completo em $RUN"; exit 1; }
T_CHK=$(sed -n 3p "$CHK/Header")
T_STOP=$(grep -oP '^\s*stop_time\s*=\s*\K[0-9.eE+-]+' "inputs.$TAG")
if awk -v a="$T_CHK" -v b="$T_STOP" 'BEGIN{exit !(a >= b)}'; then
    echo "$CHK ja esta em t = $T_CHK >= stop_time = $T_STOP: nada a fazer."; exit 0
fi

ldd "$PROB/Castro3d.gnu.MPI.ex" | grep -q "libmpi.so.40 => $MPI_LIB" || {
    echo "ABI ERRADA: libmpi nao resolve para $MPI_LIB"; exit 1; }
echo "no: $(hostname)  ranks: $SLURM_NTASKS  retomando de $CHK (t = $T_CHK de $T_STOP)"

mpirun -np "$SLURM_NTASKS" ./Castro3d.gnu.MPI.ex "inputs.$TAG" amr.restart="$CHK"
echo "=== terminou: $(ls -d plt* 2>/dev/null | tail -1) ==="
