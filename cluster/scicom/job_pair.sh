#!/bin/bash
#SBATCH --job-name=rafael-wdmag-pair
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --partition=z3-long
#SBATCH --time=10-00:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Par controlado rotacao / sem rotacao. Uso:
#   sbatch job_pair.sh norot192b
#   sbatch job_pair.sh rot192b
#
# --mem e' OBRIGATORIO. Sem ele o SLURM concede a memoria inteira do no a cada
# job, e tres jobs reivindicando 240 GB de um no de 240 GB caem no OOM killer:
# foi como 45148 morreu com SIGKILL aos 13 s, dividindo o node1 com um job de
# outro usuario. O 44987 so' sobreviveu por estar sozinho no node2.
#
# =====================================================================
# PERIGO DO CLUSTER: NUNCA rode dois jobs SEUS no mesmo no'.
#
# O /etc/slurm/slurm.epilog.clean tem uma protecao que NAO funciona aqui:
#
#     job_list=`squeue --noheader --format=%A --user=$SLURM_UID \
#                      --node=localhost`
#     for job_id in $job_list; do
#         if [ $job_id -ne $SLURM_JOB_ID ]; then exit 0; fi
#     done
#     pkill -KILL -U $SLURM_UID
#
# 'squeue --node=localhost' devolve "error: Invalid node name localhost" neste
# cluster, entao job_list sai VAZIA, o laco nunca roda, o exit 0 nunca acontece
# e o script cai direto no pkill de TODOS os seus processos no no'.
#
# Consequencia: quando qualquer job seu termina num no', todos os seus outros
# jobs naquele no' morrem com SIGKILL. Verificado duas vezes --
# 45150 morreu no instante em que 45151 foi cancelado, os dois no node1, e e'
# tambem o que matava as ondas dos arrays do 44307 (DIARIO 18, onde eu nao
# consegui fechar o diagnostico).
#
# Submeta sempre com --exclude=<nos onde voce ja' tem job>. E note que isto
# vale para TODOS os seus jobs, nao so' os deste projeto.
# =====================================================================
#
# z3-long e nao ep-short: o ep-short tem o node30 drenado e fila de varios
# usuarios com prioridade acima da nossa; o job 44958 ficou PENDING um dia
# inteiro (DIARIO 28). O z3-long tem 14 dias e costuma ter nucleos livres.

set -euo pipefail
TAG="${1:?passe norot192b ou rot192b}"
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

source /opt/ohpc/pub/apps/spack/0.22.2/share/spack/setup-env.sh
H=$(spack find -l openmpi@5.0.3 %gcc@9.4.0 2>/dev/null | grep -oE "^[a-z0-9]{7}" | head -1)
spack load "/$H"
source /opt/rh/gcc-toolset-13/enable

# Obrigatorio: as duas versoes de openmpi exportam libmpi.so.40 e o
# LD_LIBRARY_PATH do OpenHPC vence o RUNPATH do binario. Sem isto o run morre
# em MPI_Init_thread sem dizer por que.
MPI_LIB="$(readlink -f "$(dirname "$(command -v mpirun)")/../lib")"
export LD_LIBRARY_PATH="$MPI_LIB:${LD_LIBRARY_PATH:-}"

PROB="$HOME/wd-mag/Castro/Exec/science/wd_scf_stability"
RUN="$HOME/wd-mag/runs/dir_$TAG"
mkdir -p "$RUN"; cd "$RUN"

MODEL=$(grep -oP 'problem\.model_name\s*=\s*"\K[^"]+' "$PROB/inputs.$TAG")
for f in "$PROB/Castro3d.gnu.MPI.ex" "$PROB/inputs.$TAG" "$HOME/wd-mag/models/$MODEL"; do
    [ -f "$f" ] || { echo "FALTA no destino: $f"; exit 1; }
    echo "ok: $f"
done
ldd "$PROB/Castro3d.gnu.MPI.ex" | grep -q "libmpi.so.40 => $MPI_LIB" || {
    echo "ABI ERRADA: libmpi nao resolve para $MPI_LIB"; exit 1; }
ln -sf "$PROB/Castro3d.gnu.MPI.ex" .; ln -sf "$PROB/inputs.$TAG" .
ln -sf "$HOME/wd-mag/models/$MODEL" .
echo "no: $(hostname)  ranks: $SLURM_NTASKS  modelo: $MODEL"

# SEM pipe. O job_tt_probe.sh canalizava o mpirun por 'tail -40', entao nao
# havia saida de passo enquanto o run vivia -- nem dt, nem taxa de retry, nem
# avisos de densidade -- e a curva de dt do 44987 teve de ser reconstruida
# pelos timestamps dos plotfiles (DIARIO 28). O log completo vale os megabytes.
mpirun -np "$SLURM_NTASKS" ./Castro3d.gnu.MPI.ex "inputs.$TAG"
echo "=== terminou: $(ls -d plt* 2>/dev/null | tail -1) ==="
