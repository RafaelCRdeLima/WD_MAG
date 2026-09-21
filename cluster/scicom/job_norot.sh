#!/bin/bash
#SBATCH --job-name=rafael-wdmag-norot
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --partition=ep-short
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Controle sem rotacao, para a pergunta da Laura. Uso: sbatch job_norot.sh norot192
#
# Ver o cabecalho de inputs.norot192 para o que se pergunta e o que muda por
# construcao.
#
# ep-short e 7 dias porque este e' o run longo: stop_time = 60 s, o mesmo
# alcance do rot192. No lovelace isso exigia cadeia de janelas de 3 h e foi o
# que matou a campanha HZ (DIARIO 11); aqui cabe num job so'.

set -euo pipefail
TAG="${1:?passe tt192 ou tt192ctl}"
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

source /opt/ohpc/pub/apps/spack/0.22.2/share/spack/setup-env.sh
H=$(spack find -l openmpi@5.0.3 %gcc@9.4.0 2>/dev/null | grep -oE "^[a-z0-9]{7}" | head -1)
spack load "/$H"
source /opt/rh/gcc-toolset-13/enable

# OBRIGATORIO. O binario e' linkado contra o openmpi 5.0.3 do Spack, mas o
# ambiente do OpenHPC poe /opt/ohpc/pub/mpi/openmpi4-gnu9/4.1.1/lib no
# LD_LIBRARY_PATH -- e as duas versoes exportam o mesmo soname, libmpi.so.40.
# Linkers modernos gravam RUNPATH e nao RPATH, e LD_LIBRARY_PATH vence RUNPATH,
# entao o binario carregava a 4.1.1 enquanto o mpirun era a 5.0.3. O sintoma e'
# "An error occurred in MPI_Init_thread on a NULL communicator", sem mais
# detalhe nenhum. Confirmado com ldd, nao adivinhado.
MPI_LIB="$(readlink -f "$(dirname "$(command -v mpirun)")/../lib")"
export LD_LIBRARY_PATH="$MPI_LIB:${LD_LIBRARY_PATH:-}"

PROB="$HOME/wd-mag/Castro/Exec/science/wd_scf_stability"
RUN="$HOME/wd-mag/runs/dir_$TAG"
mkdir -p "$RUN"
cd "$RUN"

# Portao antes de gastar a janela (regra da DIARIO 10: verificar no DESTINO).
# A campanha HZ perdeu tres janelas porque a sincronia jurava ter entregue um
# input que nao estava la'.
MODEL=$(grep -oP 'problem\.model_name\s*=\s*"\K[^"]+' "$PROB/inputs.$TAG")
for f in "$PROB/Castro3d.gnu.MPI.ex" "$PROB/inputs.$TAG" "$HOME/wd-mag/models/$MODEL"; do
    [ -f "$f" ] || { echo "FALTA no destino: $f"; exit 1; }
    echo "ok: $f"
done
ln -sf "$PROB/Castro3d.gnu.MPI.ex" .
ln -sf "$PROB/inputs.$TAG" .
ln -sf "$HOME/wd-mag/models/$MODEL" .

# Portao de ABI: se libmpi nao vier do mpirun em uso, o run morre em
# MPI_Init_thread depois de a fila ja' ter sido gasta.
if ! ldd ./Castro3d.gnu.MPI.ex | grep -q "libmpi.so.40 => $MPI_LIB"; then
    echo "ABI ERRADA: libmpi nao resolve para $MPI_LIB"
    ldd ./Castro3d.gnu.MPI.ex | grep libmpi
    exit 1
fi
echo "libmpi ok: $MPI_LIB"
echo "no: $(hostname)  ranks: $SLURM_NTASKS  modelo: $MODEL"
mpirun -np "$SLURM_NTASKS" ./Castro3d.gnu.MPI.ex "inputs.$TAG" 2>&1 | tail -40
echo "--- como terminou ---"
ls -d plt* 2>/dev/null | tail -3
