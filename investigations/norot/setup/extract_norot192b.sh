#!/bin/bash
#SBATCH --job-name=rafael-wdmag-extr
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --partition=ep-short
#SBATCH --time=0-08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
# Extrai E_tor, E_pol, massa, raios e rotacao dos 135 plotfiles do run sem
# rotacao com reposicao de ambiente (45194, t = 0..60 s). Mesma receita do
# extract_norot.sh: corrige os cabecalhos inf antes (DIARIO 10.2), rho_cut 1e5.
# Submeter excluindo os nos com jobs seus (epilogo, DIARIO 30).
set -uo pipefail
source /opt/ohpc/pub/apps/spack/0.22.2/share/spack/setup-env.sh
H=$(spack find -l openmpi@5.0.3 %gcc@9.4.0 2>/dev/null | grep -oE "^[a-z0-9]{7}" | head -1)
spack load "/$H"
MPI_LIB="$(readlink -f "$(dirname "$(command -v mpirun)")/../lib")"
export LD_LIBRARY_PATH="$MPI_LIB:${LD_LIBRARY_PATH:-}"
cd "$HOME/wd-mag/runs/dir_norot192b"
bash "$HOME/WD_MAG/tools/patch_plotfile_inf.sh" plt[0-9]* >/dev/null 2>&1 || true
OUT="$HOME/wd-mag/runs/bt_bp_norot192b.csv"
echo "# fbtbp sobre dir_norot192b (45194, fill_ambient_bc), todos os plotfiles." > "$OUT"
echo "# rho_cut = 1.0e5 g/cm^3 -- a estrela, nao a caixa." >> "$OUT"
for p in $(ls -d plt[0-9]* | sort -V); do
    mpirun -np "$SLURM_NTASKS" "$HOME/wd-mag/toolbuild/fbtbp3d.gnu.MPI.ex" "$p" 2>/dev/null | grep -E "^ +[0-9]" >> "$OUT"
done
echo "linhas: $(grep -vc "^#" "$OUT")"
