# Ambiente do sci-com (UFES). Carregado por todo job e por qualquer sessao
# interativa. Spack para o python 3.11 (o do sistema e' 3.6.8 e nao serve),
# venv proprio para numpy/scipy.
source /opt/ohpc/pub/apps/spack/0.22.2/share/spack/setup-env.sh
spack load python@3.11.7 %gcc@11.4.0 2>/dev/null || spack load python@3.11.7
export WDMAG=${WDMAG:-$HOME/WD_MAG}
export PY=$WDMAG/.venv/bin/python3

# Uma thread por tarefa. Sem isto cada tarefa do array abre um pool de BLAS do
# tamanho do no' (16 no ry-short), dezenas de tarefas caem no mesmo no' e o OOM
# killer ceifa o array inteiro em 20 s -- foi o que aconteceu com o 44257.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
