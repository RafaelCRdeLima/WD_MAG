#!/usr/bin/env bash
# Build no sci-com (UFES). Diferencas de fundo em relacao ao lovelace:
#
#  - O gcc do sistema e' 9.4.0 e o do Spack para em 11.4. **Castro 26.07
#    inclui <format> no main.cpp, que so' existe no libstdc++ a partir do
#    GCC 13** -- o mesmo motivo que derrubou o GCC 12.2 no lovelace. A saida
#    aqui e' o gcc-toolset-13 do Rocky (/opt/rh/gcc-toolset-13), que traz
#    g++ 13.3.1 e compila <format> com -std=c++20. Verificado antes deste
#    script existir.
#  - Os openmpi do Spack foram compilados com gcc 8.5/9.4, nao com 13. O
#    wrapper mpicxx so' acrescenta flags de include/lib; o compilador de
#    verdade vem de OMPI_CXX. A biblioteca MPI e' C, entao a ABI fecha.
#    Testado: mpicxx -std=c++20 com <format> compila e roda.
#  - Sem modules tradicionais; o Spack e' quem carrega. O openmpi precisa ser
#    escolhido por hash porque ha' seis instalacoes de openmpi@5.0.3.
#
# Uso:  bash build.sh [EOS_DIR]     (default: o do GNUmakefile)
set -euo pipefail

source /opt/ohpc/pub/apps/spack/0.22.2/share/spack/setup-env.sh
H=$(spack find -l openmpi@5.0.3 %gcc@9.4.0 2>/dev/null | grep -oE "^[a-z0-9]{7}" | head -1)
[ -n "$H" ] || { echo "nao achei openmpi %gcc@9.4.0 no spack"; exit 1; }
spack load "/$H"

source /opt/rh/gcc-toolset-13/enable
export OMPI_CXX=g++ OMPI_CC=gcc OMPI_FC=gfortran
echo "g++:    $(g++ --version | head -1)"
echo "mpicxx: $(command -v mpicxx) -> $(mpicxx --showme:command)"

# Castro chama scripts com shebang "#!/usr/bin/env python"; o Rocky 8 so' tem
# python3. Mesmo shim do lovelace.
if ! command -v python >/dev/null 2>&1; then
    mkdir -p "$HOME/bin"; ln -sf "$(command -v python3)" "$HOME/bin/python"
    export PATH="$HOME/bin:$PATH"
fi
echo "python: $(python --version 2>&1)"

export CASTRO_HOME="$HOME/wd-mag/Castro"
export AMREX_HOME="$HOME/wd-mag/amrex"
export MICROPHYSICS_HOME="$HOME/wd-mag/Microphysics"
for d in "$CASTRO_HOME" "$AMREX_HOME" "$MICROPHYSICS_HOME"; do
    [ -d "$d" ] || { echo "FALTA $d"; exit 1; }
done

# Os dois patches de cluster/cenapad/castro_core.patch NAO sao opcionais.
# Sem o de Castro_io.cpp o build morre em "network_rp has not been declared"
# -- a ztwd nao tem parametros de runtime e por isso nunca inclui
# extern_parameters.H, que extern_job_info_tests.H precisa. Sem o de
# Gravity.cpp a geometria de meio-shift quebra. Reaplicar a cada clone novo.
cd "$CASTRO_HOME"
if git apply --check "$HOME/WD_MAG/cluster/scicom/castro_core.patch" 2>/dev/null; then
    git apply "$HOME/WD_MAG/cluster/scicom/castro_core.patch"
    echo "patches: aplicados agora"
elif git apply --reverse --check "$HOME/WD_MAG/cluster/scicom/castro_core.patch" 2>/dev/null; then
    echo "patches: ja aplicados"
else
    echo "ERRO: castro_core.patch nem aplica nem esta aplicado"; exit 1
fi

cd "$CASTRO_HOME/Exec/science/wd_scf_stability"
EOS_DIR_ARG="${1:-}"
MAKEARGS=(-j8)
if [ -n "$EOS_DIR_ARG" ]; then MAKEARGS+=("EOS_DIR=$EOS_DIR_ARG"); fi
echo "--- make ${MAKEARGS[*]}"
make "${MAKEARGS[@]}" 2>&1 | tail -25

BIN=$(ls -t Castro3d.*.ex 2>/dev/null | head -1)
[ -n "$BIN" ] || { echo "BUILD FALHOU: nenhum binario"; exit 1; }
if [ -n "$EOS_DIR_ARG" ]; then mv -f "$BIN" "${BIN%.ex}.${EOS_DIR_ARG}.ex"; BIN="${BIN%.ex}.${EOS_DIR_ARG}.ex"; fi
echo "BINARIO: $PWD/$BIN"
ls -la "$BIN"
