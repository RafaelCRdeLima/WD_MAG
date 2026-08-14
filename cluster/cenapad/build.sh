#!/usr/bin/env bash
# Build on Lovelace's login node. The CENAPAD guide is explicit that lovelace
# is where you compile -- the machine named "frontend" (cenapad.unicamp.br) is
# a connection server only and refuses module/qstat outright.
#
# Toolchain: OpenMPI 5.0.8 built against GCC 15.2.0, which is the module
# system's own default.
#
# GCC 12.2 was tried first, on the reasoning that GCC 15 is very new for AMReX.
# That was wrong: Castro 26.07's main.cpp includes <format>, and libstdc++ only
# gained it in GCC 13. The system gcc, 8.5.0, is far too old for anything here.
# Do not use openmpi/5.0.8-gcc-15.2.0-TESTE-NAO-USE, which the site marks as
# not-for-use.
#
# -j8 rather than -j32: this runs on a shared login node.
set -euo pipefail

module purge
module load openmpi/5.0.8-gcc-15.2.0
module list

# Castro's build scripts have "#!/usr/bin/env python" shebangs, and RHEL8 --
# which lovelace runs -- ships only python3. Without a shim the build dies at
# check_network.py with "/usr/bin/env: 'python': No such file or directory".
# A symlink on PATH is enough and is less invasive than loading anaconda,
# which would shadow other things in the environment.
if ! command -v python >/dev/null 2>&1; then
    PY3="$(command -v python3 || true)"
    [ -n "$PY3" ] || { echo "neither python nor python3 found"; exit 1; }
    mkdir -p "$HOME/bin"
    ln -sf "$PY3" "$HOME/bin/python"
    export PATH="$HOME/bin:$PATH"
    echo "shim:   $HOME/bin/python -> $PY3"
fi
echo "python: $(python --version 2>&1)"
echo "gcc:    $(gcc --version | head -1)"
echo "mpicxx: $(command -v mpicxx || echo MISSING)"
[ -n "$(command -v mpicxx)" ] || { echo "no mpicxx after module load"; exit 1; }

cd "$(dirname "${BASH_SOURCE[0]}")/../../Castro/Exec/science/wd_scf_stability" 2>/dev/null \
  || cd "${1:?pass the problem directory}"

# Optional first argument: an EOS_DIR to build instead of the one in the
# GNUmakefile. The executable is renamed to carry it, because Castro names the
# binary by dimension and compiler only -- no EOS -- so two builds collide.
#
#     ./build.sh                # ztwd, the default, -> Castro3d.gnu.MPI.ex
#     ./build.sh helmholtz      #      -> Castro3d.gnu.MPI.helm.ex
#
# The ztwd binary is moved aside first and put back afterwards, because
# `make realclean` deletes *.ex and rebuilding it costs another ten minutes.
EOS_ARG="${1:-}"
SAVED=""
if [ -n "$EOS_ARG" ] && [ -f Castro3d.gnu.MPI.ex ]; then
    SAVED="ztwd_saved.bin"          # not *.ex, so realclean leaves it alone
    mv Castro3d.gnu.MPI.ex "$SAVED"
    echo "set aside: Castro3d.gnu.MPI.ex -> $SAVED"
fi

make realclean >/dev/null 2>&1 || true
make -j8 COMP=gnu USE_MPI=TRUE ${EOS_ARG:+EOS_DIR="$EOS_ARG"} 2>&1 | tail -25

if [ -n "$EOS_ARG" ]; then
    [ -f Castro3d.gnu.MPI.ex ] || { echo "BUILD FAILED -- $SAVED left in place"; exit 1; }
    # Verify the EOS actually went in. Passing EOS_DIR that the makefile
    # ignores would produce a binary that silently runs the wrong physics,
    # which is the one failure this whole exercise cannot afford.
    if [ "$EOS_ARG" = "helmholtz" ]; then
        n=$(strings Castro3d.gnu.MPI.ex | grep -ci helm_table || true)
        echo "helm_table strings in the new binary: $n"
        [ "$n" -gt 0 ] || { echo "EOS_DIR was IGNORED -- this is a ztwd binary"; exit 1; }
        mv Castro3d.gnu.MPI.ex Castro3d.gnu.MPI.helm.ex
        echo "built: Castro3d.gnu.MPI.helm.ex"
    fi
    [ -n "$SAVED" ] && mv "$SAVED" Castro3d.gnu.MPI.ex && echo "restored: Castro3d.gnu.MPI.ex"
fi
ls -la Castro3d.gnu.MPI*.ex
ls -la ./*.ex
