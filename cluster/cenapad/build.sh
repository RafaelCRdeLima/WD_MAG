#!/usr/bin/env bash
# Build on Lovelace's login node. The CENAPAD guide is explicit that lovelace
# is where you compile -- the machine named "frontend" (cenapad.unicamp.br) is
# a connection server only and refuses module/qstat outright.
#
# Toolchain, read off `module avail` on the machine itself: OpenMPI 5.0.6 built
# against GCC 12.2.0. The system gcc is 8.5.0, which cannot build Castro (it
# needs C++20). openmpi/5.0.8-gcc-15.2.0 is the module system's default but
# GCC 15 is very new for AMReX; 12.2 is the safer pair.
#
# -j8 rather than -j32: this runs on a shared login node.
set -euo pipefail

module purge
module load openmpi/5.0.6-gcc-12.2.0
module list

echo "gcc:    $(gcc --version | head -1)"
echo "mpicxx: $(command -v mpicxx || echo MISSING)"
[ -n "$(command -v mpicxx)" ] || { echo "no mpicxx after module load"; exit 1; }

cd "$(dirname "${BASH_SOURCE[0]}")/../../Castro/Exec/science/wd_scf_stability" 2>/dev/null \
  || cd "${1:?pass the problem directory}"

make -j8 COMP=gnu USE_MPI=TRUE 2>&1 | tail -25
ls -la ./*.ex
