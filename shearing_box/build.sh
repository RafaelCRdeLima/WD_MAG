#!/usr/bin/env bash
# Build FFTW3 and SNOOPY into this directory. No root, nothing outside here.
#
#   ./build.sh [problem]      default problem: mri
#
# Why FFTW is built from source: Ubuntu ships libfftw3-double3 (runtime) but
# not libfftw3-dev, so there is no fftw3.h and no .so symlink. Installing the
# dev package needs root; building here does not, and the same tree is what we
# would need on the cluster anyway.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEPS="$ROOT/deps"
SRC="$ROOT/src/snoopy-v6.0-official"
PROBLEM="${1:-mri}"
FFTW_VER=3.3.10

# ---------------------------------------------------------------- FFTW3
if [ ! -f "$DEPS/include/fftw3.h" ]; then
    echo "== building FFTW $FFTW_VER =="
    mkdir -p "$DEPS/build"
    cd "$DEPS/build"
    [ -f "fftw-$FFTW_VER.tar.gz" ] || \
        curl -sL -o "fftw-$FFTW_VER.tar.gz" "http://www.fftw.org/fftw-$FFTW_VER.tar.gz"
    [ -d "fftw-$FFTW_VER" ] || tar xzf "fftw-$FFTW_VER.tar.gz"
    cd "fftw-$FFTW_VER"
    ./configure --prefix="$DEPS" --enable-openmp --enable-threads \
                --enable-shared --enable-sse2 --enable-avx2 >configure.log 2>&1
    make -j"$(nproc)" >make.log 2>&1
    make install >install.log 2>&1
else
    echo "== FFTW already present in deps/ =="
fi

# ---------------------------------------------------------------- SNOOPY
echo "== building SNOOPY, problem '$PROBLEM' =="
cd "$SRC"
make clean >/dev/null 2>&1 || true

# The include path goes in CFLAGS, NOT CPPFLAGS. configure detects fftw3.h via
# CPPFLAGS but its Makefile.in never substitutes CPPFLAGS, so a CPPFLAGS-only
# invocation configures cleanly and then fails every compile with
# "fftw3.h: No such file or directory".
#
# -malign-double is dropped from the stock flags. It is an i386 option that on
# x86-64 changes struct layout against the SysV ABI, and FFTW here is built
# without it; keeping it would mean linking two different ABIs.
./configure --with-problem="$PROBLEM" --enable-openmp \
    CFLAGS="-O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -I$DEPS/include" \
    LDFLAGS="-L$DEPS/lib -Wl,-rpath,$DEPS/lib" >configure.log 2>&1
make -j"$(nproc)" >make.log 2>&1

echo "== done: $SRC/snoopy =="
# Do NOT probe the binary for a version string. SNOOPY parses no arguments and
# ignores unknown ones, so `snoopy --version` does not print and exit -- it
# starts a full simulation in the current directory. Cost one wedged build.
grep -o 'The Snoopy code v[0-9.]*' "$SRC/src/snoopy.c" | head -1

cat <<'EOF'

To run:
    mkdir -p runs/<name>/data          # REQUIRED -- see README, the code
    cd runs/<name>                     # segfaults without it
    cp ../../src/snoopy-v6.0-official/src/problem/mri/snoopy.cfg .
    OMP_NUM_THREADS=8 ../../src/snoopy-v6.0-official/snoopy
EOF
