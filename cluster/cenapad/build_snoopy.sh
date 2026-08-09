#!/usr/bin/env bash
# Build SNOOPY on Lovelace's login node. Same rule as the Castro build: lovelace
# is where you compile; the machine named "frontend" (cenapad.unicamp.br) is a
# connection server and refuses module outright.
#
# Run from the directory holding snoopy_v6.0.tgz:
#     bash build_snoopy.sh
#
# Everything lands under $PWD/snoopy_build. Nothing needs root and nothing is
# installed outside that tree.
#
# Why FFTW may have to be built here too: the local workstation had
# libfftw3-double3 (runtime) without libfftw3-dev (header), and the dev package
# needs root. On the cluster the situation is usually the reverse -- a module
# exists -- so this script looks for one first and only falls back to source.
#
# -j8 rather than -j32: this is a shared login node.
set -euo pipefail

ROOT="$PWD/snoopy_build"
DEPS="$ROOT/deps"
FFTW_VER=3.3.10
PROBLEM="${1:-mri}"

module purge
module load openmpi/5.0.8-gcc-15.2.0
module list 2>&1 | head -20

mkdir -p "$ROOT"

# ---------------------------------------------------------------- FFTW3
FFTW_INC=""; FFTW_LIB=""
if module avail fftw 2>&1 | grep -qi fftw; then
    echo "== an fftw module exists; listing so you can pick one =="
    module avail fftw 2>&1 | grep -i fftw
    echo "== if one of the above is double precision, load it and re-run with"
    echo "== SNOOPY_FFTW_PREFIX=<its prefix> to skip the source build."
fi

if [ -n "${SNOOPY_FFTW_PREFIX:-}" ]; then
    FFTW_INC="$SNOOPY_FFTW_PREFIX/include"
    FFTW_LIB="$SNOOPY_FFTW_PREFIX/lib"
    echo "== using FFTW from $SNOOPY_FFTW_PREFIX =="
elif [ -f "$DEPS/include/fftw3.h" ]; then
    FFTW_INC="$DEPS/include"; FFTW_LIB="$DEPS/lib"
    echo "== FFTW already built in $DEPS =="
else
    echo "== building FFTW $FFTW_VER from source =="
    mkdir -p "$DEPS/build"; cd "$DEPS/build"
    [ -f "fftw-$FFTW_VER.tar.gz" ] || \
        curl -sL -o "fftw-$FFTW_VER.tar.gz" "http://www.fftw.org/fftw-$FFTW_VER.tar.gz"
    [ -d "fftw-$FFTW_VER" ] || tar xzf "fftw-$FFTW_VER.tar.gz"
    cd "fftw-$FFTW_VER"
    ./configure --prefix="$DEPS" --enable-openmp --enable-threads \
                --enable-shared --enable-sse2 --enable-avx2 >configure.log 2>&1
    make -j8 >make.log 2>&1
    make install >install.log 2>&1
    FFTW_INC="$DEPS/include"; FFTW_LIB="$DEPS/lib"
fi

# ---------------------------------------------------------------- SNOOPY
cd "$ROOT"
if [ ! -d snoopy6 ]; then
    [ -f ../snoopy_v6.0.tgz ] || { echo "put snoopy_v6.0.tgz one level up"; exit 1; }
    tar xzf ../snoopy_v6.0.tgz
fi
cd snoopy6
make clean >/dev/null 2>&1 || true

# The include path goes in CFLAGS, NOT CPPFLAGS. configure probes for fftw3.h
# using CPPFLAGS and reports success, but Makefile.in never substitutes
# CPPFLAGS, so a CPPFLAGS-only invocation configures cleanly then fails every
# compile with "fftw3.h: No such file or directory".
#
# -malign-double is dropped from the stock flags: it is an i386 option that on
# x86-64 changes struct layout against the SysV ABI, while FFTW here is not
# built with it.
./configure --with-problem="$PROBLEM" --enable-openmp \
    CFLAGS="-O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -I$FFTW_INC" \
    LDFLAGS="-L$FFTW_LIB -Wl,-rpath,$FFTW_LIB" >configure.log 2>&1
make -j8 >make.log 2>&1

echo "== built: $ROOT/snoopy6/snoopy =="
grep -o 'The Snoopy code v[0-9.]*' src/snoopy.c | head -1
echo
echo "Do NOT test it with 'snoopy --version'. SNOOPY parses no arguments and"
echo "ignores unknown ones, so that starts a full simulation in \$PWD."
echo "Do NOT run it on the login node at all -- submit job_box_scan.pbs."
