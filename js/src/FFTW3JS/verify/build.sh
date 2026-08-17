#!/bin/sh
# Builds ground_truth_harness. Run from this directory.
#
# Links the project's vendored, patched libfftw3.a (see
# third_party/fftw-3.3.8/kernel/trig.c and src/makeflags for why it's
# patched: -ffp-contract=off + fdlibm sin/cos instead of platform libm)
# plus the same fdlibm object files the real Fortran build links, since
# trig.c now calls fdlibm_sin/fdlibm_cos directly and expects them resolved
# at final link time.
set -e

ROOT=$(cd "$(dirname "$0")/../../../../.." && pwd)
FDLIBM_SRC="$ROOT/src/StandardMath"
FFTW_API="$ROOT/third_party/fftw-3.3.8/api"
FFTW_LIB="$ROOT/third_party/fftw-3.3.8/.libs/libfftw3.a"

gcc -c "$FDLIBM_SRC/fdlibm_sin.c" "$FDLIBM_SRC/fdlibm_cos.c" \
       "$FDLIBM_SRC/fdlibm_k_sin.c" "$FDLIBM_SRC/fdlibm_k_cos.c" \
       "$FDLIBM_SRC/fdlibm_rem_pio2_stub.c"

gcc -I"$FFTW_API" ground_truth_harness.c \
    fdlibm_sin.o fdlibm_cos.o fdlibm_k_sin.o fdlibm_k_cos.o fdlibm_rem_pio2_stub.o \
    "$FFTW_LIB" -o ground_truth_harness -lm

echo "built ./ground_truth_harness"
