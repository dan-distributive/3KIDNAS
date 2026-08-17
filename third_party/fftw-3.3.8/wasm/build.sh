#!/usr/bin/env bash
# Rebuilds fftw-module.js from a fresh emconfigure/emmake build of FFTW
# 3.3.8 (generic codelets, no SIMD, no threads -- see "Why this works" in
# README.md), fftw-driver.c's two exported entry points, and the
# project's own vendored fdlibm (../../../src/StandardMath) that
# ../../kernel/trig.c is patched to require. Produces one self-contained
# JS file (wasm binary embedded via SINGLE_FILE=1).
#
# Usage: ./build.sh
# Requires: emsdk activated (source <emsdk>/emsdk_env.sh).

set -euo pipefail
cd "$(dirname "$0")"

command -v emcc >/dev/null || { echo "emcc not found -- source emsdk_env.sh first" >&2; exit 1; }

FFTW_SRC="$(cd .. && pwd)"
FDLIBM_SRC="$(cd ../../../src/StandardMath && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# --------------------------------------------------------------------------
# FFTW's own source directory is already ./configure'd in place for a
# native build (config.h, Makefile, libfftw3.la already present -- don't
# touch it, something else may depend on that build). Autotools also
# refuses to configure a source tree out-of-tree once it's already been
# configured in-tree ("source directory already configured; run make
# distclean there first"), so the only safe path is a full copy: work on
# a disposable copy, leave the real source tree untouched.
# --------------------------------------------------------------------------
echo "== Copying FFTW source (leaving the real tree's native build untouched) =="
cp -R "$FFTW_SRC" "$WORK/fftw-src"
( cd "$WORK/fftw-src" && make distclean >/dev/null 2>&1 || true )

echo "== emconfigure: generic codelets, no SIMD, no threads, no Fortran =="
# No --enable-sse2/--enable-avx/--enable-neon/etc, no --enable-threads,
# no --enable-openmp: all default OFF, giving FFTW's portable generic-C
# codelets -- the same reasoning as never passing SoFiA's -fopenmp. This
# is also exactly the configuration this project's own native build here
# already uses (a plain `./configure`, confirmed via its config.log).
( cd "$WORK/fftw-src" && emconfigure ./configure \
    --disable-shared --enable-static \
    --disable-fortran \
    --disable-doc \
    >"$WORK/configure.log" 2>&1 ) \
  || { echo "configure failed, see $WORK/configure.log" >&2; cat "$WORK/configure.log" >&2; exit 1; }

# -ffp-contract=off: the native build's own configure picks this up
# automatically (confirmed via its Makefile: CFLAGS = -O3
# -fomit-frame-pointer -mtune=native -fstrict-aliasing -ffp-contract=off),
# but emconfigure's cross-compile probe does NOT -- confirmed directly by
# running the same ./configure invocation under emconfigure and diffing
# the resulting Makefile's CFLAGS line, which comes out missing this flag.
# That asymmetry lets the WASM build's FFTW codelets fuse multiply-add
# into a single-rounding FMA instruction while the native build's don't --
# a real, plausible source of the Fortran/JS model-cube divergence this
# was built to chase down (same class of issue as random.js's own
# documented gfortran FMA investigation, just in FFTW's codelets instead
# of ran2). Appended via sed rather than pre-set in the CFLAGS
# environment variable -- setting CFLAGS before configure SUPPRESSES its
# own auto-detection entirely (confirmed directly: doing that drops
# -O3/-fomit-frame-pointer/-mtune=native/-fstrict-aliasing too, leaving
# only the one flag), so the flag has to be added to the Makefiles
# configure already generated, not fed in ahead of time. FFTW is a
# recursive-make project -- every subdirectory (kernel/, rdft/, dft/,
# simd-support/, ...) gets its OWN independently-substituted CFLAGS line
# from the SAME configure run, confirmed directly (kernel/Makefile and
# rdft/Makefile both had the identical pre-patch CFLAGS line the top-level
# Makefile did) -- patching only the top-level Makefile would leave every
# actual codelet subdirectory still building with FMA contraction allowed.
echo "== Patching every subdirectory Makefile's CFLAGS to add -ffp-contract=off (see comment above) =="
( cd "$WORK/fftw-src" && find . -name Makefile -exec sed -i.bak 's/^CFLAGS = .*/& -ffp-contract=off/' {} \; \
  && find . -name 'Makefile.bak' -delete )
grep -rc '^CFLAGS = .*-ffp-contract=off' "$WORK/fftw-src" --include=Makefile | grep -v ':0$' | wc -l
echo "Makefiles patched (of $(find "$WORK/fftw-src" -name Makefile | wc -l) total)"

echo "== emmake: building (bench/tools failure below is expected and irrelevant -- see note) =="
# `make` also builds tools/bench, which needs its own copy of fdlibm
# linked in (via libbench2) that we don't provide here -- that's an
# unrelated CLI benchmarking utility, not part of libfftw3.a, and not
# something this port needs. Its failure is why this ignores make's exit
# code and checks for the actual artifact instead.
( cd "$WORK/fftw-src" && emmake make -j4 >"$WORK/make.log" 2>&1 ) || true
LIBFFTW3A="$WORK/fftw-src/.libs/libfftw3.a"
[ -f "$LIBFFTW3A" ] || { echo "libfftw3.a was not produced, see $WORK/make.log" >&2; tail -60 "$WORK/make.log" >&2; exit 1; }
echo "  libfftw3.a built ($(du -h "$LIBFFTW3A" | cut -f1))"

# --------------------------------------------------------------------------
# ../../kernel/trig.c is patched (see its own header comment) to call
# fdlibm_sin/fdlibm_cos instead of the platform libm, for bit-reproducible
# trig across platforms -- the same substitution already made for
# Fortran's SIN/COS elsewhere in this pipeline. libfftw3.a therefore has
# undefined fdlibm_sin/fdlibm_cos symbols by design; only the __kernel_sin/
# __kernel_cos/__ieee754_rem_pio2 dependency chain is needed to resolve
# them (__kernel_rem_pio2 -- the large-argument Payne-Hanek path -- is
# provided by a stub that aborts loudly instead, since FFTW's twiddle
# factors never approach the magnitude that path exists for).
# --------------------------------------------------------------------------
echo "== Compiling fdlibm (bit-reproducible sin/cos that trig.c requires) =="
mkdir -p "$WORK/fdlibm-obj"
( cd "$FDLIBM_SRC" && for f in fdlibm_sin.c fdlibm_cos.c fdlibm_k_sin.c fdlibm_k_cos.c fdlibm_rem_pio2_stub.c; do
    emcc -O2 -I. -c "$f" -o "$WORK/fdlibm-obj/$(basename "$f" .c).o"
  done )

echo "== Compiling driver and linking =="
emcc -O2 fftw-driver.c "$LIBFFTW3A" "$WORK"/fdlibm-obj/*.o \
  -I"$FFTW_SRC/api" \
  -s MODULARIZE=1 -s EXPORT_NAME=FftwModule \
  -s EXPORTED_FUNCTIONS='["_fftw_dft_1d_wasm","_fftw_r2c_1d_wasm","_fftw_r2c_2d_wasm","_fftw_c2r_2d_wasm","_malloc","_free"]' \
  -s EXPORTED_RUNTIME_METHODS='["ccall","cwrap","HEAPF64","getValue","setValue"]' \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s SINGLE_FILE=1 \
  -o fftw-module.js

# Same fix SoFiA's build needed and for the same reason (see
# ../../../Edequity/dcp-wasm-package/docs/patching-wasm-libraries-for-dcp.pdf):
# Emscripten's Node-environment locateFile fallback uses
# require("node:fs"), which DCP's webpack-based job.requires() bundler
# can't statically resolve even though the branch is dead in a sandboxed
# worker (no `process` global).
echo "== Patching node: scheme requires (dead code, but must still bundle) =="
node -e '
const fs = require("fs");
const path = process.argv[1];
let data = fs.readFileSync(path, "utf8");
const before = (data.match(/require\("node:fs"\)/g) || []).length
              + (data.match(/require\("node:crypto"\)/g) || []).length;
data = data.replaceAll("require(\"node:fs\")", "require(\"fs\")")
           .replaceAll("require(\"node:crypto\")", "require(\"crypto\")");
fs.writeFileSync(path, data);
console.log(`  patched ${before} occurrence(s)`);
' "$(dirname "$0")/fftw-module.js"

echo "== Done: $(dirname "$0")/fftw-module.js =="
