#!/usr/bin/env bash
# Rebuilds fdlibm-module.js: an Emscripten build of this project's
# vendored fdlibm sin/cos/atanh (../fdlibm_sin.c, ../fdlibm_cos.c, and
# their kernel/rem_pio2 dependencies) plus fdlibm-driver.c's three
# scalar exported entry points. Produces one self-contained JS file
# (wasm binary embedded via SINGLE_FILE=1), same pattern as
# ../../../third_party/fftw-3.3.8/wasm/build.sh.
#
# Motivation: standalone native benchmarking (see
# .../scratchpad/fdlibm_native_compare.js from the investigation this
# mirrors) confirmed bit-exact output (0/200000 mismatches for sin,
# cos, atanh) between this same source compiled natively and the
# verified JS port (../fdlibm.js), with native ~5.6-6.1x faster per
# call -- consistent with the JS port's software-emulated FMA (fma.js)
# being the dominant cost, vs. real hardware FMA in a compiled build.
# This wasm build is the next step: get that same native speed inside
# a DCP worker, without touching fdlibm.js's already-verified behavior
# (fdlibm.js stays as the correctness fallback / non-wasm-worker path).
#
# *** RISK BELOW: CONFIRMED REAL, NOW FIXED -- READ BEFORE CHANGING THIS ***
# The native benchmark above only proved this project's actual C build
# (clang -O on arm64, hardware FMA mandatory baseline) fuses
# fdlibm_k_sin.c/fdlibm_k_cos.c's polynomial multiply-adds, matching
# ../fma.js's software emulation exactly (0/200000 mismatches). WASM's
# scalar f64 arithmetic has NO fused multiply-add instruction in the
# base spec (relaxed-simd's fma variants only cover f32x4/f64x2 SIMD
# lanes, not plain scalars) -- so an auto-contracted `a*b+c` compiled
# for wasm legalizes back to separate fmul+fadd (double rounding), not
# a fused op. This was CONFIRMED empirically, not just theorized:
# building fdlibm_k_sin.c/fdlibm_k_cos.c as-is for wasm and bit-exact
# comparing against the verified JS port gave 1271/200000 sin
# mismatches, 1107/200000 cos mismatches, 84/200000 atanh mismatches --
# the exact same class of divergence fma.js was written to eliminate,
# just relocated from "JS vs native" to "wasm vs native".
#
# THE FIX (now applied below): fdlibm_k_sin_wasm.c/fdlibm_k_cos_wasm.c
# are wasm-only copies of the kernel files with every contracted
# expression rewritten as an explicit soft_fma() call (fdlibm_soft_fma.h,
# a straight C port of ../fma.js's Dekker/Veltkamp algorithm) instead of
# relying on compiler contraction. Verified: 0/200000 mismatches for
# sin/cos/atanh after this change, including the specific regression
# input (theta=0.18265073567382517084) that flagged the divergence.
# The canonical fdlibm_k_sin.c/fdlibm_k_cos.c are UNTOUCHED and must
# stay that way -- they're correct for the native build precisely
# because that target has real hardware FMA to contract into.
#
# Usage: ./build.sh
# Requires: emsdk activated (source <emsdk>/emsdk_env.sh).

set -euo pipefail
cd "$(dirname "$0")"

command -v emcc >/dev/null || { echo "emcc not found -- source emsdk_env.sh first" >&2; exit 1; }

FDLIBM_SRC="$(cd .. && pwd)"

echo "== Compiling fdlibm sources + driver =="
# fdlibm_k_sin.c/fdlibm_k_cos.c (the canonical, hardware-FMA-contracted
# native versions) are deliberately NOT compiled here -- see this
# script's header. Their wasm-only replacements (fdlibm_k_sin_wasm.c/
# fdlibm_k_cos_wasm.c, this directory) use explicit soft_fma() instead.
# -ffp-contract=off is defense-in-depth: makes it explicit (not just
# implicit-via-lack-of-hardware-support) that these two files must not
# have their soft_fma() call sites' surrounding arithmetic contracted.
emcc -O2 -I"$FDLIBM_SRC" -c fdlibm-driver.c -o fdlibm-driver.o
emcc -O2 -I"$FDLIBM_SRC" -c "$FDLIBM_SRC/fdlibm_sin.c" -o fdlibm_sin.o
emcc -O2 -I"$FDLIBM_SRC" -c "$FDLIBM_SRC/fdlibm_cos.c" -o fdlibm_cos.o
emcc -O2 -ffp-contract=off -I"$FDLIBM_SRC" -I. -c fdlibm_k_sin_wasm.c -o fdlibm_k_sin.o
emcc -O2 -ffp-contract=off -I"$FDLIBM_SRC" -I. -c fdlibm_k_cos_wasm.c -o fdlibm_k_cos.o
emcc -O2 -I"$FDLIBM_SRC" -c "$FDLIBM_SRC/fdlibm_atanh.c" -o fdlibm_atanh.o
# fdlibm_log1p.c (canonical) not used here either -- same hardware-FMA
# contraction issue (its Horner-chain R(z) evaluation), same fix
# (fdlibm_log1p_wasm.c, explicit soft_fma()). fdlibm_atanh.c calls
# whichever fdlibm_log1p object gets linked in -- ordinary C linkage,
# no source change needed there.
emcc -O2 -ffp-contract=off -I"$FDLIBM_SRC" -I. -c fdlibm_log1p_wasm.c -o fdlibm_log1p.o
emcc -O2 -I"$FDLIBM_SRC" -c "$FDLIBM_SRC/fdlibm_rem_pio2_stub.c" -o fdlibm_rem_pio2_stub.o

echo "== Linking =="
emcc -O2 fdlibm-driver.o fdlibm_sin.o fdlibm_cos.o fdlibm_k_sin.o fdlibm_k_cos.o \
  fdlibm_atanh.o fdlibm_log1p.o fdlibm_rem_pio2_stub.o \
  -s MODULARIZE=1 -s EXPORT_NAME=FdlibmModule \
  -s EXPORTED_FUNCTIONS='["_fd_sin_wasm","_fd_cos_wasm","_fd_atanh_wasm"]' \
  -s EXPORTED_RUNTIME_METHODS='["cwrap"]' \
  -s SINGLE_FILE=1 \
  -o fdlibm-module.js

rm -f ./*.o

# Same fix FFTW/SoFiA's builds needed and for the same reason (see
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
' "$(dirname "$0")/fdlibm-module.js"

echo "== Done: $(dirname "$0")/fdlibm-module.js =="
