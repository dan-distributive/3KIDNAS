#!/usr/bin/env bash
# Rebuilds cfitsio-module.js from a fresh emconfigure/emmake build of
# cfitsio 4.6.3 (zlib support via Emscripten's built-in zlib port, no
# curl/bzip2/threads/Fortran -- see "Why this works" in README.md) and
# cfitsio-driver.c's five exported entry points. Produces one
# self-contained JS file (wasm binary embedded via SINGLE_FILE=1).
#
# Usage: ./build.sh
# Requires: emsdk activated (source <emsdk>/emsdk_env.sh).

set -euo pipefail
cd "$(dirname "$0")"

command -v emcc >/dev/null || { echo "emcc not found -- source emsdk_env.sh first" >&2; exit 1; }

CFITSIO_SRC="$(cd .. && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# --------------------------------------------------------------------------
# cfitsio's own source directory is already ./configure'd in place for a
# native build (config.status, libcfitsio.la already present) -- don't
# touch it. Autotools also refuses to configure a source tree out-of-tree
# once it's already configured in-tree, so the only safe path is a full
# copy: work on a disposable copy, leave the real source tree untouched.
# --------------------------------------------------------------------------
echo "== Copying cfitsio source (leaving the real tree's native build untouched) =="
cp -R "$CFITSIO_SRC" "$WORK/cfitsio-src"
( cd "$WORK/cfitsio-src" && make distclean >/dev/null 2>&1 || true )

echo "== emconfigure: zlib via Emscripten's built-in port, no curl/threads =="
# --disable-curl: no real network in a sandboxed worker anyway.
# --without-zlib-check: skips configure's own AC_TRY_RUN-style zlib probe
# (can't execute a cross-compiled binary during configure); zlib itself is
# still genuinely linked in via -sUSE_ZLIB=1 below, which downloads/builds
# a real wasm zlib through Emscripten's own port system, cached after the
# first build. No --enable-reentrant: threading is opt-in and left off,
# same reasoning as SoFiA never passing -fopenmp.
( cd "$WORK/cfitsio-src" && emconfigure ./configure \
    --disable-shared --enable-static \
    --disable-curl \
    --without-zlib-check \
    CPPFLAGS="-sUSE_ZLIB=1" \
    >"$WORK/configure.log" 2>&1 ) \
  || { echo "configure failed, see $WORK/configure.log" >&2; cat "$WORK/configure.log" >&2; exit 1; }

echo "== emmake: building libcfitsio.a =="
# CFLAGS here also carries -D__arm__ (see "Why this works" in README.md):
# fitsio2.h's byte-order detection is a hardcoded #elif chain over known
# architecture macros (__x86_64__, __AARCH64EL__, __arm__, ...); wasm32
# matches none of them and silently falls through to a stale catch-all
# written for old big-endian PowerPC Macs, which is wrong here (wasm32 is
# little-endian). -D__arm__ routes it through the *already-correct*
# little-endian/32-bit-long branch instead -- confirmed via grep that
# __arm__ is checked nowhere else in cfitsio, so this has no other effect.
#
# CFLAGS also carries -UHAVE_NET_SERVICES -UHAVE_SHMEM_SERVICES
# -UHAVE_GSIFTP, undefining macros configure itself set to 1 (confirmed
# in this same command line's own DEFS, earlier on it) because it
# detected connect()/gethostbyname() as *linkable* -- true for
# Emscripten's libc stubs, but not meaningful for a sandboxed worker
# that can't do real networking. Left alone, cfileio.c's driver table
# unconditionally registers root://, shmem://, and gsiftp:// against
# drvrnet.c/drvrsmem.c/drvrgsiftp.c's real BSD-socket-using code, which
# pulls in Emscripten's Node-targeted WebSocket socket emulation --
# requiring the 'ws' npm package, which itself needs unpolyfilled Node
# core modules (zlib, stream) that broke a real deploy attempt with a
# webpack "Module not found" error, before this fix. Undefining these
# stops cfileio.c from referencing those drivers' init functions at all,
# so the linker never pulls those archive members in -- confirmed via
# `grep -c 'require("ws")'` on the output, 0 after this fix.
( cd "$WORK/cfitsio-src" && emmake make -j4 \
    CFLAGS="-g -O2 -Dg77Fortran -D__arm__ -UHAVE_NET_SERVICES -UHAVE_SHMEM_SERVICES -UHAVE_GSIFTP" \
    LDFLAGS="-sUSE_ZLIB=1" \
    libcfitsio.la >"$WORK/make.log" 2>&1 ) \
  || { echo "build failed, see $WORK/make.log" >&2; tail -60 "$WORK/make.log" >&2; exit 1; }
LIBCFITSIOA="$WORK/cfitsio-src/.libs/libcfitsio.a"
[ -f "$LIBCFITSIOA" ] || { echo "libcfitsio.a was not produced" >&2; exit 1; }
echo "  libcfitsio.a built ($(du -h "$LIBCFITSIOA" | cut -f1))"

echo "== Compiling driver and linking =="
emcc -O2 cfitsio-driver.c "$LIBCFITSIOA" \
  -sUSE_ZLIB=1 -D__arm__ \
  -I"$CFITSIO_SRC" \
  -s MODULARIZE=1 -s EXPORT_NAME=CfitsioModule \
  -s EXPORTED_FUNCTIONS='["_cfits_read_image_info_wasm","_cfits_read_image_double_wasm","_cfits_write_image_double_wasm","_cfits_read_key_str_wasm","_cfits_read_key_dbl_wasm","_cfits_status_message_wasm","_cfits_create_image_wasm","_cfits_write_key_dbl_wasm","_cfits_write_key_str_wasm","_cfits_write_image_data_wasm","_cfits_close_image_wasm","_malloc","_free"]' \
  -s EXPORTED_RUNTIME_METHODS='["ccall","cwrap","FS","HEAPF64","HEAP32","UTF8ToString","stringToUTF8"]' \
  -s FORCE_FILESYSTEM=1 \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s SINGLE_FILE=1 \
  -o cfitsio-module.js

# Same fix SoFiA's and FFTW's builds needed and for the same reason (see
# ../../../Edequity/dcp-wasm-package/docs/patching-wasm-libraries-for-dcp.pdf):
# Emscripten's Node-environment locateFile/PRNG-seed fallbacks use
# require("node:fs")/require("node:crypto"), which DCP's webpack-based
# job.requires() bundler can't statically resolve even though the branch
# is dead in a sandboxed worker (no `process` global).
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
' "$(dirname "$0")/cfitsio-module.js"

echo "== Done: $(dirname "$0")/cfitsio-module.js =="
