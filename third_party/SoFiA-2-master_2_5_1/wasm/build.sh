#!/usr/bin/env bash
# Rebuilds sofia-module.js from ../src/*.c + ../sofia.c and a fresh wcslib
# checkout. Produces one self-contained JS file (wasm binary embedded via
# SINGLE_FILE=1, no separate .wasm to fetch) suitable for job.requires()
# in a DCP worker, or plain require() under Node for local testing.
#
# Usage: ./build.sh
# Requires: emsdk activated (source <emsdk>/emsdk_env.sh) and network
# access to download wcslib source once, into ./wcslib-src-tmp (removed
# at the end).

set -euo pipefail
cd "$(dirname "$0")"

command -v emcc >/dev/null || { echo "emcc not found -- source emsdk_env.sh first" >&2; exit 1; }

WCSLIB_VERSION=8.9
WCSLIB_SHA256=82ac09ce5091b0bf06cec8f5cdeec1dabe1d06ba5dfb7ff2bdb0c1680488807b
SOFIA_SRC="$(cd .. && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "== Fetching wcslib $WCSLIB_VERSION source =="
curl -sL -o "$WORK/wcslib.tar.bz2" \
  "https://www.atnf.csiro.au/computing/software/wcs/wcslib-releases/wcslib-${WCSLIB_VERSION}.tar.bz2"
echo "$WCSLIB_SHA256  $WORK/wcslib.tar.bz2" | shasum -a 256 -c -
tar xjf "$WORK/wcslib.tar.bz2" -C "$WORK"

echo "== Assembling wcslib core sources (excluding Fortran/PGSBOX/tests/utils) =="
mkdir -p "$WORK/wcslib-c" "$WORK/include/wcslib"
cp "$WORK/wcslib-${WCSLIB_VERSION}"/C/*.c "$WORK/wcslib-${WCSLIB_VERSION}"/C/*.h "$WORK/wcslib-c/"
cp "$WORK/wcslib-${WCSLIB_VERSION}"/C/flexed/*.c "$WORK/wcslib-c/"   # pre-generated flex lexers
rm -f "$WORK/wcslib-c/getwcstab.c"    # excluded from the base library upstream too (needs cfitsio; SoFiA doesn't use -TAB)
cat > "$WORK/wcslib-c/wcsconfig.h" <<'EOF'
/* Hand-written: autoconf's AC_TRY_RUN checks can't execute a
 * cross-compiled binary, and both macros here are optional in wcslib's
 * own headers (each #ifdef-guarded with a portable fallback). */
#define HAVE_WCSLIB_VERSION
#define WCSLIB_VERSION 8.9
/* #undef HAVE_SINCOS */          /* glibc extension, not in Emscripten's libc */
#define WCSLIB_INT64 long long int
EOF
cp "$WORK"/wcslib-c/*.h "$WORK/include/wcslib/"

echo "== Compiling wcslib to wasm =="
mkdir -p "$WORK/wcslib-obj"
( cd "$WORK/wcslib-c" && for f in *.c; do
    emcc -O2 -std=gnu99 -DWCSLIB_TLS= -I. -c "$f" -o "$WORK/wcslib-obj/${f%.c}.o"
  done )
emar rcs "$WORK/libwcs-wasm.a" "$WORK"/wcslib-obj/*.o

echo "== Compiling SoFiA-2 to wasm =="
SOFIA_FILES="src/Array_dbl.c src/Array_siz.c src/Catalog.c src/common.c src/DataCube.c
  src/Flagger.c src/Header.c src/LinkerPar.c src/Map.c src/Matrix.c src/Parameter.c
  src/Path.c src/Source.c src/Stack.c src/statistics_dbl.c src/statistics_flt.c
  src/String.c src/Table.c src/WCS.c sofia.c"
mkdir -p "$WORK/sofia-obj"
( cd "$SOFIA_SRC" && for f in $SOFIA_FILES; do
    base=$(basename "$f" .c)
    emcc -O2 --std=c99 -DWCSLIB_TLS= -I"$WORK/include" -c "$f" -o "$WORK/sofia-obj/${base}.o"
  done )

echo "== Linking =="
emcc -O2 "$WORK"/sofia-obj/*.o "$WORK/libwcs-wasm.a" \
  -s MODULARIZE=1 -s EXPORT_NAME=SofiaModule \
  -s EXPORTED_RUNTIME_METHODS='["callMain","FS","ENV"]' \
  -s FORCE_FILESYSTEM=1 \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s INITIAL_MEMORY=67108864 \
  -s SINGLE_FILE=1 \
  -s EXIT_RUNTIME=1 \
  -s INVOKE_RUN=0 \
  -o "$(dirname "$0")/sofia-module.js"

# Emscripten's Node-environment fallback paths (locateFile's scriptDirectory
# lookup, and a crypto-based PRNG seed) use require("node:fs") /
# require("node:crypto"). That code is dead here -- no `process` means
# ENVIRONMENT_IS_NODE is always false in a DCP worker -- but DCP's webpack-
# based job.requires() bundler still has to statically resolve every
# require() it can see, and it only handles the legacy unprefixed form
# (require("fs")), not the node: URI scheme. Same fix duckdb-wasm needed
# (see docs/patching-wasm-libraries-for-dcp.pdf in the Edequity repo):
# drop the node: prefix so the dead branch still builds.
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
' "$(dirname "$0")/sofia-module.js"

echo "== Done: $(dirname "$0")/sofia-module.js =="
