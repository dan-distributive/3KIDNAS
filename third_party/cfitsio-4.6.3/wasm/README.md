# cfitsio on DCP (wasm build)

`cfitsio-module.js` is cfitsio 4.6.3's core C library (real,
standards-compliant FITS I/O: header keywords, image data, gzip-compressed
files), statically linked against a wasm build of zlib (via Emscripten's
own built-in zlib port), cross-compiled with Emscripten into one
self-contained JS file (`SINGLE_FILE=1`, wasm binary embedded, no separate
`.wasm` fetch). `cfitsio-wasm.js` wraps it with a plain-object-in/out API.

Like SoFiA-2 and unlike FFTW3, this is file-I/O shaped: `cfitsio-driver.c`
exports functions that operate on MEMFS paths (`fits_open_file`/
`fits_create_file` underneath), not caller-allocated buffers directly.

## Why this works

- cfitsio's `./configure` only pulls in curl (`--disable-curl` here; no
  real network in a sandboxed worker anyway), bzip2 (opt-in via
  `--with-bzip2`, never passed), or threading (`--enable-reentrant`,
  opt-in, never passed) if asked. zlib is the one hard, unconditional
  dependency (`zcompress.c` includes `zlib.h` with no `#ifdef` guard) --
  satisfied by Emscripten's built-in zlib port (`-sUSE_ZLIB=1`), which
  downloads and builds a real wasm zlib automatically, cached after the
  first build.
- The `eval.l`/`eval.y` files (cfitsio's row-filtering expression
  language) are pre-generated in the release tarball as `eval_l.c`/
  `eval_y.c`, so no flex/bison needed at build time -- the same situation
  as wcslib's pre-generated lexers in the SoFiA port.
- **A real, substantive portability bug, found and fixed**: cfitsio's own
  `fitsio2.h` decides byte order via a hardcoded `#elif` chain over known
  architecture macros (`__x86_64__`, `__AARCH64EL__`, `__arm__`,
  `__riscv`, ...). wasm32 matches none of them (Emscripten predefines
  `__wasm32__`/`__EMSCRIPTEN__`, not any of the checked names), so the
  chain silently falls through to its final `#else`: `BYTESWAPPED FALSE`,
  `CFITSIO_MACHINE NATIVE` -- a branch written for old big-endian PowerPC
  Macs ("e.g., Macs fall into this category"), wrong for little-endian
  wasm32. Confirmed directly: without the fix, cfitsio's own internal
  startup self-test caught it immediately --
  `"Byteswapping is not being done correctly on this system"` -- rather
  than silently corrupting data, which is exactly what that self-test is
  for. Fixed with `-D__arm__`: cfitsio's *already-correct*
  32-bit-little-endian branch (`grep`-confirmed to be checked nowhere
  else in the ~34MB source tree, so this has no other effect), rather
  than hand-patching the header. `build.sh` passes it at both the library
  build and driver-link steps.
- **A second real bug, found via an actual failed deploy, not anticipated
  in advance**: cfitsio's `./configure` detected `connect()` and
  `gethostbyname()` as linkable (Emscripten's libc provides real POSIX
  socket *signatures*, even though nothing behind them does real
  networking) and set `HAVE_NET_SERVICES=1`. Left alone, `cfileio.c`'s
  driver table unconditionally registers `root://`/`http://`
  (`drvrnet.c`), `shmem://` (`drvrsmem.c`, `HAVE_SHMEM_SERVICES`), and
  `gsiftp://` (`drvrgsiftp.c`, `HAVE_GSIFTP`) -- real BSD-socket-using
  code, none of it reachable or meaningful from a sandboxed worker that
  will only ever read/write local MEMFS files. Compiling it in caused
  Emscripten to link its Node-targeted WebSocket socket-emulation layer,
  which `require()`s the `ws` npm package -- which itself needs
  unpolyfilled Node core modules (`zlib`, `stream`) that DCP's
  webpack-based bundler can't resolve. This broke a real deploy attempt
  outright (`Module not found: Can't resolve 'zlib'`/`'stream'`, tracing
  straight back to `ws/lib/permessage-deflate.js` and `ws/lib/receiver.js`
  through `cfitsio-module.js`), the same *class* of dead-code-must-still-
  build problem as the `node:` scheme issue, just one dependency layer
  removed and only visible once something actually tried to bundle it.
  Fixed by undefining all three feature macros
  (`-UHAVE_NET_SERVICES -UHAVE_SHMEM_SERVICES -UHAVE_GSIFTP`) when
  building the library: `cfileio.c` stops referencing those drivers'
  init functions at all, so the linker never pulls the relevant archive
  members in, and `ws` never enters the picture. Confirmed via
  `grep -c 'require("ws")'` on the output: present before, `0` after.
- Same `node:fs`/`node:crypto` situation as the SoFiA and FFTW3 ports:
  Emscripten's generated glue has dead `ENVIRONMENT_IS_NODE`-gated
  fallbacks using the `node:` URI scheme, which DCP's webpack-based
  `job.requires()` bundler can't statically resolve even though the
  branches never execute in a sandboxed worker. `build.sh` patches them
  the same way, after linking.
- cfitsio's own source tree here (`../`) is already `./configure`'d in
  place for a native build -- `build.sh` never touches it, working on a
  disposable copy instead, the same reasoning as the FFTW3 port (autotools
  refuses to configure a source tree out-of-tree once it's configured
  in-tree).

## Verifying it

Checked against a real WALLABY HI cube
(`3KIDNASTests/TestData/WALLABY_Test_sources/WALLABY_J103554-475245/WALLABY_J103554-475245_cube.fits`),
cross-checked against two independent facts already established earlier
in this same investigation (by SoFiA-2, a completely different tool
reading the same file):

- **Dimensions**: `naxis=3, naxes=[33, 35, 53]` -- exactly matches
  SoFiA-2's own "Axis sizes: 33, 35, 53" from its own run against this
  file.
- **Pixel data**: all 61215 pixels finite, `min=-0.00697`,
  `max=0.00958` -- consistent with SoFiA-2's own "No infinite pixel
  values found" on the same file.
- **Header keyword**: `CRPIX1 = -1012` read correctly via
  `fits_read_key`.
- **Exact round-trip**: wrote a synthetic 4x3 double-precision image,
  read it back -- byte-for-byte identical data, valid FITS block
  structure (file size a multiple of 2880, FITS's fixed block size).

All of the above via `cfitsio-wasm.js`'s actual API, not just the raw
`ccall` entry points.

Confirmed on real DCP workers, not just local Node. A first deploy
attempt (`test-cfitsio-worker.js`, local `job.requires(['./cfitsio-wasm'])`)
failed at the bundling stage, before ever reaching a worker, on the
`ws`/`zlib`/`stream` issue described above; after that fix, a real
deploy (job id `MgSzZtPEqtSXwjoSRxrcZe`) ran both WALLABY cubes on two
separate remote workers and returned exact matches to the local Node
results:

```
WALLABY_J103554-475245: naxes=[33,35,53]  CRPIX1=-1012  min=-0.006969  max=0.009581  roundTrip=true  (220 ms)
WALLABY_J103458-495128: naxes=[37,37,51]  CRPIX1=-1060  min=-0.006570  max=0.007413  roundTrip=true  (205 ms)
```

Not yet published (`package/` is prepared and functionally verified via
the simulated bravojs `module.declare` contract, but `./publish` hasn't
been run for this package yet).

## Rebuilding

```bash
source <path-to-emsdk>/emsdk_env.sh
./build.sh
```

Needs network access the first time (Emscripten's zlib port downloads
zlib source once, then caches the wasm build); cfitsio's own source is
already local, copied from `../` into a disposable temp dir.

## API

```javascript
const cfitsio = require('./cfitsio-wasm');

const { naxis, naxes, bitpix } = await cfitsio.readImageInfo(fitsBytes);
const { naxes, data } = await cfitsio.readImageDouble(fitsBytes);
const value = await cfitsio.readKeyDouble(fitsBytes, 'CRPIX1');
const text = await cfitsio.readKeyString(fitsBytes, 'OBJECT');
const fitsBytes = await cfitsio.writeImageDouble(naxes, data);
```

`fitsBytes`/the return of `writeImageDouble` are plain `Uint8Array`s --
whole FITS files, ready to preload as a job's static argument or input-set
element, or to hand back as a work function's result.

## Running the worker test

```bash
node test-cfitsio-worker.js --apiKey=0x<identity> [--computeGroup=key,secret]
```

## Publishing

Packaged the same way as `../../SoFiA-2-master_2_5_1/wasm/package/`
(`sofia2wasm`) and `../../fftw-3.3.8/wasm/package/` (`fftw3wasm`),
following `docs/patching-wasm-libraries-for-dcp.pdf` in the Edequity
repo.

- **`package/build-bravojs-bundle.js`** -- flattens `cfitsio-wasm.js` +
  `cfitsio-module.js` into one file and wraps it in `module.declare`.
  Run after `build.sh`:
  ```bash
  node package/build-bravojs-bundle.js
  ```
  Output (`package/cfitsio-wasm.js`) was functionally verified -- not
  just syntax-checked -- by simulating bravojs's `module.declare(deps,
  factory)` contract in Node: read the real WALLABY cube's dimensions
  (matched `[33, 35, 53]`) and did an exact round-trip write+read, same
  as every other build variant in this directory.
- **`package/package.dcp`** -- manifest. `name` is `cfitsio4wasm` (no
  uppercase letters, confirmed free on the npm registry, matching the
  `sofia2wasm`/`fftw3wasm` naming convention); `version` is `4.6.3`,
  cfitsio's own upstream version.
- **`package/test-published-package.js`** -- run this *after* actually
  publishing, with `--apiKey=0x...`. Confirms
  `job.requires(['cfitsio4wasm/cfitsio-wasm.js'])` +
  `require('cfitsio-wasm.js')` (bare filename, no path) resolves and
  runs correctly for real, against the same WALLABY cube used throughout
  this port's verification.

To publish, from `dcp-util`'s own directory, with an absolute path (a
relative one resolves against wherever `./publish` is invoked from, not
`package.dcp`'s location):

```bash
cd /Users/dandesjardins/DCP/node_modules/dcp-util/bin/
./publish package /Users/dandesjardins/DCP/3KIDNAS/third_party/cfitsio-4.6.3/wasm/package/package.dcp
```
