# SoFiA-2 on DCP (wasm build)

## Attribution

This is a wasm build of other people's software, unmodified except for the
compiler target. All credit for SoFiA-2 and wcslib belongs to their
original authors; nothing here is original scientific or algorithmic work.
If this ever gets published as a DCP package, the listing and any
announcement of it should name both projects, not just "the wasm build":

- **SoFiA-2**: Copyright the SoFiA team; GPLv3, full text in `../LICENSE`.
  Please cite, per `../README.md`:
  - Serra, P., Westmeier, T., Giese, N., et al., 2015, MNRAS, 448, 1922,
    "SoFiA: a flexible source finder for 3D spectral line data"
  - Westmeier, T., Kitaeff, S., Pallot, D., et al., 2021, MNRAS, 506, 3962,
    "SoFiA 2 – An automated, parallel HI source finding pipeline for the
    WALLABY survey"
  - Project page: https://github.com/SoFiA-Admin/SoFiA-2
- **wcslib**: Copyright 1995-2026 Mark Calabretta, Australia Telescope
  National Facility, CSIRO. https://www.atnf.csiro.au/computing/software/wcs

`sofia-module.js` is SoFiA-2 (`../sofia.c` + `../src/*.c`, unmodified) and
wcslib 8.9 (core C library, fetched fresh by `build.sh`), cross-compiled
with Emscripten into one self-contained JS file: the wasm binary is
embedded directly in it (`SINGLE_FILE=1`), so there's no separate `.wasm`
for a sandboxed worker to fetch. `sofia-wasm.js` wraps it with the
buffer-in/buffer-out API a DCP work function actually wants.

## Why this works

- SoFiA-2 itself is clean C99: no `pthread`/`fork`/`mmap`/`dlopen`, no
  `getenv`, no SIMD intrinsics beyond a couple of portable
  `__builtin_bswap*` calls. All file access is plain `fopen`/`fread`/
  `fwrite`, which maps directly onto Emscripten's in-memory filesystem
  (MEMFS) -- preload input bytes, run, read output bytes back out.
- It uses OpenMP (`#pragma omp` throughout `DataCube.c`/`LinkerPar.c`)
  for intra-node multi-threading, but every pragma is additive: without
  `-fopenmp` on the compiler invocation, they're silently ignored and
  the code runs sequentially. `build.sh` never passes `-fopenmp`, so
  there's no threads-inside-a-sandboxed-worker problem to solve. DCP's
  own sharding (see `test-sofia-worker.js`) replaces the parallelism this
  would otherwise have bought on one machine.
- The one external dependency, wcslib, is itself portable C with no
  threading or unusual syscalls. `build.sh` fetches its source fresh
  (checksum-verified) and compiles the ~20 core files SoFiA actually
  needs (`wcs.c`, `wcshdr.c`, `wcsfix.c`, and their dependencies) into a
  wasm static archive, skipping the Fortran bindings, PGSBOX, and the
  test/utility programs. `wcsconfig.h` is hand-written rather than run
  through wcslib's own `./configure`, because autoconf's `AC_TRY_RUN`
  feature checks can't execute a cross-compiled binary; both macros it
  would otherwise set are `#ifdef`-guarded with portable fallbacks in
  wcslib's headers, so this is a small, safe substitute, not a risky one.
- Emscripten's own generated glue includes two Node-only fallback
  branches (a `locateFile` helper and a crypto-based PRNG seed), each
  gated behind `ENVIRONMENT_IS_NODE` and each using `require("node:fs")`
  / `require("node:crypto")`. Both are dead code in a DCP worker (no
  `process` global means that check is always false), but DCP's
  webpack-based `job.requires()` bundler still has to statically resolve
  every `require()` it can see, and it only handles the legacy
  unprefixed form (`require("fs")`), not the `node:` URI scheme --
  confirmed by a real bundling failure (`UnhandledSchemeError`) the
  first time this job was deployed. `build.sh` patches both call sites
  to drop the `node:` prefix after linking, the same fix
  `docs/patching-wasm-libraries-for-dcp.pdf` in the Edequity repo
  documents for the identical problem in duckdb-wasm.

## Verifying it

`build.sh` was run end to end and the resulting module was exercised
against a real WALLABY HI data cube (`3KIDNASTests/TestData/
WALLABY_Test_sources/WALLABY_J103554-475245/WALLABY_J103554-475245_cube.fits`):
the S+C finder, linker, and WCS parameterisation all ran correctly, and
the one detected source's fitted sky position (RA/Dec computed via the
wasm build of wcslib) matched the source name already embedded in that
test file's own filename -- real coordinate math, not a stub.

Confirmed on a real DCP worker, not just local Node. A first deploy
attempt (`test-sofia-worker.js`) failed at the bundling stage, before
ever reaching a worker, on the `node:fs`/`node:crypto` issue described
above; after that patch, a real deploy (job id `8cG9uHpYtEGFDcyeo6hzvM`)
ran both WALLABY cubes on two separate remote workers and returned
`exit=0`, one detected source each, matching the local Node runs
exactly:

```
WALLABY_J103554-475245: exit=0  sources=1  files=[result_cat.txt, result_cat.xml]  (249 ms)
WALLABY_J103458-495128: exit=0  sources=1  files=[result_cat.txt, result_cat.xml]  (301 ms)
```

So: wasm bytes decode and instantiate correctly in the sandbox, MEMFS
preload/`callMain`/read-back works, and wcslib's WCS math runs correctly
remotely, all with no code changes beyond the `node:` scheme patch.

## Rebuilding

```bash
source <path-to-emsdk>/emsdk_env.sh
./build.sh
```

Fetches wcslib source fresh each time (checksum-pinned to 8.9), so it
needs network access; everything else is local.

## API

```javascript
const sofia = require('./sofia-wasm');

const { exitCode, log, files } = await sofia.run({
  cube: cubeBytes,   // Uint8Array, a FITS data cube
  par: parText,       // parameter file text, see ../template_par_file.par
  extraFiles: {},      // optional: weights/mask/noise cubes, keyed by MEMFS path
});
// files: Map<filename, Uint8Array> of whatever output.write* products
// the par text requested (catalogue, mask, moment maps, ...).
```

`par` must set `input.data = /work/cube.fits` and should set
`output.directory = /work/out` (created for you); everything else is the
same parameter grammar as the real `sofia` CLI, see
`../template_par_file.par` for the full list.

See `test-sofia-worker.js` for a full DCP job built on this via local
(unpublished) `job.requires(['./sofia-wasm'])`.

## Publishing

`package/` holds everything needed to publish this as a reusable DCP
package -- `sofia` -- so other job authors get it with one
`job.requires()` line and never see any of the Emscripten/MEMFS detail
above. Not yet published; this is prepared, not done.

Local `job.requires()` tolerates plain CommonJS because it goes through
a real webpack build at deploy time (confirmed: that's the same build
that needed the `node:` scheme patch above). `publish` does not -- it
looks like a raw file upload, no build step -- so the published artifact
has to already be in bravojs's module format
(`module.declare([], function(require, exports, module) {...})`), per
`docs/patching-wasm-libraries-for-dcp.pdf` in the Edequity repo (the
same requirement duckdb-wasm hit).

- **`package/build-bravojs-bundle.js`** -- flattens `sofia-wasm.js` +
  `sofia-module.js` into one file and wraps it in `module.declare`.
  Doesn't need a real bundler (unlike duckdb-wasm's 129-module Apache
  Arrow tree): this is two local files with one `require()` between
  them and zero external npm dependencies, so a plain closure
  composition is the whole job. Run after `build.sh`:
  ```bash
  node package/build-bravojs-bundle.js
  ```
  Output (`package/sofia-wasm.js`) was functionally verified -- not just
  syntax-checked -- by simulating bravojs's `module.declare(deps,
  factory)` contract in Node and running the real WALLABY test cube
  through it; produced the same correct catalogue as every other build
  variant in this directory.
- **`package/package.dcp`** -- the manifest. `name` has no uppercase
  letters (`sofia`); `version` is `2.5.1+wcslib.8.9`, both upstream
  versions actually bundled inside (SoFiA-2's own version, then wcslib's
  as semver build metadata) rather than an arbitrary counter, so the
  published version stays meaningfully tied to what's inside it -- bump
  both halves whenever either upstream source is updated.
- **`package/test-published-package.js`** -- run this *after* actually
  publishing, with `--apiKey=0x...`. Confirms
  `job.requires(['sofia2wasm/sofia-wasm.js'])` +
  `require('sofia-wasm.js')` (bare filename, no path -- a different
  convention from the local `./sofia-wasm` form) resolves and runs
  correctly for real, the same way `dcp-wasm-package/package/
  test-published-package.js` confirmed it for `duckdbwasm` in the
  Edequity repo. `test-sofia-worker.js` already confirmed the
  *unpublished* build works on a real worker (job id
  `8cG9uHpYtEGFDcyeo6hzvM`); this is a narrower, separate check that the
  bravojs wrap specifically survives the publish step.

To actually publish, from `dcp-util`'s own directory, with an absolute
path (a relative one resolves against wherever `./publish` is invoked
from, not `package.dcp`'s location -- confirmed gotcha in the companion
doc):

```bash
cd node_modules/dcp-util/bin/
./publish package /absolute/path/to/wasm/package/package.dcp
```
