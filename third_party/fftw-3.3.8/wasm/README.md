# FFTW3 on DCP (wasm build)

`fftw-module.js` is FFTW 3.3.8's core C library, generic (no-SIMD)
double-precision codelets, statically linked against this pipeline's own
vendored fdlibm (`../../../src/StandardMath`, which `../kernel/trig.c` is
already patched to require -- see "Why this works" below), cross-compiled
with Emscripten into one self-contained JS file (`SINGLE_FILE=1`, wasm
binary embedded, no separate `.wasm` fetch). `fftw-wasm.js` wraps it with
a plain-array-in/plain-array-out API.

This is a different shape from the SoFiA-2 port
(`../../SoFiA-2-master_2_5_1/wasm/`): FFTW is a *library*, not a CLI tool,
so there's no `callMain()`/file-I/O story here -- `fftw-driver.c` exports
two C functions directly (`fftw_dft_1d_wasm`, `fftw_r2c_1d_wasm`), called
from JS via `ccall`/`cwrap` against caller-allocated heap buffers.

## Why this works

- FFTW's `./configure` only turns on SIMD (SSE/AVX/NEON/AltiVec/...),
  threads, or OpenMP if you explicitly pass `--enable-*` for them --  all
  default off. `build.sh` passes none of them, the same reasoning as
  SoFiA's build never passing `-fopenmp`: no
  threads/SIMD-intrinsics-in-a-sandboxed-worker problem to solve, just
  FFTW's portable generic-C codelets. This is also exactly the
  configuration this project's *own* native build of FFTW here already
  uses -- confirmed via its `config.log`, which shows a plain
  `./configure` with no flags.
- Codelets are pre-generated in the release tarball (the `genfft`
  OCaml-based generator is only needed building from the git repo, per
  FFTW's own `README`), so no exotic build tooling is needed, the same
  situation as wcslib's pre-generated flex lexers in the SoFiA port.
- `../kernel/trig.c` is patched (see its own header comment, tagged
  "3KIDNAS FFTW3-port project") to call `fdlibm_sin`/`fdlibm_cos` --
  this pipeline's own vendored, bit-reproducible trig implementation --
  instead of the platform's `libm`, for the same reason Fortran's
  `SIN`/`COS` intrinsics get the same substitution elsewhere in this
  pipeline: different platforms' `libm`s can differ from each other and
  from fdlibm by up to ~1 ULP, and this removes that risk from FFTW's
  internal twiddle-factor computation. That patch leaves `libfftw3.a`
  with undefined `fdlibm_sin`/`fdlibm_cos` symbols by design, resolved
  by compiling five files from `../../../src/StandardMath` alongside it:
  `fdlibm_sin.c`, `fdlibm_cos.c`, their `fdlibm_k_sin.c`/`fdlibm_k_cos.c`
  kernels, and `fdlibm_rem_pio2_stub.c` -- a deliberate stub for the
  large-argument Payne-Hanek reduction path (`__kernel_rem_pio2`, only
  reached for `|x| >~ 2^20 * pi/2`) that aborts loudly instead of
  silently doing the wrong thing, since FFTW's twiddle-factor arguments
  (`2*pi*m/n` for practical transform sizes) never come close.
- FFTW's own source tree here (`../`) is already `./configure`'d in
  place for a native build (`config.h`, `libfftw3.la` already present) --
  `build.sh` never touches it. Autotools also refuses to configure a
  source tree out-of-tree once it's configured in-tree, so the only safe
  path is a full disposable copy in a temp dir; the real source tree is
  read-only from `build.sh`'s point of view.
- `make` also builds `tools/bench`, which needs its own separate copy of
  fdlibm linked in via `libbench2` that this build doesn't provide (an
  unrelated CLI benchmarking utility, not part of `libfftw3.a`) -- its
  failure is expected and irrelevant; `build.sh` checks for the actual
  `libfftw3.a` artifact rather than `make`'s exit code.
- Same `node:fs`/`node:crypto` situation as the SoFiA port: Emscripten's
  generated glue has a dead `ENVIRONMENT_IS_NODE`-gated `locateFile`
  fallback using the `node:` URI scheme, which DCP's webpack-based
  `job.requires()` bundler can't statically resolve even though the
  branch never executes in a sandboxed worker. `build.sh` patches it the
  same way, after linking.

## Verifying it

Checked against an independent naive O(n²) DFT (no relation to FFTW's
own code), computed in the same test scripts, not trusted from memory:

- **Local, both entry points**: a unit cosine at bin 3, N=16 -- matched
  the exact closed-form theoretical value (amplitude N/2 = 8.000 at bins
  3 and 13) as well as the naive-DFT reference to ~1e-14. A real-to-complex
  transform at N=17 (odd, to exercise the general non-power-of-2 case)
  matched to the same precision.
- **Via `fftw-wasm.js`'s actual API** (`dft1d`/`r2c1d`, not the raw
  `ccall` calls): four signals of deliberately awkward sizes (16, 17 odd,
  64, 100 composite), matched to ~1e-13.
- **Via the exact work-function logic in `test-fftw-worker.js`**, run
  locally with `progress()` stubbed before ever spending real funds
  deploying it: all four signals pass at the job's own 1e-8 tolerance
  (actual diffs ~1e-13 to 1e-14).

Not yet verified: an actual DCP worker. `test-fftw-worker.js` is ready
but hasn't been deployed yet -- see "Running the worker test" below.

## Rebuilding

```bash
source <path-to-emsdk>/emsdk_env.sh
./build.sh
```

No network access needed (unlike SoFiA's wcslib fetch): FFTW's source is
already local, copied from `../` into a disposable temp dir.

## API

```javascript
const fftw = require('./fftw-wasm');

const { re, im } = await fftw.dft1d(reArray, imArray);  // complex -> complex, length n
const { re, im } = await fftw.r2c1d(realArray);          // real -> complex, n/2+1 bins
```

`dft1d` takes `{ inverse: true }` as a third argument for `FFTW_BACKWARD`
instead of `FFTW_FORWARD`. `r2c1d` is forward-only (r2c is inherently one
direction); the complex spectrum comes back as the non-redundant half,
FFTW's standard r2c output size.

## Running the worker test

```bash
node test-fftw-worker.js --apiKey=0x<identity> [--computeGroup=key,secret]
```

Ships four signals (input set) and a shared numerical tolerance (static
argument) to real DCP workers via local (unpublished)
`job.requires(['./fftw-wasm'])`; each worker runs both DFT entry points
and checks its own result against an independent naive DFT computed in
the same work function. Prints PASS/FAIL per slice and a summary; exits
1 if anything failed.

## Publishing

Confirmed on a real DCP worker (via `test-fftw-worker.js`), then
packaged for publishing the same way as
`../../SoFiA-2-master_2_5_1/wasm/package/` (`sofia2wasm`), following
`docs/patching-wasm-libraries-for-dcp.pdf` in the Edequity repo.

- **`package/build-bravojs-bundle.js`** -- flattens `fftw-wasm.js` +
  `fftw-module.js` into one file and wraps it in `module.declare`. Run
  after `build.sh`:
  ```bash
  node package/build-bravojs-bundle.js
  ```
  Output (`package/fftw-wasm.js`) was functionally verified -- not just
  syntax-checked -- by simulating bravojs's `module.declare(deps,
  factory)` contract in Node and running a signal through both DFT entry
  points; matched an independent naive DFT to ~1e-14, same as every
  other build variant in this directory.
- **`package/package.dcp`** -- manifest. `name` is `fftw3wasm` (no
  uppercase letters, confirmed free on the npm registry -- DCP's
  publisher rejects names that collide with a real npm package, which
  bit the SoFiA port's first attempt); `version` is `3.3.8`, FFTW's own
  upstream version (unlike the SoFiA port, there's no separate
  independently-versioned library to encode as build metadata -- fdlibm
  here is this project's own vendored source, not a distinct upstream
  dependency).
- **`package/test-published-package.js`** -- run this *after* actually
  publishing, with `--apiKey=0x...`. Confirms
  `job.requires(['fftw3wasm/fftw-wasm.js'])` +
  `require('fftw-wasm.js')` (bare filename, no path) resolves and runs
  correctly for real.

To publish, from `dcp-util`'s own directory, with an absolute path (a
relative one resolves against wherever `./publish` is invoked from, not
`package.dcp`'s location):

```bash
cd /Users/dandesjardins/DCP/node_modules/dcp-util/bin/
./publish package /Users/dandesjardins/DCP/3KIDNAS/third_party/fftw-3.3.8/wasm/package/package.dcp
```
