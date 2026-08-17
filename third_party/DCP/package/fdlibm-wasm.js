module.declare(["./fdlibm-wasm","./fdlibm-module"], function (require, exports, module) {
/**
 * @file        fdlibm-wasm.js
 * @description Single entry point for calling this project's vendored
 *              fdlibm sin/cos/atanh (compiled natively, no software FMA
 *              emulation -- see build.sh's header for the benchmark this
 *              is based on) inside a DCP worker. require()s
 *              fdlibm-module.js, the Emscripten build produced by
 *              build.sh, with the wasm binary embedded directly in the
 *              JS (SINGLE_FILE=1).
 *
 *              fd_sin_wasm/fd_cos_wasm/fd_atanh_wasm are plain scalar
 *              double-in/double-out C functions (fdlibm-driver.c) --
 *              no caller-allocated buffers, no malloc/HEAPF64
 *              bookkeeping, unlike fftw-wasm.js's dft1d/r2c1d.
 *
 * @usage       const fdlibm = require('./fdlibm-wasm');
 *              await fdlibm.warmUp();           // once, forces async wasm instantiation
 *              const s = fdlibm.fdSinWasm(x);    // sync from here on
 *              const c = fdlibm.fdCosWasm(x);
 *              const a = fdlibm.fdAtanhWasm(x);
 *
 *              This module is NOT yet wired into any production call
 *              site (buildTiltedRingModel/ring_ParticleGeneration still
 *              use ../fdlibm.js's fdSin/fdCos). It exists so the wasm
 *              build can be benchmarked in-worker (JS<->wasm scalar
 *              call overhead wasn't measurable without emcc available
 *              in the dev environment that wrote this file) before
 *              deciding whether to swap the hot path over.
 */
'use strict';

// .default because fdlibm-module.js exports via `exports.default = ...`
// (property assignment on the given exports object), not `module.exports =
// ...` (whole-object reassignment) -- the latter silently breaks under the
// real DCP sandbox's module.declare() wrapping (confirmed via a real
// dispatch: "Cannot set properties of undefined (setting 'default')" from
// fdlibm-module.js's own raw Emscripten UMD tail, which used exactly that
// pattern). See fftw-wasm.js's published bundle for the same underlying
// issue solved a different way (full inlining instead of a require()
// boundary) -- this file still crosses a real module boundary, so it needs
// the property-assignment-only export style every other file in this
// package already uses.
const { default: createFdlibmModule } = require('./fdlibm-module');

let modulePromise = null;
let resolvedModule = null;
let fdSinFn = null;
let fdCosFn = null;
let fdAtanhFn = null;

function getModule() {
  if (!modulePromise) {
    modulePromise = createFdlibmModule().then((m) => {
      resolvedModule = m;
      fdSinFn = m.cwrap('fd_sin_wasm', 'number', ['number']);
      fdCosFn = m.cwrap('fd_cos_wasm', 'number', ['number']);
      fdAtanhFn = m.cwrap('fd_atanh_wasm', 'number', ['number']);
      return m;
    });
  }
  return modulePromise;
}

/** Forces the one-time async wasm instantiation. Call once before using fdSinWasm/fdCosWasm/fdAtanhWasm. */
async function warmUp() {
  await getModule();
}

function requireWarm() {
  if (!resolvedModule) {
    throw new Error('fdlibm-wasm: module not warmed up -- call `await fdlibm.warmUp()` once first');
  }
}

function fdSinWasm(x) { requireWarm(); return fdSinFn(x); }
function fdCosWasm(x) { requireWarm(); return fdCosFn(x); }
function fdAtanhWasm(x) { requireWarm(); return fdAtanhFn(x); }

module.exports = { warmUp, fdSinWasm, fdCosWasm, fdAtanhWasm };

});
