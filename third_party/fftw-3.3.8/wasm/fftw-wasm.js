/**
 * @file        fftw-wasm.js
 * @description Single entry point for running FFTW 3.3.8 (double-precision,
 *              generic/no-SIMD codelets, see build.sh) inside a DCP worker.
 *              job.requires(['./fftw-wasm']) is enough -- this file
 *              require()s fftw-module.js, the Emscripten build of libfftw3
 *              (statically linked against this pipeline's own vendored
 *              fdlibm, which ../kernel/trig.c is patched to require for
 *              bit-reproducible twiddle factors), with the wasm binary
 *              embedded directly in the JS (SINGLE_FILE=1).
 *
 *              Unlike SoFiA-2 (a CLI tool driven by callMain() + file I/O),
 *              FFTW is a library: fftw-driver.c exports two plain C
 *              functions operating on caller-allocated heap buffers, called
 *              here via ccall/cwrap. dft1d/r2c1d below hide the
 *              malloc/HEAPF64/free bookkeeping that implies.
 *
 * @usage       const fftw = require('./fftw-wasm');
 *              const { re, im } = await fftw.dft1d(reArray, imArray);        // complex -> complex
 *              const { re, im } = await fftw.r2c1d(realArray);               // real -> complex (n/2+1 bins)
 *
 *              SYNC API -- for hot paths that can't await (e.g. a synchronous
 *              optimizer objective function): call `await fftw.warmUp()` ONCE
 *              (forces the one genuinely async step, wasm instantiation),
 *              then use dft1dSync/r2c1dSync freely -- same cwrap'd C calls as
 *              the async versions, just skipping the per-call Promise/await
 *              once the module is already resident. Throws if called before
 *              warmUp() resolves.
 */
'use strict';

const createFftwModule = require('./fftw-module');

let modulePromise = null;
let resolvedModule = null;
function getModule() {
  if (!modulePromise) {
    modulePromise = createFftwModule().then((m) => { resolvedModule = m; return m; });
  }
  return modulePromise;
}

/** Forces the one-time async wasm instantiation. Call once before using the *Sync API. */
async function warmUp() {
  await getModule();
}

/** Returns the already-instantiated Module synchronously. Throws if warmUp() hasn't resolved yet. */
function getModuleSync() {
  if (!resolvedModule) {
    throw new Error('fftw-wasm: module not warmed up -- call `await fftw.warmUp()` once before using the sync API');
  }
  return resolvedModule;
}

function allocDoubles(Module, arr) {
  const p = Module._malloc(arr.length * 8);
  Module.HEAPF64.set(arr, p / 8);
  return p;
}
function readDoubles(Module, p, n) {
  return Array.from(Module.HEAPF64.subarray(p / 8, p / 8 + n));
}

const FFTW_FORWARD = -1;
const FFTW_BACKWARD = 1;

// cwrap'd function cache -- see dft1dCore/r2c1dCore's comment for why this
// exists: found via a real DCP worker dispatch where a single 2D convolve
// call (composed of ~84 of these 1D calls -- N0 rows + NC columns, see
// FFTW3WasmRank2.js's decomposition) was measured at 677ms in the real
// sandboxed worker vs an estimated ~2-3ms locally in plain Node. Re-cwrap'ing
// on every call was invisible in local V8 (cwrap is cheap there) but is the
// leading suspect for that gap in whatever engine the real sandbox uses --
// caching the wrapped function once, like getModule()/getModuleSync() above
// already cache the module instance, removes the redundant work regardless
// of whether it's the actual root cause.
let cachedDftFn = null;
let cachedR2cFn = null;

// ---------------------------------------------------------------------------
// Core implementations -- pure functions of an already-resolved Module.
// Both the async (dft1d/r2c1d, await getModule() first) and sync
// (dft1dSync/r2c1dSync, getModuleSync() first) public entry points below
// call these same two functions, so there's exactly one implementation of
// the actual malloc/cwrap/free sequence to keep correct.
// ---------------------------------------------------------------------------
function dft1dCore(Module, reIn, imIn, opts) {
  if (reIn.length !== imIn.length) throw new Error('dft1d: reIn/imIn length mismatch');
  const n = reIn.length;
  if (!cachedDftFn) {
    cachedDftFn = Module.cwrap('fftw_dft_1d_wasm', 'number', ['number', 'number', 'number', 'number', 'number', 'number']);
  }
  const fn = cachedDftFn;

  const reInP = allocDoubles(Module, reIn);
  const imInP = allocDoubles(Module, imIn);
  const reOutP = Module._malloc(n * 8);
  const imOutP = Module._malloc(n * 8);
  try {
    const rc = fn(reInP, imInP, reOutP, imOutP, n, opts.inverse ? FFTW_BACKWARD : FFTW_FORWARD);
    if (rc !== 0) throw new Error(`fftw_dft_1d_wasm failed, code ${rc}`);
    return { re: readDoubles(Module, reOutP, n), im: readDoubles(Module, imOutP, n) };
  } finally {
    [reInP, imInP, reOutP, imOutP].forEach((p) => Module._free(p));
  }
}

function r2c1dCore(Module, reIn) {
  const n = reIn.length;
  const nc = Math.floor(n / 2) + 1;
  if (!cachedR2cFn) {
    cachedR2cFn = Module.cwrap('fftw_r2c_1d_wasm', 'number', ['number', 'number', 'number', 'number']);
  }
  const fn = cachedR2cFn;

  const reInP = allocDoubles(Module, reIn);
  const reOutP = Module._malloc(nc * 8);
  const imOutP = Module._malloc(nc * 8);
  try {
    const rc = fn(reInP, reOutP, imOutP, n);
    if (rc !== 0) throw new Error(`fftw_r2c_1d_wasm failed, code ${rc}`);
    return { re: readDoubles(Module, reOutP, nc), im: readDoubles(Module, imOutP, nc) };
  } finally {
    [reInP, reOutP, imOutP].forEach((p) => Module._free(p));
  }
}

let cachedR2c2dFn = null;
let cachedC2r2dFn = null;

// Native 2D r2c/c2r -- calls FFTW's own 2D planner directly (ONE
// fftw_execute per transform) instead of composing from N0+NC separate
// 1D calls the way FFTW3WasmRank2.js's rdft2R2cSync/rdft2C2rSync do.
// That composition was the dominant cost of the whole objective-function
// evaluation (~85-90% of eval time, found via direct per-eval
// instrumentation) -- these two functions are the real fix.
function r2c2dCore(Module, n0, n1, input) {
  const nc = Math.floor(n1 / 2) + 1;
  if (!cachedR2c2dFn) {
    cachedR2c2dFn = Module.cwrap('fftw_r2c_2d_wasm', 'number', ['number', 'number', 'number', 'number', 'number']);
  }
  const fn = cachedR2c2dFn;

  const inP = allocDoubles(Module, input);
  const reOutP = Module._malloc(n0 * nc * 8);
  const imOutP = Module._malloc(n0 * nc * 8);
  try {
    const rc = fn(inP, reOutP, imOutP, n0, n1);
    if (rc !== 0) throw new Error(`fftw_r2c_2d_wasm failed, code ${rc}`);
    return { re: readDoubles(Module, reOutP, n0 * nc), im: readDoubles(Module, imOutP, n0 * nc) };
  } finally {
    [inP, reOutP, imOutP].forEach((p) => Module._free(p));
  }
}

function c2r2dCore(Module, n0, n1, reIn, imIn) {
  const nc = Math.floor(n1 / 2) + 1;
  if (!cachedC2r2dFn) {
    cachedC2r2dFn = Module.cwrap('fftw_c2r_2d_wasm', 'number', ['number', 'number', 'number', 'number', 'number']);
  }
  const fn = cachedC2r2dFn;

  const reInP = allocDoubles(Module, reIn);
  const imInP = allocDoubles(Module, imIn);
  const outP = Module._malloc(n0 * n1 * 8);
  try {
    const rc = fn(reInP, imInP, outP, n0, n1);
    if (rc !== 0) throw new Error(`fftw_c2r_2d_wasm failed, code ${rc}`);
    return readDoubles(Module, outP, n0 * n1);
  } finally {
    [reInP, imInP, outP].forEach((p) => Module._free(p));
  }
}

/** Synchronous native 2D real->complex DFT. requires warmUp() first. */
function r2c2dSync(n0, n1, input) {
  return r2c2dCore(getModuleSync(), n0, n1, input);
}

/** Synchronous native 2D complex->real (inverse) DFT. requires warmUp() first. */
function c2r2dSync(n0, n1, reIn, imIn) {
  return c2r2dCore(getModuleSync(), n0, n1, reIn, imIn);
}

/**
 * Complex-to-complex 1D DFT.
 * @param {number[]} reIn - real parts, length n
 * @param {number[]} imIn - imaginary parts, length n
 * @param {Object} [opts]
 * @param {boolean} [opts.inverse] - FFTW_BACKWARD instead of FFTW_FORWARD (default false)
 * @returns {Promise<{re: number[], im: number[]}>}
 */
async function dft1d(reIn, imIn, opts = {}) {
  const Module = await getModule();
  return dft1dCore(Module, reIn, imIn, opts);
}

/**
 * Real-to-complex 1D DFT (forward only). Returns the non-redundant half
 * of the spectrum: n/2+1 bins, FFTW's standard r2c output size.
 * @param {number[]} reIn - real samples, length n
 * @returns {Promise<{re: number[], im: number[]}>}
 */
async function r2c1d(reIn) {
  const Module = await getModule();
  return r2c1dCore(Module, reIn);
}

/** Synchronous twin of dft1d -- requires warmUp() to have already resolved. */
function dft1dSync(reIn, imIn, opts = {}) {
  return dft1dCore(getModuleSync(), reIn, imIn, opts);
}

/** Synchronous twin of r2c1d -- requires warmUp() to have already resolved. */
function r2c1dSync(reIn) {
  return r2c1dCore(getModuleSync(), reIn);
}

module.exports = { dft1d, r2c1d, warmUp, dft1dSync, r2c1dSync, r2c2dSync, c2r2dSync };
