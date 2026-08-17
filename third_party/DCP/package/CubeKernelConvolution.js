module.declare(["./CalculateBeamKernel.js","./DataCube.js","./FFTW3WasmRank2.js"], function (require, exports, module) {
'use strict';

// =============================================================================
// CubeKernelConvolution.js
// High-fidelity port of src/ConvolveCube/CubeKernelConvolution.f
// (CubeKernelConvolutionMod)
//
// PORTING NOTES
// -------------
// FFT convolution via FFTW3WasmRank2.js -- the REAL compiled FFTW3 (through
// third_party/fftw-3.3.8/wasm/fftw-wasm.js's synchronous 1D primitives),
// not FFTW3JS's hand-ported engine (retired from this file -- see git
// history / FFTW3JS/ for the prior version, kept only as size-limited
// reference/verification fixtures). Since this is the actual library, not a
// curated port, it handles ANY N0xN1 natively -- no isTransformSupported
// gate or ndarray-fft fallback needed anymore.
//
// SYNC, NOT ASYNC: FFTW3WasmRank2.js exposes rdft2R2cSync/rdft2C2rSync, both
// synchronous once the wasm module has been warmed up once (see warmUp()
// re-exported below) -- this keeps buildComplexKernel/convolve2DChannel/
// cubeBeamConvolution synchronous, so amoeba/amotry/tiltedRingModelComparison
// need no changes to call them, matching this project's already bit-verified
// (against Fortran) optimizer as-is. Callers MUST await warmUp() once before
// the optimizer loop starts (typically at the top of the launcher's async
// work function) -- calling the sync functions before that throws.
//
// FFTW halfcomplex format (output of RFFT for size n):
//   out[k]   = re[k]   for k = 0..n/2
//   out[n-k] = im[k]   for k = 1..n/2-1
//   im[0] = 0, im[n/2] = 0 (DC and Nyquist are purely real)
//
// 2D FFT is performed as row-FFTs then column-FFTs, matching how FFTW
// implements r2c_2d internally (see FFTW3WasmRank2.js's header).
//
// PRECISION
// ---------
// Real compiled FFTW3, double precision throughout the transform. Verified
// against FFTW3JS (the old bit-exact-vs-real-FFTW3 hand-port) at the real
// production size (57x53): agrees to ~1e-14, far below the float32
// rounding everything here is subjected to anyway (see
// FFTW3WasmRank2.js's self-test).
// =============================================================================

const f32  = Math.fround;

const { calculate2DBeamKernel, makeWrappedArray } = require('./CalculateBeamKernel.js');

const { flatIndxCalc } = require('./DataCube.js');
// Native (rdft2R2cSyncNative/rdft2C2rSyncNative), not composed -- see
// FFTW3WasmRank2.js's header on those for why: the composed 1D-row-then-
// column path was found to be ~85-90% of the whole objective-function
// eval's cost (150-170 individual 1D FFTW calls per channel). The native
// versions call FFTW's own 2D planner directly (one fftw_execute per
// transform) and are verified round-trip-correct + matching the composed
// path exactly (see FFTW3WasmRank2.js's self-test).
const { rdft2R2cSyncNative: rdft2R2cSync, rdft2C2rSyncNative: rdft2C2rSync, warmUp } = require('./FFTW3WasmRank2.js');


// ---------------------------------------------------------------------------
// buildComplexKernel
// Builds the beam kernel's 2D FFT via the real compiled FFTW3 (see module
// header). Stored as the same { re, im, ps0, ps1 } planar-Float64Array
// shape convolve2DChannel's pointwise-multiply code expects.
// ---------------------------------------------------------------------------
function buildComplexKernel(b) {
  const ps0 = b.paddedSize[0];
  const ps1 = b.paddedSize[1];
  const n   = b.nRadialCells;
  const kSz = 2 * n + 1;

  // Pad kernel into ps0 x ps1
  const paddedKernel = new Float64Array(ps0 * ps1);
  for (let i = 1; i <= kSz; i++) {
    for (let j = 1; j <= kSz; j++) {
      const ki  = i - n - 1;
      const kj  = j - n - 1;
      paddedKernel[(i - 1) * ps1 + (j - 1)] = b.kernel[(ki + n) * kSz + (kj + n)];
    }
  }

  // Wrap
  const centKernel    = [Math.trunc(kSz / 2), Math.trunc(kSz / 2)];
  const wrappedKernel = new Float64Array(ps0 * ps1);
  makeWrappedArray([ps0, ps1], centKernel, paddedKernel, wrappedKernel);

  // rdft2R2cSync returns a single interleaved complex array (Re,Im pairs);
  // de-interleave into the planar {re,im} shape convolve2DChannel expects.
  const NC = Math.floor(ps1 / 2) + 1;
  const interleaved = rdft2R2cSync(ps0, ps1, wrappedKernel);
  const nBins = ps0 * NC;
  const re = new Float64Array(nBins);
  const im = new Float64Array(nBins);
  for (let idx = 0; idx < nBins; idx++) {
    re[idx] = interleaved[2 * idx];
    im[idx] = interleaved[2 * idx + 1];
  }
  b.complexKernel = { re, im, ps0, ps1 };

  b.complexKernelCreated = true;
}


// ---------------------------------------------------------------------------
// getConvScratch
// convolve2DChannel runs once per channel per objective-function evaluation
// (~100 channels x hundreds of evaluations per fit). ps0/ps1 are fixed for
// the whole fit, so allocate its padded scratch buffer once and cache it on
// the beam object instead of `new Float64Array(...)` on every call -- cuts
// both allocation count and GC pressure. Safe to cache on `b`: it's a fresh
// per-realization object built inside bootstrap(), never shared across
// DCP slice invocations.
// ---------------------------------------------------------------------------
function getConvScratch(b) {
  const ps0 = b.complexKernel.ps0;
  const ps1 = b.complexKernel.ps1;
  const sc = b._convScratch;
  if (sc && sc.ps0 === ps0 && sc.ps1 === ps1) return sc;

  const fresh = { ps0, ps1, padded: new Float64Array(ps0 * ps1) };
  b._convScratch = fresh;
  return fresh;
}


// ---------------------------------------------------------------------------
// convolve2DChannel
// FFT-based 2D convolution via the real compiled FFTW3 (FFTW3WasmRank2.js).
// rdft2C2rSync is UNNORMALIZED -- Fortran's own Convolve2D applies the
// normalization explicitly and sequentially (RealConvolve=RealConvolve/
// SizePad(1)/SizePad(2), double precision throughout, TwoDConvolution.f),
// so this replicates that exact two-step division rather than one combined
// /(ps0*ps1), with a single f32 rounding only at the final assignment
// (matching RealConvolve (double) -> ConvolvedArray (real) in Fortran).
// ---------------------------------------------------------------------------
function convolve2DChannel(sliceIn, nPixels, b, sliceOut) {
  const nx = nPixels[0];
  const ny = nPixels[1];
  const { re: kRe, im: kIm, ps0, ps1 } = b.complexKernel;
  const { padded } = getConvScratch(b);

  padded.fill(0);
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < ny; j++)
      padded[i * ps1 + j] = sliceIn[i * ny + j];

  const complex = rdft2R2cSync(ps0, ps1, padded);
  const NC = Math.floor(ps1 / 2) + 1;
  const nBins = ps0 * NC;
  for (let idx = 0; idx < nBins; idx++) {
    const re1 = complex[2 * idx],     im1 = complex[2 * idx + 1];
    const re2 = kRe[idx],             im2 = kIm[idx];
    complex[2 * idx]     = re1 * re2 - im1 * im2;
    complex[2 * idx + 1] = re1 * im2 + im1 * re2;
  }

  const back = rdft2C2rSync(ps0, ps1, complex);
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < ny; j++)
      sliceOut[i * ny + j] = f32(back[i * ps1 + j] / ps0 / ps1);
}


// ---------------------------------------------------------------------------
// Convolution timing accumulator -- diagnoses whether round-2 DCP worker
// time is actually going into FFT convolution (called once per
// objective-function evaluation, ~100 channels each, hundreds of
// evaluations per fit -- see getConvScratch's comment above) as opposed to
// somewhere else in the optimizer. Module-level, reset at the start of each
// bootstrap() call in bootstrap-fit-launcher.js, read at the end.
// ---------------------------------------------------------------------------
let convolveMs = 0;
let convolveCalls = 0;

function resetConvolveStats() {
  convolveMs = 0;
  convolveCalls = 0;
}

function getConvolveStats() {
  return { convolveMs, convolveCalls };
}


// ---------------------------------------------------------------------------
// cubeBeamConvolution
// ---------------------------------------------------------------------------
function cubeBeamConvolution(dc, b) {
  const t0 = Date.now();
  if (!b.complexKernelCreated) buildComplexKernel(b);

  const dh  = dc.dh;
  const nx  = dh.nPixels[0];
  const ny  = dh.nPixels[1];
  const nch = dh.nChannels;
  const sliceIn  = new Float32Array(nx * ny);
  const sliceOut = new Float32Array(nx * ny);

  for (let ch = 0; ch < nch; ch++) {
    for (let i = 0; i < nx; i++)
      for (let j = 0; j < ny; j++)
        sliceIn[i * ny + j] = dc.flux[flatIndxCalc(i, j, ch, dh)];

    convolve2DChannel(sliceIn, dh.nPixels, b, sliceOut);

    for (let i = 0; i < nx; i++)
      for (let j = 0; j < ny; j++)
        dc.flux[flatIndxCalc(i, j, ch, dh)] = sliceOut[i * ny + j];
  }
  convolveMs += Date.now() - t0;
  convolveCalls += 1;
}

module.exports = {
  cubeBeamConvolution, buildComplexKernel, convolve2DChannel,
  resetConvolveStats, getConvolveStats, warmUp,
};

if (require.main === module) {
  const { Beam2D, allocate_Beam2D } = require('./Beam.js');
  const { DataCube, allocateDataCube } = require('./DataCube.js');

  (async () => {
    // Real fftw3wasm's one-time async instantiation -- must resolve before
    // any of the (now synchronous) convolution calls below.
    await warmUp();

    const dc = new DataCube();
    const dh = dc.dh;
    dh.nPixels[0] = 32; dh.nPixels[1] = 32; dh.nChannels = 4;
    dh.pixelSize[0] = f32(-1.0); dh.pixelSize[1] = f32(1.0);
    dh.channelSize = f32(10.0);
    dh.refLocation[0] = f32(16.0); dh.refLocation[1] = f32(16.0);
    dh.refLocation[2] = f32(2.0);
    dh.refVal[0] = f32(0.0); dh.refVal[1] = f32(0.0); dh.refVal[2] = f32(1000.0);
    allocateDataCube(dc);

    dc.flux[flatIndxCalc(16, 16, 1, dh)] = f32(1.0);
    const totalBefore = dc.flux.reduce((a, v) => a + v, 0);

    const b = new Beam2D();
    b.beamFWHM = f32(3.0); b.sigmaLengths = f32(3.0);
    b.pixelSize[0] = f32(-1.0); b.pixelSize[1] = f32(1.0);
    allocate_Beam2D(b, new Int32Array([32, 32]));
    calculate2DBeamKernel(b, b.pixelSize);
    cubeBeamConvolution(dc, b);

    const totalAfter = dc.flux.reduce((a, v) => f32(a + v), f32(0));
    const centre     = dc.flux[flatIndxCalc(16, 16, 1, dh)];
    const neighbour  = dc.flux[flatIndxCalc(17, 16, 1, dh)];
    const empty      = dc.flux[flatIndxCalc(16, 16, 0, dh)];

    // NOTE on precision: the "(Fortran: ...)" reference values below are only
    // recorded to 6 decimals (their original provenance predates this
    // comment and wasn't preserved at higher precision) -- toFixed(6)
    // matching here is a coarse sanity check only, NOT a measurement of this
    // engine's actual precision ceiling. For a rigorous, bit-level answer
    // (checked directly against real compiled FFTW3, not just a 6-decimal
    // string), see FFTW3WasmRank2.js's own self-test: delta-function inputs
    // are bit-exact; general/random inputs agree to within a few tens of
    // ULP in double precision (~1e-14 relative), far below the float32
    // rounding this whole pipeline applies to every result anyway.
    console.log('=== cubeBeamConvolution (real fftw3wasm engine) ===');
    console.log('total flux before:', totalBefore.toFixed(6));
    console.log('total flux after: ', totalAfter.toFixed(6), '(expect ~1.0)');
    console.log('centre pixel:     ', centre.toFixed(6), '(Fortran: 0.098133, only recorded to 6dp -- see note above)');
    console.log('neighbour pixel:  ', neighbour.toFixed(6), '(Fortran: 0.072111, only recorded to 6dp -- see note above)');
    console.log('empty channel:    ', empty.toFixed(6), '(expect 0)');
    console.log('flux conserved:   ', Math.abs(totalAfter - totalBefore) < 0.01 ? 'OK' : 'FAIL');
    console.log('centre is max:    ', centre > neighbour ? 'OK' : 'FAIL');
    console.log('channels isolated:', empty === 0 ? 'OK' : 'FAIL');
  })();
}

});
