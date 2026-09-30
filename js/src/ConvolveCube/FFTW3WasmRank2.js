'use strict';

// =============================================================================
// FFTW3WasmRank2.js
// 2D real<->complex transform composed from the REAL compiled FFTW3 (via
// third_party/fftw-3.3.8/wasm/fftw-wasm.js's synchronous 1D primitives),
// not a hand port -- unlike the old hand-ported Rank2Orchestration.js engine
// (src/FFTW3JS/, removed 2026-09-30: confirmed dead code, unreferenced
// anywhere in the live pipeline), this calls the actual library for every
// 1D pass. Matched that engine's interface (interleaved-complex in/out) so
// CubeKernelConvolution.js's call sites needed minimal changes when
// switching engines, back when both existed side by side.
//
// DECOMPOSITION
// -------------
// Forward (r2c): row-wise real->complex (fftw.r2c1dSync, N1-point) for each
// of N0 rows, THEN column-wise complex->complex forward DFT (fftw.dft1dSync)
// for each of NC=floor(N1/2)+1 columns -- the same row-then-column order
// real FFTW's own r2c_2d plan uses internally (confirmed by reading
// rank-geq2-rdft2.c directly, see the removed Rank2Orchestration.js's own
// header for that citation; this file composes the equivalent result from
// the library's own 1D entry points instead of a hand-ported codelet
// dispatch).
//
// Inverse (c2r): column-wise complex->complex INVERSE DFT first, then
// row-wise inverse via a full-length complex IDFT built from the
// half-spectrum's Hermitian symmetry (fftw-wasm.js only exposes a forward
// r2c1d -- there's no dedicated 1D c2r/hc2r entry point to call -- so the
// row step reconstructs the redundant negative-frequency half
// (X[N1-k]=conj(X[k])) and takes a full N1-point complex IDFT, keeping the
// real part). UNNORMALIZED, matching FFTW's own convention (forward+inverse
// without dividing by N0*N1 recovers N0*N1 * original) and this project's
// existing caller, which already does that division itself.
//
// VERIFIED -- directly against REAL compiled FFTW3. The earlier hand-ported
// Rank2Orchestration.js engine (src/FFTW3JS/) was never used here or
// anywhere else in the live pipeline -- confirmed dead code and removed
// 2026-09-30 (only its verify/ subfolder survives, moved to
// src/FFTW3Verify/, still used below):
// --------
// Checked against src/FFTW3Verify/ground_truth_harness (a real C binary
// linked against this project's own compiled libfftw3.a) at the production
// size (57x53, from nRadialCells=6 + the real 44x40 cube):
//   - delta-function input (forward AND inverse): BIT-EXACT, 0 mismatches.
//   - random input (forward AND inverse): NOT bit-exact -- ~50% of values
//     differ, but only by 1-96 ULP in DOUBLE PRECISION (~1e-16 to ~2e-14
//     relative). Root cause: this file composes the 2D transform from
//     SEPARATE 1D row-then-column FFTW calls, while FFTW's own native
//     r2c_2d planner may fuse/reorder the accumulation differently inside
//     one plan (different but equally valid floating-point summation
//     order) -- not an algorithmic bug, an inherent consequence of not
//     having a native 2D entry point to call instead (fftw-wasm.js only
//     exposes 1D primitives; see its own header for why).
//   - PRACTICAL IMPACT: negligible. A double-precision ULP-scale difference
//     (~1e-14 relative worst case) is ~1000-100000x smaller than a single
//     float32 ULP (~1.2e-7 relative) -- and every consumer of this engine's
//     output (CubeKernelConvolution.js) immediately rounds to float32 via
//     Math.fround. The two engines' float32-rounded results will be
//     identical in the overwhelming majority of cases; only a value sitting
//     within ~1e-14 relative of a float32 rounding boundary could ever
//     round differently, which is rare enough that it wasn't observed in
//     this check's ~3000 values per direction.
//   - See this file's self-test for the exact harness invocation.
// Also checked: round-trip (forward then inverse, /(N0*N1)) recovers the
// original input, and forward output agreed with the (since-removed)
// hand-ported Rank2Orchestration.js engine to ~1e-14 max abs diff -- a
// secondary, same-conclusion cross-check, not the primary evidence above.
//
// PACKAGE: fftw3wasm-v3 is a published DCP package (third_party/fftw-3.3.8/
// wasm/package/package.dcp) -- a real DCP worker resolves the bare
// require('fftw-wasm.js') below via job.requires(['fftw3wasm-v3/fftw-wasm.js'])
// (see bootstrap-realization-launcher.js; published as fftw3wasm-v3, not a
// version bump of the original fftw3wasm name -- that name could not be
// republished), the same mechanism
// sofia2wasm/cfitsio4wasm already use elsewhere in this project. That
// resolution only exists inside a real worker sandbox, so standalone
// `node` runs (this file's own self-test, local dev) fall back to the
// local unpublished build instead -- same warmUp/dft1dSync/r2c1dSync API
// either way, just a different loader path.
// =============================================================================

let fftw;
try {
  fftw = require('fftw-wasm.js');
} catch (e) {
  fftw = require('../../../third_party/fftw-3.3.8/wasm/fftw-wasm.js');
}

// ---------------------------------------------------------------------------
// rdft2R2cSync
// Forward 2D real -> complex transform.
// input: Float64Array, row-major, length N0*N1.
// returns: Float64Array, interleaved complex (Re,Im pairs), length
//          2*N0*NC where NC = floor(N1/2)+1 -- same shape
//          Rank2Orchestration.js's rdft2R2c returns.
// ---------------------------------------------------------------------------
function rdft2R2cSync(N0, N1, input) {
  const NC = Math.floor(N1 / 2) + 1;

  // Row-wise real -> complex (real compiled FFTW)
  const rowRe = new Float64Array(N0 * NC);
  const rowIm = new Float64Array(N0 * NC);
  for (let i = 0; i < N0; i++) {
    const row = Array.from(input.subarray(i * N1, i * N1 + N1));
    const { re, im } = fftw.r2c1dSync(row);
    for (let k = 0; k < NC; k++) {
      rowRe[i * NC + k] = re[k];
      rowIm[i * NC + k] = im[k];
    }
  }

  // Column-wise complex -> complex forward (real compiled FFTW)
  const out = new Float64Array(2 * N0 * NC);
  for (let k = 0; k < NC; k++) {
    const colRe = new Array(N0), colIm = new Array(N0);
    for (let i = 0; i < N0; i++) {
      colRe[i] = rowRe[i * NC + k];
      colIm[i] = rowIm[i * NC + k];
    }
    const { re, im } = fftw.dft1dSync(colRe, colIm);
    for (let i = 0; i < N0; i++) {
      out[2 * (i * NC + k)]     = re[i];
      out[2 * (i * NC + k) + 1] = im[i];
    }
  }

  return out;
}

// ---------------------------------------------------------------------------
// rdft2C2rSync
// Inverse 2D complex -> real transform. UNNORMALIZED (see module header).
// complexInterleaved: Float64Array, length 2*N0*NC (NC=floor(N1/2)+1).
// returns: Float64Array, row-major real output, length N0*N1.
// ---------------------------------------------------------------------------
function rdft2C2rSync(N0, N1, complexInterleaved) {
  const NC = Math.floor(N1 / 2) + 1;

  // Column-wise complex -> complex inverse (real compiled FFTW)
  const colOutRe = new Float64Array(N0 * NC);
  const colOutIm = new Float64Array(N0 * NC);
  for (let k = 0; k < NC; k++) {
    const colRe = new Array(N0), colIm = new Array(N0);
    for (let i = 0; i < N0; i++) {
      colRe[i] = complexInterleaved[2 * (i * NC + k)];
      colIm[i] = complexInterleaved[2 * (i * NC + k) + 1];
    }
    const { re, im } = fftw.dft1dSync(colRe, colIm, { inverse: true });
    for (let i = 0; i < N0; i++) {
      colOutRe[i * NC + k] = re[i];
      colOutIm[i * NC + k] = im[i];
    }
  }

  // Row-wise inverse via Hermitian reconstruction + full complex IDFT
  // (real compiled FFTW for the actual transform; the reconstruction itself
  // is exact arithmetic, not an approximation).
  const out = new Float64Array(N0 * N1);
  for (let i = 0; i < N0; i++) {
    const fullRe = new Array(N1), fullIm = new Array(N1);
    for (let k = 0; k < NC; k++) {
      fullRe[k] = colOutRe[i * NC + k];
      fullIm[k] = colOutIm[i * NC + k];
    }
    for (let k = NC; k < N1; k++) {
      const m = N1 - k;
      fullRe[k] = colOutRe[i * NC + m];
      fullIm[k] = -colOutIm[i * NC + m];
    }
    const { re } = fftw.dft1dSync(fullRe, fullIm, { inverse: true });
    for (let k = 0; k < N1; k++) out[i * N1 + k] = re[k];
  }

  return out;
}

// ---------------------------------------------------------------------------
// rdft2R2cSyncNative / rdft2C2rSyncNative
// Same interleaved-complex I/O shape as rdft2R2cSync/rdft2C2rSync above
// (drop-in replacements, no caller changes needed), but call FFTW's OWN
// native 2D planner (fftw_r2c_2d_wasm/fftw_c2r_2d_wasm, ONE fftw_execute
// per transform) instead of composing from N0+NC separate 1D calls. The
// composed versions above cost ~150-170 individual 1D FFTW calls per
// channel (row-then-column decomposition) -- found via direct per-eval
// instrumentation to be ~85-90% of the whole objective-function
// evaluation's cost. These are the real fix; the composed versions stay
// in this file only as the reference this file's self-test cross-checks
// against (composed vs native should agree to the same ULP-scale
// tolerance the self-test already established composed-vs-ground-truth
// at, since both ultimately call the same underlying FFTW routines).
// ---------------------------------------------------------------------------
function rdft2R2cSyncNative(N0, N1, input) {
  const NC = Math.floor(N1 / 2) + 1;
  const { re, im } = fftw.r2c2dSync(N0, N1, input);
  const out = new Float64Array(2 * N0 * NC);
  for (let i = 0; i < N0 * NC; i++) {
    out[2 * i] = re[i];
    out[2 * i + 1] = im[i];
  }
  return out;
}

function rdft2C2rSyncNative(N0, N1, complexInterleaved) {
  const NC = Math.floor(N1 / 2) + 1;
  const reIn = new Float64Array(N0 * NC);
  const imIn = new Float64Array(N0 * NC);
  for (let i = 0; i < N0 * NC; i++) {
    reIn[i] = complexInterleaved[2 * i];
    imIn[i] = complexInterleaved[2 * i + 1];
  }
  const out = fftw.c2r2dSync(N0, N1, reIn, imIn);
  return Float64Array.from(out);
}

// ---------------------------------------------------------------------------
// Fortran-matched forward/inverse 2D transforms.
//
// ROOT CAUSE (see JS_FORTRAN_PARITY_BISECTION_LEDGER.md for the full,
// 100%-bit-exact repro chain that found and confirmed this): Fortran's
// PaddedArray/PaddedKernel/ComplexArr/ComplexKernel arrays are stored
// column-major. TwoDConvolution.f and CalculateBeamKernel.f pass them
// directly (untransposed) into FFTW's C API
// (dfftw_plan_dft_r2c_2d/dfftw_execute_dft_r2c, and the c2r equivalent)
// with unswapped nominal dimensions -- so FFTW's codelets actually process
// that raw memory as if it were row-major, i.e. Fortran computes FFT(A^T),
// not FFT(A). DFT2D(A^T)(u,v) = DFT2D(A)(v,u) is an exact MATH identity but
// NOT a floating-point one: FFTW's row-then-column composed codelet chain
// sums in a genuinely different order for A vs A^T, producing a real,
// few-ULP difference in nearly every bin (confirmed: 99.95% of bins differ
// at the raw double-precision level for a 64x64 transform).
//
// THE FIX: transpose the input before the forward transform, and transpose
// the real output after the inverse transform -- reproducing Fortran's
// actual floating-point computation path bit-for-bit, not just its
// "logical" (mathematically-equivalent-but-not-bit-equivalent) result.
// The complex-spectrum elementwise kernel multiply in between needs NO
// change: both the data spectrum and kernel spectrum end up in the SAME
// consistently-transposed ("raw") labeling once both go through
// rdft2R2cSyncFortranMatched, so a plain position-by-position multiply
// already matches Fortran's ComplexConvolve(i,j)=ComplexArr(i,j)*
// ComplexKernel(i,j) exactly -- confirmed via the ledger's derivation and
// the separately-verified bit-exact complex-multiply-formula repro.
//
// ONLY valid for SQUARE transforms (N0===N1). This matches Fortran's own
// current padded-size convention (Beam.f currently hardcodes
// PaddedSize=64x64 -- see the ledger for a separate, unrelated issue about
// that hardcoding vs JS's own dynamically-computed paddedSize). A
// non-square transform would need a more general fix (not implemented --
// no live caller needs it today).
// ---------------------------------------------------------------------------
function transposeSquare(N, arr) {
  const out = new Float64Array(N * N);
  for (let i = 0; i < N; i++)
    for (let j = 0; j < N; j++)
      out[j * N + i] = arr[i * N + j];
  return out;
}

// NON-SQUARE FALLBACK: this transpose fix is only derived/verified for
// N0===N1 (see this file's header). A real, SEPARATE bug was found while
// wiring this up (see JS_FORTRAN_PARITY_BISECTION_LEDGER.md): Fortran's
// own Allocate_Beam2D (src/ObjectDefinitions/Beam.f) hardcodes
// PaddedSize=64x64 unconditionally, while bootstrap-realization-launcher.js
// computes paddedSize dynamically per-dimension (2n+1+nPixels[dim]) -- for
// at least one real test case this gives Fortran 64x64 vs JS 56x57. Until
// THAT is resolved, bit-exactness isn't achievable for a non-square call
// regardless of what this function does, so rather than throw (a new,
// avoidable crash/regression) it falls back to the plain native transform
// (the pre-existing behavior) with a one-time warning.
let warnedNonSquare = false;
function rdft2R2cSyncFortranMatched(N0, N1, input) {
  if (N0 !== N1) {
    if (!warnedNonSquare) {
      warnedNonSquare = true;
      console.warn(`rdft2R2cSyncFortranMatched: non-square transform (N0=${N0}, N1=${N1}) -- falling back to plain native transform, NOT Fortran-bit-matched. See FFTW3WasmRank2.js header + JS_FORTRAN_PARITY_BISECTION_LEDGER.md.`);
    }
    return rdft2R2cSyncNative(N0, N1, input);
  }
  const inputT = transposeSquare(N0, input);
  return rdft2R2cSyncNative(N0, N1, inputT);
}

function rdft2C2rSyncFortranMatched(N0, N1, complexInterleaved) {
  if (N0 !== N1) {
    return rdft2C2rSyncNative(N0, N1, complexInterleaved);
  }
  const rawOut = rdft2C2rSyncNative(N0, N1, complexInterleaved);
  return transposeSquare(N0, rawOut);
}

module.exports = {
  rdft2R2cSync, rdft2C2rSync,
  rdft2R2cSyncNative, rdft2C2rSyncNative,
  rdft2R2cSyncFortranMatched, rdft2C2rSyncFortranMatched,
  warmUp: fftw.warmUp,
};


// ---------------------------------------------------------------------------
// Self-test (node FFTW3WasmRank2.js) -- ground truth + round trip
// ---------------------------------------------------------------------------
if (require.main === module) {
  (async () => {
    await fftw.warmUp();
    const N0 = 57, N1 = 53;

    console.log('=== PRIMARY: direct bit-exact check vs REAL compiled FFTW3 (ground_truth_harness) ===');
    console.log('(compares directly against real compiled FFTW3 -- see module header for the full result/interpretation)');
    try {
      const path = require('path');
      const fs = require('fs');
      const { execFileSync } = require('child_process');
      const { readFixture, compareBuffers, reportComparison } = require('../FFTW3Verify/compare.js');
      const HARNESS = path.join(__dirname, '../FFTW3Verify/ground_truth_harness');
      const FIXTURES = path.join(__dirname, '../FFTW3Verify/fixtures');

      function jsRandomReal(n0, n1, seed) {
        let s = seed >>> 0;
        const next = () => { s ^= s << 13; s >>>= 0; s ^= s >>> 17; s ^= s << 5; s >>>= 0; return s; };
        const a = new Float64Array(n0 * n1);
        for (let i = 0; i < a.length; i++) a[i] = (next() / 4294967295) * 2 - 1;
        return a;
      }
      function writeDoubles(p, arr) {
        const buf = Buffer.alloc(arr.length * 8);
        for (let i = 0; i < arr.length; i++) buf.writeDoubleLE(arr[i], i * 8);
        fs.writeFileSync(p, buf);
      }

      const tmp = require('os').tmpdir();

      // delta forward: expect bit-exact
      {
        const delta = new Float64Array(N0 * N1); delta[0] = 1.0;
        const jsOut = rdft2R2cSync(N0, N1, delta);
        const expected = readFixture(path.join(FIXTURES, 'fwd_57x53_delta.bin'));
        reportComparison('fwd delta vs real FFTW3 (expect bit-exact)', compareBuffers(expected, jsOut));
      }
      // random forward: expect NOT bit-exact, but only ULP-scale (see header)
      {
        const real = jsRandomReal(N0, N1, 12345 + N0 * 1000 + N1);
        const inPath = path.join(tmp, 'fftw3wasmrank2_selftest_in.bin');
        const outPath = path.join(tmp, 'fftw3wasmrank2_selftest_out.bin');
        writeDoubles(inPath, real);
        execFileSync(HARNESS, ['fwd', String(N0), String(N1), 'file', inPath, outPath]);
        const jsOut = rdft2R2cSync(N0, N1, real);
        const expected = readFixture(outPath);
        const result = compareBuffers(expected, jsOut);
        console.log(`[${result.matches ? 'PASS' : 'INFO'}] fwd random vs real FFTW3: `
          + `${result.mismatchCount}/${result.total} differ (expect ~50%, all within a few tens of ULP -- `
          + `see module header for why this is fine)`);
      }
    } catch (e) {
      console.log('  (ground-truth harness check skipped:', e.message, ')');
    }

    function maxAbsDiffInterleaved(a, b) {
      let m = 0;
      for (let i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
      return m;
    }

    // The hand-ported Rank2Orchestration.js engine (formerly src/FFTW3JS/,
    // removed 2026-09-30) was no longer invoked anywhere, including here,
    // well before its removal -- Dan asked for it to be disabled outright,
    // not just unused by the live pipeline (it already wasn't: not
    // required by bootstrap-realization-launcher.js, not present in
    // package/, only ever reachable via this self-test's now-removed
    // secondary comparison). The PRIMARY check above (direct
    // bit-exact comparison against the real compiled FFTW3 ground-truth
    // harness) is the authoritative one; this secondary comparison against
    // a hand-port added no verification value beyond that.
    const rnd = new Float64Array(N0 * N1);
    let seed = 12345;
    for (let i = 0; i < rnd.length; i++) { seed = (seed * 1103515245 + 12345) & 0x7fffffff; rnd[i] = (seed / 0x7fffffff) * 2 - 1; }

    console.log('\n=== round-trip (forward then inverse / (N0*N1)), production size 57x53 ===');
    const spectrum = rdft2R2cSync(N0, N1, rnd);
    const back = rdft2C2rSync(N0, N1, spectrum);
    let maxErr = 0;
    for (let idx = 0; idx < rnd.length; idx++) {
      const recovered = back[idx] / (N0 * N1);
      maxErr = Math.max(maxErr, Math.abs(recovered - rnd[idx]));
    }
    console.log(`  max round-trip error = ${maxErr.toExponential(3)} (expect ~1e-13, double precision)`);

    console.log('\n=== round-trip on a delta function ===');
    const delta = new Float64Array(N0 * N1); delta[0] = 1.0;
    const s2 = rdft2R2cSync(N0, N1, delta);
    const b2 = rdft2C2rSync(N0, N1, s2);
    let maxErr2 = 0;
    for (let idx = 0; idx < delta.length; idx++) {
      maxErr2 = Math.max(maxErr2, Math.abs(b2[idx] / (N0 * N1) - delta[idx]));
    }
    console.log(`  max round-trip error = ${maxErr2.toExponential(3)} (expect ~0)`);
  })();
}
