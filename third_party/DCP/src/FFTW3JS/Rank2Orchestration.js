'use strict';

// =============================================================================
// Rank2Orchestration.js
// Faithful JS port of FFTW3's 2D r2c/c2r orchestration -- how the row (1D
// r2hc/hc2r) and column (1D full complex DFT) sub-transforms from
// GenericSolver1D.js combine into a full 2D real<->complex transform.
// Ports rdft/rank-geq2-rdft2.c's apply_r2hc / apply_hc2r (the actual
// execution order and re/im wiring FFTW3 uses -- confirmed via reading
// that file directly, not assumed) plus rdft/rdft2-rdft.c's hc2c/c2hc
// bridging (already fused into GenericSolver1D.js's rdftGenericR2hc /
// rdftGenericHc2r, see that file's header for why fusing them is safe:
// hc2c/c2hc do zero arithmetic).
//
// EXECUTION ORDER (confirmed by reading rank-geq2-rdft2.c's apply_r2hc/
// apply_hc2r, NOT assumed from the printed plan text, which always shows
// cldr before cldc regardless of actual execution order):
//   forward (r2c): row transform (r2hc) FIRST, then column transform
//                  (full complex DFT) SECOND, in place on the row output.
//   inverse (c2r): column transform (full complex DFT) FIRST -- with real
//                  and imaginary arrays SWAPPED relative to forward (the
//                  standard "get an IDFT by swapping re/im roles" trick;
//                  rank-geq2-rdft2.c's own comment: "HC2R must swap re/im
//                  parts to get IDFT") -- then row transform (hc2r) SECOND.
//
// DATA LAYOUT (this port's own convention, matching
// verify/ground_truth_harness.c's direct C-API fftw_plan_dft_r2c_2d(N0, N1,
// ...) usage -- NOT necessarily Fortran's SizePad(1)/SizePad(2) order,
// which is reversed by FFTW's legacy Fortran interface; that mapping is
// resolved separately when this engine gets wired into
// CubeKernelConvolution.js):
//   real array:    row-major, in[i*N1 + j],       i=0..N0-1, j=0..N1-1
//   complex array: row-major, out[i*NC + k] (complex index),
//                  interleaved as out2[2*(i*NC+k)]=Re, +1=Im,
//                  where NC = floor(N1/2)+1, i=0..N0-1, k=0..NC-1
//
// "Row" = transform along the N1 axis (FFTW's r2hc/hc2r, repeated N0
// times); "column" = transform along the N0 axis (FFTW's full complex
// dft-generic, repeated NC times) -- matching the plan output structure
// confirmed via fftw_print_plan (e.g. verify/fixtures/plan_47x43.txt:
// "rdft-generic-r2hc-43 x47" = size-43 row transform done 47 times,
// "dft-generic-47 x22" = size-47 column transform done 22=NC times).
//
// Source: third_party/fftw-3.3.8/rdft/rank-geq2-rdft2.c, FFTW3 (c) 2003,
// 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const { rdftGenericR2hc, rdftGenericHc2r, dftGeneric } = require('./GenericSolver1D.js');
const { dftRader37 } = require('./RaderSolver.js');
const { dftRadix3CT } = require('./CompositeSolver1D.js');
const { isPrime } = require('./PlanTable.js');
const { dftComposite, isFullyPortedComplex } = require('./Composite1D.js');
const { dftR2HCInterleaved, isFullyPortedR2HC, dftHC2RFromInterleaved, isFullyPortedHC2R } = require('./RealEngine1D.js');

// Column (full complex, size N0) dispatch: FFTW3's planner picks
// dft-rader-37 over dft-generic-37 specifically for N=37 (36=6*6 factors
// favorably -- see RaderSolver.js's header), dft-ct-dit/3 (radix-3 CT) for
// N=3*p where p is an odd prime (e.g. 57=3*19 -- see CompositeSolver1D.js's
// header), and dft-generic for every other odd prime this pipeline's real
// galaxy sizes produce (43, 47, ...). This dispatch must match FFTW3's
// actual per-size solver choice, not just "any correct DFT", since the
// whole point is bit-exactness against the real compiled library.
//
// isFullyPortedComplex() is checked FIRST, ahead of the two pre-existing
// special cases below -- found the hard way: dft1D used to check
// `n%3===0 && isPrime(n/3)` before anything else, which was written to
// catch n=57 (3*19) but ALSO wrongly matches n=9 (3*3, since isPrime(3) is
// true) and routes it into dftRadix3CT/dftGeneric -- even though real
// FFTW3's actual choice for a size-9 column is the DIRECT n1_9 codelet
// (confirmed via fftw_print_plan), which Composite1D already predicts
// correctly. isFullyPortedComplex(57) is false (57's m=19 sub-transform has
// no ported codelet, so it correctly defers to the old, fixture-verified
// dftRadix3CT path) so this reordering doesn't touch the one case the old
// special case actually needs to keep handling.
//
// The two special cases below predate Composite1D.js's general engine and
// stay exactly as they are otherwise -- proven, production-verified,
// untouched -- for whatever isFullyPortedComplex() doesn't yet cover.
function dft1D(n, ri, riOff, ii, iiOff, is, ro, roOff, io, ioOff, os) {
  if (isFullyPortedComplex(n)) {
    return dftComposite(n, ri, riOff, ii, iiOff, is, ro, roOff, io, ioOff, os);
  }
  if (n === 37) return dftRader37(n, ri, riOff, ii, iiOff, is, ro, roOff, io, ioOff, os);
  if (n % 3 === 0 && n / 3 > 1 && isPrime(n / 3)) {
    const [reOut, imOut] = dftRadix3CT(n, sliceStrided(ri, riOff, is, n), sliceStrided(ii, iiOff, is, n));
    for (let i = 0; i < n; i++) {
      ro[roOff + i * os] = reOut[i];
      io[ioOff + i * os] = imOut[i];
    }
    return;
  }
  return dftGeneric(n, ri, riOff, ii, iiOff, is, ro, roOff, io, ioOff, os);
}

function sliceStrided(arr, off, stride, n) {
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) out[i] = arr[off + i * stride];
  return out;
}

// ---------------------------------------------------------------------------
// rdft2R2c -- forward real-to-complex 2D transform.
// real: Float64Array, length N0*N1, row-major (real[i*N1+j]).
// Returns: Float64Array, length N0*NC*2 (interleaved complex), NC=floor(N1/2)+1.
// Unnormalized, matching dfftw_execute_dft_r2c (Fortran applies its own
// explicit normalization after -- see CubeKernelConvolution.js).
// ---------------------------------------------------------------------------
function rdft2R2c(N0, N1, real) {
  const NC = Math.floor(N1 / 2) + 1;
  const complex = new Float64Array(N0 * NC * 2);

  // Row pass (apply_r2hc's cldr step): for each of the N0 rows, a size-N1
  // r2hc transform, real row (stride 1) -> complex row (interleaved,
  // complex-stride 1 meaning 2 reals per bin -- rdftGenericR2hc's `os`
  // parameter is in units of "reals per unit increment of bin index", so
  // os=2 for tightly-packed interleaved complex output).
  //
  // RealEngine1D's dftR2HCInterleaved is used ONLY when isFullyPortedR2HC
  // confirms every codelet the real FFTW3 planner would pick for N1 is
  // ported (R2HC/forward direction only so far -- see RealEngine1D.js's
  // header); rdftGenericR2hc (already independently verified) remains the
  // default otherwise.
  const useEngine = isFullyPortedR2HC(N1);
  for (let i = 0; i < N0; i++) {
    if (useEngine) {
      dftR2HCInterleaved(N1, real, i * N1, 1, complex, i * NC * 2, 2);
    } else {
      rdftGenericR2hc(N1, real, i * N1, 1, complex, i * NC * 2, 2);
    }
  }

  // Column pass (apply_r2hc's cldc step): for each of the NC frequency
  // columns, a size-N0 full complex DFT, in place, on the (strided) column
  // of the row-transformed array. Column entries are NC complex bins apart
  // in row-major layout, i.e. stride NC*2 reals per row increment.
  const ri = new Float64Array(N0), ii = new Float64Array(N0);
  const ro = new Float64Array(N0), io = new Float64Array(N0);
  for (let k = 0; k < NC; k++) {
    for (let i = 0; i < N0; i++) {
      const off = (i * NC + k) * 2;
      ri[i] = complex[off];
      ii[i] = complex[off + 1];
    }
    dft1D(N0, ri, 0, ii, 0, 1, ro, 0, io, 0, 1);
    for (let i = 0; i < N0; i++) {
      const off = (i * NC + k) * 2;
      complex[off] = ro[i];
      complex[off + 1] = io[i];
    }
  }

  return complex;
}

// ---------------------------------------------------------------------------
// rdft2C2r -- inverse complex-to-real 2D transform.
// complex: Float64Array, length N0*NC*2 (interleaved), NC=floor(N1/2)+1.
// Returns: Float64Array, length N0*N1, row-major.
// Unnormalized, matching dfftw_execute_dft_c2r (Fortran applies its own
// explicit /SizePad(1)/SizePad(2) normalization after -- replicate that at
// the call site, not here -- see CubeKernelConvolution.js).
//
// NOTE: per FFTW3 semantics, a 2D c2r transform destroys its complex input
// array (see ground_truth_harness.c's header comment on FFTW_PRESERVE_INPUT
// -- Fortran's own pipeline has always effectively used DESTROY_INPUT
// semantics here too, confirmed by tracing the legacy Fortran interface).
// This port copies rather than destroying the caller's array (JS has no
// "in place on caller's buffer" concept worth fighting for here), which is
// numerically equivalent -- the caller just shouldn't assume `complex`
// survives if it later switches to a true in-place variant for performance.
// ---------------------------------------------------------------------------
function rdft2C2r(N0, N1, complexIn) {
  const NC = Math.floor(N1 / 2) + 1;
  const complex = Float64Array.from(complexIn); // working copy

  // Column pass FIRST (apply_hc2r's cldc step), with re/im SWAPPED to turn
  // the forward complex-DFT codelet into an inverse one (rank-geq2-rdft2.c:
  // "HC2R must swap re/im parts to get IDFT" -- cldc->apply(cldc, ci, cr,
  // ci, cr), i.e. what's normally the imaginary array is passed as "real"
  // and vice versa, for both input AND output).
  const ri = new Float64Array(N0), ii = new Float64Array(N0);
  const ro = new Float64Array(N0), io = new Float64Array(N0);
  for (let k = 0; k < NC; k++) {
    for (let i = 0; i < N0; i++) {
      const off = (i * NC + k) * 2;
      ri[i] = complex[off + 1]; // swapped: "real" input = actual imaginary
      ii[i] = complex[off];     // swapped: "imag" input = actual real
    }
    dft1D(N0, ri, 0, ii, 0, 1, ro, 0, io, 0, 1);
    for (let i = 0; i < N0; i++) {
      const off = (i * NC + k) * 2;
      complex[off + 1] = ro[i]; // swapped back on write, same convention
      complex[off] = io[i];
    }
  }

  // Row pass SECOND (apply_hc2r's cldr step): size-N1 hc2r transform per
  // row, complex row (interleaved, stride 2) -> real row (stride 1).
  //
  // RealEngine1D's dftHC2RFromInterleaved is used ONLY when
  // isFullyPortedHC2R confirms every codelet the real FFTW3 planner would
  // pick for N1 is ported; rdftGenericHc2r (already independently verified)
  // remains the default otherwise -- same gating pattern as rdft2R2c's row
  // pass.
  const real = new Float64Array(N0 * N1);
  const useEngineHC2R = isFullyPortedHC2R(N1);
  for (let i = 0; i < N0; i++) {
    if (useEngineHC2R) {
      dftHC2RFromInterleaved(N1, complex, i * NC * 2, 2, real, i * N1, 1);
    } else {
      rdftGenericHc2r(N1, complex, i * NC * 2, 2, real, i * N1, 1);
    }
  }

  return real;
}

module.exports = { rdft2R2c, rdft2C2r };

// ---------------------------------------------------------------------------
// Self-test (node Rank2Orchestration.js) -- smoke test only, NOT a
// substitute for the real bit-exact comparison against
// verify/ground_truth_harness.c fixtures (see verify/run_phase1_tests.js).
// Round-trips a small 2D real array and checks flux conservation (DC bin
// should equal the array sum) plus that the round trip returns N0*N1 times
// the original.
// ---------------------------------------------------------------------------
if (require.main === module) {
  const N0 = 5, N1 = 7;
  const real = new Float64Array(N0 * N1);
  for (let i = 0; i < real.length; i++) real[i] = i + 1;

  console.log('=== rdft2R2c / rdft2C2r round-trip (5x7) ===');
  const complex = rdft2R2c(N0, N1, real);
  const sum = real.reduce((a, b) => a + b, 0);
  console.log('DC bin (complex[0]):', complex[0], '(expect sum =', sum, ')');

  const back = rdft2C2r(N0, N1, complex);
  const scaled = Array.from(back).map((v) => v / (N0 * N1));
  const maxDiff = Math.max(...scaled.map((v, idx) => Math.abs(v - real[idx])));
  console.log('round-trip max |diff| after /N0N1 scaling:', maxDiff, '(expect ~0)');
}
