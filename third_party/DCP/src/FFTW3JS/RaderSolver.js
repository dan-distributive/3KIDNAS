'use strict';

// =============================================================================
// RaderSolver.js
// Faithful JS port of FFTW3's Rader's-algorithm solver, for the ONE concrete
// case this pipeline's real galaxy sizes actually need: N=37. Confirmed via
// fftw_print_plan (verify/ground_truth_harness.c) that FFTW3's planner picks
// `dft-rader-37` (not the plain `dft-generic-37` GenericSolver1D.js already
// covers) for the column (full complex) sub-transform whenever N0=37 --
// because 36 = N-1 factors favorably (6*6) into small radices, which makes
// Rader's algorithm cost-competitive under FFTW3's planner in a way it isn't
// for 43 (42=2*3*7) or 47 (46=2*23).
//
// Rader's algorithm turns a prime-N DFT into a cyclic convolution of size
// N-1, computed as: permute input by a primitive-root generator, forward-DFT
// the permuted sequence (size N-1), multiply pointwise by a precomputed
// "omega" array (itself the DFT of the convolution kernel), forward-DFT
// again (FFTW reuses one forward-DFT plan twice rather than building a
// separate inverse -- see the "fft(x*)* = ifft(x)" identity noted in
// dft/rader.c's header comment), then un-permute by the inverse generator.
//
// For N=37, N-1=36=6*6 is itself computed via genuine (Cooley-Tukey) radix-6
// decimation-in-time, using FFTW3's own hand-tuned n1_6 (no-twiddle) and
// t1_6 (twiddle) codelets -- confirmed via fftw_print_plan's dump:
//   (dft-rader-37 (dft-ct-dit/6 (dftw-direct-6/10 "t1_6") (dft-direct-6-x6 "n1_6")) ...)
// appearing three times (cld1, cld2, cld_omega in dft/rader.c's mkP -- all
// three are structurally identical size-36 DFTs, just wired to different
// buffers at the call site).
//
// FMA note: unlike fdlibm.js, NO software-FMA emulation is needed here.
// n1_6.c/t1_6.c's FMA/FNMS macros (kernel/ifftw.h) are plain C expressions
// (`(a*b)+c`), and this project's vendored FFTW3 build has `-ffp-contract=off`
// on every generated Makefile (see src/makeflags's header comment) --
// confirmed via `-S` disassembly of these exact codelets under this
// project's actual build flags: zero fmadd/fmsub/fnmadd/fnmsub instructions.
// So plain, unfused JS arithmetic is the byte-exact match here.
//
// Sources: third_party/fftw-3.3.8/{dft/rader.c, kernel/rader.c, kernel/
// primes.c, kernel/twiddle.c, dft/ct.c, dft/dftw-direct.c,
// dft/scalar/codelets/{n1_6,t1_6}.c}, FFTW3 (c) 2003, 2007-14 Matteo Frigo
// and MIT, GPLv2+.
// =============================================================================

const { realCexp } = require('./Trig.js');
const { n1_6: radix6NoTwiddle } = require('./Codelets/complex/n1_6.js');
const { t1_6: radix6Twiddle } = require('./Codelets/complex/t1_6.js');

// ---------------------------------------------------------------------------
// kernel/primes.c: MULMOD, power_mod, find_generator -- plain Number
// arithmetic is exact here (all operands stay < 37, nowhere near
// Number.MAX_SAFE_INTEGER), so safe_mulmod's overflow-avoidance path is
// never needed.
// ---------------------------------------------------------------------------
function mulmod(x, y, p) {
  return (x * y) % p;
}

function powerMod(n, m, p) {
  if (m === 0) return 1;
  if (m % 2 === 0) {
    const x = powerMod(n, m / 2, p);
    return mulmod(x, x, p);
  }
  return mulmod(n, powerMod(n, m - 1, p), p);
}

function getPrimeFactors(n) {
  const primef = [];
  primef.push(2);
  while ((n & 1) === 0) n >>= 1;
  if (n === 1) return primef;
  for (let i = 3; i * i <= n; i += 2) {
    if (n % i === 0) {
      primef.push(i);
      while (n % i === 0) n /= i;
    }
  }
  if (n === 1) return primef;
  primef.push(n);
  return primef;
}

function findGenerator(p) {
  if (p === 2) return 1;
  const pm1 = p - 1;
  const primef = getPrimeFactors(pm1);
  let n = 2;
  for (let i = 0; i < primef.length; i++) {
    if (powerMod(n, pm1 / primef[i], p) === 1) {
      i = -1;
      n++;
    }
  }
  return n;
}

// ---------------------------------------------------------------------------
// Size-36 complex DFT (FFT_SIGN=-1), via dft/ct.c's apply_dit composed with
// the two codelets above: r=6, m=6 decimation-in-time.
//
// Stage 1 (n1_6, "dft-direct-6-x6"): decimate input by phase p=0..5 (element
// p+k*6 for k=0..5), radix-6 DFT each phase, store phase p's 6 outputs at
// physical positions [p*6 .. p*6+5] (block-major).
// Stage 2 (t1_6, "dftw-direct-6/10"): for each column m=0..5, gather the 6
// block values at position m (physical [p*6+m] for p=0..5), twiddle-multiply
// by exp(-2*pi*i*m*p/36) for p=1..5, radix-6-combine, write back in place.
// By construction physical position [p*6+m] after stage 2 holds X[m+p*6] --
// which equals the physical index itself, so no final permutation is needed
// (see dft/ct.c's apply_dit / this file's header comment for the derivation).
// ---------------------------------------------------------------------------
function dftSize36(reIn, imIn) {
  const workRe = new Float64Array(36), workIm = new Float64Array(36);

  const riP = new Float64Array(6), iiP = new Float64Array(6);
  for (let p = 0; p < 6; p++) {
    for (let k = 0; k < 6; k++) {
      riP[k] = reIn[p + k * 6];
      iiP[k] = imIn[p + k * 6];
    }
    const [roP, ioP] = radix6NoTwiddle(riP, iiP);
    for (let k = 0; k < 6; k++) {
      workRe[p * 6 + k] = roP[k];
      workIm[p * 6 + k] = ioP[k];
    }
  }

  const br = new Float64Array(6), bi = new Float64Array(6);
  const Wc = new Float64Array(6), Ws = new Float64Array(6);
  for (let m = 0; m < 6; m++) {
    for (let p = 0; p < 6; p++) {
      br[p] = workRe[p * 6 + m];
      bi[p] = workIm[p * 6 + m];
    }
    for (let p = 1; p < 6; p++) {
      const [c, s] = realCexp((m * p) % 36, 36);
      Wc[p] = c; Ws[p] = s;
    }
    const [outR, outI] = radix6Twiddle(br, bi, Wc, Ws);
    for (let p = 0; p < 6; p++) {
      workRe[p * 6 + m] = outR[p];
      workIm[p * 6 + m] = outI[p];
    }
  }

  return [workRe, workIm];
}

// ---------------------------------------------------------------------------
// dft/rader.c: mkomega -- precompute the (N-1)-length "omega" array (the
// DFT of the Rader convolution kernel b[i] = exp(-2*pi*i*g^-i/N)/(N-1)).
// Note the twiddle factors here use N=37 (the OUTER prime size, mkomega's
// own `n` parameter), NOT N-1=36 -- only the final "FFT the omega array"
// step uses the size-36 sub-transform.
// ---------------------------------------------------------------------------
function computeOmega(n, ginv) {
  const scale = n - 1;
  const omRe = new Float64Array(n - 1), omIm = new Float64Array(n - 1);
  let gpower = 1;
  for (let i = 0; i < n - 1; i++) {
    const [c, s] = realCexp(gpower, n);
    omRe[i] = c / scale;
    omIm[i] = (-1 * s) / scale; // FFT_SIGN * w[1] / scale, FFT_SIGN=-1
    gpower = mulmod(gpower, ginv, n);
  }
  // gpower === 1 here (g^-(n-1) == g^0 mod n)

  const [outRe, outIm] = dftSize36(omRe, omIm);
  const omega = new Float64Array((n - 1) * 2);
  for (let i = 0; i < n - 1; i++) {
    omega[2 * i] = outRe[i];
    omega[2 * i + 1] = outIm[i];
  }
  return omega;
}

const RADER_N = 37;
const RADER_G = findGenerator(RADER_N);
const RADER_GINV = powerMod(RADER_G, RADER_N - 2, RADER_N);
let _omegaCache = null;
function getOmega() {
  if (!_omegaCache) _omegaCache = computeOmega(RADER_N, RADER_GINV);
  return _omegaCache;
}

// ---------------------------------------------------------------------------
// dft/rader.c: apply() -- full Rader's-algorithm complex DFT of length 37.
// Same call signature as GenericSolver1D.js's dftGeneric, so it's a drop-in
// replacement at the Rank2Orchestration.js dispatch site.
// ---------------------------------------------------------------------------
function dftRader37(n, ri, riOff, ii, iiOff, is, ro, roOff, io, ioOff, os) {
  if (n !== RADER_N) throw new Error(`dftRader37: only n=37 is supported, got n=${n}`);

  const r = RADER_N, g = RADER_G, ginv = RADER_GINV;
  const omega = getOmega();

  const r0 = ri[riOff], i0 = ii[iiOff];

  const bufRe = new Float64Array(r - 1), bufIm = new Float64Array(r - 1);
  let gpower = 1;
  for (let k = 0; k < r - 1; k++) {
    bufRe[k] = ri[riOff + gpower * is];
    bufIm[k] = ii[iiOff + gpower * is];
    gpower = mulmod(gpower, g, r);
  }
  // gpower === 1 here (g^(r-1) == 1 mod r)

  const [dftRe, dftIm] = dftSize36(bufRe, bufIm);
  for (let k = 0; k < r - 1; k++) {
    ro[roOff + (k + 1) * os] = dftRe[k];
    io[ioOff + (k + 1) * os] = dftIm[k];
  }

  ro[roOff] = r0 + ro[roOff + os];
  io[ioOff] = i0 + io[ioOff + os];

  for (let k = 0; k < r - 1; k++) {
    const rW = omega[2 * k], iW = omega[2 * k + 1];
    const rB = ro[roOff + (k + 1) * os], iB = io[ioOff + (k + 1) * os];
    ro[roOff + (k + 1) * os] = rW * rB - iW * iB;
    io[ioOff + (k + 1) * os] = -(rW * iB + iW * rB);
  }

  ro[roOff + os] += r0;
  io[ioOff + os] -= i0;

  const inRe = new Float64Array(r - 1), inIm = new Float64Array(r - 1);
  for (let k = 0; k < r - 1; k++) {
    inRe[k] = ro[roOff + (k + 1) * os];
    inIm[k] = io[ioOff + (k + 1) * os];
  }
  const [dft2Re, dft2Im] = dftSize36(inRe, inIm);

  gpower = 1;
  for (let k = 0; k < r - 1; k++) {
    ro[roOff + gpower * os] = dft2Re[k];
    io[ioOff + gpower * os] = -dft2Im[k];
    gpower = mulmod(gpower, ginv, r);
  }
  // gpower === 1 here
}

// =============================================================================
// General Rader's-algorithm path -- same algorithm as dftRader37 above, but
// parameterized over n instead of hardcoded to 37, using Composite1D.js's
// dftComposite(n-1, ...) for the two "FFT the size-(n-1) buffer" steps
// instead of the hand-inlined size-36 special case. Kept as a SEPARATE
// function from dftRader37 (not a refactor of it) so the original,
// long-verified n=37 path stays byte-for-byte untouched -- this is purely
// additive, matching the "existing files are a stable back-compat facade"
// principle (see FFTW3JS/README.md).
//
// Requiring Composite1D.js here (rather than at module top) avoids a
// load-order issue: Composite1D.js's own require of Codelets/complex/index
// happens before RaderSolver.js finishes defining its exports, but since
// dftRader() is only ever CALLED (not required) after both modules have
// fully loaded, a top-level require is actually fine too -- done at the top
// for clarity, no lazy-require needed, now that n1_6/t1_6 moved out (see
// this file's earlier comment block).
// =============================================================================
const { dftComposite, isFullyPortedComplex, isSingleLevelOrDirect } = require('./Composite1D.js');

// True iff dftRader(n, ...) is safe to trust as FFTW3's exact bit path:
// applicability matches dft/rader.c's own gate (n prime, n>RADER_MAX_SLOW,
// n-1 factors into {2,3,5} -- see Planner/chooseDecomposition.js's identical
// condition) AND the n-1 sub-transform itself is fully ported, where "fully
// ported" means isFullyPortedComplexEmbedded(n-1, undefined) -- the
// genuinely in-place-aware predictor (Composite1D.js), which models BOTH
// q1 (dft/dftw-directsq.c) AND dft-buffered (dft/buffered.c) as competing
// nested-cofactor solvers, cost-compared against each other exactly like
// real FFTW's own LIFO search (see Planner/chooseDecomposition.js's
// bufferedCost header for the buffered derivation).
//
// HISTORY (this session, in order -- kept because each step's lesson still
// applies to future changes here):
//  1. An EARLIER attempt to relax the original isSingleLevelOrDirect-based
//     gate was based on a broken test harness (a parsing bug silently
//     treated every comparison as a false-positive "bit-exact" match --
//     NaN vs NaN compares false, so a max-tracking loop that only updates
//     on `d > m` never updates at all). Once fixed, n=113 (bluestein) came
//     back genuinely wrong (~1.2e-14), confirming the restriction was
//     needed.
//  2. A first attempt at a genuinely in-place-aware predictor
//     (isFullyPortedComplexEmbedded modeling ONLY q1, no buffered) got the
//     q1 math itself right (independently verified bit-exact for radixes
//     2, 3, 5, 8 -- n=24/27/125/128) but was UNSAFE as a GENERAL
//     nested-cofactor predictor: real FFTW's in-place planner very often
//     picks dft-buffered instead of either q1 or the naive out-of-place
//     structure, and a cost model that doesn't know buffered exists will
//     confidently predict the wrong structure. A broad safety sweep found
//     ~1300+ mismatches once this started trusting q1-only structures for
//     arbitrary nested cofactors -- reverted to isSingleLevelOrDirect at
//     that point.
//  3. dft-buffered was then properly modeled (bufferedCost in
//     chooseDecomposition.js, ported from dft/buffered.c + kernel/
//     buffered.c + rdft/dft-r2hc.c + rdft/rank0.c) and wired in as a
//     genuine cost-compared competitor to q1, not an unconditional
//     fallback -- this file's isFullyPortedRader was updated to use
//     isFullyPortedComplexEmbedded again, THIS time with buffered in the
//     model. See Composite1D.js's isFullyPortedComplexEmbedded header for
//     the current state.
//
// HARD RULE (established and reconfirmed multiple times this session): any
// change here must be confirmed via a harness proven to detect a KNOWN-BAD
// case before trusting a broad "all clear" sweep result -- and the sweep
// must cover the FULL range of newly-trusted sizes, not just a handful of
// hand-picked spot checks, no matter how solid the underlying math seemed
// in isolation.
function isFullyPortedRader(n) {
  const RADER_MAX_SLOW = 32;
  if (!(n > RADER_MAX_SLOW)) return false;
  let m = n - 1;
  for (const p of [2, 3, 5]) { while (m % p === 0) m /= p; }
  if (m !== 1) return false;
  return isFullyPortedComplex(n - 1) && isSingleLevelOrDirect(n - 1);
}

const _generatorCache = new Map(); // n -> {g, ginv}
function getGeneratorFor(n) {
  let e = _generatorCache.get(n);
  if (!e) {
    const g = findGenerator(n);
    const ginv = powerMod(g, n - 2, n);
    e = { g, ginv };
    _generatorCache.set(n, e);
  }
  return e;
}

// General version of computeOmega (see that function's header) -- uses
// dftComposite(n-1, ...) instead of the hardcoded dftSize36.
function computeOmegaGeneral(n, ginv) {
  const scale = n - 1;
  const omRe = new Float64Array(n - 1), omIm = new Float64Array(n - 1);
  let gpower = 1;
  for (let i = 0; i < n - 1; i++) {
    const [c, s] = realCexp(gpower, n);
    omRe[i] = c / scale;
    omIm[i] = (-1 * s) / scale;
    gpower = mulmod(gpower, ginv, n);
  }

  const outRe = new Float64Array(n - 1), outIm = new Float64Array(n - 1);
  dftComposite(n - 1, omRe, 0, omIm, 0, 1, outRe, 0, outIm, 0, 1);
  const omega = new Float64Array((n - 1) * 2);
  for (let i = 0; i < n - 1; i++) {
    omega[2 * i] = outRe[i];
    omega[2 * i + 1] = outIm[i];
  }
  return omega;
}

const _omegaCacheGeneral = new Map(); // n -> Float64Array
function getOmegaGeneral(n, ginv) {
  let om = _omegaCacheGeneral.get(n);
  if (!om) {
    om = computeOmegaGeneral(n, ginv);
    _omegaCacheGeneral.set(n, om);
  }
  return om;
}

// dft/rader.c: apply() -- general Rader's-algorithm complex DFT, any prime
// n. Same call signature/structure as dftRader37, parameterized by n.
function dftRader(n, ri, riOff, ii, iiOff, is, ro, roOff, io, ioOff, os) {
  const { g, ginv } = getGeneratorFor(n);
  const omega = getOmegaGeneral(n, ginv);

  const r0 = ri[riOff], i0 = ii[iiOff];

  const bufRe = new Float64Array(n - 1), bufIm = new Float64Array(n - 1);
  let gpower = 1;
  for (let k = 0; k < n - 1; k++) {
    bufRe[k] = ri[riOff + gpower * is];
    bufIm[k] = ii[iiOff + gpower * is];
    gpower = mulmod(gpower, g, n);
  }

  const dftRe = new Float64Array(n - 1), dftIm = new Float64Array(n - 1);
  dftComposite(n - 1, bufRe, 0, bufIm, 0, 1, dftRe, 0, dftIm, 0, 1);
  for (let k = 0; k < n - 1; k++) {
    ro[roOff + (k + 1) * os] = dftRe[k];
    io[ioOff + (k + 1) * os] = dftIm[k];
  }

  ro[roOff] = r0 + ro[roOff + os];
  io[ioOff] = i0 + io[ioOff + os];

  for (let k = 0; k < n - 1; k++) {
    const rW = omega[2 * k], iW = omega[2 * k + 1];
    const rB = ro[roOff + (k + 1) * os], iB = io[ioOff + (k + 1) * os];
    ro[roOff + (k + 1) * os] = rW * rB - iW * iB;
    io[ioOff + (k + 1) * os] = -(rW * iB + iW * rB);
  }

  ro[roOff + os] += r0;
  io[ioOff + os] -= i0;

  const inRe = new Float64Array(n - 1), inIm = new Float64Array(n - 1);
  for (let k = 0; k < n - 1; k++) {
    inRe[k] = ro[roOff + (k + 1) * os];
    inIm[k] = io[ioOff + (k + 1) * os];
  }
  const dft2Re = new Float64Array(n - 1), dft2Im = new Float64Array(n - 1);
  dftComposite(n - 1, inRe, 0, inIm, 0, 1, dft2Re, 0, dft2Im, 0, 1);

  gpower = 1;
  for (let k = 0; k < n - 1; k++) {
    ro[roOff + gpower * os] = dft2Re[k];
    io[ioOff + gpower * os] = -dft2Im[k];
    gpower = mulmod(gpower, ginv, n);
  }
}

module.exports = {
  dftRader37,
  mulmod, powerMod, findGenerator,
  radix6NoTwiddle, radix6Twiddle, dftSize36,
  computeOmega,
  RADER_N, RADER_G, RADER_GINV,
  dftRader, isFullyPortedRader,
};

// ---------------------------------------------------------------------------
// Self-test (node RaderSolver.js) -- smoke test only, NOT a substitute for
// the real cross-language verification (see verify/run_phase1_tests.js,
// which exercises this via Rank2Orchestration.js's 37x37/37x43 fixtures).
// Checks: found generator/ginv match the hand-derived values (g=2, ginv=19
// for N=37), and a DC-term sanity check (X[0] must equal the plain sum).
// ---------------------------------------------------------------------------
if (require.main === module) {
  console.log('=== RaderSolver self-test (n=37) ===');
  console.log(`g=${RADER_G} (expect 2), ginv=${RADER_GINV} (expect 19)`);
  console.log(`mulmod(g,ginv,37)=${mulmod(RADER_G, RADER_GINV, 37)} (expect 1)`);

  const n = 37;
  const ri = new Float64Array(n), ii = new Float64Array(n);
  for (let i = 0; i < n; i++) ri[i] = i + 1;
  const ro = new Float64Array(n), io = new Float64Array(n);
  dftRader37(n, ri, 0, ii, 0, 1, ro, 0, io, 0, 1);
  const sum = Array.from(ri).reduce((a, b) => a + b, 0);
  console.log(`DC (ro[0]): ${ro[0]} (expect sum = ${sum})`);
}
