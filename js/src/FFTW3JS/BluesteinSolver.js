'use strict';

// =============================================================================
// BluesteinSolver.js
// Faithful JS port of dft/bluestein.c's chirp-z (Bluestein's algorithm)
// solver, generalized over n from the start (unlike RaderSolver.js's
// original N=37-only path before it was generalized) -- this file only ever
// had a general form. Applicable to any prime n>24 (BLUESTEIN_MAX_SLOW),
// used by real FFTW3 whenever Rader's algorithm doesn't apply (n-1 doesn't
// factor into {2,3,5}) -- see Planner/chooseDecomposition.js's identical
// applicability condition.
//
// ALGORITHM (transcribed literally from dft/bluestein.c's apply(), variable
// names kept close to the source so it's directly diffable against it):
//   1. nb = smallest 5-smooth number >= 2n-1 (choose_transform_size).
//   2. w[k] = exp(2*pi*i*k^2/(2n)) for k=0..n-1 (bluestein_sequence),
//      computed via an incremental k^2 update (ksq += 2k+1, wrapped mod 2n)
//      to avoid overflow, matching the C source exactly.
//   3. W = DFT(conjugate-symmetric zero-padded w)/nb, precomputed once
//      (mktwiddle) -- W[i]=W[nb-i]=w[i]/nb for i=1..n-1, W[0]=w[0]/nb.
//   4. apply: b[i] = x[i]*conj(w[i]) for i<n, zero-padded to nb; b=DFT(b);
//      b = (pointwise combine with W, exact formula transcribed as-is, see
//      below -- this step plus the next DFT together realize the "multiply
//      by conj(w), FFT, pointwise-multiply by W, inverse-FFT-via-forward-FFT
//      trick" chirp-z convolution); b=DFT(b) again; output[i] = combine
//      b[i] (read with re/im SWAPPED, matching the source's "IFFT via FFT
//      with re/im input/output swapped" comment) with w[i].
//
// The two inner "DFT of size nb" calls reuse the SAME dftInterleavedInPlace
// helper (itself built on Composite1D.dftComposite) -- matching dft/
// bluestein.c's own reuse of one cldf plan for both calls.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const { realCexp } = require('./Trig.js');
const { dftComposite, isFullyPortedComplex, isSingleLevelOrDirect } = require('./Composite1D.js');

const BLUESTEIN_MAX_SLOW = 24;

function isPrime(n) {
  if (n <= 1) return false;
  if (n % 2 === 0) return n === 2;
  for (let i = 3; i * i <= n; i += 2) if (n % i === 0) return false;
  return true;
}

function factorsIntoSmallPrimes(n) {
  for (const p of [2, 3, 5]) { while (n % p === 0) n /= p; }
  return n === 1;
}

function chooseTransformSize(minsz) {
  let n = minsz;
  while (!factorsIntoSmallPrimes(n)) n++;
  return n;
}

// True iff dftBluestein(n, ...) is safe to trust as FFTW3's exact bit path:
// applicability matches dft/bluestein.c's own gate (n prime, n>16, n>
// BLUESTEIN_MAX_SLOW -- see Planner/chooseDecomposition.js's identical
// condition) AND the size-nb sub-transform itself is fully ported, where
// "fully ported" means isFullyPortedComplexEmbedded(nb, undefined) -- the
// genuinely in-place-aware predictor (Composite1D.js), which models BOTH
// q1 (dft/dftw-directsq.c) AND dft-buffered (dft/buffered.c) as competing
// nested-cofactor solvers, cost-compared exactly like real FFTW's own LIFO
// search (see Planner/chooseDecomposition.js's bufferedCost header).
//
// HISTORY (this session, in order -- kept because each step's lesson still
// applies to future changes here): found the hard way, N=113 (nb=225, a
// two-level nested CT decomposition) verified WRONG (~1.2e-14 off real
// FFTW) despite the ORIGINAL isFullyPortedComplex(225) reporting true,
// root-caused via verify/stride_test.c and verify/bluestein_ref.c --
// replaced with isFullyPortedComplex(nb) && isSingleLevelOrDirect(nb). A
// first attempt at a genuinely in-place-aware replacement (q1-only, no
// buffered) got the q1 math right (independently verified bit-exact for
// radixes 2, 3, 5, 8 -- n=24/27/125/128) but was UNSAFE as a GENERAL
// predictor: real FFTW's in-place planner very often picks dft-buffered
// instead, and a cost model that doesn't know it exists confidently
// predicts the wrong structure -- a broad safety sweep found ~1300+
// mismatches (e.g. n=131, nb=270: real FFTW uses dft-buffered-27-x10(t2_10,
// vrank-x10(t1_3, n1_9)); the q1-only predictor wrongly claimed
// dft-ct-dit/3(t1_3, dft-ct-dif/3(q1_3, ...)) instead). dft-buffered was
// then properly modeled (bufferedCost in chooseDecomposition.js, ported
// from dft/buffered.c + kernel/buffered.c + rdft/dft-r2hc.c +
// rdft/rank0.c) and wired in as a genuine cost-compared competitor to q1,
// restoring this gate to isFullyPortedComplexEmbedded.
//
// HARD RULE (established and reconfirmed multiple times this session): any
// change here must be confirmed via a harness proven to detect a KNOWN-BAD
// case (e.g. n=113) before trusting a broad "all clear" sweep result -- and
// the sweep must cover the FULL range of newly-trusted sizes, not just a
// handful of hand-picked spot checks, no matter how solid the underlying
// math seemed in isolation.
function isFullyPortedBluestein(n) {
  if (!isPrime(n) || !(n > 16) || !(n > BLUESTEIN_MAX_SLOW)) return false;
  const nb = chooseTransformSize(2 * n - 1);
  return isFullyPortedComplex(nb) && isSingleLevelOrDirect(nb);
}

// DFT a size-n complex array IN PLACE, interleaved (buf[2i]=Re, buf[2i+1]=Im).
function dftInterleavedInPlace(n, buf) {
  const ro = new Float64Array(n), io = new Float64Array(n);
  dftComposite(n, buf, 0, buf, 1, 2, ro, 0, io, 0, 1);
  for (let i = 0; i < n; i++) { buf[2 * i] = ro[i]; buf[2 * i + 1] = io[i]; }
}

const _chirpCache = new Map(); // n -> { w, W, nb }
function getChirp(n) {
  let e = _chirpCache.get(n);
  if (e) return e;

  const nb = chooseTransformSize(2 * n - 1);
  const n2 = 2 * n;

  // bluestein_sequence: w[k] = exp(2*pi*i*k^2/(2n)), incremental ksq.
  const w = new Float64Array(2 * n);
  let ksq = 0;
  for (let k = 0; k < n; k++) {
    const [c, s] = realCexp(ksq, n2);
    w[2 * k] = c; w[2 * k + 1] = s;
    ksq += 2 * k + 1;
    while (ksq > n2) ksq -= n2;
  }

  // mktwiddle: W = zero-padded conjugate-symmetric copy of w/nb, then DFT.
  const W = new Float64Array(2 * nb);
  W[0] = w[0] / nb; W[1] = w[1] / nb;
  for (let i = 1; i < n; i++) {
    W[2 * i] = W[2 * (nb - i)] = w[2 * i] / nb;
    W[2 * i + 1] = W[2 * (nb - i) + 1] = w[2 * i + 1] / nb;
  }
  dftInterleavedInPlace(nb, W);

  e = { w, W, nb };
  _chirpCache.set(n, e);
  return e;
}

// dft/bluestein.c: apply() -- Bluestein's-algorithm complex DFT of any
// prime n. Same call signature as GenericSolver1D.js's dftGeneric.
function dftBluestein(n, ri, riOff, ii, iiOff, is, ro, roOff, io, ioOff, os) {
  const { w, W, nb } = getChirp(n);

  const b = new Float64Array(2 * nb); // zero-padded by construction
  for (let i = 0; i < n; i++) {
    const xr = ri[riOff + i * is], xi = ii[iiOff + i * is];
    const wr = w[2 * i], wi = w[2 * i + 1];
    b[2 * i] = xr * wr + xi * wi;
    b[2 * i + 1] = xi * wr - xr * wi;
  }

  dftInterleavedInPlace(nb, b);

  for (let i = 0; i < nb; i++) {
    const xr = b[2 * i], xi = b[2 * i + 1];
    const wr = W[2 * i], wi = W[2 * i + 1];
    b[2 * i] = xi * wr + xr * wi;
    b[2 * i + 1] = xr * wr - xi * wi;
  }

  dftInterleavedInPlace(nb, b);

  for (let i = 0; i < n; i++) {
    const xi = b[2 * i], xr = b[2 * i + 1]; // swapped read, per source comment
    const wr = w[2 * i], wi = w[2 * i + 1];
    ro[roOff + i * os] = xr * wr + xi * wi;
    io[ioOff + i * os] = xi * wr - xr * wi;
  }
}

module.exports = {
  dftBluestein, isFullyPortedBluestein, chooseTransformSize, factorsIntoSmallPrimes,
};

// ---------------------------------------------------------------------------
// Self-test (node BluesteinSolver.js) -- smoke test only, NOT a substitute
// for real cross-language verification.
// ---------------------------------------------------------------------------
if (require.main === module) {
  function handDFT(re, im) {
    const n = re.length;
    const outRe = new Array(n).fill(0), outIm = new Array(n).fill(0);
    for (let k = 0; k < n; k++) {
      for (let j = 0; j < n; j++) {
        const ang = -2 * Math.PI * k * j / n;
        outRe[k] += re[j] * Math.cos(ang) - im[j] * Math.sin(ang);
        outIm[k] += re[j] * Math.sin(ang) + im[j] * Math.cos(ang);
      }
    }
    return [outRe, outIm];
  }
  for (const n of [17, 19, 23, 29, 31, 43, 47, 53, 103]) {
    const ri = Array.from({ length: n }, (_, i) => Math.sin(i * 1.7 + 0.3));
    const ii = Array.from({ length: n }, (_, i) => Math.cos(i * 0.9 + 1.1));
    const ro = new Float64Array(n), io = new Float64Array(n);
    try {
      dftBluestein(n, new Float64Array(ri), 0, new Float64Array(ii), 0, 1, ro, 0, io, 0, 1);
      const [hre, him] = handDFT(ri, ii);
      let maxDiff = 0;
      for (let i = 0; i < n; i++) maxDiff = Math.max(maxDiff, Math.abs(ro[i] - hre[i]), Math.abs(io[i] - him[i]));
      console.log(`n=${n}: maxDiff vs hand-DFT=${maxDiff} fullyPorted=${isFullyPortedBluestein(n)}`);
    } catch (e) {
      console.log(`n=${n}: THROWN (${e.message}) fullyPorted=${isFullyPortedBluestein(n)}`);
    }
  }
}
