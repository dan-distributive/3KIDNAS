'use strict';

// =============================================================================
// DhtRaderSolver.js
// Faithful JS port of FFTW3's REAL-side prime-N solver chain:
//   rdft/rdft-dht.c  (R2HC/HC2R <-> DHT adapter, "r2hc-dht"/"hc2r-dht")
//   rdft/dht-rader.c (Rader's algorithm for a prime-size DHT, "dht-rader")
//
// WHY THIS EXISTS: real FFTW never applies the COMPLEX-domain Rader/
// Bluestein (RaderSolver.js/BluesteinSolver.js) to a real-to-halfcomplex
// (R2HC) or halfcomplex-to-real (HC2R) problem directly -- confirmed via
// fftw_fprint_plan for many prime N (193, 257, 1009, ...): the ACTUAL
// structure is always "(r2hc-dht (dht-rader-N/(N-1) ...))", a genuinely
// different algorithm family that stays in the REAL domain throughout.
// Before this file, RealEngine1D.js's ONLY candidate for a prime N (or a
// prime FACTOR inside a composite decomposition) beyond GENERIC_SAFE_MAX
// was the O(n^2) generic fallback -- this is what actually unlocks large
// primes at O(n log n) cost, matching real FFTW's own choice.
//
// THE KEY SIMPLIFICATION THAT MAKES THIS TRACTABLE: dht-rader.c's Rader
// convolution is computed via TWO R2HC transforms of size npad=N-1 (NOT
// complex DFTs -- see R2HC_ONLY_CONV in the source), both at UNIT STRIDE
// (mktensor_1d(npad,1,1)) -- i.e. the EXACT SAME chooseReal('R2HC',npad)/
// executePlanR2HC(...) machinery already built and extensively verified
// elsewhere in this file, with NO stride-mismatch risk (contrast
// RaderSolver.js/BluesteinSolver.js, whose complex sub-plans use a
// NON-unit stride and therefore need the separate isSingleLevelOrDirect
// caution -- not needed here).
//
// FFTW's planner detail worth recording (kernel/planner.c's search(),
// relax_tab): NO_SLOW is only a soft preference -- if no plan is found
// with it active, the search RETRIES with NO_SLOW (then NO_UGLY) relaxed.
// This is why dht-rader-1009/1008 appears in ground truth even though
// 1008=2^4*3^2*7 does NOT satisfy dht-rader.c's own "n-1 factors into
// {2,3,5}" NO_SLOW-gated preference -- once NO_SLOW is relaxed (nothing
// else applies for a bare prime), dht-rader is still cost-selected over
// generic (O(n log n)-ish vs O(n^2)). This port does not replicate the
// full progressive-relaxation search; instead it exploits the fact that
// for a PRIME n, dht-rader (once its own n-1 sub-transform is fully
// ported) essentially ALWAYS beats generic under real FFTW's actual cost
// model for any n large enough to matter, and confirms this by direct
// sweep against fftw_fprint_plan rather than by re-deriving relax_tab.
//
// Only the UNPADDED variant (npad = n-1, FFTW's "pad=0" solver instance)
// is implemented -- matches every ground-truth case checked so far
// (n=193,257,1009,...); the "pad=1" zero-padded-convolution variant is
// not yet needed and is a natural future extension if a prime is found
// where n-1's own R2HC decomposition can't be trusted but some 5-smooth
// padded size's can.
//
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const { realCexp } = require('./Trig.js');
const { mulmod, powerMod, findGenerator } = require('./RaderSolver.js');
const { chooseReal } = require('./Planner/chooseDecomposition.js');

// Lazy require -- RealEngine1D.js requires THIS file to execute the
// 'dht-rader' plan node, so a top-level require here would cycle. Safe
// because executePlanR2HC is only ever CALLED (not required) after both
// modules have finished their initial module.exports population -- same
// pattern/reasoning as RaderSolver.js's own Composite1D.js require note.
function execR2HC(plan, n, x) {
  return require('./RealEngine1D.js').executePlanR2HC(plan, n, x);
}

const RADER_MAX_SLOW = 32; // kernel/ifftw.h -- shared with RaderSolver.js

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

// rdft/dht-rader.c: mkomega -- precompute the (n-1)-length "omega" array,
// the R2HC transform of the real Hartley-domain convolution kernel
// b[i] = cas(2*pi*g^-i/n)/(n-1), cas(x)=cos(x)+sin(x). Uses ginv (NOT g)
// for its own permutation -- transcribed exactly as the C source orders
// it (a DIFFERENT generator-power walk than apply()'s initial permute).
function computeOmegaDht(n, npad, ginv) {
  const scale = npad;
  const omega = new Float64Array(npad);
  let gpower = 1;
  for (let i = 0; i < n - 1; i++) {
    const [c, s] = realCexp(gpower, n);
    omega[i] = (c + s) / scale;
    gpower = mulmod(gpower, ginv, n);
  }
  const plan = chooseReal('R2HC', npad).plan;
  return execR2HC(plan, npad, omega);
}

const _omegaCache = new Map(); // n -> Float64Array
function getOmega(n, ginv) {
  let om = _omegaCache.get(n);
  if (!om) {
    om = computeOmegaDht(n, n - 1, ginv);
    _omegaCache.set(n, om);
  }
  return om;
}

// rdft/dht-rader.c: apply() (pad=0, R2HC_ONLY_CONV=1 branch) -- Rader's
// algorithm DHT of prime length n. x: real array length n -> real array
// length n (the DHT, NOT packed-halfcomplex -- see applyR2HCViaDht/
// applyHC2RViaDht below for the R2HC/HC2R conversion layer).
function dhtRader(n, x) {
  const { g, ginv } = getGeneratorFor(n);
  const npad = n - 1;
  const omega = getOmega(n, ginv);

  const buf0 = new Float64Array(npad);
  let gpower = 1;
  for (let k = 0; k < n - 1; k++) {
    buf0[k] = x[gpower];
    gpower = mulmod(gpower, g, n);
  }
  // gpower === 1 here (g^(n-1) == 1 mod n)

  const plan = chooseReal('R2HC', npad).plan;
  const buf1 = execR2HC(plan, npad, buf0);

  const O = new Float64Array(n);
  const r0 = x[0];
  O[0] = r0 + buf1[0];

  const buf2 = new Float64Array(npad);
  buf2[0] = buf1[0] * omega[0];
  const kNy = npad / 2;
  for (let k = 1; k < kNy; k++) {
    const rW = omega[k], iW = omega[npad - k];
    const rB = buf1[k], iB = buf1[npad - k];
    const a = rW * rB - iW * iB;
    const b = rW * iB + iW * rB;
    buf2[k] = a + b;
    buf2[npad - k] = a - b;
  }
  buf2[kNy] = buf1[kNy] * omega[kNy]; // Nyquist (npad is always even, n odd prime)

  buf2[0] += r0;

  const buf3 = execR2HC(plan, npad, buf2);

  O[1] = buf3[0];
  gpower = ginv;
  let k = 1;
  for (; k < kNy; k++) {
    O[gpower] = buf3[k] + buf3[npad - k];
    gpower = mulmod(gpower, ginv, n);
  }
  O[gpower] = buf3[k]; // k === kNy here
  k++;
  gpower = mulmod(gpower, ginv, n);
  for (; k < npad; k++) {
    O[gpower] = buf3[npad - k] - buf3[k];
    gpower = mulmod(gpower, ginv, n);
  }
  // gpower === 1 here

  return O;
}

// rdft/rdft-dht.c: apply_r2hc -- R2HC(x) via DHT(x) post-processing
// (FFT_SIGN===-1 branch, matching this project's convention throughout).
function applyR2HCViaDht(n, x) {
  const O = dhtRader(n, x);
  for (let i = 1; i < n - i; i++) {
    const a = 0.5 * O[i];
    const b = 0.5 * O[n - i];
    O[i] = a + b;
    O[n - i] = b - a;
  }
  return O;
}

// rdft/rdft-dht.c: apply_hc2r_save -- HC2R(packedHC) via DHT, WITHOUT
// destroying the input (matches this port's existing genericFallbackHC2R
// convention of never mutating its input array).
function applyHC2RViaDht(n, packedHC) {
  // n is always an odd prime here (n > RADER_MAX_SLOW), so i===n-i never
  // triggers -- the C source's even-n Nyquist branch is dead code in this
  // context and is correctly omitted (matches dht-rader.c's own reasoning,
  // not a simplification).
  const buf = new Float64Array(n);
  buf[0] = packedHC[0];
  for (let i = 1; i < n - i; i++) {
    const a = packedHC[i];
    const b = packedHC[n - i];
    buf[i] = a - b;
    buf[n - i] = a + b;
  }
  return dhtRader(n, buf);
}

// True iff applyR2HCViaDht(n,...)/applyHC2RViaDht(n,...) are safe to trust
// as FFTW3's exact bit path: n prime, n > RADER_MAX_SLOW (matches
// dht-rader.c's own hard gate), and the n-1 R2HC sub-transform this
// solver's cld1/cld2 literally ARE is itself fully ported -- reusing
// isFullyPortedR2HC directly (not re-deriving a separate condition) is
// correct here specifically because, unlike RaderSolver.js/
// BluesteinSolver.js's complex sub-plans, dht-rader's own sub-plans use
// UNIT stride, identical to how isFullyPortedR2HC(n-1) is itself always
// evaluated -- no stride-mismatch risk.
function isFullyPortedDhtRader(n, isFullyPortedR2HC) {
  if (!(n > RADER_MAX_SLOW)) return false;
  return isFullyPortedR2HC(n - 1);
}

module.exports = {
  dhtRader, applyR2HCViaDht, applyHC2RViaDht, isFullyPortedDhtRader,
  computeOmegaDht, getGeneratorFor, RADER_MAX_SLOW,
};

// ---------------------------------------------------------------------------
// Self-test (node DhtRaderSolver.js) -- smoke test only, NOT a substitute
// for real cross-language verification (see verify/ ground-truth sweeps).
// Checks applyR2HCViaDht against a hand-rolled O(n^2) reference R2HC for
// several primes with a nicely-factoring n-1.
// ---------------------------------------------------------------------------
if (require.main === module) {
  function handR2HC(n, x) {
    const O = new Float64Array(n);
    const half = Math.floor(n / 2);
    for (let k = 0; k <= half; k++) {
      let re = 0;
      for (let j = 0; j < n; j++) re += x[j] * Math.cos((2 * Math.PI * k * j) / n);
      O[k] = re;
    }
    for (let k = 1; k < Math.ceil(n / 2); k++) {
      let im = 0;
      for (let j = 0; j < n; j++) im += x[j] * Math.sin((2 * Math.PI * k * j) / n);
      O[n - k] = -im;
    }
    return O;
  }
  for (const n of [37, 41, 43, 61, 73, 97, 109, 127, 151, 181, 193, 211, 241]) {
    const x = Float64Array.from({ length: n }, (_, i) => Math.sin(i * 1.7 + 0.3));
    let ok = true, maxDiff = 0;
    try {
      const O = applyR2HCViaDht(n, x);
      const ref = handR2HC(n, x);
      for (let i = 0; i < n; i++) maxDiff = Math.max(maxDiff, Math.abs(O[i] - ref[i]));
      ok = maxDiff < 1e-8;
    } catch (e) {
      ok = false;
      console.log(`n=${n}: THREW ${e.message}`);
      continue;
    }
    console.log(`n=${n}: maxDiff vs hand-R2HC=${maxDiff} ${ok ? 'OK' : 'FAIL'}`);
  }
}
