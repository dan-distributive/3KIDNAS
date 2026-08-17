'use strict';

// =============================================================================
// Trig.js
// Faithful JS port of FFTW3's real_cexp (third_party/fftw-3.3.8/kernel/trig.c),
// part of the hand-port of FFTW3's actual r2c/c2r algorithm to JS (see
// third_party/DCP/src/FFTW3JS/ and the project plan for why: Fortran calls
// real FFTW3 directly for beam convolution; the JS port used a different,
// independent FFT library (ndarray-fft), and that difference -- not any
// remaining RNG/rounding-order bug -- was proven to be the root cause of the
// whole Fortran/JS optimizer basin divergence chased earlier this project).
//
// SCOPE
// -----
// This build of FFTW3 (see third_party/fftw-3.3.8/config.h: FFTW_LDOUBLE and
// FFTW_QUAD are both undefined) uses `trigreal == double` throughout, so
// there is exactly one twiddle-factor code path to port: real_cexp's
// octant-based argument reduction followed by a call to cos/sin. FFTW3's
// OWN `real_cexp` already calls this project's own vendored fdlibm
// (third_party/fftw-3.3.8/kernel/trig.c was patched to call fdlibm_cos/
// fdlibm_sin instead of the platform libm, for exactly the same
// cross-platform-determinism reason fdlibm.js exists) -- so this JS port
// calls the SAME fdlibm.js already used elsewhere in this codebase, giving
// both sides an identical trig algorithm, not just "close enough" native
// Math.cos/sin.
//
// real_cexp computes cos(2*pi*m/n) and sin(2*pi*m/n) for integer m, n, but
// NOT by naively calling cos/sin on the raw angle -- it first folds m/n
// into the first octant (0 <= theta <= pi/4) via three conditional swaps
// (each exploiting a symmetry of sin/cos), improving accuracy, then
// un-folds the result the same way. This exact reduction (not just "some
// argument reduction") needs replicating bit-for-bit, since the whole
// point of this port is matching FFTW3's specific arithmetic, not just
// producing *a* correct twiddle factor.
//
// Source: third_party/fftw-3.3.8/kernel/trig.c (real_cexp, by2pi, K2PI),
// FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+. See that file's
// header for the full license notice -- a derivative port of it carries
// the same obligations.
// =============================================================================

const { fdCos, fdSin } = require('../StandardMath/fdlibm.js');

// FFTW3: static const trigreal K2PI = KTRIG(6.2831853071795864769252867665590057683943388);
// Written as the identical decimal literal (not 2*Math.PI) so both languages'
// compile-time/parse-time rounding to the nearest double are provably the
// same operation on the same digit string, not just "should be equal".
const K2PI = 6.2831853071795864769252867665590057683943388;

// FFTW3: #define by2pi(m, n) ((K2PI * (m)) / (n))
function by2pi(m, n) {
  return (K2PI * m) / n;
}

// ---------------------------------------------------------------------------
// realCexp -- port of real_cexp(INT m, INT n, trigreal *out)
//
// Returns [c, s] = [cos(2*pi*m/n), sin(2*pi*m/n)].
//
// Line-for-line transcription of the C control flow, including the exact
// order of the three octant-reduction conditionals and the exact order of
// the three un-fold conditionals afterward -- this order is not arbitrary
// (each step assumes the previous ones already ran), so it is preserved
// exactly rather than reordered for "clarity".
// ---------------------------------------------------------------------------
function realCexp(m, n) {
  let octant = 0;
  const quarter_n = n;

  n = n * 4;
  m = m * 4;

  if (m < 0) m = m + n;
  if (m > n - m) { m = n - m; octant |= 4; }
  if (m - quarter_n > 0) { m = m - quarter_n; octant |= 2; }
  if (m > quarter_n - m) { m = quarter_n - m; octant |= 1; }

  const theta = by2pi(m, n);
  let c = fdCos(theta);
  let s = fdSin(theta);

  if (octant & 1) { const t = c; c = s; s = t; }
  if (octant & 2) { const t = c; c = -s; s = t; }
  if (octant & 4) { s = -s; }

  return [c, s];
}

module.exports = { realCexp, K2PI, by2pi };

// ---------------------------------------------------------------------------
// Self-test (node Trig.js) -- smoke test only, NOT a substitute for real
// cross-language verification against the compiled ground-truth harness
// (see verify/ground_truth_harness.c). Spot-checks realCexp against
// Math.cos/Math.sin(2*pi*m/n) computed the "obvious" (unreduced) way, for a
// spread of m/n including values that exercise all three octant folds.
// ---------------------------------------------------------------------------
if (require.main === module) {
  const cases = [
    [0, 47], [1, 47], [5, 47], [11, 47], [23, 47],
    [1, 37], [18, 37], [36, 37],
    [1, 4], [1, 8], [3, 8],
    [-1, 47], [-23, 47],
  ];
  console.log('=== realCexp vs naive Math.cos/sin(2*pi*m/n) ===');
  for (const [m, n] of cases) {
    const [c, s] = realCexp(m, n);
    const theta = (2 * Math.PI * m) / n;
    const c2 = Math.cos(theta);
    const s2 = Math.sin(theta);
    console.log(
      `m=${m} n=${n}  c=${c} (naive ${c2}, diff ${c - c2})  `
      + `s=${s} (naive ${s2}, diff ${s - s2})`
    );
  }
}
