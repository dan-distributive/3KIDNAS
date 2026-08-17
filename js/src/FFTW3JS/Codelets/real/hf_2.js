'use strict';

// =============================================================================
// hf_2.js -- faithful JS port of rdft/scalar/r2cf/hf_2.c (non-FMA branch),
// FFTW3's radix-2 Cooley-Tukey twiddle/combine codelet for the R2HC
// (real-input, DIT) recursive real engine.
//
// Structurally distinct from the complex side's t1_r: because the input is
// real, the frequency-domain data has Hermitian symmetry, so this combines
// PACKED (Re,Im) pairs per phase rather than raw complex values, and the
// output role assignment (which slot becomes Re vs Im, and sign) is
// genuinely different from t1_2's combine -- transcribed literally
// variable-for-variable from the C source rather than derived by analogy,
// to avoid introducing a sign error.
//
// cr[p]/ci[p] (p=0..1): phase p's LOCAL (Re,Im) pair at this column (from
// its own size-m sub-transform's packed-halfcomplex block). Wc[1]/Ws[1]:
// twiddle factor for phase 1 (phase 0 is always untwiddled).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

function hf_2(cr, ci, Wc, Ws) {
  const T1 = cr[0], T8 = ci[0];
  const T3 = cr[1], T5 = ci[1];
  const T6 = Wc[1] * T3 + Ws[1] * T5;
  const T7 = Wc[1] * T5 - Ws[1] * T3;

  const outCr = new Float64Array(2), outCi = new Float64Array(2);
  outCi[0] = T1 - T6;
  outCr[0] = T1 + T6;
  outCr[1] = T7 - T8;
  outCi[1] = T7 + T8;
  return [outCr, outCi];
}

module.exports = { hf_2 };
