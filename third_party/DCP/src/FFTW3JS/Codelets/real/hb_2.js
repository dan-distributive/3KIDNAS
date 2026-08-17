'use strict';

// =============================================================================
// hb_2.js -- faithful JS port of rdft/scalar/r2cb/hb_2.c (non-FMA branch),
// FFTW3's radix-2 Cooley-Tukey twiddle/combine codelet for the HC2R
// (halfcomplex-input, DIF) recursive real engine -- the inverse-direction
// counterpart to hf_2.js. Same (cr,ci,Wc,Ws) -> (outCr,outCi) physical
// position convention as hf_2.js; transcribed literally, not derived by
// analogy (see hf_2.js's header for why).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

function hb_2(cr, ci, Wc, Ws) {
  const T1 = cr[0], T2 = ci[0];
  const T6 = T1 - T2;
  const T3 = ci[1], T4 = cr[1];
  const T8 = T3 + T4;

  const outCr = new Float64Array(2), outCi = new Float64Array(2);
  outCr[0] = T1 + T2;
  outCi[0] = T3 - T4;
  const T5 = Wc[1], T7 = Ws[1];
  outCr[1] = T5 * T6 - T7 * T8;
  outCi[1] = T7 * T6 + T5 * T8;

  return [outCr, outCi];
}

module.exports = { hb_2 };
