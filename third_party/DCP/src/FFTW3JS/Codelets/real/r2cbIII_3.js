'use strict';

// =============================================================================
// r2cbIII_3.js -- faithful JS port of rdft/scalar/r2cb/r2cbIII_3.c
// (non-FMA branch). HC2RIII ("shifted") radix-3 direct codelet -- the
// BACKWARD counterpart of r2cfII_3.js. INPUT in[0..2] uses r2cfII_3.js's
// OUTPUT convention (in[0]=Cr[0], in[1]=Cr[1], in[2]=Ci[0]); OUTPUT
// out[0..2] uses r2cfII_3.js's INPUT convention (out[2k]=R0[k],
// out[2k+1]=R1[k]).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP2_000000000 = 2.000000000000000000000000000000000000000000000;
const KP1_732050807 = 1.732050807568877293527446341505872366942805254;

function r2cbIII_3(inArr) {
  const T4 = inArr[2];
  const T5 = KP1_732050807 * T4;
  const T1 = inArr[1];
  const T2 = inArr[0];
  const T3 = T2 - T1;

  const out = new Float64Array(3);
  out[0] = KP2_000000000 * T2 + T1; // FMA
  out[2] = -(T3 + T5);
  out[1] = T3 - T5;
  return out;
}

module.exports = { r2cbIII_3 };
