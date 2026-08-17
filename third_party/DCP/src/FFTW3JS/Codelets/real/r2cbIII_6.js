'use strict';

// =============================================================================
// r2cbIII_6.js -- faithful JS port of rdft/scalar/r2cb/r2cbIII_6.c
// (non-FMA branch). HC2RIII ("shifted") radix-6 direct codelet -- the
// BACKWARD counterpart of r2cfII_6.js. INPUT in[0..5] uses r2cfII_6.js's
// OUTPUT convention (in[k]=Cr[k] for k=0..2, in[5-k]=Ci[k] for k=0..2);
// OUTPUT out[0..5] uses r2cfII_6.js's INPUT convention (out[2k]=R0[k],
// out[2k+1]=R1[k]).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP2_000000000 = 2.000000000000000000000000000000000000000000000;
const KP1_732050807 = 1.732050807568877293527446341505872366942805254;

function r2cbIII_6(inArr) {
  const T1 = inArr[1];
  const T6 = inArr[4];
  const T2 = inArr[2];
  const T3 = inArr[0];
  const T4 = T2 + T3;
  const T5 = KP1_732050807 * (T2 - T3);
  const T7 = inArr[3];
  const T8 = inArr[5];
  const T9 = T7 + T8;
  const Tb = KP1_732050807 * (T7 - T8);

  const out = new Float64Array(6);
  out[0] = KP2_000000000 * (T1 + T4);
  const Ta = KP2_000000000 * T6 + T9; // FMA
  out[1] = -(T5 + Ta);
  const Tc = KP2_000000000 * T1 - T4; // FMS
  out[2] = Tb - Tc;
  out[3] = KP2_000000000 * (T6 - T9);
  out[4] = Tc + Tb;
  out[5] = T5 - Ta;
  return out;
}

module.exports = { r2cbIII_6 };
