'use strict';

// =============================================================================
// r2cbIII_5.js -- faithful JS port of rdft/scalar/r2cb/r2cbIII_5.c
// (non-FMA branch). HC2RIII ("shifted") radix-5 direct codelet -- the
// BACKWARD counterpart of r2cfII_5.js. INPUT in[0..4] uses r2cfII_5.js's
// OUTPUT convention (in[0..2]=Cr, in[3]=Ci[1], in[4]=Ci[0]); OUTPUT
// out[0..4] uses r2cfII_5.js's INPUT convention (out[2k]=R0[k],
// out[2k+1]=R1[k]).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP2_000000000 = 2.000000000000000000000000000000000000000000000;
const KP1_118033988 = 1.118033988749894848204586834365638117720309180;
const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP1_175570504 = 1.175570504584946258337411909278145537195304875;
const KP1_902113032 = 1.902113032590307144232878666758764286811397268;

function r2cbIII_5(inArr) {
  const T8 = inArr[3];
  const T9 = inArr[4];
  const Ta = KP1_902113032 * T8 + KP1_175570504 * T9; // FMA
  const Tc = KP1_175570504 * T8 - KP1_902113032 * T9; // FNMS
  const T1 = inArr[2];
  const T2 = inArr[1];
  const T3 = inArr[0];
  const T4 = T2 + T3;
  const T5 = KP500000000 * T4 - T1; // FMS
  const T6 = KP1_118033988 * (T3 - T2);

  const out = new Float64Array(5);
  out[0] = KP2_000000000 * T4 + T1; // FMA
  const Tb = T6 - T5;
  out[2] = Tb + Tc;
  out[3] = Tc - Tb;
  const T7 = T5 + T6;
  out[1] = T7 - Ta;
  out[4] = -(T7 + Ta);

  return out;
}

module.exports = { r2cbIII_5 };
