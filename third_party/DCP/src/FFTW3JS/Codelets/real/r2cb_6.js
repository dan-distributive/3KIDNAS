'use strict';

// =============================================================================
// r2cb_6.js -- faithful JS port of rdft/scalar/r2cb/r2cb_6.c (non-FMA branch).
// O[0..5] packed halfcomplex (O[0..2]=Re0..Re2, O[3..5]=Im2..Im0 via
// O[n-k]=Im_k) -> x[0..5] real, UNNORMALIZED.
//
// Output uses the same R0/R1 (stride-2) convention as r2cb_5.js: R0[k]
// writes x[2k], R1[k] writes x[2k+1] -- see r2cf_6.js's header for the
// direct-r2c.c derivation (same stride construction, roles swapped).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP2_000000000 = 2.000000000000000000000000000000000000000000000;
const KP1_732050807 = 1.732050807568877293527446341505872366942805254;

function r2cb_6(O) {
  const T1 = O[0];
  const T2 = O[3];
  const T3 = T1 - T2;
  const T7 = T1 + T2;
  const Ta = O[4];
  const Tb = O[5];
  const Tc = KP1_732050807 * (Ta - Tb);
  const Te = KP1_732050807 * (Ta + Tb);
  const T4 = O[2];
  const T5 = O[1];
  const T6 = T4 - T5;
  const T8 = T4 + T5;

  const x = new Float64Array(6);
  x[3] = KP2_000000000 * T6 + T3; // R1[WS(rs,1)]
  x[0] = KP2_000000000 * T8 + T7; // R0[0]
  const T9 = T7 - T8;
  x[4] = T9 - Tc; // R0[WS(rs,2)]
  x[2] = T9 + Tc; // R0[WS(rs,1)]
  const Td = T3 - T6;
  x[1] = Td - Te; // R1[0]
  x[5] = Td + Te; // R1[WS(rs,2)]

  return x;
}

module.exports = { r2cb_6 };
