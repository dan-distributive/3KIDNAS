'use strict';

// =============================================================================
// r2cfII_6.js -- faithful JS port of rdft/scalar/r2cf/r2cfII_6.c (non-FMA
// branch). R2HCII ("shifted") radix-6 direct codelet -- the "cldm"
// middle-column combine for an EVEN outer radix (see r2cfII_8.js's header
// for the general even-radix derivation: clean 3-and-3 Cr/Ci split,
// out[5-k]=Ci[k]).
// INPUT: ph[0..5], R0[k]=ph[2k] (k=0..2), R1[k]=ph[2k+1] (k=0..2).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function r2cfII_6(ph) {
  const Ta = ph[3];
  const T5 = ph[5];
  const T6 = ph[1];
  const T7 = KP866025403 * (T5 - T6);
  const T9 = T5 + T6;
  const T1 = ph[0];
  const T3 = ph[2];
  const T2 = ph[4];
  const T8 = KP866025403 * (T2 + T3);
  const T4 = KP500000000 * (T3 - T2) + T1; // FMA
  const Tb = KP500000000 * T9 + Ta; // FMA

  const out = new Float64Array(6);
  out[0] = T4 - T7;
  out[1] = T1 + T2 - T3;
  out[2] = T4 + T7;
  out[3] = T8 - Tb;
  out[4] = Ta - T9;
  out[5] = -(T8 + Tb);
  return out;
}

module.exports = { r2cfII_6 };
