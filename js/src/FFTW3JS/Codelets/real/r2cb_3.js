'use strict';

// =============================================================================
// r2cb_3.js -- faithful JS port of rdft/scalar/r2cb/r2cb_3.c (non-FMA branch).
// O[0..2] packed halfcomplex (O[0]=Re0,O[1]=Re1,O[2]=Im1) -> x[0..2] real.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP2_000000000 = 2.0;
const KP1_732050807 = 1.732050807568877293527446341505872366942805254;

function r2cb_3(O) {
  const T4 = O[2];
  const T5 = KP1_732050807 * T4;
  const T1 = O[0], T2 = O[1];
  const T3 = T1 - T2;

  const x = new Float64Array(3);
  x[0] = KP2_000000000 * T2 + T1;
  x[2] = T3 + T5; // R0[WS(rs,1)] -- even index 2
  x[1] = T3 - T5; // R1[0] -- odd index 1
  return x;
}

module.exports = { r2cb_3 };
