'use strict';

// =============================================================================
// r2cb_4.js -- faithful JS port of rdft/scalar/r2cb/r2cb_4.c (non-FMA branch).
// O[0..3] packed halfcomplex (O[0]=Re0,O[1]=Re1,O[2]=Re2,O[3]=Im1) ->
// x[0..3] real.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP2_000000000 = 2.0;

function r2cb_4(O) {
  const T4 = O[1];
  const T5 = KP2_000000000 * T4;
  const T7 = O[3];
  const T8 = KP2_000000000 * T7;
  const T1 = O[0], T2 = O[2];
  const T3 = T1 + T2;
  const T6 = T1 - T2;

  const x = new Float64Array(4);
  x[2] = T3 - T5; // R0[WS(rs,1)] -- even index 2
  x[3] = T6 + T8; // R1[WS(rs,1)] -- odd index 3
  x[0] = T3 + T5; // R0[0]
  x[1] = T6 - T8; // R1[0]
  return x;
}

module.exports = { r2cb_4 };
