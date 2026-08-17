'use strict';

// =============================================================================
// r2cb_5.js -- faithful JS port of rdft/scalar/r2cb/r2cb_5.c (non-FMA branch).
// O[0..4] packed halfcomplex (O[0]=Re0,O[1]=Re1,O[2]=Re2,O[3]=Im2,O[4]=Im1)
// -> x[0..4] real.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP2_000000000 = 2.0;
const KP1_118033988 = 1.118033988749894848204586834365638117720309180;
const KP500000000 = 0.5;
const KP1_902113032 = 1.902113032590307144232878666758764286811397268;
const KP1_175570504 = 1.175570504584946258337411909278145537195304875;

function r2cb_5(O) {
  const T8 = O[4], T9 = O[3];
  const Ta = KP1_175570504 * T8 - KP1_902113032 * T9;
  const Tc = KP1_902113032 * T8 + KP1_175570504 * T9;
  const T1 = O[0], T2 = O[1], T3 = O[2];
  const T4 = T2 + T3;
  const T5 = T1 - KP500000000 * T4;
  const T6 = KP1_118033988 * (T2 - T3);

  // R0 walks even indices in order {0,2,4}: R0[0]->x[0], R0[WS(rs,1)]->x[2],
  // R0[WS(rs,2)]->x[4] (confirmed against r2cf_5.js's identical mapping).
  const x = new Float64Array(5);
  x[0] = KP2_000000000 * T4 + T1; // R0[0]
  const Tb = T6 + T5;
  x[1] = Tb - Tc; // R1[0]
  x[4] = Tb + Tc; // R0[WS(rs,2)] -- even index 4
  const T7 = T5 - T6;
  x[2] = T7 - Ta; // R0[WS(rs,1)] -- even index 2
  x[3] = T7 + Ta; // R1[WS(rs,1)] -- odd index 3

  return x;
}

module.exports = { r2cb_5 };
