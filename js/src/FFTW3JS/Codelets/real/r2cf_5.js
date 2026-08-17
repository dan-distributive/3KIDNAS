'use strict';

// =============================================================================
// r2cf_5.js -- faithful JS port of rdft/scalar/r2cf/r2cf_5.c (non-FMA branch).
// x[0..4] (real) -> O[0..4] packed halfcomplex: O[0]=Re0,O[1]=Re1,O[2]=Re2,
// O[3]=Im2,O[4]=Im1.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;

function r2cf_5(x) {
  const Ta = x[0];
  const T1 = x[4], T2 = x[1];
  const T7 = T2 + T1;
  const T4 = x[2], T5 = x[3];
  const T8 = T4 + T5;
  const T3 = T1 - T2;
  const Tb = T7 + T8;
  const T6 = T4 - T5;

  const O = new Float64Array(5);
  O[4] = KP951056516 * T3 - KP587785252 * T6;
  O[0] = Ta + Tb;
  O[3] = KP951056516 * T6 + KP587785252 * T3;
  const T9 = KP559016994 * (T7 - T8);
  const Tc = Ta - KP250000000 * Tb;
  O[1] = T9 + Tc;
  O[2] = Tc - T9;
  return O;
}

module.exports = { r2cf_5 };
