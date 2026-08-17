'use strict';

// =============================================================================
// r2cf_4.js -- faithful JS port of rdft/scalar/r2cf/r2cf_4.c (non-FMA branch).
// x[0..3] (real) -> O[0..3] packed halfcomplex: O[0]=Re0,O[1]=Re1,O[2]=Re2
// (Nyquist), O[3]=Im1. Zero multiplications, so no FMA concerns.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

function r2cf_4(x) {
  const T1 = x[0], T2 = x[2];
  const T3 = T1 + T2;
  const T4 = x[1], T5 = x[3];
  const T6 = T4 + T5;
  const O = new Float64Array(4);
  O[1] = T1 - T2;
  O[3] = T5 - T4;
  O[2] = T3 - T6;
  O[0] = T3 + T6;
  return O;
}

module.exports = { r2cf_4 };
