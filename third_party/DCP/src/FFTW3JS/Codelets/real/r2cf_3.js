'use strict';

// =============================================================================
// r2cf_3.js -- faithful JS port of rdft/scalar/r2cf/r2cf_3.c (non-FMA branch).
// x[0..2] (real) -> O[0..2] packed halfcomplex: O[0]=Re0, O[1]=Re1, O[2]=Im1.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP500000000 = 0.5;

function r2cf_3(x) {
  const T1 = x[0], T2 = x[1], T3 = x[2];
  const T4 = T2 + T3;
  const O = new Float64Array(3);
  O[1] = T1 - KP500000000 * T4;
  O[2] = KP866025403 * (T3 - T2);
  O[0] = T1 + T4;
  return O;
}

module.exports = { r2cf_3 };
