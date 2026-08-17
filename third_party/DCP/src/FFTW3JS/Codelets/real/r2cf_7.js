'use strict';

// =============================================================================
// r2cf_7.js -- faithful JS port of rdft/scalar/r2cf/r2cf_7.c (non-FMA branch).
// x[0..6] (real) -> O[0..6] packed halfcomplex: O[0..3]=Re0..Re3,
// O[4..6]=Im3..Im1 (O[n-k]=Im_k, same convention as r2cf_5.js).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP222520933 = 0.222520933956314404288902564496794759466355569;
const KP900968867 = 0.900968867902419126236102319507445051165919162;
const KP623489801 = 0.623489801858733530525004884004239810632274731;
const KP433883739 = 0.433883739117558120475768332848358754609990728;
const KP781831482 = 0.781831482468029808708444526674057750232334519;
const KP974927912 = 0.974927912181823607018131682993931217232785801;

function r2cf_7(x) {
  const T1 = x[0];
  const T8 = x[1], T9 = x[6];
  const Ta = T8 + T9;
  const Tb = T9 - T8;

  const T2 = x[2], T3 = x[5];
  const T4 = T2 + T3;
  const Td = T3 - T2;
  const T5 = x[3], T6 = x[4];
  const T7 = T5 + T6;
  const Tc = T6 - T5;

  const O = new Float64Array(7);
  const M5 = KP974927912 * Tb - KP781831482 * Tc;
  O[5] = M5 - KP433883739 * Td;
  const M6 = KP781831482 * Tb + KP974927912 * Td;
  O[6] = M6 + KP433883739 * Tc;
  const M2a = KP623489801 * T7 + T1, M2b = KP900968867 * T4 + KP222520933 * Ta;
  O[2] = M2a - M2b;
  const M4 = KP433883739 * Tb + KP974927912 * Tc;
  O[4] = M4 - KP781831482 * Td;
  const M3a = KP623489801 * T4 + T1, M3b = KP222520933 * T7 + KP900968867 * Ta;
  O[3] = M3a - M3b;
  const M1a = KP623489801 * Ta + T1, M1b = KP900968867 * T7 + KP222520933 * T4;
  O[1] = M1a - M1b;
  O[0] = T1 + Ta + T4 + T7;
  return O;
}

module.exports = { r2cf_7 };
