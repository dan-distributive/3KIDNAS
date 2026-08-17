'use strict';

// =============================================================================
// r2cf_6.js -- faithful JS port of rdft/scalar/r2cf/r2cf_6.c (non-FMA branch).
// x[0..5] (real) -> O[0..5] packed halfcomplex: O[0..3]=Re0..Re3,
// O[4..5]=Im2..Im1 (O[n-k]=Im_k, same convention as r2cf_5/7/9.js).
//
// Unlike the odd-radix codelets (r2cf_5/7/9.js), this EVEN-radix direct
// codelet's C source takes R0/R1 (not a flat array): rdft/direct-r2c.c's
// apply_r2hc calls ego->k(I, I+rs0, ...) with the codelet's own "rs" stride
// object set to 2*rs0 (X(mkstride)(n, 2*rs) in mkplan) -- so R0[k] reads
// x[2k] and R1[k] reads x[2k+1] (confirmed by walking mkplan's stride
// construction directly, same derivation as r2cfII_7.js's header). Every
// R0[WS(rs,k)]/R1[WS(rs,k)] below is transcribed as x[2*k]/x[2*k+1].
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function r2cf_6(x) {
  const T1 = x[0];
  const T2 = x[3];
  const T3 = T1 - T2;
  const Td = T1 + T2;
  const T7 = x[4];
  const T8 = x[1];
  const T9 = T7 - T8;
  const Tc = T7 + T8;
  const T4 = x[2];
  const T5 = x[5];
  const T6 = T4 - T5;
  const Tb = T4 + T5;

  const O = new Float64Array(6);
  O[5] = KP866025403 * (T9 - T6);
  const Ta = T6 + T9;
  O[1] = T3 - KP500000000 * Ta;
  O[3] = T3 + Ta;
  O[4] = KP866025403 * (Tb - Tc);
  const Te = Tb + Tc;
  O[2] = Td - KP500000000 * Te;
  O[0] = Td + Te;

  return O;
}

module.exports = { r2cf_6 };
