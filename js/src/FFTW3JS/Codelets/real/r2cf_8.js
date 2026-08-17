'use strict';

// =============================================================================
// r2cf_8.js -- faithful JS port of rdft/scalar/r2cf/r2cf_8.c (non-FMA branch).
// x[0..7] (real) -> O[0..7] packed halfcomplex: O[0..4]=Re0..Re4,
// O[5..7]=Im3..Im1 (O[n-k]=Im_k). R0/R1 (stride-2) input convention -- see
// r2cf_6.js's header for the direct-r2c.c derivation.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP707106781 = 0.707106781186547524400844362104849039284835938;

function r2cf_8(x) {
  const T1 = x[0];
  const T2 = x[4];
  const T3 = T1 + T2;
  const T7 = T1 - T2;
  const Tb = x[7];
  const Tc = x[3];
  const Td = Tb - Tc;
  const Tj = Tb + Tc;
  const T4 = x[2];
  const T5 = x[6];
  const T6 = T4 + T5;
  const Tg = T4 - T5;
  const T8 = x[1];
  const T9 = x[5];
  const Ta = T8 - T9;
  const Ti = T8 + T9;

  const O = new Float64Array(8);
  O[2] = T3 - T6;
  O[6] = Tj - Ti;
  const Te = KP707106781 * (Ta + Td);
  O[3] = T7 - Te;
  O[1] = T7 + Te;
  const Tf = KP707106781 * (Td - Ta);
  O[7] = Tf - Tg;
  O[5] = Tg + Tf;
  const Th = T3 + T6;
  const Tk = Ti + Tj;
  O[4] = Th - Tk;
  O[0] = Th + Tk;

  return O;
}

module.exports = { r2cf_8 };
