'use strict';

// =============================================================================
// hb_3.js -- faithful JS port of rdft/scalar/r2cb/hb_3.c (non-FMA branch).
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hb_2.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP500000000 = 0.5;

function hb_3(cr, ci, Wc, Ws) {
  const T1 = cr[0], T3 = ci[0];
  const T2 = cr[1];
  const T4 = T2 + T3;
  const Ta = T1 - KP500000000 * T4;
  const Te = KP866025403 * (T2 - T3);

  const T6 = ci[1], T7 = cr[2];
  const T5 = ci[2];
  const T8 = T6 - T7;
  const Tb = KP866025403 * (T6 + T7);
  const Tf = T5 - KP500000000 * T8;

  const outCr = new Float64Array(3), outCi = new Float64Array(3);
  outCr[0] = T1 + T4;
  outCi[0] = T5 + T8;

  const Tc = Ta - Tb;
  const Tg = Te + Tf;
  const T9 = Wc[1], Td = Ws[1];
  outCr[1] = T9 * Tc - Td * Tg;
  outCi[1] = T9 * Tg + Td * Tc;

  const Ti = Ta + Tb;
  const Tk = Tf - Te;
  const Th = Wc[2], Tj = Ws[2];
  outCr[2] = Th * Ti - Tj * Tk;
  outCi[2] = Th * Tk + Tj * Ti;

  return [outCr, outCi];
}

module.exports = { hb_3 };
