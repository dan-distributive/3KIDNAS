'use strict';

// =============================================================================
// hb_4.js -- faithful JS port of rdft/scalar/r2cb/hb_4.c (non-FMA branch).
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hb_2/3.js. Note the
// input pairing here genuinely differs from hf_4/hb_2/hb_3 (cr[0] combines
// with ci[1], not ci[0]) -- transcribed literally by physical index exactly
// as the C source reads them, not "cleaned up" to a uniform pattern.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

function hb_4(cr, ci, Wc, Ws) {
  const T1 = cr[0], T2 = ci[1];
  const T3 = T1 + T2;
  const Ti = T1 - T2;
  const T4 = cr[1], T5 = ci[0];
  const T6 = T4 + T5;
  const Tm = T4 - T5;
  const Ta = ci[3], Tb = cr[2];
  const Tc = Ta - Tb;
  const Tn = Ta + Tb;
  const Td = ci[2], Te = cr[3];
  const Tf = Td - Te;
  const Tj = Td + Te;

  const outCr = new Float64Array(4), outCi = new Float64Array(4);
  outCr[0] = T3 + T6;
  outCi[0] = Tc + Tf;

  const T8 = T3 - T6, Tg = Tc - Tf;
  const T7 = Wc[2], T9 = Ws[2];
  outCr[2] = T7 * T8 - T9 * Tg;
  outCi[2] = T9 * T8 + T7 * Tg;

  const Tk = Ti - Tj, To = Tm + Tn;
  const Th = Wc[1], Tl = Ws[1];
  outCr[1] = Th * Tk - Tl * To;
  outCi[1] = Th * To + Tl * Tk;

  const Tq = Ti + Tj, Ts = Tn - Tm;
  const Tp = Wc[3], Tr = Ws[3];
  outCr[3] = Tp * Tq - Tr * Ts;
  outCi[3] = Tp * Ts + Tr * Tq;

  return [outCr, outCi];
}

module.exports = { hb_4 };
