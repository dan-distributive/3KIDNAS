'use strict';

// =============================================================================
// hf_4.js -- faithful JS port of rdft/scalar/r2cf/hf_4.c (non-FMA branch).
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hf_2.js/hf_3.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

function hf_4(cr, ci, Wc, Ws) {
  const T1 = cr[0], Tp = ci[0];

  const T3 = cr[2], T5 = ci[2];
  const T6 = Wc[2] * T3 + Ws[2] * T5;
  const To = Wc[2] * T5 - Ws[2] * T3;

  const T9 = cr[1], Tb = ci[1];
  const Tc = Wc[1] * T9 + Ws[1] * Tb;
  const Tk = Wc[1] * Tb - Ws[1] * T9;

  const Te = cr[3], Tg = ci[3];
  const Th = Wc[3] * Te + Ws[3] * Tg;
  const Tl = Wc[3] * Tg - Ws[3] * Te;

  const outCr = new Float64Array(4), outCi = new Float64Array(4);

  const T7 = T1 + T6;
  const Ti = Tc + Th;
  outCi[1] = T7 - Ti;
  outCr[0] = T7 + Ti;
  const Tj = T1 - T6;
  const Tm = Tk - Tl;
  outCi[0] = Tj - Tm;
  outCr[1] = Tj + Tm;

  const Tn = Tk + Tl;
  const Tq = To + Tp;
  outCr[2] = Tn - Tq;
  outCi[3] = Tn + Tq;
  const Tr = Th - Tc;
  const Ts = Tp - To;
  outCr[3] = Tr - Ts;
  outCi[2] = Tr + Ts;

  return [outCr, outCi];
}

module.exports = { hf_4 };
