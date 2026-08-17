'use strict';

// =============================================================================
// hf_3.js -- faithful JS port of rdft/scalar/r2cf/hf_3.c (non-FMA branch).
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hf_2.js -- see that
// file's header for why this is transcribed literally rather than derived
// by analogy with t1_3.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP500000000 = 0.5;

function hf_3(cr, ci, Wc, Ws) {
  const T1 = cr[0], Ti = ci[0];

  const T3 = cr[1], T5 = ci[1];
  const T6 = Wc[1] * T3 + Ws[1] * T5;
  const Te = Wc[1] * T5 - Ws[1] * T3;

  const T8 = cr[2], Ta = ci[2];
  const Tb = Wc[2] * T8 + Ws[2] * Ta;
  const Tf = Wc[2] * Ta - Ws[2] * T8;

  const Tc = T6 + Tb;
  const Tj = Te + Tf;

  const outCr = new Float64Array(3), outCi = new Float64Array(3);
  outCr[0] = T1 + Tc;
  const Td = T1 - KP500000000 * Tc;
  const Tg = KP866025403 * (Te - Tf);
  outCi[0] = Td - Tg;
  outCr[1] = Td + Tg;
  outCi[2] = Tj + Ti;
  const Th = KP866025403 * (Tb - T6);
  const Tk = Ti - KP500000000 * Tj;
  outCr[2] = Th - Tk;
  outCi[1] = Th + Tk;

  return [outCr, outCi];
}

module.exports = { hf_3 };
