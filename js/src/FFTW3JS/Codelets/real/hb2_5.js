'use strict';

// =============================================================================
// hb2_5.js -- faithful JS port of rdft/scalar/r2cb/hb2_5.c (non-FMA branch).
// Alternate-codegen ("twiddle-log3/precompute-twiddles") sibling of
// hb_5.js -- same math, different rounding (see hf2_5.js's header for the
// forward-direction analogue). twinstr only trig-generates W^1, W^3; W^2
// and W^4 are DERIVED via complex products of those two. Same
// (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hb_5.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.25;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP559016994 = 0.559016994374947424102293417182819058860154590;

function hb2_5(cr, ci, Wc, Ws) {
  const Th = Wc[1], Tk = Ws[1], Ti = Wc[3], Tl = Ws[3];
  const Tj = Th * Ti, Tw = Tk * Ti, Tm = Tk * Tl, Tv = Th * Tl;
  const Tn = Tj + Tm, TP = Tv + Tw, Tx = Tv - Tw, TN = Tj - Tm;

  const T1 = cr[0];
  const T2 = cr[1], T3 = ci[0];
  const T4 = T2 + T3, Ty = T2 - T3;
  const T5 = cr[2], T6 = ci[1];
  const T7 = T5 + T6, Tz = T5 - T6;

  const Tp = KP559016994 * (T4 - T7);
  const TK = KP951056516 * Ty + KP587785252 * Tz;
  const TA = KP587785252 * Ty - KP951056516 * Tz;
  const T8 = T4 + T7;
  const To = T1 - KP250000000 * T8;

  const T9 = ci[4];
  const Ta = ci[3], Tb = cr[4];
  const Tc = Ta - Tb, Tr = Ta + Tb;
  const Td = ci[2], Te = cr[3];
  const Tf = Td - Te, Ts = Td + Te;

  const Tt = KP587785252 * Tr - KP951056516 * Ts;
  const TI = KP951056516 * Tr + KP587785252 * Ts;
  const TC = KP559016994 * (Tc - Tf);
  const Tg = Tc + Tf;
  const TB = T9 - KP250000000 * Tg;

  const outCr = new Float64Array(5), outCi = new Float64Array(5);
  outCr[0] = T1 + T8;
  outCi[0] = T9 + Tg;

  const Tq = To - Tp;
  const Tu = Tq - Tt, TF = Tq + Tt;
  const TD = TB - TC;
  const TE = TA + TD, TG = TD - TA;

  outCr[2] = Tn * Tu - Tx * TE;
  outCi[2] = Tn * TE + Tx * Tu;
  outCr[3] = Ti * TF - Tl * TG;
  outCi[3] = Ti * TG + Tl * TF;

  const TH = Tp + To;
  const TJ = TH - TI, TO = TH + TI;
  const TL = TC + TB;
  const TM = TK + TL, TQ = TL - TK;

  outCr[1] = Th * TJ - Tk * TM;
  outCi[1] = Th * TM + Tk * TJ;
  outCr[4] = TN * TO - TP * TQ;
  outCi[4] = TN * TQ + TP * TO;

  return [outCr, outCi];
}

module.exports = { hb2_5 };
