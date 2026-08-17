'use strict';

// =============================================================================
// hf2_5.js -- faithful JS port of rdft/scalar/r2cf/hf2_5.c (non-FMA branch).
// Alternate-codegen ("twiddle-log3/precompute-twiddles") sibling of hf_5.js
// -- same math, different rounding, matching complex-side t2_r vs t1_r
// (see t2_5.js's header). twinstr only trig-generates W^1, W^3 (two raw
// pairs); the other needed multiple (W^2, W^4) is DERIVED via complex
// products of those two, same trick as t2_5.js. Same (cr,ci,Wc,Ws) ->
// (outCr,outCi) convention as hf_5.js -- only indices 1, 3 are actually
// read from Wc/Ws; every other index is deliberately ignored in favor of
// recomputing the equivalent product chain, to match FFTW's exact rounding.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;

function hf2_5(cr, ci, Wc, Ws) {
  const T2 = Wc[1], T4 = Ws[1], T7 = Wc[3], T9 = Ws[3];
  const T8 = T2 * T7, Te = T4 * T7, Ta = T4 * T9, Td = T2 * T9;
  const Tb = T8 - Ta, Tl = Td - Te, Tf = Td + Te, Tj = T8 + Ta;

  const T1 = cr[0], TI = ci[0];

  const T3 = cr[1], T5 = ci[1];
  const T6 = T2 * T3 + T4 * T5;
  const Tw = T2 * T5 - T4 * T3;
  const To = cr[3], Tp = ci[3];
  const Tq = T7 * To + T9 * Tp;
  const TA = T7 * Tp - T9 * To;

  const Tc = cr[4], Tg = ci[4];
  const Th = Tb * Tc + Tf * Tg;
  const Tx = Tb * Tg - Tf * Tc;
  const Tk = cr[2], Tm = ci[2];
  const Tn = Tj * Tk + Tl * Tm;
  const Tz = Tj * Tm - Tl * Tk;

  const Ty = Tw - Tx;
  const TB = Tz - TA;
  const TG = Tn - Tq;
  const TF = Th - T6;
  const TJ = Tw + Tx;
  const TK = Tz + TA;
  const TL = TJ + TK;
  const Ti = T6 + Th;
  const Tr = Tn + Tq;
  const Ts = Ti + Tr;

  const outCr = new Float64Array(5), outCi = new Float64Array(5);
  outCr[0] = T1 + Ts;

  const TC = KP951056516 * Ty + KP587785252 * TB;
  const TE = KP951056516 * TB - KP587785252 * Ty;
  const Tt = KP559016994 * (Ti - Tr);
  const Tu = T1 - KP250000000 * Ts;
  const Tv = Tt + Tu, TD = Tu - Tt;
  outCi[0] = Tv - TC;
  outCi[1] = TD + TE;
  outCr[1] = Tv + TC;
  outCr[2] = TD - TE;

  outCi[4] = TL + TI;

  const TH = KP587785252 * TF + KP951056516 * TG;
  const TP = KP951056516 * TF - KP587785252 * TG;
  const TM = TI - KP250000000 * TL;
  const TN = KP559016994 * (TJ - TK);
  const TO = TM - TN, TQ = TN + TM;
  outCr[3] = TH - TO;
  outCi[3] = TP + TQ;
  outCi[2] = TH + TO;
  outCr[4] = TP - TQ;

  return [outCr, outCi];
}

module.exports = { hf2_5 };
