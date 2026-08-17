'use strict';

// =============================================================================
// hf_5.js -- faithful JS port of rdft/scalar/r2cf/hf_5.c (non-FMA branch).
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hf_2/3/4.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;

function hf_5(cr, ci, Wc, Ws) {
  const T1 = cr[0], TE = ci[0];

  const T3 = cr[1], T5 = ci[1];
  const T6 = Wc[1] * T3 + Ws[1] * T5;
  const Ts = Wc[1] * T5 - Ws[1] * T3;

  const Tj = cr[3], Tl = ci[3];
  const Tm = Wc[3] * Tj + Ws[3] * Tl;
  const Tw = Wc[3] * Tl - Ws[3] * Tj;

  const T8 = cr[4], Ta = ci[4];
  const Tb = Wc[4] * T8 + Ws[4] * Ta;
  const Tt = Wc[4] * Ta - Ws[4] * T8;

  const Te = cr[2], Tg = ci[2];
  const Th = Wc[2] * Te + Ws[2] * Tg;
  const Tv = Wc[2] * Tg - Ws[2] * Te;

  const Tu = Ts - Tt;
  const Tx = Tv - Tw;
  const TC = Th - Tm;
  const TB = Tb - T6;
  const TF = Ts + Tt;
  const TG = Tv + Tw;
  const TH = TF + TG;
  const Tc = T6 + Tb;
  const Tn = Th + Tm;
  const To = Tc + Tn;

  const outCr = new Float64Array(5), outCi = new Float64Array(5);
  outCr[0] = T1 + To;

  const Ty = KP951056516 * Tu + KP587785252 * Tx;
  const TA = KP951056516 * Tx - KP587785252 * Tu;
  const Tp = KP559016994 * (Tc - Tn);
  const Tq = T1 - KP250000000 * To;
  const Tr = Tp + Tq;
  const Tz = Tq - Tp;
  outCi[0] = Tr - Ty;
  outCi[1] = Tz + TA;
  outCr[1] = Tr + Ty;
  outCr[2] = Tz - TA;

  outCi[4] = TH + TE;

  const TD = KP587785252 * TB + KP951056516 * TC;
  const TL = KP951056516 * TB - KP587785252 * TC;
  const TI = TE - KP250000000 * TH;
  const TJ = KP559016994 * (TF - TG);
  const TK = TI - TJ;
  const TM = TJ + TI;
  outCr[3] = TD - TK;
  outCi[3] = TL + TM;
  outCi[2] = TD + TK;
  outCr[4] = TL - TM;

  return [outCr, outCi];
}

module.exports = { hf_5 };
