'use strict';

// =============================================================================
// hb_5.js -- faithful JS port of rdft/scalar/r2cb/hb_5.c (non-FMA branch).
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hb_2/3/4.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.25;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP559016994 = 0.559016994374947424102293417182819058860154590;

function hb_5(cr, ci, Wc, Ws) {
  const T1 = cr[0];

  const T2 = cr[1], T3 = ci[0];
  const T4 = T2 + T3;
  const Tq = T2 - T3;
  const T5 = cr[2], T6 = ci[1];
  const T7 = T5 + T6;
  const Tr = T5 - T6;

  const Tj = KP559016994 * (T4 - T7);
  const TG = KP951056516 * Tq + KP587785252 * Tr;
  const Ts = KP587785252 * Tq - KP951056516 * Tr;
  const T8 = T4 + T7;
  const Ti = T1 - KP250000000 * T8;

  const T9 = ci[4];
  const Ta = ci[3], Tb = cr[4];
  const Tc = Ta - Tb;
  const Tl = Ta + Tb;
  const Td = ci[2], Te = cr[3];
  const Tf = Td - Te;
  const Tm = Td + Te;

  const Tn = KP587785252 * Tl - KP951056516 * Tm;
  const TD = KP951056516 * Tl + KP587785252 * Tm;
  const Tu = KP559016994 * (Tc - Tf);
  const Tg = Tc + Tf;
  const Tt = T9 - KP250000000 * Tg;

  const outCr = new Float64Array(5), outCi = new Float64Array(5);
  outCr[0] = T1 + T8;
  outCi[0] = T9 + Tg;

  const Tk = Ti - Tj;
  const To = Tk - Tn;
  const Ty = Tk + Tn;
  const Tv = Tt - Tu;
  const Tw = Ts + Tv;
  const TA = Tv - Ts;

  const Th = Wc[2], Tp = Ws[2];
  outCr[2] = Th * To - Tp * Tw;
  outCi[2] = Th * Tw + Tp * To;
  const Tx = Wc[3], Tz = Ws[3];
  outCr[3] = Tx * Ty - Tz * TA;
  outCi[3] = Tx * TA + Tz * Ty;

  const TC = Tj + Ti;
  const TE = TC - TD;
  const TK = TC + TD;
  const TH = Tu + Tt;
  const TI = TG + TH;
  const TM = TH - TG;

  const TB = Wc[1], TF = Ws[1];
  outCr[1] = TB * TE - TF * TI;
  outCi[1] = TB * TI + TF * TE;
  const TJ = Wc[4], TL = Ws[4];
  outCr[4] = TJ * TK - TL * TM;
  outCi[4] = TJ * TM + TL * TK;

  return [outCr, outCi];
}

module.exports = { hb_5 };
