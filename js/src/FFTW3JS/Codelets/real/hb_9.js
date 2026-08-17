'use strict';

// =============================================================================
// hb_9.js -- faithful JS port of rdft/scalar/r2cb/hb_9.c (non-FMA branch).
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hb_2/3/4/5/7.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP984807753 = 0.984807753012208059366743024589523013670643252;
const KP173648177 = 0.173648177666930348851716626769314796000375677;
const KP342020143 = 0.342020143325668733044099614682259580763083368;
const KP939692620 = 0.939692620785908384054109277324731469936208134;
const KP642787609 = 0.642787609686539326322643409907263432907559884;
const KP766044443 = 0.766044443118978035202392650555416673935832457;
const KP500000000 = 0.5;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function hb_9(cr, ci, Wc, Ws) {
  const T1 = cr[0];
  const Th = ci[8];
  const T2 = cr[3], T3 = ci[2];
  const T4 = T2 + T3;
  const T14 = KP866025403 * (T2 - T3);
  const Ti = ci[5], Tj = cr[6];
  const Tk = Ti - Tj;
  const TL = KP866025403 * (Ti + Tj);

  const T5 = T1 + T4;
  const Tl = Th + Tk;
  const TK = T1 - KP500000000 * T4;
  const TM = TK - TL, T1o = TK + TL;
  const T15 = Th - KP500000000 * Tk;
  const T16 = T14 + T15, T1y = T15 - T14;

  const T6 = cr[1];
  const T7 = cr[4], T8 = ci[1];
  const T9 = T7 + T8;
  const TN = T6 - KP500000000 * T9;
  const TQ = KP866025403 * (T7 - T8);
  const Tm = ci[7], Tn = ci[4], To = cr[7];
  const Tp = Tn - To;
  const TO = KP866025403 * (Tn + To);
  const TR = Tm - KP500000000 * Tp;

  const Tb = cr[2], Tc = ci[3], Td = ci[0];
  const Te = Tc + Td;
  const TU = Tb - KP500000000 * Te;
  const TX = KP866025403 * (Tc - Td);
  const Tr = ci[6], Ts = cr[5], Tt = cr[8];
  const Tu = Ts + Tt;
  const TV = KP866025403 * (Ts - Tt);
  const TY = KP500000000 * Tu + Tr;

  const Ta = T6 + T9;
  const Tf = Tb + Te;
  const Tg = Ta + Tf;
  const Tq = Tm + Tp;
  const Tv = Tr - Tu;
  const Tw = Tq + Tv;
  const TP = TN - TO;
  const TS = TQ + TR;
  const TT = KP766044443 * TP - KP642787609 * TS;
  const T17 = KP766044443 * TS + KP642787609 * TP;
  const T1s = TU - TV;
  const T1t = TY - TX;
  const T1u = KP939692620 * T1s + KP342020143 * T1t;
  const T1A = KP342020143 * T1s - KP939692620 * T1t;

  const T1p = TN + TO;
  const T1q = TR - TQ;
  const T1r = KP173648177 * T1p - KP984807753 * T1q;
  const T1z = KP173648177 * T1q + KP984807753 * T1p;
  const TW = TU + TV;
  const TZ = TX + TY;
  const T10 = KP173648177 * TW - KP984807753 * TZ;
  const T18 = KP984807753 * TW + KP173648177 * TZ;

  const outCr = new Float64Array(9), outCi = new Float64Array(9);
  outCr[0] = T5 + Tg;
  outCi[0] = Tl + Tw;

  const Ty = T5 - KP500000000 * Tg;
  const Tz = KP866025403 * (Tv - Tq);
  const TA = Ty - Tz, TG = Ty + Tz;
  const TC = Tl - KP500000000 * Tw;
  const TD = KP866025403 * (Ta - Tf);
  const TE = TC - TD, TI = TD + TC;

  const Tx = Wc[6], TB = Ws[6];
  outCr[6] = Tx * TA - TB * TE;
  outCi[6] = Tx * TE + TB * TA;
  const TF = Wc[3], TH = Ws[3];
  outCr[3] = TF * TG - TH * TI;
  outCi[3] = TF * TI + TH * TG;

  const T1d = KP866025403 * (T18 - T17);
  const T1h = KP866025403 * (TT - T10);
  const T11 = TT + T10;
  const T12 = TM + T11;
  const T1c = TM - KP500000000 * T11;
  const T19 = T17 + T18;
  const T1a = T16 + T19;
  const T1g = T16 - KP500000000 * T19;
  const TJ = Wc[1], T13 = Ws[1];
  outCr[1] = TJ * T12 - T13 * T1a;
  outCi[1] = TJ * T1a + T13 * T12;

  const T1k = T1c + T1d;
  const T1m = T1h + T1g;
  const T1j = Wc[4], T1l = Ws[4];
  outCr[4] = T1j * T1k - T1l * T1m;
  outCi[4] = T1j * T1m + T1l * T1k;

  const T1e = T1c - T1d;
  const T1i = T1g - T1h;
  const T1b = Wc[7], T1f = Ws[7];
  outCr[7] = T1b * T1e - T1f * T1i;
  outCi[7] = T1b * T1i + T1f * T1e;

  const T1F = KP866025403 * (T1A - T1z);
  const T1J = KP866025403 * (T1r + T1u);
  const T1v = T1r - T1u;
  const T1w = T1o + T1v;
  const T1E = T1o - KP500000000 * T1v;
  const T1B = T1z + T1A;
  const T1C = T1y + T1B;
  const T1I = T1y - KP500000000 * T1B;
  const T1n = Wc[2], T1x = Ws[2];
  outCr[2] = T1n * T1w - T1x * T1C;
  outCi[2] = T1n * T1C + T1x * T1w;

  const T1M = T1F + T1E;
  const T1O = T1I + T1J;
  const T1L = Wc[5], T1N = Ws[5];
  outCr[5] = T1L * T1M - T1N * T1O;
  outCi[5] = T1L * T1O + T1N * T1M;

  const T1G = T1E - T1F;
  const T1K = T1I - T1J;
  const T1D = Wc[8], T1H = Ws[8];
  outCr[8] = T1D * T1G - T1H * T1K;
  outCi[8] = T1D * T1K + T1H * T1G;

  return [outCr, outCi];
}

module.exports = { hb_9 };
