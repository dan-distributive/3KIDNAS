'use strict';

// =============================================================================
// hf_9.js -- faithful JS port of rdft/scalar/r2cf/hf_9.c (non-FMA branch).
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hf_2/3/4/5/7.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP642787609 = 0.642787609686539326322643409907263432907559884;
const KP766044443 = 0.766044443118978035202392650555416673935832457;
const KP939692620 = 0.939692620785908384054109277324731469936208134;
const KP342020143 = 0.342020143325668733044099614682259580763083368;
const KP984807753 = 0.984807753012208059366743024589523013670643252;
const KP173648177 = 0.173648177666930348851716626769314796000375677;
const KP500000000 = 0.5;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function hf_9(cr, ci, Wc, Ws) {
  const T1 = cr[0], T1B = ci[0];

  const T3 = cr[3], T5 = ci[3];
  const T6 = Wc[3] * T3 + Ws[3] * T5;
  const TO = Wc[3] * T5 - Ws[3] * T3;
  const T8 = cr[6], Ta = ci[6];
  const Tb = Wc[6] * T8 + Ws[6] * Ta;
  const TP = Wc[6] * Ta - Ws[6] * T8;

  const TQ = KP866025403 * (TO - TP);
  const T1A = KP866025403 * (Tb - T6);
  const Tc = T6 + Tb;
  const TN = T1 - KP500000000 * Tc;
  const T1C = TO + TP;
  const T1D = T1B - KP500000000 * T1C;

  const Tw = cr[2], Ty = ci[2];
  const Tz = Wc[2] * Tw + Ws[2] * Ty;
  const T13 = Wc[2] * Ty - Ws[2] * Tw;
  const TB = cr[5], TD = ci[5];
  const TE = Wc[5] * TB + Ws[5] * TD;
  const T14 = Wc[5] * TD - Ws[5] * TB;
  const TG = cr[8], TI = ci[8];
  const TJ = Wc[8] * TG + Ws[8] * TI;
  const T15 = Wc[8] * TI - Ws[8] * TG;

  const TK = TE + TJ;
  const T16 = T14 + T15;
  const TL = Tz + TK;
  const T1x = T13 + T16;
  const T17 = T13 - KP500000000 * T16;
  const T18 = KP866025403 * (TJ - TE);
  const T19 = T17 - T18, T1o = T18 + T17;
  const T1a = Tz - KP500000000 * TK;
  const T1b = KP866025403 * (T14 - T15);
  const T1c = T1a - T1b, T1n = T1a + T1b;

  const Tf = cr[1], Th = ci[1];
  const Ti = Wc[1] * Tf + Ws[1] * Th;
  const TX = Wc[1] * Th - Ws[1] * Tf;
  const Tk = cr[4], Tm = ci[4];
  const Tn = Wc[4] * Tk + Ws[4] * Tm;
  const TT = Wc[4] * Tm - Ws[4] * Tk;
  const Tp = cr[7], Tr = ci[7];
  const Ts = Wc[7] * Tp + Ws[7] * Tr;
  const TU = Wc[7] * Tr - Ws[7] * Tp;

  const Tt = Tn + Ts;
  const TY = TT + TU;
  const Tu = Ti + Tt;
  const T1w = TX + TY;
  const TS = Ti - KP500000000 * Tt;
  const TV = KP866025403 * (TT - TU);
  const TW = TS - TV, T1k = TS + TV;
  const TZ = TX - KP500000000 * TY;
  const T10 = KP866025403 * (Ts - Tn);
  const T11 = TZ - T10, T1l = T10 + TZ;

  const outCr = new Float64Array(9), outCi = new Float64Array(9);

  const T1y = KP866025403 * (T1w - T1x);
  const Td = T1 + Tc;
  const TM = Tu + TL;
  const T1v = Td - KP500000000 * TM;
  outCr[0] = Td + TM;
  outCr[3] = T1v + T1y;
  outCi[2] = T1v - T1y;

  const TR = TN - TQ;
  const T1I = T1D - T1A;
  const T12 = KP173648177 * TW + KP984807753 * T11;
  const T1d = KP342020143 * T19 - KP939692620 * T1c;
  const T1e = T12 + T1d;
  const T1K = KP866025403 * (T1d - T12);
  const T1g = KP173648177 * T11 - KP984807753 * TW;
  const T1h = KP342020143 * T1c + KP939692620 * T19;
  const T1i = KP866025403 * (T1g + T1h);
  const T1H = T1g - T1h;

  outCr[2] = TR + T1e;
  outCi[6] = T1H + T1I;
  const T1f = TR - KP500000000 * T1e;
  outCi[0] = T1f - T1i;
  outCi[3] = T1f + T1i;
  const T1J = KP500000000 * T1H - T1I;
  outCr[5] = T1J - T1K;
  outCr[8] = T1K + T1J;

  const T1L = KP866025403 * (TL - Tu);
  const T1M = T1C + T1B;
  const T1N = T1w + T1x;
  const T1O = T1M - KP500000000 * T1N;
  outCr[6] = T1L - T1O;
  outCi[8] = T1N + T1M;
  outCi[5] = T1L + T1O;

  const T1j = TN + TQ;
  const T1E = T1A + T1D;
  const T1m = KP766044443 * T1k + KP642787609 * T1l;
  const T1p = KP173648177 * T1n + KP984807753 * T1o;
  const T1q = T1m + T1p;
  const T1z = KP866025403 * (T1p - T1m);
  const T1s = KP766044443 * T1l - KP642787609 * T1k;
  const T1t = KP173648177 * T1o - KP984807753 * T1n;
  const T1u = KP866025403 * (T1s - T1t);
  const T1F = T1s + T1t;

  outCr[1] = T1j + T1q;
  const T1r = T1j - KP500000000 * T1q;
  outCi[1] = T1r - T1u;
  outCr[4] = T1r + T1u;
  outCi[7] = T1F + T1E;
  const T1G = T1E - KP500000000 * T1F;
  outCr[7] = T1z - T1G;
  outCi[4] = T1z + T1G;

  return [outCr, outCi];
}

module.exports = { hf_9 };
