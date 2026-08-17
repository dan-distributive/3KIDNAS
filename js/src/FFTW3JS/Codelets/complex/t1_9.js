'use strict';

// =============================================================================
// t1_9.js -- faithful JS port of dft/scalar/codelets/t1_9.c (non-FMA
// branch), FFTW3's plain (TW_FULL) radix-9 Cooley-Tukey twiddle/combine
// codelet. Unlike the t2_* family, this reads ALL r-1=8 raw twiddle pairs
// directly from Wc[1..8]/Ws[1..8] (same convention as t1_5/t1_6/t1_7/t1_8),
// no precompute-twiddles derivation needed.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP939692620 = 0.939692620785908384054109277324731469936208134;
const KP342020143 = 0.342020143325668733044099614682259580763083368;
const KP984807753 = 0.984807753012208059366743024589523013670643252;
const KP173648177 = 0.173648177666930348851716626769314796000375677;
const KP642787609 = 0.642787609686539326322643409907263432907559884;
const KP766044443 = 0.766044443118978035202392650555416673935832457;
const KP500000000 = 0.5;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function t1_9(br, bi, Wc, Ws) {
  const T1 = br[0], T1B = bi[0];

  const T3 = br[3], T5 = bi[3], T2 = Wc[3], T4 = Ws[3];
  const T6 = T2 * T3 + T4 * T5;
  const TO = T2 * T5 - T4 * T3;

  const T8 = br[6], Ta = bi[6], T7 = Wc[6], T9 = Ws[6];
  const Tb = T7 * T8 + T9 * Ta;
  const TP = T7 * Ta - T9 * T8;

  const TQ = KP866025403 * (TO - TP);
  const T1G = KP866025403 * (Tb - T6);
  const Tc = T6 + Tb;
  const TN = T1 - KP500000000 * Tc;
  const T1A = TO + TP;
  const T1H = T1B - KP500000000 * T1A;

  const Tw = br[2], Ty = bi[2], Tv = Wc[2], Tx = Ws[2];
  const Tz = Tv * Tw + Tx * Ty;
  const T19 = Tv * Ty - Tx * Tw;

  const TB = br[5], TD = bi[5], TA = Wc[5], TC = Ws[5];
  const TE = TA * TB + TC * TD;
  const T14 = TA * TD - TC * TB;

  const TG = br[8], TI = bi[8], TF = Wc[8], TH = Ws[8];
  const TJ = TF * TG + TH * TI;
  const T15 = TF * TI - TH * TG;

  const TK = TE + TJ;
  const T1a = T14 + T15;
  const TL = Tz + TK;
  const T1x = T19 + T1a;

  const T13 = Tz - KP500000000 * TK;
  const T16 = KP866025403 * (T14 - T15);
  const T17 = T13 + T16;
  const T1o = T13 - T16;
  const T18 = KP866025403 * (TJ - TE);
  const T1b = T19 - KP500000000 * T1a;
  const T1c = T18 + T1b;
  const T1n = T1b - T18;

  const Tf = br[1], Th = bi[1], Te = Wc[1], Tg = Ws[1];
  const Ti = Te * Tf + Tg * Th;
  const TY = Te * Th - Tg * Tf;

  const Tk = br[4], Tm = bi[4], Tj = Wc[4], Tl = Ws[4];
  const Tn = Tj * Tk + Tl * Tm;
  const TT = Tj * Tm - Tl * Tk;

  const Tp = br[7], Tr = bi[7], To = Wc[7], Tq = Ws[7];
  const Ts = To * Tp + Tq * Tr;
  const TU = To * Tr - Tq * Tp;

  const Tt = Tn + Ts;
  const TZ = TT + TU;
  const Tu = Ti + Tt;
  const T1w = TY + TZ;

  const TS = Ti - KP500000000 * Tt;
  const TV = KP866025403 * (TT - TU);
  const TW = TS + TV;
  const T1k = TS - TV;
  const TX = KP866025403 * (Ts - Tn);
  const T10 = TY - KP500000000 * TZ;
  const T11 = TX + T10;
  const T1l = T10 - TX;

  const outR = new Float64Array(9), outI = new Float64Array(9);

  const T1y = KP866025403 * (T1w - T1x);
  const Td = T1 + Tc;
  const TM = Tu + TL;
  const T1v = Td - KP500000000 * TM;
  outR[0] = Td + TM;
  outR[3] = T1v + T1y;
  outR[6] = T1v - T1y;

  const T1D = KP866025403 * (TL - Tu);
  const T1z = T1w + T1x;
  const T1C = T1A + T1B;
  const T1E = T1C - KP500000000 * T1z;
  outI[0] = T1z + T1C;
  outI[6] = T1E - T1D;
  outI[3] = T1D + T1E;

  const TR = TN + TQ;
  const T1I = T1G + T1H;
  const T12 = KP766044443 * TW + KP642787609 * T11;
  const T1d = KP173648177 * T17 + KP984807753 * T1c;
  const T1e = T12 + T1d;
  const T1J = KP866025403 * (T1d - T12);
  const T1g = KP766044443 * T11 - KP642787609 * TW;
  const T1h = KP173648177 * T1c - KP984807753 * T17;
  const T1i = KP866025403 * (T1g - T1h);
  const T1F = T1g + T1h;

  outR[1] = TR + T1e;
  outI[1] = T1F + T1I;
  const T1f = TR - KP500000000 * T1e;
  outR[7] = T1f - T1i;
  outR[4] = T1f + T1i;
  const T1K = T1I - KP500000000 * T1F;
  outI[4] = T1J + T1K;
  outI[7] = T1K - T1J;

  const T1j = TN - TQ;
  const T1M = T1H - T1G;
  const T1m = KP173648177 * T1k + KP984807753 * T1l;
  const T1p = KP342020143 * T1n - KP939692620 * T1o;
  const T1q = T1m + T1p;
  const T1N = KP866025403 * (T1p - T1m);
  const T1s = KP173648177 * T1l - KP984807753 * T1k;
  const T1t = KP342020143 * T1o + KP939692620 * T1n;
  const T1u = KP866025403 * (T1s + T1t);
  const T1L = T1s - T1t;

  outR[2] = T1j + T1q;
  outI[2] = T1L + T1M;
  const T1r = T1j - KP500000000 * T1q;
  outR[8] = T1r - T1u;
  outR[5] = T1r + T1u;
  const T1O = T1M - KP500000000 * T1L;
  outI[5] = T1N + T1O;
  outI[8] = T1O - T1N;

  return [outR, outI];
}

module.exports = { t1_9 };
