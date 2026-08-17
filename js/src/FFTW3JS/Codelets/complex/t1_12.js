'use strict';

// =============================================================================
// t1_12.js -- faithful JS port of dft/scalar/codelets/t1_12.c (non-FMA
// branch), FFTW3's plain (TW_FULL) radix-12 Cooley-Tukey twiddle/combine
// codelet. Reads all r-1=11 raw twiddle pairs directly from Wc[1..11]/
// Ws[1..11] (same convention as t1_5/t1_6/t1_7/t1_8/t1_9).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP500000000 = 0.5;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function t1_12(br, bi, Wc, Ws) {
  const T1 = br[0], T1W = bi[0];

  const T3 = br[4], T5 = bi[4], T2 = Wc[4], T4 = Ws[4];
  const T6 = T2 * T3 + T4 * T5;
  const T16 = T2 * T5 - T4 * T3;

  const T8 = br[8], Ta = bi[8], T7 = Wc[8], T9 = Ws[8];
  const Tb = T7 * T8 + T9 * Ta;
  const T17 = T7 * Ta - T9 * T8;

  const T18 = KP866025403 * (T16 - T17);
  const T21 = KP866025403 * (Tb - T6);
  const Tc = T6 + Tb;
  const T15 = T1 - KP500000000 * Tc;
  const T1V = T16 + T17;
  const T22 = T1W - KP500000000 * T1V;

  const TO = br[9], TQ = bi[9], TN = Wc[9], TP = Ws[9];
  const TR = TN * TO + TP * TQ;
  const T1E = TN * TQ - TP * TO;

  const TY = br[5], T10 = bi[5], TX = Wc[5], TZ = Ws[5];
  const T11 = TX * TY + TZ * T10;
  const T1n = TX * T10 - TZ * TY;

  const TT = br[1], TV = bi[1], TS = Wc[1], TU = Ws[1];
  const TW = TS * TT + TU * TV;
  const T1m = TS * TV - TU * TT;

  const T1o = KP866025403 * (T1m - T1n);
  const T1D = KP866025403 * (T11 - TW);
  const T12 = TW + T11;
  const T1l = TR - KP500000000 * T12;
  const T1F = T1m + T1n;
  const T1G = T1E - KP500000000 * T1F;

  const Tf = br[6], Th = bi[6], Te = Wc[6], Tg = Ws[6];
  const Ti = Te * Tf + Tg * Th;
  const T1S = Te * Th - Tg * Tf;

  const Tp = br[2], Tr = bi[2], To = Wc[2], Tq = Ws[2];
  const Ts = To * Tp + Tq * Tr;
  const T1c = To * Tr - Tq * Tp;

  const Tk = br[10], Tm = bi[10], Tj = Wc[10], Tl = Ws[10];
  const Tn = Tj * Tk + Tl * Tm;
  const T1b = Tj * Tm - Tl * Tk;

  const T1d = KP866025403 * (T1b - T1c);
  const T24 = KP866025403 * (Ts - Tn);
  const Tt = Tn + Ts;
  const T1a = Ti - KP500000000 * Tt;
  const T1T = T1b + T1c;
  const T25 = T1S - KP500000000 * T1T;

  const Tx = br[3], Tz = bi[3], Tw = Wc[3], Ty = Ws[3];
  const TA = Tw * Tx + Ty * Tz;
  const T1z = Tw * Tz - Ty * Tx;

  const TH = br[11], TJ = bi[11], TG = Wc[11], TI = Ws[11];
  const TK = TG * TH + TI * TJ;
  const T1i = TG * TJ - TI * TH;

  const TC = br[7], TE = bi[7], TB = Wc[7], TD = Ws[7];
  const TF = TB * TC + TD * TE;
  const T1h = TB * TE - TD * TC;

  const T1j = KP866025403 * (T1h - T1i);
  const T1y = KP866025403 * (TK - TF);
  const TL = TF + TK;
  const T1g = TA - KP500000000 * TL;
  const T1A = T1h + T1i;
  const T1B = T1z - KP500000000 * T1A;

  const outR = new Float64Array(12), outI = new Float64Array(12);

  const Td = T1 + Tc;
  const Tu = Ti + Tt;
  const Tv = Td + Tu;
  const T1N = Td - Tu;
  const T1U = T1S + T1T;
  const T1X = T1V + T1W;
  const T1Y = T1U + T1X;
  const T20 = T1X - T1U;

  const TM = TA + TL;
  const T13 = TR + T12;
  const T14 = TM + T13;
  const T1Z = TM - T13;
  const T1O = T1z + T1A;
  const T1P = T1E + T1F;
  const T1Q = T1O - T1P;
  const T1R = T1O + T1P;

  outR[6] = Tv - T14;
  outI[6] = T1Y - T1R;
  outR[0] = Tv + T14;
  outI[0] = T1R + T1Y;
  outR[3] = T1N - T1Q;
  outI[3] = T1Z + T20;
  outR[9] = T1N + T1Q;
  outI[9] = T20 - T1Z;

  const T1r = T15 + T18;
  const T1s = T1a + T1d;
  const T1t = T1r + T1s;
  const T1x = T1r - T1s;
  const T23 = T21 + T22;
  const T26 = T24 + T25;
  const T27 = T23 - T26;
  const T2a = T26 + T23;

  const T1u = T1g + T1j;
  const T1v = T1l + T1o;
  const T1w = T1u + T1v;
  const T28 = T1u - T1v;
  const T1C = T1y + T1B;
  const T1H = T1D + T1G;
  const T1I = T1C - T1H;
  const T29 = T1C + T1H;

  outR[10] = T1t - T1w;
  outI[10] = T2a - T29;
  outR[4] = T1t + T1w;
  outI[4] = T29 + T2a;
  outR[7] = T1x - T1I;
  outI[7] = T28 + T27;
  outR[1] = T1x + T1I;
  outI[1] = T27 - T28;

  const T19 = T15 - T18;
  const T1e = T1a - T1d;
  const T1f = T19 + T1e;
  const T1J = T19 - T1e;
  const T2b = T25 - T24;
  const T2c = T22 - T21;
  const T2d = T2b + T2c;
  const T2f = T2c - T2b;

  const T1k = T1g - T1j;
  const T1p = T1l - T1o;
  const T1q = T1k + T1p;
  const T2g = T1k - T1p;
  const T1K = T1B - T1y;
  const T1L = T1G - T1D;
  const T1M = T1K - T1L;
  const T2e = T1K + T1L;

  outR[2] = T1f - T1q;
  outI[2] = T2d - T2e;
  outR[8] = T1f + T1q;
  outI[8] = T2e + T2d;
  outR[11] = T1J - T1M;
  outI[11] = T2g + T2f;
  outR[5] = T1J + T1M;
  outI[5] = T2f - T2g;

  return [outR, outI];
}

module.exports = { t1_12 };
