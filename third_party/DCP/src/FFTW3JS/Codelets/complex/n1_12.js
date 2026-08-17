'use strict';

// =============================================================================
// n1_12.js -- faithful JS port of dft/scalar/codelets/n1_12.c (non-FMA
// branch), FFTW3's direct (base-case) radix-12 complex DFT codelet.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP500000000 = 0.5;

function n1_12(ri, ii) {
  const T1 = ri[0], T2 = ri[4], T3 = ri[8];
  const T4 = T2 + T3;
  const T5 = T1 + T4;
  const TR = T1 - KP500000000 * T4;
  const TA = KP866025403 * (T3 - T2);

  const To = ii[0], Tp = ii[4], Tq = ii[8];
  const Tr = Tp + Tq;
  const Ts = To + Tr;
  const TS = KP866025403 * (Tp - Tq);
  const Tz = To - KP500000000 * Tr;

  const T6 = ri[6], T7 = ri[10], T8 = ri[2];
  const T9 = T7 + T8;
  const Ta = T6 + T9;
  const TU = T6 - KP500000000 * T9;
  const TD = KP866025403 * (T8 - T7);

  const Tt = ii[6], Tu = ii[10], Tv = ii[2];
  const Tw = Tu + Tv;
  const Tx = Tt + Tw;
  const TV = KP866025403 * (Tu - Tv);
  const TC = Tt - KP500000000 * Tw;

  const Tc = ri[3], Td = ri[7], Te = ri[11];
  const Tf = Td + Te;
  const Tg = Tc + Tf;
  const T1a = KP866025403 * (Te - Td);
  const TG = Tc - KP500000000 * Tf;

  const T1b = ii[3], TH = ii[7], TI = ii[11];
  const T1c = TH + TI;
  const TJ = KP866025403 * (TH - TI);
  const T1u = T1b + T1c;
  const T1d = T1b - KP500000000 * T1c;

  const Th = ri[9], Ti = ri[1], Tj = ri[5];
  const Tk = Ti + Tj;
  const Tl = Th + Tk;
  const T1f = KP866025403 * (Tj - Ti);
  const TL = Th - KP500000000 * Tk;

  const T1g = ii[9], TM = ii[1], TN = ii[5];
  const T1h = TM + TN;
  const TO = KP866025403 * (TM - TN);
  const T1v = T1g + T1h;
  const T1i = T1g - KP500000000 * T1h;

  const ro = new Float64Array(12), io = new Float64Array(12);

  const Tb = T5 + Ta, Tm = Tg + Tl;
  ro[6] = Tb - Tm;
  ro[0] = Tb + Tm;

  const T1x = Ts + Tx, T1y = T1u + T1v;
  io[6] = T1x - T1y;
  io[0] = T1x + T1y;
  const Tn = Tg - Tl, Ty = Ts - Tx;
  io[3] = Tn + Ty;
  io[9] = Ty - Tn;

  const T1t = T5 - Ta, T1w = T1u - T1v;
  ro[3] = T1t - T1w;
  ro[9] = T1t + T1w;

  const TZ = TA + Tz, T10 = TD + TC;
  const T11 = TZ - T10, T1l = TZ + T10;
  const T1e = T1a + T1d, T1j = T1f + T1i;
  const T1k = T1e - T1j, T1m = T1e + T1j;

  const T12 = TG + TJ, T13 = TL + TO;
  const T14 = T12 - T13, T18 = T12 + T13;
  const T15 = TR + TS, T16 = TU + TV;
  const T17 = T15 + T16, T19 = T15 - T16;

  io[1] = T11 - T14;
  ro[1] = T19 + T1k;
  io[7] = T11 + T14;
  ro[7] = T19 - T1k;
  ro[10] = T17 - T18;
  io[10] = T1l - T1m;
  ro[4] = T17 + T18;
  io[4] = T1l + T1m;

  const TB = Tz - TA, TE = TC - TD;
  const TF = TB - TE, T1r = TB + TE;
  const T1o = T1d - T1a, T1p = T1i - T1f;
  const T1q = T1o - T1p, T1s = T1o + T1p;

  const TK = TG - TJ, TP = TL - TO;
  const TQ = TK - TP, TY = TK + TP;
  const TT = TR - TS, TW = TU - TV;
  const TX = TT + TW, T1n = TT - TW;

  io[5] = TF - TQ;
  ro[5] = T1n + T1q;
  io[11] = TF + TQ;
  ro[11] = T1n - T1q;
  ro[2] = TX - TY;
  io[2] = T1r - T1s;
  ro[8] = TX + TY;
  io[8] = T1r + T1s;

  return [ro, io];
}

module.exports = { n1_12 };
