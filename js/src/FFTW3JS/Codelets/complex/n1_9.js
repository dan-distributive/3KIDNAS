'use strict';

// =============================================================================
// n1_9.js -- faithful JS port of dft/scalar/codelets/n1_9.c (non-FMA branch),
// FFTW3's direct (base-case) radix-9 complex DFT codelet.
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

function n1_9(ri, ii) {
  const T1 = ri[0], T2 = ri[3], T3 = ri[6];
  const T4 = T2 + T3;
  const T5 = T1 + T4;
  const TO = KP866025403 * (T3 - T2);
  const Th = T1 - KP500000000 * T4;

  const TP = ii[0], Ti = ii[3], Tj = ii[6];
  const TQ = Ti + Tj;
  const Tk = KP866025403 * (Ti - Tj);
  const T1g = TP + TQ;
  const TR = TP - KP500000000 * TQ;

  const T6 = ri[1], Ts = ii[1];
  const T7 = ri[4], T8 = ri[7];
  const T9 = T7 + T8;
  const Tr = KP866025403 * (T8 - T7);
  const Tn = ii[4], To = ii[7];
  const Tp = KP866025403 * (Tn - To);
  const Tt = Tn + To;
  const Ta = T6 + T9;
  const T1c = Ts + Tt;
  const Tm = T6 - KP500000000 * T9;
  const Tq = Tm + Tp;
  const TW = Tm - Tp;
  const Tu = Ts - KP500000000 * Tt;
  const Tv = Tr + Tu;
  const TX = Tu - Tr;

  const Tb = ri[2], TD = ii[2];
  const Tc = ri[5], Td = ri[8];
  const Te = Tc + Td;
  const TC = KP866025403 * (Td - Tc);
  const Ty = ii[5], Tz = ii[8];
  const TA = KP866025403 * (Ty - Tz);
  const TE = Ty + Tz;
  const Tf = Tb + Te;
  const T1d = TD + TE;
  const Tx = Tb - KP500000000 * Te;
  const TB = Tx + TA;
  const T10 = Tx - TA;
  const TF = TD - KP500000000 * TE;
  const TG = TC + TF;
  const TZ = TF - TC;

  const ro = new Float64Array(9), io = new Float64Array(9);

  const T1e = KP866025403 * (T1c - T1d);
  const Tg = Ta + Tf;
  const T1b = T5 - KP500000000 * Tg;
  ro[0] = T5 + Tg;
  ro[3] = T1b + T1e;
  ro[6] = T1b - T1e;
  const T1f = KP866025403 * (Tf - Ta);
  const T1h = T1c + T1d;
  const T1i = T1g - KP500000000 * T1h;
  io[3] = T1f + T1i;
  io[0] = T1g + T1h;
  io[6] = T1i - T1f;

  const Tl = Th + Tk;
  const TS = TO + TR;
  const Tw = KP766044443 * Tq + KP642787609 * Tv;
  const TH = KP173648177 * TB + KP984807753 * TG;
  const TI = Tw + TH;
  const TN = KP866025403 * (TH - Tw);
  const TK = KP766044443 * Tv - KP642787609 * Tq;
  const TL = KP173648177 * TG - KP984807753 * TB;
  const TM = KP866025403 * (TK - TL);
  const TT = TK + TL;
  ro[1] = Tl + TI;
  io[1] = TS + TT;
  const TJ = Tl - KP500000000 * TI;
  ro[7] = TJ - TM;
  ro[4] = TJ + TM;
  const TU = TS - KP500000000 * TT;
  io[4] = TN + TU;
  io[7] = TU - TN;

  const TV = Th - Tk;
  const T14 = TR - TO;
  const TY = KP173648177 * TW + KP984807753 * TX;
  const T11 = KP342020143 * TZ - KP939692620 * T10;
  const T12 = TY + T11;
  const T13 = KP866025403 * (T11 - TY);
  const T15 = KP173648177 * TX - KP984807753 * TW;
  const T16 = KP342020143 * T10 + KP939692620 * TZ;
  const T17 = T15 - T16;
  const T1a = KP866025403 * (T15 + T16);
  ro[2] = TV + T12;
  io[2] = T14 + T17;
  const T18 = T14 - KP500000000 * T17;
  io[5] = T13 + T18;
  io[8] = T18 - T13;
  const T19 = TV - KP500000000 * T12;
  ro[8] = T19 - T1a;
  ro[5] = T19 + T1a;

  return [ro, io];
}

module.exports = { n1_9 };
