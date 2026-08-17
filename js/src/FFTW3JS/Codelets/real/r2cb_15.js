'use strict';

// =============================================================================
// r2cb_15.js -- faithful JS port of rdft/scalar/r2cb/r2cb_15.c (non-FMA
// branch). O[0..14] packed halfcomplex (O[0..7]=Re0..Re7, O[8..14]=Im7..
// Im1, same convention as r2cf_15.js) -> x[0..14] real.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP1_118033988 = 1.118033988749894848204586834365638117720309180;
const KP1_902113032 = 1.902113032590307144232878666758764286811397268;
const KP1_175570504 = 1.175570504584946258337411909278145537195304875;
const KP500000000 = 0.5;
const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP2_000000000 = 2.0;
const KP1_732050807 = 1.732050807568877293527446341505872366942805254;

function r2cb_15(O) {
  const Tg = O[10];
  const Th = KP1_732050807 * Tg;
  const T1 = O[0], T2 = O[5];
  const Tf = T1 - T2;
  const T3 = KP2_000000000 * T2 + T1;
  const Tu = Tf - Th, Ti = Tf + Th;

  const T4 = O[3], TD = O[12], T9 = O[6], TI = O[9];
  const T5 = O[7], T6 = O[2];
  const T7 = T5 + T6;
  const Ta = O[4], Tb = O[1];
  const Tc = Ta + Tb;

  const Tp = O[11], Tq = O[14];
  const Tr = KP866025403 * (Tp + Tq);
  const TH = Tp - Tq;
  const Tk = O[8], Tl = O[13];
  const Tm = KP866025403 * (Tk - Tl);
  const TC = Tk + Tl;

  const TB = KP866025403 * (T5 - T6);
  const TZ = TD - TC;
  const T10 = TI - TH;
  const TE = KP500000000 * TC + TD;
  const TG = KP866025403 * (Ta - Tb);
  const TJ = KP500000000 * TH + TI;
  const Tj = T4 - KP500000000 * T7;
  const Tn = Tj - Tm, Tv = Tj + Tm;
  const To = T9 - KP500000000 * Tc;
  const Ts = To - Tr, Tw = To + Tr;
  const T8 = T4 + T7;
  const Td = T9 + Tc;
  const Te = T8 + Td;

  const x = new Float64Array(15);
  x[0] = KP2_000000000 * Te + T3;

  const T11 = KP1_175570504 * TZ - KP1_902113032 * T10;
  const T13 = KP1_902113032 * TZ + KP1_175570504 * T10;
  const TW = T3 - KP500000000 * Te;
  const TX = KP1_118033988 * (T8 - Td);
  const TY = TW - TX;
  const T12 = TX + TW;
  x[12] = TY - T11;
  x[9] = T12 + T13;
  x[3] = TY + T11;
  x[6] = T12 - T13;

  const TP = KP1_118033988 * (Tn - Ts);
  const Tt = Tn + Ts;
  const TO = Ti - KP500000000 * Tt;
  const TR = TE - TB;
  const TS = TJ - TG;
  const TT = KP1_175570504 * TR - KP1_902113032 * TS;
  const TV = KP1_902113032 * TR + KP1_175570504 * TS;
  x[5] = KP2_000000000 * Tt + Ti;
  const TU = TP + TO;
  x[11] = TU - TV;
  x[14] = TU + TV;
  const TQ = TO - TP;
  x[2] = TQ - TT;
  x[8] = TQ + TT;

  const Tz = KP1_118033988 * (Tv - Tw);
  const Tx = Tv + Tw;
  const Ty = Tu - KP500000000 * Tx;
  const TF = TB + TE;
  const TK = TG + TJ;
  const TL = KP1_175570504 * TF - KP1_902113032 * TK;
  const TN = KP1_902113032 * TF + KP1_175570504 * TK;
  x[10] = KP2_000000000 * Tx + Tu;
  const TM = Tz + Ty;
  x[1] = TM - TN;
  x[4] = TM + TN;
  const TA = Ty - Tz;
  x[7] = TA - TL;
  x[13] = TA + TL;

  return x;
}

module.exports = { r2cb_15 };
