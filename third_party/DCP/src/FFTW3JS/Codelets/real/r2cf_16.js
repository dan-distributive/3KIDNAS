'use strict';

// =============================================================================
// r2cf_16.js -- faithful JS port of rdft/scalar/r2cf/r2cf_16.c (non-FMA
// branch). x[0..15] (real) -> O[0..15] packed halfcomplex: O[0..8]=Re0..Re8,
// O[9..15]=Im7..Im1 (O[n-k]=Im_k). R0/R1 (stride-2) input convention -- see
// r2cf_6.js's header.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function r2cf_16(x) {
  const T1 = x[0];
  const T2 = x[8];
  const T3 = T1 + T2;
  const T4 = x[4];
  const T5 = x[12];
  const T6 = T4 + T5;
  const T7 = T3 + T6;
  const Tz = T1 - T2;
  const Ti = T4 - T5;

  const T8 = x[2];
  const T9 = x[10];
  const Ta = T8 + T9;
  const Tg = T8 - T9;
  const Tb = x[14];
  const Tc = x[6];
  const Td = Tb + Tc;
  const Tf = Tb - Tc;

  const Te = Ta + Td;
  const TA = KP707106781 * (Tg + Tf);
  const Th = KP707106781 * (Tf - Tg);

  const Tk = x[15];
  const Tl = x[7];
  const Tm = Tk - Tl;
  const TN = Tk + Tl;
  const Tn = x[3];
  const To = x[11];
  const Tp = Tn - To;
  const TO = Tn + To;

  const Tq = KP382683432 * Tm - KP923879532 * Tp;
  const TV = TN + TO;
  const TF = KP923879532 * Tm + KP382683432 * Tp;
  const TP = TN - TO;

  const Tr = x[1];
  const Ts = x[9];
  const Tt = Tr - Ts;
  const TK = Tr + Ts;
  const Tu = x[5];
  const Tv = x[13];
  const Tw = Tu - Tv;
  const TL = Tu + Tv;

  const Tx = KP382683432 * Tt + KP923879532 * Tw;
  const TU = TK + TL;
  const TE = KP923879532 * Tt - KP382683432 * Tw;
  const TM = TK - TL;

  const O = new Float64Array(16);
  O[4] = T7 - Te;
  O[12] = TV - TU;

  const Tj = Th - Ti;
  const Ty = Tq - Tx;
  O[15] = Tj + Ty;
  O[9] = Ty - Tj;
  const TD = Tz + TA;
  const TG = TE + TF;
  O[7] = TD - TG;
  O[1] = TD + TG;

  const TB = Tz - TA;
  const TC = Tx + Tq;
  O[5] = TB - TC;
  O[3] = TB + TC;
  const TH = Ti + Th;
  const TI = TF - TE;
  O[13] = TH + TI;
  O[11] = TI - TH;

  const TJ = T3 - T6;
  const TQ = KP707106781 * (TM + TP);
  O[6] = TJ - TQ;
  O[2] = TJ + TQ;

  const TR = Td - Ta;
  const TS = KP707106781 * (TP - TM);
  O[14] = TR + TS;
  O[10] = TS - TR;
  const TT = T7 + Te;
  const TW = TU + TV;
  O[8] = TT - TW;
  O[0] = TT + TW;

  return O;
}

module.exports = { r2cf_16 };
