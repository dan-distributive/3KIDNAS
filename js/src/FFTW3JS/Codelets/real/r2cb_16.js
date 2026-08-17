'use strict';

// =============================================================================
// r2cb_16.js -- faithful JS port of rdft/scalar/r2cb/r2cb_16.c (non-FMA
// branch). O[0..15] packed halfcomplex -> x[0..15] real, UNNORMALIZED.
// Output uses the R0/R1 (stride-2) convention -- see r2cb_6.js's header.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP1_847759065 = 1.847759065022573512256366378793576573644833252;
const KP765366864 = 0.765366864730179543456919968060797733522689125;
const KP1_414213562 = 1.414213562373095048801688724209698078569671875;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;

function r2cb_16(O) {
  const T7 = O[2];
  const T8 = O[6];
  const TE = T7 - T8;
  const Tj = O[14];
  const Tk = O[10];
  const TF = Tj + Tk;
  const T9 = KP2_000000000 * (T7 + T8);
  const TS = KP1_414213562 * (TE + TF);
  const Tl = KP2_000000000 * (Tj - Tk);
  const TG = KP1_414213562 * (TE - TF);

  const T4 = O[4];
  const T5 = KP2_000000000 * T4;
  const TB = O[12];
  const TC = KP2_000000000 * TB;
  const T1 = O[0];
  const T2 = O[8];
  const T3 = T1 + T2;
  const TA = T1 - T2;

  const T6 = T3 + T5;
  const TR = TA + TC;
  const Ti = T3 - T5;
  const TD = TA - TC;

  const Tb = O[1];
  const Tc = O[7];
  const Td = Tb + Tc;
  const TI = Tb - Tc;
  const To = O[15];
  const Tp = O[9];
  const Tq = To - Tp;
  const TM = To + Tp;

  const Te = O[5];
  const Tf = O[3];
  const Tg = Te + Tf;
  const TL = Te - Tf;
  const Tr = O[11];
  const Ts = O[13];
  const Tt = Tr - Ts;
  const TJ = Tr + Ts;

  const Tn = Td - Tg;
  const Tu = Tq - Tt;
  const TV = TM - TL;
  const TU = TI + TJ;
  const TN = TL + TM;
  const TK = TI - TJ;

  const x = new Float64Array(16);
  const Ta = T6 + T9;
  const Th = KP2_000000000 * (Td + Tg);
  x[8] = Ta - Th; // R0[WS(rs,4)]
  x[0] = Ta + Th; // R0[0]
  const TT = TR - TS;
  const TW = KP765366864 * TU - KP1_847759065 * TV; // FNMS
  x[11] = TT - TW; // R1[WS(rs,5)]
  x[3] = TT + TW; // R1[WS(rs,1)]

  const TX = TR + TS;
  const TY = KP1_847759065 * TU + KP765366864 * TV; // FMA
  x[7] = TX - TY; // R1[WS(rs,3)]
  x[15] = TX + TY; // R1[WS(rs,7)]
  const Tm = Ti - Tl;
  const Tv = KP1_414213562 * (Tn - Tu);
  x[10] = Tm - Tv; // R0[WS(rs,5)]
  x[2] = Tm + Tv; // R0[WS(rs,1)]

  const Tw = Ti + Tl;
  const Tx = KP1_414213562 * (Tn + Tu);
  x[6] = Tw - Tx; // R0[WS(rs,3)]
  x[14] = Tw + Tx; // R0[WS(rs,7)]
  const TH = TD + TG;
  const TO = KP1_847759065 * TK - KP765366864 * TN; // FNMS
  x[9] = TH - TO; // R1[WS(rs,4)]
  x[1] = TH + TO; // R1[0]

  const TP = TD - TG;
  const TQ = KP765366864 * TK + KP1_847759065 * TN; // FMA
  x[5] = TP - TQ; // R1[WS(rs,2)]
  x[13] = TP + TQ; // R1[WS(rs,6)]
  const Ty = T6 - T9;
  const Tz = KP2_000000000 * (Tt + Tq);
  x[4] = Ty - Tz; // R0[WS(rs,2)]
  x[12] = Ty + Tz; // R0[WS(rs,6)]

  return x;
}

module.exports = { r2cb_16 };
