'use strict';

// =============================================================================
// r2cb_20.js -- faithful JS port of rdft/scalar/r2cb/r2cb_20.c (non-FMA
// branch). O[0..19] packed halfcomplex -> x[0..19] real, UNNORMALIZED.
// Output uses the R0/R1 (stride-2) convention -- see r2cb_6.js's header.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP1_118033988 = 1.118033988749894848204586834365638117720309180;
const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP1_902113032 = 1.902113032590307144232878666758764286811397268;
const KP1_175570504 = 1.175570504584946258337411909278145537195304875;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;

function r2cb_20(O) {
  const T4 = O[5];
  const T5 = KP2_000000000 * T4;
  const Tr = O[15];
  const Ts = KP2_000000000 * Tr;
  const T1 = O[0];
  const T2 = O[10];
  const T3 = T1 + T2;
  const Tq = T1 - T2;
  const T6 = T3 - T5;
  const TF = Tq - Ts;
  const Tm = T3 + T5;
  const Tt = Tq + Ts;

  const T7 = O[4];
  const T8 = O[6];
  const T9 = T7 + T8;
  const Tu = T7 - T8;
  const TM = O[16];
  const TN = O[14];
  const TO = TM - TN;
  const T1b = TM + TN;

  const Ta = O[9];
  const Tb = O[1];
  const Tc = Ta + Tb;
  const T1a = Ta - Tb;
  const Tv = O[11];
  const Tw = O[19];
  const Tx = Tv + Tw;
  const TP = Tv - Tw;

  const Te = O[8];
  const Tf = O[2];
  const Tg = Te + Tf;
  const Tz = Te - Tf;
  const TR = O[12];
  const TS = O[18];
  const TT = TR - TS;
  const T1e = TR + TS;

  const Th = O[7];
  const Ti = O[3];
  const Tj = Th + Ti;
  const T1d = Th - Ti;
  const TA = O[13];
  const TB = O[17];
  const TC = TA + TB;
  const TU = TB - TA;

  const TQ = TO - TP;
  const T1n = T1e - T1d;
  const T1f = T1d + T1e;
  const T12 = TP + TO;
  const T1m = T1b - T1a;
  const TV = TT - TU;
  const T13 = TU + TT;
  const T1c = T1a + T1b;
  const Td = T9 - Tc;
  const Tk = Tg - Tj;
  const Tl = Td + Tk;
  const Ty = Tu + Tx;
  const TD = Tz - TC;
  const TE = Ty + TD;
  const Tn = T9 + Tc;
  const To = Tg + Tj;
  const Tp = Tn + To;
  const TG = Tu - Tx;
  const TH = Tz + TC;
  const TI = TG + TH;

  const x = new Float64Array(20);
  x[10] = KP2_000000000 * Tl + T6; // R0[WS(rs,5)]
  x[15] = KP2_000000000 * TE + Tt; // R1[WS(rs,7)]
  x[5] = KP2_000000000 * TI + TF; // R1[WS(rs,2)]
  x[0] = KP2_000000000 * Tp + Tm; // R0[0]

  {
    const TW = KP1_175570504 * TQ - KP1_902113032 * TV; // FNMS
    const TY = KP1_902113032 * TQ + KP1_175570504 * TV; // FMA
    const TJ = T6 - KP500000000 * Tl; // FNMS
    const TK = KP1_118033988 * (Td - Tk);
    const TL = TJ - TK;
    const TX = TK + TJ;
    x[2] = TL - TW; // R0[WS(rs,1)]
    x[14] = TX + TY; // R0[WS(rs,7)]
    x[18] = TL + TW; // R0[WS(rs,9)]
    x[6] = TX - TY; // R0[WS(rs,3)]
  }
  {
    const T1g = KP1_175570504 * T1c - KP1_902113032 * T1f; // FNMS
    const T1i = KP1_902113032 * T1c + KP1_175570504 * T1f; // FMA
    const T17 = TF - KP500000000 * TI; // FNMS
    const T18 = KP1_118033988 * (TG - TH);
    const T19 = T17 - T18;
    const T1h = T18 + T17;
    x[17] = T19 - T1g; // R1[WS(rs,8)]
    x[9] = T1h + T1i; // R1[WS(rs,4)]
    x[13] = T19 + T1g; // R1[WS(rs,6)]
    x[1] = T1h - T1i; // R1[0]
  }
  {
    const T1o = KP1_175570504 * T1m - KP1_902113032 * T1n; // FNMS
    const T1q = KP1_902113032 * T1m + KP1_175570504 * T1n; // FMA
    const T1j = Tt - KP500000000 * TE; // FNMS
    const T1k = KP1_118033988 * (Ty - TD);
    const T1l = T1j - T1k;
    const T1p = T1k + T1j;
    x[7] = T1l - T1o; // R1[WS(rs,3)]
    x[19] = T1p + T1q; // R1[WS(rs,9)]
    x[3] = T1l + T1o; // R1[WS(rs,1)]
    x[11] = T1p - T1q; // R1[WS(rs,5)]
  }
  {
    const T14 = KP1_175570504 * T12 - KP1_902113032 * T13; // FNMS
    const T16 = KP1_902113032 * T12 + KP1_175570504 * T13; // FMA
    const TZ = Tm - KP500000000 * Tp; // FNMS
    const T10 = KP1_118033988 * (Tn - To);
    const T11 = TZ - T10;
    const T15 = T10 + TZ;
    x[12] = T11 - T14; // R0[WS(rs,6)]
    x[4] = T15 + T16; // R0[WS(rs,2)]
    x[8] = T11 + T14; // R0[WS(rs,4)]
    x[16] = T15 - T16; // R0[WS(rs,8)]
  }

  return x;
}

module.exports = { r2cb_20 };
