'use strict';

// =============================================================================
// r2cb_32.js -- faithful JS port of rdft/scalar/r2cb/r2cb_32.c (non-FMA
// branch). O[0..31] packed halfcomplex -> x[0..31] real, UNNORMALIZED.
// Output uses the R0/R1 (stride-2) convention -- see r2cb_6.js's header.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP1_662939224 = 1.662939224605090474157576755235811513477121624;
const KP1_111140466 = 1.111140466039204449485661627897065748749874382;
const KP1_961570560 = 1.961570560806460898252364472268478073947867462;
const KP390180644 = 0.390180644032256535696569736954044481855383236;
const KP765366864 = 0.765366864730179543456919968060797733522689125;
const KP1_847759065 = 1.847759065022573512256366378793576573644833252;
const KP707106781 = 0.707106781186547524400844362104849039284835938;
const KP1_414213562 = 1.414213562373095048801688724209698078569671875;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;

function r2cb_32(O) {
  const T7 = O[4];
  const T8 = O[12];
  const T1w = T7 - T8;
  const Tz = O[28];
  const TA = O[20];
  const T1x = Tz + TA;
  const T9 = KP2_000000000 * (T7 + T8);
  const T2c = KP1_414213562 * (T1w + T1x);
  const TB = KP2_000000000 * (Tz - TA);
  const T1y = KP1_414213562 * (T1w - T1x);

  const T4 = O[8];
  const T5 = KP2_000000000 * T4;
  const T1t = O[24];
  const T1u = KP2_000000000 * T1t;
  const T1 = O[0];
  const T2 = O[16];
  const T3 = T1 + T2;
  const T1s = T1 - T2;
  const T6 = T3 + T5;
  const T2b = T1s + T1u;
  const Ty = T3 - T5;
  const T1v = T1s - T1u;

  const Tb = O[2];
  const Tc = O[14];
  const Td = Tb + Tc;
  const T1A = Tb - Tc;
  const TE = O[30];
  const TF = O[18];
  const TG = TE - TF;
  const T1E = TE + TF;

  const Te = O[10];
  const Tf = O[6];
  const Tg = Te + Tf;
  const T1D = Te - Tf;
  const TH = O[22];
  const TI = O[26];
  const TJ = TH - TI;
  const T1B = TH + TI;

  const Th = KP2_000000000 * (Td + Tg);
  const T2e = T1A + T1B;
  const T2f = T1E - T1D;
  const TD = Td - Tg;
  const TK = TG - TJ;
  const T1C = T1A - T1B;
  const T1F = T1D + T1E;
  const T1h = KP2_000000000 * (TJ + TG);

  const Tj = O[1];
  const Tk = O[15];
  const Tl = Tj + Tk;
  const T1I = Tj - Tk;
  const TX = O[31];
  const TY = O[17];
  const TZ = TX - TY;
  const T1X = TX + TY;

  const Tm = O[9];
  const Tn = O[7];
  const To = Tm + Tn;
  const T1W = Tm - Tn;
  const T10 = O[23];
  const T11 = O[25];
  const T12 = T10 - T11;
  const T1J = T10 + T11;

  const Tp = Tl + To;
  const T2i = T1I + T1J;
  const T2m = T1X - T1W;
  const TN = Tl - To;
  const T13 = TZ - T12;
  const T1K = T1I - T1J;
  const T1Y = T1W + T1X;
  const T1k = T12 + TZ;

  const Tq = O[5];
  const Tr = O[11];
  const Ts = Tq + Tr;
  const T1L = Tq - Tr;
  const TR = O[27];
  const TS = O[21];
  const TT = TR - TS;
  const T1M = TR + TS;

  const Tt = O[3];
  const Tu = O[13];
  const Tv = Tt + Tu;
  const T1O = Tt - Tu;
  const TO = O[19];
  const TP = O[29];
  const TQ = TO - TP;
  const T1P = TP + TO;

  const Tw = Ts + Tv;
  const TU = TQ - TT;
  const T1l = TT + TQ;
  const TW = Ts - Tv;

  const T1T = T1L + T1M;
  const T1U = T1O + T1P;
  const T1V = KP707106781 * (T1T - T1U);
  const T2j = KP707106781 * (T1T + T1U);
  const T1N = T1L - T1M;
  const T1Q = T1O - T1P;
  const T1R = KP707106781 * (T1N + T1Q);
  const T2l = KP707106781 * (T1N - T1Q);

  const x = new Float64Array(32);

  {
    const Tx = KP2_000000000 * (Tp + Tw);
    const T1r = KP2_000000000 * (T1l + T1k);
    const Ta = T6 + T9;
    const Ti = Ta + Th;
    const T1q = Ta - Th;
    x[16] = Ti - Tx; // R0[WS(rs,8)]
    x[24] = T1q + T1r; // R0[WS(rs,12)]
    x[0] = Ti + Tx; // R0[0]
    x[8] = T1q - T1r; // R0[WS(rs,4)]
  }
  {
    const T1g = T6 - T9;
    const T1i = T1g - T1h;
    const T1o = T1g + T1h;
    const T1j = Tp - Tw;
    const T1m = T1k - T1l;
    const T1n = KP1_414213562 * (T1j - T1m);
    const T1p = KP1_414213562 * (T1j + T1m);
    x[20] = T1i - T1n; // R0[WS(rs,10)]
    x[28] = T1o + T1p; // R0[WS(rs,14)]
    x[4] = T1i + T1n; // R0[WS(rs,2)]
    x[12] = T1o - T1p; // R0[WS(rs,6)]
  }
  {
    const TC = Ty - TB;
    const TL = KP1_414213562 * (TD - TK);
    const TM = TC + TL;
    const T16 = TC - TL;
    const TV = TN + TU;
    const T14 = TW + T13;
    const T15 = KP1_847759065 * TV - KP765366864 * T14; // FNMS
    const T17 = KP765366864 * TV + KP1_847759065 * T14; // FMA
    x[18] = TM - T15; // R0[WS(rs,9)]
    x[26] = T16 + T17; // R0[WS(rs,13)]
    x[2] = TM + T15; // R0[WS(rs,1)]
    x[10] = T16 - T17; // R0[WS(rs,5)]
  }
  {
    const T2r = T2b + T2c;
    const T2s = KP1_847759065 * T2e + KP765366864 * T2f; // FMA
    const T2t = T2r - T2s;
    const T2x = T2r + T2s;
    const T2u = T2i + T2j;
    const T2v = T2m - T2l;
    const T2w = KP390180644 * T2u - KP1_961570560 * T2v; // FNMS
    const T2y = KP1_961570560 * T2u + KP390180644 * T2v; // FMA
    x[23] = T2t - T2w; // R1[WS(rs,11)]
    x[31] = T2x + T2y; // R1[WS(rs,15)]
    x[7] = T2t + T2w; // R1[WS(rs,3)]
    x[15] = T2x - T2y; // R1[WS(rs,7)]
  }
  {
    const T18 = Ty + TB;
    const T19 = KP1_414213562 * (TD + TK);
    const T1a = T18 - T19;
    const T1e = T18 + T19;
    const T1b = TN - TU;
    const T1c = T13 - TW;
    const T1d = KP765366864 * T1b - KP1_847759065 * T1c; // FNMS
    const T1f = KP1_847759065 * T1b + KP765366864 * T1c; // FMA
    x[22] = T1a - T1d; // R0[WS(rs,11)]
    x[30] = T1e + T1f; // R0[WS(rs,15)]
    x[6] = T1a + T1d; // R0[WS(rs,3)]
    x[14] = T1e - T1f; // R0[WS(rs,7)]
  }
  {
    const T23 = T1v - T1y;
    const T24 = KP765366864 * T1C + KP1_847759065 * T1F; // FMA
    const T25 = T23 - T24;
    const T29 = T23 + T24;
    const T26 = T1K - T1R;
    const T27 = T1Y - T1V;
    const T28 = KP1_111140466 * T26 - KP1_662939224 * T27; // FNMS
    const T2a = KP1_662939224 * T26 + KP1_111140466 * T27; // FMA
    x[21] = T25 - T28; // R1[WS(rs,10)]
    x[29] = T29 + T2a; // R1[WS(rs,14)]
    x[5] = T25 + T28; // R1[WS(rs,2)]
    x[13] = T29 - T2a; // R1[WS(rs,6)]
  }
  {
    const T2d = T2b - T2c;
    const T2g = KP765366864 * T2e - KP1_847759065 * T2f; // FNMS
    const T2h = T2d + T2g;
    const T2p = T2d - T2g;
    const T2k = T2i - T2j;
    const T2n = T2l + T2m;
    const T2o = KP1_662939224 * T2k - KP1_111140466 * T2n; // FNMS
    const T2q = KP1_111140466 * T2k + KP1_662939224 * T2n; // FMA
    x[19] = T2h - T2o; // R1[WS(rs,9)]
    x[27] = T2p + T2q; // R1[WS(rs,13)]
    x[3] = T2h + T2o; // R1[WS(rs,1)]
    x[11] = T2p - T2q; // R1[WS(rs,5)]
  }
  {
    const T1z = T1v + T1y;
    const T1G = KP1_847759065 * T1C - KP765366864 * T1F; // FNMS
    const T1H = T1z + T1G;
    const T21 = T1z - T1G;
    const T1S = T1K + T1R;
    const T1Z = T1V + T1Y;
    const T20 = KP1_961570560 * T1S - KP390180644 * T1Z; // FNMS
    const T22 = KP390180644 * T1S + KP1_961570560 * T1Z; // FMA
    x[17] = T1H - T20; // R1[WS(rs,8)]
    x[25] = T21 + T22; // R1[WS(rs,12)]
    x[1] = T1H + T20; // R1[0]
    x[9] = T21 - T22; // R1[WS(rs,4)]
  }

  return x;
}

module.exports = { r2cb_32 };
