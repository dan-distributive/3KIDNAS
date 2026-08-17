'use strict';

// =============================================================================
// r2cbIII_16.js -- faithful JS port of rdft/scalar/r2cb/r2cbIII_16.c
// (non-FMA branch). HC2RIII ("shifted") radix-16 direct codelet -- the
// BACKWARD counterpart of r2cfII_16.js (see that file's header). INPUT
// in[0..15] uses r2cfII_16.js's OUTPUT convention (in[k]=Cr[k] for k=0..7,
// in[15-k]=Ci[k] for k=0..7); OUTPUT out[0..15] (phase p) uses
// r2cfII_16.js's INPUT convention (out[2k]=R0[k], out[2k+1]=R1[k]).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP1_961570560 = 1.961570560806460898252364472268478073947867462;
const KP390180644 = 0.390180644032256535696569736954044481855383236;
const KP1_111140466 = 1.111140466039204449485661627897065748749874382;
const KP1_662939224 = 1.662939224605090474157576755235811513477121624;
const KP707106781 = 0.707106781186547524400844362104849039284835938;
const KP1_414213562 = 1.414213562373095048801688724209698078569671875;
const KP765366864 = 0.765366864730179543456919968060797733522689125;
const KP1_847759065 = 1.847759065022573512256366378793576573644833252;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;

function r2cbIII_16(inArr) {
  const T1 = inArr[0];
  const T2 = inArr[7];
  const T3 = T1 + T2;
  const Tf = T1 - T2;
  const TA = inArr[15];
  const TB = inArr[8];
  const TC = TA + TB;
  const TV = TB - TA;

  const T4 = inArr[4];
  const T5 = inArr[3];
  const T6 = T4 + T5;
  const Tz = T4 - T5;
  const Tg = inArr[11];
  const Th = inArr[12];
  const Ti = Tg + Th;
  const TU = Tg - Th;

  const T7 = T3 + T6;
  const TW = TU + TV;
  const T13 = TV - TU;
  const Tj = Tf - Ti;
  const TD = Tz + TC;
  const TK = Tz - TC;
  const TP = T3 - T6;
  const TH = Tf + Ti;

  const T8 = inArr[2];
  const T9 = inArr[5];
  const Ta = T8 + T9;
  const Tk = T8 - T9;
  const Tl = inArr[13];
  const Tm = inArr[10];
  const Tn = Tl + Tm;
  const TR = Tl - Tm;

  const Tb = inArr[1];
  const Tc = inArr[6];
  const Td = Tb + Tc;
  const Tp = Tb - Tc;
  const Tq = inArr[14];
  const Tr = inArr[9];
  const Ts = Tq + Tr;
  const TQ = Tr - Tq;

  const Te = Ta + Td;
  const TX = Ta - Td;
  const T12 = TR + TQ;
  const To = Tk - Tn;
  const Tt = Tp - Ts;
  const Tx = Tp + Ts;
  const TS = TQ - TR;
  const Tw = Tk + Tn;

  const out = new Float64Array(16);
  out[0] = KP2_000000000 * (T7 + Te);
  out[8] = KP2_000000000 * (T13 - T12);
  const TT = TP + TS;
  const TY = TW - TX;
  out[2] = KP1_847759065 * TT + KP765366864 * TY; // FMA
  out[10] = KP1_847759065 * TY - KP765366864 * TT; // FNMS

  const T11 = T7 - Te;
  const T14 = T12 + T13;
  out[4] = KP1_414213562 * (T11 + T14);
  out[12] = KP1_414213562 * (T14 - T11);
  const TZ = TP - TS;
  const T10 = TX + TW;
  out[6] = KP765366864 * TZ + KP1_847759065 * T10; // FMA
  out[14] = KP765366864 * T10 - KP1_847759065 * TZ; // FNMS

  const TI = KP707106781 * (Tw + Tx);
  const TJ = TH - TI;
  const TN = TH + TI;
  const TL = KP707106781 * (To - Tt);
  const TM = TK - TL;
  const TO = TL + TK;
  out[3] = KP1_662939224 * TJ + KP1_111140466 * TM; // FMA
  out[15] = KP390180644 * TO - KP1_961570560 * TN; // FNMS
  out[11] = KP1_662939224 * TM - KP1_111140466 * TJ; // FNMS
  out[7] = KP390180644 * TN + KP1_961570560 * TO; // FMA

  const Tu = KP707106781 * (To + Tt);
  const Tv = Tj + Tu;
  const TF = Tj - Tu;
  const Ty = KP707106781 * (Tw - Tx);
  const TE = Ty + TD;
  const TG = Ty - TD;
  out[1] = KP1_961570560 * Tv - KP390180644 * TE; // FNMS
  out[13] = KP1_111140466 * TG - KP1_662939224 * TF; // FNMS
  out[9] = -(KP390180644 * Tv + KP1_961570560 * TE); // -FMA
  out[5] = KP1_111140466 * TF + KP1_662939224 * TG; // FMA

  return out;
}

module.exports = { r2cbIII_16 };
