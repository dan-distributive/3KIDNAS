'use strict';

// =============================================================================
// hb2_8.js -- faithful JS port of rdft/scalar/r2cb/hb2_8.c (non-FMA branch).
// Alternate-codegen ("twiddle-log3/precompute-twiddles") EVEN-radix (r=8)
// backward twiddle codelet -- same role as hb_8.js, different rounding
// (see hf2_8.js's header for the forward-direction analogue). twinstr
// trig-generates only W^1, W^3, W^7; W^2, W^4, W^6 are derived. Same
// (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hb2_5.js -- note the
// cross-indexed cr[k+1]/ci[k] pairing pattern (an established property of
// this codelet family's packed convention, transcribed literally from the
// C source's WS(rs,k) indices, not re-derived).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP707106781 = 0.707106781186547524400844362104849039284835938;

function hb2_8(cr, ci, Wc, Ws) {
  const Tf = Wc[1], Ti = Ws[1], Tg = Wc[3], Tj = Ws[3];
  const Th = Tf * Tg, To = Ti * Tg, Tk = Ti * Tj, Tn = Tf * Tj;
  const Tl = Th - Tk, Tp = Tn + To, TP = Th + Tk, TR = Tn - To;
  const TF = Wc[7], TG = Ws[7];
  const TH = Tf * TF + Ti * TG; // FMA
  const T15 = TP * TG - TR * TF; // FNMS
  const TL = Tf * TG - Ti * TF; // FNMS
  const TT = TP * TF + TR * TG; // FMA

  const T1 = cr[0];
  const T2 = ci[3];
  const T3 = T1 + T2;
  const TU = T1 - T2;
  const Tt = ci[5];
  const Tu = cr[6];
  const Tv = Tt - Tu;
  const TV = Tt + Tu;

  const T4 = cr[2];
  const T5 = ci[1];
  const T6 = T4 + T5;
  const T16 = T4 - T5;
  const Tq = ci[7];
  const Tr = cr[4];
  const Ts = Tq - Tr;
  const T17 = Tq + Tr;

  const T7 = T3 + T6;
  const T1f = TU + TV;
  const T1i = T17 - T16;
  const Tw = Ts + Tv;
  const TI = T3 - T6;
  const TW = TU - TV;
  const T18 = T16 + T17;
  const TM = Ts - Tv;

  const T8 = cr[1];
  const T9 = ci[2];
  const Ta = T8 + T9;
  const TX = T8 - T9;
  const TA = ci[4];
  const TB = cr[7];
  const TC = TA - TB;
  const T11 = TA + TB;

  const Tb = ci[0];
  const Tc = cr[3];
  const Td = Tb + Tc;
  const T10 = Tb - Tc;
  const Tx = ci[6];
  const Ty = cr[5];
  const Tz = Tx - Ty;
  const TY = Tx + Ty;

  const Te = Ta + Td;
  const T19 = TX + TY;
  const T1a = T10 + T11;
  const TD = Tz + TC;
  const TJ = TC - Tz;
  const TZ = TX - TY;
  const T12 = T10 - T11;
  const TN = Ta - Td;

  const outCr = new Float64Array(8), outCi = new Float64Array(8);
  outCr[0] = T7 + Te;
  outCi[0] = Tw + TD;
  const Tm = T7 - Te;
  const TE = Tw - TD;
  outCr[4] = Tl * Tm - Tp * TE; // FNMS
  outCi[4] = Tp * Tm + Tl * TE; // FMA

  const TQ = TI + TJ;
  const TS = TN + TM;
  outCr[2] = TP * TQ - TR * TS; // FNMS
  outCi[2] = TP * TS + TR * TQ; // FMA
  const TK = TI - TJ;
  const TO = TM - TN;
  outCr[6] = TH * TK - TL * TO; // FNMS
  outCi[6] = TH * TO + TL * TK; // FMA

  const T1g = KP707106781 * (T19 + T1a);
  const T1h = T1f - T1g;
  const T1l = T1f + T1g;
  const T1j = KP707106781 * (TZ - T12);
  const T1k = T1i + T1j;
  const T1m = T1i - T1j;
  outCr[3] = Tg * T1h - Tj * T1k; // FNMS
  outCi[3] = Tj * T1h + Tg * T1k; // FMA
  outCr[7] = TF * T1l - TG * T1m; // FNMS
  outCi[7] = TG * T1l + TF * T1m; // FMA

  const T13 = KP707106781 * (TZ + T12);
  const T14 = TW - T13;
  const T1d = TW + T13;
  const T1b = KP707106781 * (T19 - T1a);
  const T1c = T18 - T1b;
  const T1e = T18 + T1b;
  outCr[5] = TT * T14 - T15 * T1c; // FNMS
  outCi[5] = T15 * T14 + TT * T1c; // FMA
  outCr[1] = Tf * T1d - Ti * T1e; // FNMS
  outCi[1] = Ti * T1d + Tf * T1e; // FMA

  return [outCr, outCi];
}

module.exports = { hb2_8 };
