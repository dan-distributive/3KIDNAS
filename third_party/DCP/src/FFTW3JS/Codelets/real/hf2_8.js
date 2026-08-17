'use strict';

// =============================================================================
// hf2_8.js -- faithful JS port of rdft/scalar/r2cf/hf2_8.c (non-FMA branch).
// Alternate-codegen ("twiddle-log3/precompute-twiddles") EVEN-radix (r=8)
// twiddle codelet -- same role as hf_8.js (an outer hc2hc-direct combine),
// different rounding. twinstr trig-generates only W^1, W^3, W^7 (three raw
// pairs, matching hf2_5.js's "only some powers stored, rest derived via
// complex products" pattern) -- W^2, W^4, W^6 are derived here via complex
// products of those three, exactly matching the C source's operation order.
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hf2_5.js -- only
// indices 1, 3, 7 are actually read from Wc/Ws.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP707106781 = 0.707106781186547524400844362104849039284835938;

function hf2_8(cr, ci, Wc, Ws) {
  const T2 = Wc[1], T5 = Ws[1], T3 = Wc[3], T6 = Ws[3];
  const T4 = T2 * T3, Tb = T5 * T3, T7 = T5 * T6, Ta = T2 * T6;
  const T8 = T4 - T7, Tc = Ta + Tb, Tg = T4 + T7, Ti = Ta - Tb;
  const Tl = Wc[7], Tm = Ws[7];
  const Tn = T2 * Tl + T5 * Tm; // FMA
  const Tz = Tg * Tm - Ti * Tl; // FNMS
  const Tp = T2 * Tm - T5 * Tl; // FNMS
  const Tx = Tg * Tl + Ti * Tm; // FMA

  const T1 = cr[0], T1c = ci[0];
  const T9 = cr[4], Td = ci[4];
  const Te = T8 * T9 + Tc * Td; // FMA
  const T1b = T8 * Td - Tc * T9; // FNMS
  const Tf = T1 + Te;
  const T1j = T1c - T1b;
  const TL = T1 - Te;
  const T1d = T1b + T1c;

  const TD = cr[7], TE = ci[7];
  const TF = Tl * TD + Tm * TE; // FMA
  const TW = Tl * TE - Tm * TD; // FNMS
  const TG = cr[3], TH = ci[3];
  const TI = T3 * TG + T6 * TH; // FMA
  const TX = T3 * TH - T6 * TG; // FNMS
  const TJ = TF + TI;
  const T16 = TW + TX;
  const TV = TF - TI;
  const TY = TW - TX;

  const Th = cr[2], Tj = ci[2];
  const Tk = Tg * Th + Ti * Tj; // FMA
  const TM = Tg * Tj - Ti * Th; // FNMS
  const To = cr[6], Tq = ci[6];
  const Tr = Tn * To + Tp * Tq; // FMA
  const TN = Tn * Tq - Tp * To; // FNMS
  const Ts = Tk + Tr;
  const T1i = Tk - Tr;
  const TO = TM - TN;
  const T1a = TM + TN;

  const Tu = cr[1], Tv = ci[1];
  const Tw = T2 * Tu + T5 * Tv; // FMA
  const TR = T2 * Tv - T5 * Tu; // FNMS
  const Ty = cr[5], TA = ci[5];
  const TB = Tx * Ty + Tz * TA; // FMA
  const TS = Tx * TA - Tz * Ty; // FNMS
  const TC = Tw + TB;
  const T17 = TR + TS;
  const TQ = Tw - TB;
  const TT = TR - TS;

  const outCr = new Float64Array(8), outCi = new Float64Array(8);

  const Tt = Tf + Ts;
  const TK = TC + TJ;
  outCi[3] = Tt - TK;
  outCr[0] = Tt + TK;
  const T1f = TJ - TC;
  const T1g = T1d - T1a;
  outCr[6] = T1f - T1g;
  outCi[5] = T1f + T1g;

  const T11 = TL - TO;
  const T1m = T1j - T1i;
  const T12 = TQ - TT;
  const T13 = TV + TY;
  const T14 = KP707106781 * (T12 + T13);
  const T1l = KP707106781 * (T13 - T12);
  outCr[3] = T11 - T14;
  outCi[6] = T1l + T1m;
  outCi[0] = T11 + T14;
  outCr[5] = T1l - T1m;

  const T19 = T17 + T16;
  const T1e = T1a + T1d;
  outCr[4] = T19 - T1e;
  outCi[7] = T19 + T1e;
  const T15 = Tf - Ts;
  const T18 = T16 - T17;
  outCr[2] = T15 - T18;
  outCi[1] = T15 + T18;

  const TP = TL + TO;
  const T1k = T1i + T1j;
  const TU = TQ + TT;
  const TZ = TV - TY;
  const T10 = KP707106781 * (TU + TZ);
  const T1h = KP707106781 * (TZ - TU);
  outCi[2] = TP - T10;
  outCi[4] = T1h + T1k;
  outCr[1] = TP + T10;
  outCr[7] = T1h - T1k;

  return [outCr, outCi];
}

module.exports = { hf2_8 };
