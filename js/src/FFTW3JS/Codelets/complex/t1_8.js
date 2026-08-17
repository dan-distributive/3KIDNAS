'use strict';

// =============================================================================
// t1_8.js -- faithful JS port of dft/scalar/codelets/t1_8.c (non-FMA branch),
// FFTW3's radix-8 Cooley-Tukey twiddle/combine codelet (DIT).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP707106781 = 0.707106781186547524400844362104849039284835938;

function t1_8(br, bi, Wc, Ws) {
  const T1 = br[0], T18 = bi[0];
  const T6 = Wc[4] * br[4] + Ws[4] * bi[4];
  const T17 = Wc[4] * bi[4] - Ws[4] * br[4];
  const T7 = T1 + T6, T1e = T18 - T17, TH = T1 - T6, T19 = T17 + T18;

  const Tz = Wc[7] * br[7] + Ws[7] * bi[7];
  const TS = Wc[7] * bi[7] - Ws[7] * br[7];
  const TE = Wc[3] * br[3] + Ws[3] * bi[3];
  const TT = Wc[3] * bi[3] - Ws[3] * br[3];
  const TF = Tz + TE, T13 = TS + TT, TR = Tz - TE, TU = TS - TT;

  const Tc = Wc[2] * br[2] + Ws[2] * bi[2];
  const TI = Wc[2] * bi[2] - Ws[2] * br[2];
  const Th = Wc[6] * br[6] + Ws[6] * bi[6];
  const TJ = Wc[6] * bi[6] - Ws[6] * br[6];
  const Ti = Tc + Th, T1f = Tc - Th, TK = TI - TJ, T16 = TI + TJ;

  const To = Wc[1] * br[1] + Ws[1] * bi[1];
  const TN = Wc[1] * bi[1] - Ws[1] * br[1];
  const Tt = Wc[5] * br[5] + Ws[5] * bi[5];
  const TO = Wc[5] * bi[5] - Ws[5] * br[5];
  const Tu = To + Tt, T12 = TN + TO, TM = To - Tt, TP = TN - TO;

  const outR = new Float64Array(8), outI = new Float64Array(8);

  const Tj = T7 + Ti, TG = Tu + TF;
  outR[4] = Tj - TG;
  outR[0] = Tj + TG;

  const T15 = T12 + T13, T1a = T16 + T19;
  outI[0] = T15 + T1a;
  outI[4] = T1a - T15;
  const T11 = T7 - Ti, T14 = T12 - T13;
  outR[6] = T11 - T14;
  outR[2] = T11 + T14;

  const T1b = TF - Tu, T1c = T19 - T16;
  outI[2] = T1b + T1c;
  outI[6] = T1c - T1b;

  const TX = TH - TK, T1g = T1e - T1f, TY = TP - TM, TZ = TR + TU;
  const T10 = KP707106781 * (TY - TZ);
  const T1d = KP707106781 * (TY + TZ);
  outR[7] = TX - T10;
  outI[5] = T1g - T1d;
  outR[3] = TX + T10;
  outI[1] = T1d + T1g;

  const TL = TH + TK, T1i = T1f + T1e, TQ = TM + TP, TV = TR - TU;
  const TW = KP707106781 * (TQ + TV);
  const T1h = KP707106781 * (TV - TQ);
  outR[5] = TL - TW;
  outI[7] = T1i - T1h;
  outR[1] = TL + TW;
  outI[3] = T1h + T1i;

  return [outR, outI];
}

module.exports = { t1_8 };
