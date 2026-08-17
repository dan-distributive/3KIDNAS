'use strict';

// =============================================================================
// t1_6.js -- faithful JS port of dft/scalar/codelets/t1_6.c (non-FMA branch),
// FFTW3's radix-6 Cooley-Tukey twiddle/combine codelet (DIT).
// Moved here from RaderSolver.js -- see n1_6.js's header for why.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP500000000 = 0.5;

function t1_6(br, bi, Wc, Ws) {
  const T1 = br[0], TN = bi[0];

  const T3 = br[3], T5 = bi[3];
  const T6 = Wc[3] * T3 + Ws[3] * T5;
  const TM = Wc[3] * T5 - Ws[3] * T3;
  const T7 = T1 - T6, TS = TN - TM, Tv = T1 + T6, TO = TM + TN;

  const Tk = br[4], Tm = bi[4];
  const Tn = Wc[4] * Tk + Ws[4] * Tm;
  const TD = Wc[4] * Tm - Ws[4] * Tk;
  const Tp = br[1], Tr = bi[1];
  const Ts = Wc[1] * Tp + Ws[1] * Tr;
  const TE = Wc[1] * Tr - Ws[1] * Tp;
  const Tt = Tn - Ts, TJ = TD + TE, Tx = Tn + Ts, TF = TD - TE;

  const T9 = br[2], Tb = bi[2];
  const Tc = Wc[2] * T9 + Ws[2] * Tb;
  const TA = Wc[2] * Tb - Ws[2] * T9;
  const Te = br[5], Tg = bi[5];
  const Th = Wc[5] * Te + Ws[5] * Tg;
  const TB = Wc[5] * Tg - Ws[5] * Te;
  const Ti = Tc - Th, TI = TA + TB, Tw = Tc + Th, TC = TA - TB;

  const outR = new Float64Array(6), outI = new Float64Array(6);

  const TG = KP866025403 * (TC - TF);
  const Tu = Ti + Tt;
  const Tz = T7 - KP500000000 * Tu;
  outR[3] = T7 + Tu;
  outR[1] = Tz + TG;
  outR[5] = Tz - TG;
  const TR = KP866025403 * (Tt - Ti);
  const TT = TC + TF;
  const TU = TS - KP500000000 * TT;
  outI[1] = TR + TU;
  outI[3] = TT + TS;
  outI[5] = TU - TR;

  const TK = KP866025403 * (TI - TJ);
  const Ty = Tw + Tx;
  const TH = Tv - KP500000000 * Ty;
  outR[0] = Tv + Ty;
  outR[4] = TH + TK;
  outR[2] = TH - TK;
  const TQ = KP866025403 * (Tx - Tw);
  const TL = TI + TJ;
  const TP = TO - KP500000000 * TL;
  outI[0] = TL + TO;
  outI[4] = TQ + TP;
  outI[2] = TP - TQ;

  return [outR, outI];
}

module.exports = { t1_6 };
