'use strict';

// =============================================================================
// t1_7.js -- faithful JS port of dft/scalar/codelets/t1_7.c (non-FMA branch),
// FFTW3's radix-7 Cooley-Tukey twiddle/combine codelet (DIT).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP222520933 = 0.222520933956314404288902564496794759466355569;
const KP900968867 = 0.900968867902419126236102319507445051165919162;
const KP623489801 = 0.623489801858733530525004884004239810632274731;
const KP433883739 = 0.433883739117558120475768332848358754609990728;
const KP781831482 = 0.781831482468029808708444526674057750232334519;
const KP974927912 = 0.974927912181823607018131682993931217232785801;

function t1_7(br, bi, Wc, Ws) {
  const T1 = br[0], TR = bi[0];

  const T6 = Wc[1] * br[1] + Ws[1] * bi[1];
  const TA = Wc[1] * bi[1] - Ws[1] * br[1];
  const Tb = Wc[6] * br[6] + Ws[6] * bi[6];
  const TB = Wc[6] * bi[6] - Ws[6] * br[6];
  const Tc = T6 + Tb, TS = Tb - T6, TC = TA - TB, TO = TA + TB;

  const Th = Wc[2] * br[2] + Ws[2] * bi[2];
  const TG = Wc[2] * bi[2] - Ws[2] * br[2];
  const Tm = Wc[5] * br[5] + Ws[5] * bi[5];
  const TH = Wc[5] * bi[5] - Ws[5] * br[5];
  const Tn = Th + Tm, TT = Tm - Th, TI = TG - TH, TP = TG + TH;

  const Ts = Wc[3] * br[3] + Ws[3] * bi[3];
  const TD = Wc[3] * bi[3] - Ws[3] * br[3];
  const Tx = Wc[4] * br[4] + Ws[4] * bi[4];
  const TE = Wc[4] * bi[4] - Ws[4] * br[4];
  const Ty = Ts + Tx, TU = Tx - Ts, TF = TD - TE, TQ = TD + TE;

  const outR = new Float64Array(7), outI = new Float64Array(7);
  outR[0] = T1 + Tc + Tn + Ty;
  outI[0] = TO + TP + TQ + TR;

  const TJ = (KP974927912 * TC - KP781831482 * TF) - KP433883739 * TI;
  const Tz = (KP623489801 * Ty + T1) - (KP900968867 * Tn + KP222520933 * Tc);
  outR[5] = Tz - TJ;
  outR[2] = Tz + TJ;
  const TX = (KP974927912 * TS - KP781831482 * TU) - KP433883739 * TT;
  const TY = (KP623489801 * TQ + TR) - (KP900968867 * TP + KP222520933 * TO);
  outI[2] = TX + TY;
  outI[5] = TY - TX;

  const TL = (KP781831482 * TC + KP974927912 * TI) + KP433883739 * TF;
  const TK = (KP623489801 * Tc + T1) - (KP900968867 * Ty + KP222520933 * Tn);
  outR[6] = TK - TL;
  outR[1] = TK + TL;
  const TV = (KP781831482 * TS + KP974927912 * TT) + KP433883739 * TU;
  const TW = (KP623489801 * TO + TR) - (KP900968867 * TQ + KP222520933 * TP);
  outI[1] = TV + TW;
  outI[6] = TW - TV;

  const TN = (KP433883739 * TC + KP974927912 * TF) - KP781831482 * TI;
  const TM = (KP623489801 * Tn + T1) - (KP222520933 * Ty + KP900968867 * Tc);
  outR[4] = TM - TN;
  outR[3] = TM + TN;
  const TZ = (KP433883739 * TS + KP974927912 * TU) - KP781831482 * TT;
  const T10 = (KP623489801 * TP + TR) - (KP222520933 * TQ + KP900968867 * TO);
  outI[3] = TZ + T10;
  outI[4] = T10 - TZ;

  return [outR, outI];
}

module.exports = { t1_7 };
