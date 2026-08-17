'use strict';

// =============================================================================
// hf_7.js -- faithful JS port of rdft/scalar/r2cf/hf_7.c (non-FMA branch).
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hf_2/3/4/5.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP222520933 = 0.222520933956314404288902564496794759466355569;
const KP900968867 = 0.900968867902419126236102319507445051165919162;
const KP623489801 = 0.623489801858733530525004884004239810632274731;
const KP433883739 = 0.433883739117558120475768332848358754609990728;
const KP974927912 = 0.974927912181823607018131682993931217232785801;
const KP781831482 = 0.781831482468029808708444526674057750232334519;

function hf_7(cr, ci, Wc, Ws) {
  const T1 = cr[0], TT = ci[0];

  const T3 = cr[1], T5 = ci[1];
  const T6 = Wc[1] * T3 + Ws[1] * T5;
  const TA = Wc[1] * T5 - Ws[1] * T3;
  const T8 = cr[6], Ta = ci[6];
  const Tb = Wc[6] * T8 + Ws[6] * Ta;
  const TB = Wc[6] * Ta - Ws[6] * T8;
  const Tc = T6 + Tb, TV = TA + TB, TC = TA - TB, TO = Tb - T6;

  const Te = cr[2], Tg = ci[2];
  const Th = Wc[2] * Te + Ws[2] * Tg;
  const TG = Wc[2] * Tg - Ws[2] * Te;
  const Tj = cr[5], Tl = ci[5];
  const Tm = Wc[5] * Tj + Ws[5] * Tl;
  const TH = Wc[5] * Tl - Ws[5] * Tj;
  const Tn = Th + Tm, TS = TG + TH, TI = TG - TH, TP = Th - Tm;

  const Tp = cr[3], Tr = ci[3];
  const Ts = Wc[3] * Tp + Ws[3] * Tr;
  const TD = Wc[3] * Tr - Ws[3] * Tp;
  const Tu = cr[4], Tw = ci[4];
  const Tx = Wc[4] * Tu + Ws[4] * Tw;
  const TE = Wc[4] * Tw - Ws[4] * Tu;
  const Ty = Ts + Tx, TU = TD + TE, TF = TD - TE, TQ = Tx - Ts;

  const outCr = new Float64Array(7), outCi = new Float64Array(7);
  outCr[0] = T1 + Tc + Tn + Ty;

  const ML = KP781831482 * TC + KP974927912 * TI;
  const TL = ML + KP433883739 * TF;
  const MKa = KP623489801 * Tc + T1, MKb = KP900968867 * Ty + KP222520933 * Tn;
  const TK = MKa - MKb;
  outCi[0] = TK - TL;
  outCr[1] = TK + TL;
  outCi[6] = TV + TS + TU + TT;
  const MZ = KP781831482 * TO + KP433883739 * TQ;
  const TZ = MZ - KP974927912 * TP;
  const M10a = KP623489801 * TV + TT, M10b = KP900968867 * TU + KP222520933 * TS;
  const T10 = M10a - M10b;
  outCr[6] = TZ - T10;
  outCi[5] = TZ + T10;

  const MX = KP974927912 * TO + KP433883739 * TP;
  const TX = MX - KP781831482 * TQ;
  const MYa = KP623489801 * TU + TT, MYb = KP900968867 * TS + KP222520933 * TV;
  const TY = MYa - MYb;
  outCr[5] = TX - TY;
  outCi[4] = TX + TY;
  const MR = KP433883739 * TO + KP781831482 * TP;
  const TR = MR + KP974927912 * TQ;
  const MWa = KP623489801 * TS + TT, MWb = KP222520933 * TU + KP900968867 * TV;
  const TW = MWa - MWb;
  outCr[4] = TR - TW;
  outCi[3] = TR + TW;

  const MN = KP433883739 * TC + KP974927912 * TF;
  const TN = MN - KP781831482 * TI;
  const MMa = KP623489801 * Tn + T1, MMb = KP222520933 * Ty + KP900968867 * Tc;
  const TM = MMa - MMb;
  outCi[2] = TM - TN;
  outCr[3] = TM + TN;
  const MJ = KP974927912 * TC - KP781831482 * TF;
  const TJ = MJ - KP433883739 * TI;
  const Mza = KP623489801 * Ty + T1, Mzb = KP900968867 * Tn + KP222520933 * Tc;
  const Tz = Mza - Mzb;
  outCi[1] = Tz - TJ;
  outCr[2] = Tz + TJ;

  return [outCr, outCi];
}

module.exports = { hf_7 };
