'use strict';

// =============================================================================
// hb_7.js -- faithful JS port of rdft/scalar/r2cb/hb_7.c (non-FMA branch).
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hb_2/3/4/5.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP222520933 = 0.222520933956314404288902564496794759466355569;
const KP900968867 = 0.900968867902419126236102319507445051165919162;
const KP623489801 = 0.623489801858733530525004884004239810632274731;
const KP781831482 = 0.781831482468029808708444526674057750232334519;
const KP974927912 = 0.974927912181823607018131682993931217232785801;
const KP433883739 = 0.433883739117558120475768332848358754609990728;

function hb_7(cr, ci, Wc, Ws) {
  const T1 = cr[0];
  const T2 = cr[1], T3 = ci[0];
  const T4 = T2 + T3, Tu = T2 - T3;
  const T5 = cr[2], T6 = ci[1];
  const T7 = T5 + T6, Tw = T5 - T6;
  const T8 = cr[3], T9 = ci[2];
  const Ta = T8 + T9, Tv = T8 - T9;

  const Mx = KP433883739 * Tu + KP974927912 * Tv;
  const Tx = Mx - KP781831482 * Tw;
  const MI = KP781831482 * Tu + KP974927912 * Tw;
  const TI = MI + KP433883739 * Tv;
  const MV = KP974927912 * Tu - KP781831482 * Tv;
  const TV = MV - KP433883739 * Tw;
  const MQa = KP623489801 * Ta + T1, MQb = KP900968867 * T7 + KP222520933 * T4;
  const TQ = MQa - MQb;
  const MEa = KP623489801 * T4 + T1, MEb = KP900968867 * Ta + KP222520933 * T7;
  const TE = MEa - MEb;
  const Mma = KP623489801 * T7 + T1, Mmb = KP222520933 * Ta + KP900968867 * T4;
  const Tm = Mma - Mmb;

  const Tb = ci[6];
  const Tc = ci[5], Td = cr[6];
  const Te = Tc - Td, Tp = Tc + Td;
  const Tf = ci[4], Tg = cr[5];
  const Th = Tf - Tg, Tn = Tf + Tg;
  const Ti = ci[3], Tj = cr[4];
  const Tk = Ti - Tj, To = Ti + Tj;

  const Mq = KP781831482 * Tn - KP974927912 * To;
  const Tq = Mq - KP433883739 * Tp;
  const MF = KP781831482 * Tp + KP974927912 * Tn;
  const TF = MF + KP433883739 * To;
  const MR = KP433883739 * Tn + KP781831482 * To;
  const TR = MR - KP974927912 * Tp;
  const MUa = KP623489801 * Tk + Tb, MUb = KP900968867 * Th + KP222520933 * Te;
  const TU = MUa - MUb;
  const MJa = KP623489801 * Te + Tb, MJb = KP900968867 * Tk + KP222520933 * Th;
  const TJ = MJa - MJb;
  const Mta = KP623489801 * Th + Tb, Mtb = KP222520933 * Tk + KP900968867 * Te;
  const Tt = Mta - Mtb;

  const outCr = new Float64Array(7), outCi = new Float64Array(7);
  outCr[0] = T1 + T4 + T7 + Ta;
  outCi[0] = Tb + Te + Th + Tk;

  const Tr = Tm - Tq, Ty = Tt - Tx;
  const Tl = Wc[4], Ts = Ws[4];
  outCr[4] = Tl * Tr - Ts * Ty;
  outCi[4] = Tl * Ty + Ts * Tr;

  const TY = TQ + TR, T10 = TV + TU;
  const TX = Wc[2], TZ = Ws[2];
  outCr[2] = TX * TY - TZ * T10;
  outCi[2] = TX * T10 + TZ * TY;

  const TA = Tm + Tq, TC = Tx + Tt;
  const Tz = Wc[3], TB = Ws[3];
  outCr[3] = Tz * TA - TB * TC;
  outCi[3] = Tz * TC + TB * TA;

  const TM = TE + TF, TO = TJ - TI;
  const TL = Wc[6], TN = Ws[6];
  outCr[6] = TL * TM - TN * TO;
  outCi[6] = TL * TO + TN * TM;

  const TS = TQ - TR, TW = TU - TV;
  const TP = Wc[5], TT = Ws[5];
  outCr[5] = TP * TS - TT * TW;
  outCi[5] = TP * TW + TT * TS;

  const TG = TE - TF, TK = TI + TJ;
  const TD = Wc[1], TH = Ws[1];
  outCr[1] = TD * TG - TH * TK;
  outCi[1] = TD * TK + TH * TG;

  return [outCr, outCi];
}

module.exports = { hb_7 };
