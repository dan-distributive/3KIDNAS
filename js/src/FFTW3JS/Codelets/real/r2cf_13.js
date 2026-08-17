'use strict';

// =============================================================================
// r2cf_13.js -- faithful JS port of rdft/scalar/r2cf/r2cf_13.c (non-FMA
// branch). x[0..12] (real) -> O[0..12] packed halfcomplex: O[0..6]=Re0..
// Re6, O[7..12]=Im6..Im1 (O[n-k]=Im_k, same convention as r2cf_5/7/9/11.js).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP083333333 = 0.083333333333333333333333333333333333333333333;
const KP075902986 = 0.075902986037193865983102897245103540356428373;
const KP251768516 = 0.251768516431883313623436926934233488546674281;
const KP503537032 = 0.503537032863766627246873853868466977093348562;
const KP113854479 = 0.113854479055790798974654345867655310534642560;
const KP265966249 = 0.265966249214837287587521063842185948798330267;
const KP387390585 = 0.387390585467617292130675966426762851778775217;
const KP300462606 = 0.300462606288665774426601772289207995520941381;
const KP132983124 = 0.132983124607418643793760531921092974399165133;
const KP258260390 = 0.258260390311744861420450644284508567852516811;
const KP2_000000000 = 2.0;
const KP1_732050807 = 1.732050807568877293527446341505872366942805254;
const KP300238635 = 0.300238635966332641462884626667381504676006424;
const KP011599105 = 0.011599105605768290721655456654083252189827041;
const KP156891391 = 0.156891391051584611046832726756003269660212636;
const KP256247671 = 0.256247671582936600958684654061725059144125175;
const KP174138601 = 0.174138601152135905005660794929264742616964676;
const KP575140729 = 0.575140729474003121368385547455453388461001608;
const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP500000000 = 0.5;

function r2cf_13(x) {
  const T13 = x[0];

  const Tc = x[8], Td = x[5];
  const Te = Tc - Td, TO = Tc + Td;

  const T6 = x[1], T7 = x[3], T8 = x[9];
  const T9 = T7 + T8;
  const Ta = T6 + T9;
  const Tv = T7 - T8;
  const To = T6 - KP500000000 * T9;

  const T1 = x[12], T2 = x[10], T3 = x[4];
  const T4 = T2 + T3;
  const T5 = T1 + T4;
  const Tw = T2 - T3;
  const Tp = T1 - KP500000000 * T4;

  const Tf = x[11], Tg = x[6];
  const Th = Tf - Tg, Tr = Tf + Tg;
  const Ti = x[7], Tj = x[2];
  const Tk = Ti - Tj, Ts = Ti + Tj;

  const Tl = Th + Tk;
  const TP = Tr + Ts;
  const Tb = T5 - Ta;
  const Tm = Te + Tl;
  const TW = Ta + T5;
  const TX = TO + TP;
  const T14 = TW + TX;

  const TS = Tv + Tw;
  const TT = Th - Tk;
  const TU = TS - TT, T10 = TS + TT;
  const Tx = KP866025403 * (Tv - Tw);
  const Ty = Te - KP500000000 * Tl;
  const Tz = Tx + Ty, TB = Ty - Tx;

  const Tq = To - Tp;
  const Tt = KP866025403 * (Tr - Ts);
  const Tu = Tq - Tt, TC = Tq + Tt;
  const TN = To + Tp;
  const TQ = TO - KP500000000 * TP;
  const TR = TN - TQ, T11 = TN + TQ;

  const O = new Float64Array(13);
  O[0] = T13 + T14;

  const Tn = KP575140729 * Tb - KP174138601 * Tm;
  const TG = KP174138601 * Tb + KP575140729 * Tm;

  const TA = KP256247671 * Tu - KP156891391 * Tz;
  const TD = KP011599105 * TB - KP300238635 * TC;
  const TE = TA + TD;
  const TF = KP1_732050807 * (TD - TA);
  const TH = KP300238635 * TB + KP011599105 * TC;
  const TI = KP256247671 * Tz + KP156891391 * Tu;
  const TJ = TH - TI;
  const TM = KP1_732050807 * (TI + TH);

  O[8] = KP2_000000000 * TE + Tn;
  O[12] = KP2_000000000 * TJ + TG;
  const TK = TG - TJ;
  O[9] = TF - TK;
  O[10] = TF + TK;
  const TL = Tn - TE;
  O[11] = TL - TM;
  O[7] = TL + TM;

  const TV = KP258260390 * TR - KP132983124 * TU;
  const TY = KP300462606 * (TW - TX);
  const TZ = KP2_000000000 * TV + TY;
  const T1b = TY - TV;

  const T17 = KP387390585 * TU + KP265966249 * TR;
  const T18 = KP113854479 * T10 - KP503537032 * T11;
  const T19 = T17 - T18, T1e = T17 + T18;
  const T12 = KP251768516 * T10 + KP075902986 * T11;
  const T15 = T13 - KP083333333 * T14;
  const T16 = KP2_000000000 * T12 + T15;
  const T1a = T15 - T12;

  O[1] = TZ + T16;
  O[5] = T16 - TZ;
  const T1c = T1a - T1b;
  O[2] = T19 + T1c;
  O[6] = T1c - T19;
  const T1d = T1b + T1a;
  O[3] = T1d - T1e;
  O[4] = T1e + T1d;

  return O;
}

module.exports = { r2cf_13 };
