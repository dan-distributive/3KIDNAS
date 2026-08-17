'use strict';

// =============================================================================
// r2cb_13.js -- faithful JS port of rdft/scalar/r2cb/r2cb_13.c (non-FMA
// branch). O[0..12] packed halfcomplex (O[0..6]=Re0..Re6, O[7..12]=Im6..
// Im1, same convention as r2cf_13.js) -> x[0..12] real.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP1_007074065 = 1.007074065727533254493747707736933954186697125;
const KP227708958 = 0.227708958111581597949308691735310621069285120;
const KP531932498 = 0.531932498429674575175042127684371897596660533;
const KP774781170 = 0.774781170935234584261351932853525703557550433;
const KP265966249 = 0.265966249214837287587521063842185948798330267;
const KP516520780 = 0.516520780623489722840901288569017135705033622;
const KP151805972 = 0.151805972074387731966205794490207080712856746;
const KP503537032 = 0.503537032863766627246873853868466977093348562;
const KP166666666 = 0.166666666666666666666666666666666666666666667;
const KP600925212 = 0.600925212577331548853203544578415991041882762;
const KP500000000 = 0.5;
const KP256247671 = 0.256247671582936600958684654061725059144125175;
const KP156891391 = 0.156891391051584611046832726756003269660212636;
const KP348277202 = 0.348277202304271810011321589858529485233929352;
const KP1_150281458 = 1.150281458948006242736771094910906776922003215;
const KP300238635 = 0.300238635966332641462884626667381504676006424;
const KP011599105 = 0.011599105605768290721655456654083252189827041;
const KP1_732050807 = 1.732050807568877293527446341505872366942805254;
const KP2_000000000 = 2.0;

function r2cb_13(O) {
  const Ts = O[12], Tt = O[10], Tu = O[9];
  const Tv = Tt - Tu;
  const Tw = KP2_000000000 * Ts - Tv;
  const TE = KP1_732050807 * (Tt + Tu);
  const TC = O[8], Tx = O[7], Ty = O[11];
  const TB = Tx + Ty;
  const Tz = KP1_732050807 * (Tx - Ty);
  const TD = TB - KP2_000000000 * TC;

  const TA = Tw + Tz;
  const TF = TD - TE;
  const TG = KP011599105 * TA + KP300238635 * TF;
  const TS = KP300238635 * TA - KP011599105 * TF;

  const TP = Ts + Tv;
  const TQ = TB + TC;
  const TR = KP1_150281458 * TP - KP348277202 * TQ;
  const T15 = KP348277202 * TP + KP1_150281458 * TQ;
  const TH = Tw - Tz;
  const TI = TE + TD;
  const TJ = KP156891391 * TH + KP256247671 * TI;
  const TT = KP156891391 * TI - KP256247671 * TH;

  const T1 = O[0];
  const T7 = O[5], T8 = O[2], T9 = O[6];
  const Ta = T8 + T9;
  const Tb = T7 + Ta;
  const Ti = T7 - KP500000000 * Ta;
  const Tf = T8 - T9;

  const T2 = O[1], T3 = O[3], T4 = O[4];
  const T5 = T3 + T4;
  const T6 = T2 + T5;
  const Th = T2 - KP500000000 * T5;
  const Te = T3 - T4;

  const Tm = KP600925212 * (T6 - Tb);
  const Tc = T6 + Tb;
  const Td = T1 - KP166666666 * Tc;
  const Tg = Te + Tf;
  const Tj = Th + Ti;
  const Tk = KP503537032 * Tg + KP151805972 * Tj;
  const Tn = Th - Ti;
  const To = Te - Tf;
  const Tp = KP516520780 * Tn - KP265966249 * To;

  const x = new Float64Array(13);
  x[0] = KP2_000000000 * Tc + T1;

  const TK = KP1_732050807 * (TG + TJ);
  const T1b = KP1_732050807 * (TS - TT);
  const TU = TS + TT;
  const TV = TR - TU;
  const T12 = KP2_000000000 * TU + TR;
  const T14 = TG - TJ;
  const T16 = KP2_000000000 * T14 - T15;
  const T18 = T14 + T15;
  const TM = KP774781170 * To + KP531932498 * Tn;
  const TN = KP227708958 * Tg - KP1_007074065 * Tj;
  const TO = TM - TN;
  const T1a = TM + TN;

  const Tl = Td - Tk;
  const Tq = Tm - Tp;
  const Tr = Tl - Tq;
  const T17 = Tq + Tl;
  const TZ = KP2_000000000 * Tk + Td;
  const T10 = KP2_000000000 * Tp + Tm;
  const T11 = TZ - T10;
  const T13 = T10 + TZ;

  x[5] = T11 - T12;
  x[12] = T13 - T16;
  x[1] = T13 + T16;
  x[8] = T11 + T12;

  const TL = Tr - TK;
  const TW = TO - TV;
  x[7] = TL - TW;
  x[2] = TL + TW;
  const T19 = T17 - T18;
  const T1c = T1a + T1b;
  x[3] = T19 - T1c;
  x[9] = T1c + T19;

  const T1d = T1a - T1b;
  const T1e = T17 + T18;
  x[4] = T1d + T1e;
  x[10] = T1e - T1d;
  const TX = Tr + TK;
  const TY = TO + TV;
  x[6] = TX - TY;
  x[11] = TX + TY;

  return x;
}

module.exports = { r2cb_13 };
