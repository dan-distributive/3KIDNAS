'use strict';

// =============================================================================
// r2cb_9.js -- faithful JS port of rdft/scalar/r2cb/r2cb_9.c (non-FMA branch).
// O[0..8] packed halfcomplex (O[0..4]=Re0..Re4, O[5..8]=Im4..Im1, same
// convention as r2cf_9.js) -> x[0..8] real.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP984807753 = 0.984807753012208059366743024589523013670643252;
const KP173648177 = 0.173648177666930348851716626769314796000375677;
const KP300767466 = 0.300767466360870593278543795225003852144476517;
const KP1_705737063 = 1.705737063904886419256501927880148143872040591;
const KP642787609 = 0.642787609686539326322643409907263432907559884;
const KP766044443 = 0.766044443118978035202392650555416673935832457;
const KP1_326827896 = 1.326827896337876792410842639271782594433726619;
const KP1_113340798 = 1.113340798452838732905825904094046265936583811;
const KP500000000 = 0.5;
const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP2_000000000 = 2.0;
const KP1_732050807 = 1.732050807568877293527446341505872366942805254;

function r2cb_9(O) {
  const Ta = O[6];
  const Tb = KP1_732050807 * Ta;
  const T1 = O[0], T2 = O[3];
  const T9 = T1 - T2;
  const T3 = KP2_000000000 * T2 + T1;
  const Tq = T9 + Tb, Tc = T9 - Tb;

  const T4 = O[1], Tk = O[8];
  const T5 = O[4], T6 = O[2];
  const T7 = T5 + T6;
  const Ti = KP866025403 * (T5 - T6);
  const Te = O[5], Tf = O[7];
  const Tg = KP866025403 * (Te + Tf);
  const Tj = Tf - Te;

  const T8 = T4 + T7;
  const Tl = KP500000000 * Tj + Tk;
  const Tm = Ti + Tl, Ts = Tl - Ti;
  const Td = T4 - KP500000000 * T7;
  const Th = Td - Tg, Tr = Td + Tg;

  const x = new Float64Array(9);
  x[0] = KP2_000000000 * T8 + T3;
  const Tw = T3 - T8;
  const Tx = KP1_732050807 * (Tk - Tj);
  x[3] = Tw - Tx;
  x[6] = Tw + Tx;

  const Tp = KP1_113340798 * Th + KP1_326827896 * Tm;
  const Tn = KP766044443 * Th - KP642787609 * Tm;
  const To = Tc - Tn;
  x[1] = KP2_000000000 * Tn + Tc;
  x[7] = To + Tp;
  x[4] = To - Tp;
  const Tv = KP1_705737063 * Tr + KP300767466 * Ts;
  const Tt = KP173648177 * Tr - KP984807753 * Ts;
  const Tu = Tq - Tt;
  x[2] = KP2_000000000 * Tt + Tq;
  x[8] = Tu + Tv;
  x[5] = Tu - Tv;

  return x;
}

module.exports = { r2cb_9 };
