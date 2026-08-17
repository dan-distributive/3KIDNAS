'use strict';

// =============================================================================
// r2cb_12.js -- faithful JS port of rdft/scalar/r2cb/r2cb_12.c (non-FMA
// branch). O[0..11] packed halfcomplex -> x[0..11] real, UNNORMALIZED.
// Output uses the R0/R1 (stride-2) convention -- see r2cb_6.js's header.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP1_732050807 = 1.732050807568877293527446341505872366942805254;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;

function r2cb_12(O) {
  const T8 = O[3];
  const T9 = O[5];
  const Ta = O[1];
  const Tb = T9 + Ta;
  const Tm = KP2_000000000 * T8 - Tb;
  const TA = KP1_732050807 * (T9 - Ta);
  const Tw = O[9];
  const Tn = O[7];
  const To = O[11];
  const Tx = Tn + To;
  const Tp = KP1_732050807 * (Tn - To);
  const TB = KP2_000000000 * Tw + Tx;

  const Te = O[8];
  const Tf = KP1_732050807 * Te;
  const T1 = O[0];
  const T2 = O[4];
  const Td = T1 - T2;
  const T3 = KP2_000000000 * T2 + T1;
  const Tr = Td - Tf;
  const Tg = Td + Tf;

  const Ti = O[10];
  const Tj = KP1_732050807 * Ti;
  const T4 = O[6];
  const T5 = O[2];
  const Th = T4 - T5;
  const T6 = KP2_000000000 * T5 + T4;
  const Ts = Th + Tj;
  const Tk = Th - Tj;

  const x = new Float64Array(12);
  const T7 = T3 + T6;
  const Tc = KP2_000000000 * (T8 + Tb);
  x[6] = T7 - Tc; // R0[WS(rs,3)]
  x[0] = T7 + Tc; // R0[0]

  const Tl = Tg + Tk;
  const Tq = Tm - Tp;
  x[2] = Tl - Tq; // R0[WS(rs,1)]
  x[8] = Tl + Tq; // R0[WS(rs,4)]
  const TD = Tg - Tk;
  const TE = TB - TA;
  x[5] = TD - TE; // R1[WS(rs,2)]
  x[11] = TD + TE; // R1[WS(rs,5)]

  const Tz = Tr - Ts;
  const TC = TA + TB;
  x[1] = Tz - TC; // R1[0]
  x[7] = Tz + TC; // R1[WS(rs,3)]

  const Tv = T3 - T6;
  const Ty = KP2_000000000 * (Tw - Tx);
  x[9] = Tv - Ty; // R1[WS(rs,4)]
  x[3] = Tv + Ty; // R1[WS(rs,1)]
  const Tt = Tr + Ts;
  const Tu = Tm + Tp;
  x[10] = Tt - Tu; // R0[WS(rs,5)]
  x[4] = Tt + Tu; // R0[WS(rs,2)]

  return x;
}

module.exports = { r2cb_12 };
