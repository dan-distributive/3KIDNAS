'use strict';

// =============================================================================
// r2cf_12.js -- faithful JS port of rdft/scalar/r2cf/r2cf_12.c (non-FMA
// branch). x[0..11] (real) -> O[0..11] packed halfcomplex: O[0..6]=Re0..Re6,
// O[7..11]=Im5..Im1 (O[n-k]=Im_k). R0/R1 (stride-2) input convention -- see
// r2cf_6.js's header.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP500000000 = 0.500000000000000000000000000000000000000000000;

function r2cf_12(x) {
  const T1 = x[0];
  const T2 = x[4];
  const T3 = x[8];
  const T4 = T2 + T3;
  const T5 = T1 + T4;
  const Tp = T1 - KP500000000 * T4;
  const Tb = T3 - T2;

  const Tj = x[3];
  const Tk = x[7];
  const Tl = x[11];
  const Tm = Tk + Tl;
  const Tn = Tj - KP500000000 * Tm;
  const Ty = Tl - Tk;
  const Tt = Tj + Tm;

  const T6 = x[6];
  const T7 = x[10];
  const T8 = x[2];
  const T9 = T7 + T8;
  const Ta = T6 + T9;
  const Tq = T6 - KP500000000 * T9;
  const Tc = T8 - T7;

  const Te = x[9];
  const Tf = x[1];
  const Tg = x[5];
  const Th = Tf + Tg;
  const Ti = Te - KP500000000 * Th;
  const Tz = Tg - Tf;
  const Tu = Te + Th;

  const O = new Float64Array(12);
  O[3] = T5 - Ta;
  O[9] = Tt - Tu;
  const Td = KP866025403 * (Tb - Tc);
  const To = Ti - Tn;
  O[11] = Td + To;
  O[7] = To - Td;

  const Tx = Tp - Tq;
  const TA = KP866025403 * (Ty - Tz);
  O[5] = Tx - TA;
  O[1] = Tx + TA;
  const Tv = T5 + Ta;
  const Tw = Tt + Tu;
  O[6] = Tv - Tw;
  O[0] = Tv + Tw;

  const Tr = Tp + Tq;
  const Ts = Tn + Ti;
  O[2] = Tr - Ts;
  O[4] = Tr + Ts;
  const TB = Ty + Tz;
  const TC = Tb + Tc;
  O[10] = KP866025403 * (TB - TC);
  O[8] = KP866025403 * (TC + TB);

  return O;
}

module.exports = { r2cf_12 };
