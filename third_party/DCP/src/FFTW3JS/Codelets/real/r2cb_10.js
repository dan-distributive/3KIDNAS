'use strict';

// =============================================================================
// r2cb_10.js -- faithful JS port of rdft/scalar/r2cb/r2cb_10.c (non-FMA
// branch). O[0..9] packed halfcomplex -> x[0..9] real, UNNORMALIZED. Output
// uses the R0/R1 (stride-2) convention -- see r2cb_6.js's header.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP1_902113032 = 1.902113032590307144232878666758764286811397268;
const KP1_175570504 = 1.175570504584946258337411909278145537195304875;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;
const KP1_118033988 = 1.118033988749894848204586834365638117720309180;

function r2cb_10(O) {
  const T1 = O[0];
  const T2 = O[5];
  const T3 = T1 - T2;
  const Tb = T1 + T2;
  const Tl = O[6];
  const Tm = O[9];
  const Tn = Tl - Tm;
  const Tv = Tl + Tm;
  const Ti = O[8];
  const Tj = O[7];
  const Tk = Ti - Tj;
  const Tu = Ti + Tj;
  const T4 = O[2];
  const T5 = O[3];
  const T6 = T4 - T5;
  const Tc = T4 + T5;
  const T7 = O[4];
  const T8 = O[1];
  const T9 = T7 - T8;
  const Td = T7 + T8;
  const Ta = T6 + T9;
  const Ts = KP1_118033988 * (T6 - T9);
  const Te = Tc + Td;
  const Tg = KP1_118033988 * (Tc - Td);

  const x = new Float64Array(10);
  x[5] = T3 + KP2_000000000 * Ta; // R1[WS(rs,2)]
  x[0] = Tb + KP2_000000000 * Te; // R0[0]

  const To = KP1_175570504 * Tk - KP1_902113032 * Tn;
  const Tq = KP1_175570504 * Tn + KP1_902113032 * Tk;
  const Tf = Tb - KP500000000 * Te;
  const Th = Tf - Tg;
  const Tp = Tg + Tf;
  x[2] = Th - To; // R0[WS(rs,1)]
  x[4] = Tp + Tq; // R0[WS(rs,2)]
  x[8] = Th + To; // R0[WS(rs,4)]
  x[6] = Tp - Tq; // R0[WS(rs,3)]

  const Tw = KP1_175570504 * Tu - KP1_902113032 * Tv;
  const Ty = KP1_175570504 * Tv + KP1_902113032 * Tu;
  const Tr = T3 - KP500000000 * Ta;
  const Tt = Tr - Ts;
  const Tx = Ts + Tr;
  x[7] = Tt - Tw; // R1[WS(rs,3)]
  x[9] = Tx + Ty; // R1[WS(rs,4)]
  x[3] = Tt + Tw; // R1[WS(rs,1)]
  x[1] = Tx - Ty; // R1[0]

  return x;
}

module.exports = { r2cb_10 };
