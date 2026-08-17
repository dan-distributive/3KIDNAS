'use strict';

// =============================================================================
// r2cb_8.js -- faithful JS port of rdft/scalar/r2cb/r2cb_8.c (non-FMA branch).
// O[0..7] packed halfcomplex -> x[0..7] real, UNNORMALIZED. Output uses the
// R0/R1 (stride-2) convention -- see r2cb_6.js's header.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP1_414213562 = 1.414213562373095048801688724209698078569671875;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;

function r2cb_8(O) {
  const T4 = O[2];
  const T5 = KP2_000000000 * T4;
  const Tf = O[6];
  const Tg = KP2_000000000 * Tf;
  const T1 = O[0];
  const T2 = O[4];
  const T3 = T1 + T2;
  const Te = T1 - T2;
  const T7 = O[1];
  const T8 = O[3];
  const T9 = KP2_000000000 * (T7 + T8);
  const Ti = T7 - T8;
  const Tb = O[7];
  const Tc = O[5];
  const Td = KP2_000000000 * (Tb - Tc);
  const Tj = Tb + Tc;

  const x = new Float64Array(8);
  const T6 = T3 + T5;
  x[4] = T6 - T9; // R0[WS(rs,2)]
  x[0] = T6 + T9; // R0[0]
  const Ta = T3 - T5;
  x[2] = Ta - Td; // R0[WS(rs,1)]
  x[6] = Ta + Td; // R0[WS(rs,3)]
  const Th = Te - Tg;
  const Tk = KP1_414213562 * (Ti - Tj);
  x[5] = Th - Tk; // R1[WS(rs,2)]
  x[1] = Th + Tk; // R1[0]
  const Tl = Te + Tg;
  const Tm = KP1_414213562 * (Ti + Tj);
  x[3] = Tl - Tm; // R1[WS(rs,1)]
  x[7] = Tl + Tm; // R1[WS(rs,3)]

  return x;
}

module.exports = { r2cb_8 };
