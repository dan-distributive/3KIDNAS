'use strict';

// =============================================================================
// n1_4.js -- faithful JS port of dft/scalar/codelets/n1_4.c (non-FMA branch),
// FFTW3's direct (base-case) radix-4 complex DFT codelet. Zero multiplications
// (desc = {16, "n1_4", {16,0,0,0}, ...}), so no FMA concerns at all.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

function n1_4(ri, ii) {
  const T1 = ri[0], T2 = ri[2];
  const T3 = T1 + T2;
  const Tb = T1 - T2;
  const T7 = ii[0], T8 = ii[2];
  const T9 = T7 - T8;
  const Tf = T7 + T8;

  const T4 = ri[1], T5 = ri[3];
  const T6 = T4 + T5;
  const Ta = T4 - T5;
  const Tc = ii[1], Td = ii[3];
  const Te = Tc - Td;
  const Tg = Tc + Td;

  const ro = new Float64Array(4), io = new Float64Array(4);
  ro[2] = T3 - T6;
  io[2] = Tf - Tg;
  ro[0] = T3 + T6;
  io[0] = Tf + Tg;
  io[1] = T9 - Ta;
  ro[1] = Tb + Te;
  io[3] = Ta + T9;
  ro[3] = Tb - Te;
  return [ro, io];
}

module.exports = { n1_4 };
