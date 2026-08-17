'use strict';

// =============================================================================
// n1_3.js -- faithful JS port of dft/scalar/codelets/n1_3.c (non-FMA branch),
// FFTW3's direct (base-case) radix-3 complex DFT codelet.
// FMA(a,b,c)=a*b+c, FNMS(a,b,c)=c-a*b macro-expand to plain two-step
// arithmetic here (ARCH_PREFERS_FMA undefined, -ffp-contract=off) --
// confirmed zero fmadd/fnmadd/etc in the compiled object file.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP500000000 = 0.5;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function n1_3(ri, ii) {
  const T1 = ri[0], Ta = ii[0];
  const T2 = ri[1], T3 = ri[2];
  const T4 = T2 + T3;
  const T9 = KP866025403 * (T3 - T2);
  const T6 = ii[1], T7 = ii[2];
  const T8 = KP866025403 * (T6 - T7);
  const Tb = T6 + T7;

  const ro = new Float64Array(3), io = new Float64Array(3);
  ro[0] = T1 + T4;
  io[0] = Ta + Tb;
  const T5 = T1 - KP500000000 * T4;
  ro[2] = T5 - T8;
  ro[1] = T5 + T8;
  const Tc = Ta - KP500000000 * Tb;
  io[1] = T9 + Tc;
  io[2] = Tc - T9;
  return [ro, io];
}

module.exports = { n1_3 };
