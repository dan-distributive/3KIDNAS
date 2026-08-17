'use strict';

// =============================================================================
// n1_2.js -- faithful JS port of dft/scalar/codelets/n1_2.c (non-FMA branch),
// FFTW3's direct (base-case, no recursion) radix-2 complex DFT codelet.
// Ports the literal arithmetic; no FMA fusion needed since n1_2 has zero
// multiplications (desc = {4, "n1_2", {4,0,0,0}, ...}).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

function n1_2(ri, ii) {
  const T1 = ri[0], T2 = ri[1];
  const T3 = ii[0], T4 = ii[1];
  const ro = new Float64Array(2), io = new Float64Array(2);
  ro[1] = T1 - T2;
  ro[0] = T1 + T2;
  io[1] = T3 - T4;
  io[0] = T3 + T4;
  return [ro, io];
}

module.exports = { n1_2 };
