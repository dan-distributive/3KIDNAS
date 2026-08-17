'use strict';

// =============================================================================
// t1_2.js -- faithful JS port of dft/scalar/codelets/t1_2.c (non-FMA branch),
// FFTW3's radix-2 Cooley-Tukey twiddle/combine codelet (decimation-in-time).
// br/bi[0] is the untwiddled DC block; br/bi[1] is phase 1, twiddled by
// Wc[1]/Ws[1] = cos/sin(2*pi*col*1/n) before combining -- same convention as
// CompositeSolver1D.js's radix3Twiddle.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

function t1_2(br, bi, Wc, Ws) {
  const T1 = br[0], T8 = bi[0];
  const T3 = br[1], T5 = bi[1];
  const T2 = Wc[1], T4 = Ws[1];
  const T6 = T2 * T3 + T4 * T5;
  const T7 = T2 * T5 - T4 * T3;

  const outR = new Float64Array(2), outI = new Float64Array(2);
  outR[1] = T1 - T6;
  outI[1] = T8 - T7;
  outR[0] = T1 + T6;
  outI[0] = T7 + T8;
  return [outR, outI];
}

module.exports = { t1_2 };
