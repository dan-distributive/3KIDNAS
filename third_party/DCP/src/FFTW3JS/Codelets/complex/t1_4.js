'use strict';

// =============================================================================
// t1_4.js -- faithful JS port of dft/scalar/codelets/t1_4.c (non-FMA branch),
// FFTW3's radix-4 Cooley-Tukey twiddle/combine codelet (decimation-in-time).
// br/bi[0] is the untwiddled DC block; phases 1,2,3 are twiddled by
// Wc[p]/Ws[p] = cos/sin(2*pi*col*p/n) before combining.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

function t1_4(br, bi, Wc, Ws) {
  const T1 = br[0], Tp = bi[0];

  const T6 = Wc[2] * br[2] + Ws[2] * bi[2];
  const To = Wc[2] * bi[2] - Ws[2] * br[2];

  const Tc = Wc[1] * br[1] + Ws[1] * bi[1];
  const Tk = Wc[1] * bi[1] - Ws[1] * br[1];

  const Th = Wc[3] * br[3] + Ws[3] * bi[3];
  const Tl = Wc[3] * bi[3] - Ws[3] * br[3];

  const outR = new Float64Array(4), outI = new Float64Array(4);

  const T7 = T1 + T6;
  const Ti = Tc + Th;
  outR[2] = T7 - Ti;
  outR[0] = T7 + Ti;
  const Tn = Tk + Tl;
  const Tq = To + Tp;
  outI[0] = Tn + Tq;
  outI[2] = Tq - Tn;

  const Tj = T1 - T6;
  const Tm = Tk - Tl;
  outR[3] = Tj - Tm;
  outR[1] = Tj + Tm;
  const Tr = Tp - To;
  const Ts = Tc - Th;
  outI[1] = Tr - Ts;
  outI[3] = Ts + Tr;

  return [outR, outI];
}

module.exports = { t1_4 };
