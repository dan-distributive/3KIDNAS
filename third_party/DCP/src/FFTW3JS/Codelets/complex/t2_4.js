'use strict';

// =============================================================================
// t2_4.js -- faithful JS port of dft/scalar/codelets/t2_4.c (non-FMA
// branch), FFTW3's "twiddle-log3 / precompute-twiddles" radix-4 twiddle
// codelet. twinstr only trig-generates W^1, W^3; W^2 is derived via one
// complex multiply: (T6,T8) = conj(W^1)*W^3 = W^2.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

function t2_4(br, bi, Wc, Ws) {
  const T2 = Wc[1], T4 = Ws[1], T3 = Wc[3], T5 = Ws[3];
  const T6 = T2 * T3 + T4 * T5;
  const T8 = T2 * T5 - T4 * T3;

  const T1 = br[0], Tp = bi[0];
  const T7 = br[2], T9 = bi[2];
  const Ta = T6 * T7 + T8 * T9;
  const To = T6 * T9 - T8 * T7;

  const Tc = br[1], Td = bi[1];
  const Te = T2 * Tc + T4 * Td;
  const Tk = T2 * Td - T4 * Tc;
  const Tf = br[3], Tg = bi[3];
  const Th = T3 * Tf + T5 * Tg;
  const Tl = T3 * Tg - T5 * Tf;

  const outR = new Float64Array(4), outI = new Float64Array(4);

  const Tb = T1 + Ta;
  const Ti = Te + Th;
  outR[2] = Tb - Ti;
  outR[0] = Tb + Ti;
  const Tn = Tk + Tl;
  const Tq = To + Tp;
  outI[0] = Tn + Tq;
  outI[2] = Tq - Tn;

  const Tj = T1 - Ta;
  const Tm = Tk - Tl;
  outR[3] = Tj - Tm;
  outR[1] = Tj + Tm;
  const Tr = Tp - To;
  const Ts = Te - Th;
  outI[1] = Tr - Ts;
  outI[3] = Ts + Tr;

  return [outR, outI];
}

module.exports = { t2_4 };
