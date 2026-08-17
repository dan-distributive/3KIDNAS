'use strict';

// =============================================================================
// n1_6.js -- faithful JS port of dft/scalar/codelets/n1_6.c (non-FMA branch),
// FFTW3's direct (base-case) radix-6 complex DFT codelet.
//
// Moved here from RaderSolver.js (which originally hand-inlined this codelet
// for its N=37 special case, before RaderSolver.js was generalized to use
// Composite1D.js for arbitrary N -- keeping it inline there would have
// created a circular require (RaderSolver -> Composite1D -> this registry ->
// RaderSolver)). RaderSolver.js now imports it from here, same as every
// other radix.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP500000000 = 0.5;

function n1_6(ri, ii) {
  const T1 = ri[0], T2 = ri[3];
  const T3 = T1 - T2, Tb = T1 + T2;
  const To = ii[0], Tp = ii[3];
  const Tq = To - Tp, Tx = To + Tp;
  const T4 = ri[2], T5 = ri[5];
  const T6 = T4 - T5, Tc = T4 + T5;
  const T7 = ri[4], T8 = ri[1];
  const T9 = T7 - T8, Td = T7 + T8;
  const Ta = T6 + T9, Te = Tc + Td;
  const Tg = ii[2], Th = ii[5];
  const Ti = Tg - Th, Tu = Tg + Th;
  const Tj = ii[4], Tk = ii[1];
  const Tl = Tj - Tk, Tv = Tj + Tk;
  const Tr = Ti + Tl, Ty = Tu + Tv;

  const ro = new Float64Array(6), io = new Float64Array(6);
  ro[3] = T3 + Ta; io[3] = Tq + Tr;
  ro[0] = Tb + Te; io[0] = Tx + Ty;

  const Tf = T3 - KP500000000 * Ta;
  const Tm = KP866025403 * (Ti - Tl);
  ro[5] = Tf - Tm; ro[1] = Tf + Tm;
  const Tn = KP866025403 * (T9 - T6);
  const Ts = Tq - KP500000000 * Tr;
  io[1] = Tn + Ts; io[5] = Ts - Tn;

  const Tt = Tb - KP500000000 * Te;
  const Tw = KP866025403 * (Tu - Tv);
  ro[2] = Tt - Tw; ro[4] = Tt + Tw;
  const Tz = Tx - KP500000000 * Ty;
  const TA = KP866025403 * (Td - Tc);
  io[2] = Tz - TA; io[4] = TA + Tz;

  return [ro, io];
}

module.exports = { n1_6 };
