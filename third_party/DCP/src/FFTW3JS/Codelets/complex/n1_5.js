'use strict';

// =============================================================================
// n1_5.js -- faithful JS port of dft/scalar/codelets/n1_5.c (non-FMA branch),
// FFTW3's direct (base-case) radix-5 complex DFT codelet.
// FMA(a,b,c)=a*b+c, FNMS(a,b,c)=c-a*b macro-expand to plain two-step
// arithmetic (ARCH_PREFERS_FMA undefined, -ffp-contract=off) -- confirmed
// zero fmadd/fnmadd/etc in the compiled n1_5.o object file.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.25;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP559016994 = 0.559016994374947424102293417182819058860154590;

function n1_5(ri, ii) {
  const T1 = ri[0], To = ii[0];

  const T2 = ri[1], T3 = ri[4];
  const T4 = T2 + T3;
  const T5 = ri[2], T6 = ri[3];
  const T7 = T5 + T6;
  const T8 = T4 + T7;
  const Tt = T5 - T6;
  const T9 = KP559016994 * (T4 - T7);
  const Ts = T2 - T3;

  const Tc = ii[1], Td = ii[4];
  const Tl = Tc + Td;
  const Tf = ii[2], Tg = ii[3];
  const Tm = Tf + Tg;
  const Te = Tc - Td;
  const Tp = Tl + Tm;
  const Th = Tf - Tg;
  const Tn = KP559016994 * (Tl - Tm);

  const ro = new Float64Array(5), io = new Float64Array(5);
  ro[0] = T1 + T8;
  io[0] = To + Tp;

  const Ti = KP951056516 * Te + KP587785252 * Th;
  const Tk = KP951056516 * Th - KP587785252 * Te;
  const Ta = T1 - KP250000000 * T8;
  const Tb = T9 + Ta;
  const Tj = Ta - T9;
  ro[4] = Tb - Ti;
  ro[3] = Tj + Tk;
  ro[1] = Tb + Ti;
  ro[2] = Tj - Tk;

  const Tu = KP951056516 * Ts + KP587785252 * Tt;
  const Tv = KP951056516 * Tt - KP587785252 * Ts;
  const Tq = To - KP250000000 * Tp;
  const Tr = Tn + Tq;
  const Tw = Tq - Tn;
  io[1] = Tr - Tu;
  io[3] = Tw - Tv;
  io[4] = Tu + Tr;
  io[2] = Tv + Tw;

  return [ro, io];
}

module.exports = { n1_5 };
