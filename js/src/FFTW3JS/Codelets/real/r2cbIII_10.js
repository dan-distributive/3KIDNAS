'use strict';

// =============================================================================
// r2cbIII_10.js -- faithful JS port of rdft/scalar/r2cb/r2cbIII_10.c
// (non-FMA branch). HC2RIII ("shifted") radix-10 direct codelet -- the
// BACKWARD counterpart of r2cfII_10.js. INPUT in[0..9] uses r2cfII_10.js's
// OUTPUT convention (in[k]=Cr[k] for k=0..4, in[9-k]=Ci[k] for k=0..4);
// OUTPUT out[0..9] uses r2cfII_10.js's INPUT convention (out[2k]=R0[k],
// out[2k+1]=R1[k]).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP1_902113032 = 1.902113032590307144232878666758764286811397268;
const KP1_175570504 = 1.175570504584946258337411909278145537195304875;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;
const KP1_118033988 = 1.118033988749894848204586834365638117720309180;

function r2cbIII_10(inArr) {
  const T1 = inArr[2];
  const To = inArr[7];

  const T2 = inArr[4];
  const T3 = inArr[0];
  const T4 = T2 + T3;
  const T5 = inArr[3];
  const T6 = inArr[1];
  const T7 = T5 + T6;
  const T8 = T4 + T7;
  const Tq = T5 - T6;
  const Ta = KP1_118033988 * (T7 - T4);
  const Tp = T2 - T3;

  const Tc = inArr[5];
  const Td = inArr[9];
  const Tm = Tc + Td;
  const Tf = inArr[8];
  const Tg = inArr[6];
  const Tl = Tg + Tf;
  const Te = Tc - Td;
  const Ts = KP1_118033988 * (Tl + Tm);
  const Th = Tf - Tg;
  const Tn = Tl - Tm;

  const out = new Float64Array(10);
  out[0] = KP2_000000000 * (T1 + T8);
  out[5] = KP2_000000000 * (Tn - To);

  {
    const Ti = KP1_175570504 * Te - KP1_902113032 * Th; // FNMS
    const Tj = KP1_902113032 * Te + KP1_175570504 * Th; // FMA
    const T9 = KP500000000 * T8 - KP2_000000000 * T1; // FNMS
    const Tb = T9 - Ta;
    const Tk = T9 + Ta;
    out[2] = Tb + Ti;
    out[6] = Tk + Tj;
    out[8] = Ti - Tb;
    out[4] = Tj - Tk;
  }
  {
    const Tr = KP1_902113032 * Tp + KP1_175570504 * Tq; // FMA
    const Tv = KP1_902113032 * Tq - KP1_175570504 * Tp; // FNMS
    const Tt = KP500000000 * Tn + KP2_000000000 * To; // FMA
    const Tu = Ts + Tt;
    const Tw = Tt - Ts;
    out[1] = -(Tr + Tu);
    out[7] = Tw - Tv;
    out[9] = Tr - Tu;
    out[3] = Tv + Tw;
  }

  return out;
}

module.exports = { r2cbIII_10 };
