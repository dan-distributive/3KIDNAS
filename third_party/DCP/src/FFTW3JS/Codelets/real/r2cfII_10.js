'use strict';

// =============================================================================
// r2cfII_10.js -- faithful JS port of rdft/scalar/r2cf/r2cfII_10.c
// (non-FMA branch). R2HCII ("shifted") radix-10 direct codelet -- the
// "cldm" middle-column combine for an EVEN outer radix (see r2cfII_8.js's
// header for the general even-radix derivation: clean 5-and-5 Cr/Ci
// split, out[9-k]=Ci[k]).
// INPUT: ph[0..9], R0[k]=ph[2k] (k=0..4), R1[k]=ph[2k+1] (k=0..4).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.250000000000000000000000000000000000000000000;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP559016994 = 0.559016994374947424102293417182819058860154590;

function r2cfII_10(ph) {
  const T1 = ph[0];
  const To = ph[5];

  const T2 = ph[4];
  const T3 = ph[6];
  const T4 = T2 - T3;
  const T5 = ph[8];
  const T6 = ph[2];
  const T7 = T5 - T6;
  const T8 = T4 + T7;
  const Tq = T5 + T6;
  const T9 = KP559016994 * (T4 - T7);
  const Tp = T2 + T3;

  const Tc = ph[1];
  const Td = ph[9];
  const Tm = Tc + Td;
  const Tf = ph[3];
  const Tg = ph[7];
  const Tl = Tf + Tg;
  const Te = Tc - Td;
  const Ts = KP559016994 * (Tm + Tl);
  const Th = Tf - Tg;
  const Tn = Tl - Tm;

  const out = new Float64Array(10);
  out[2] = T1 + T8;
  out[7] = Tn - To;

  {
    const Ti = KP951056516 * Te + KP587785252 * Th; // FMA
    const Tk = KP951056516 * Th - KP587785252 * Te; // FNMS
    const Ta = T1 - KP250000000 * T8; // FNMS
    const Tb = T9 + Ta;
    const Tj = Ta - T9;
    out[4] = Tb - Ti;
    out[3] = Tj + Tk;
    out[0] = Tb + Ti;
    out[1] = Tj - Tk;
  }
  {
    const Tr = KP951056516 * Tp + KP587785252 * Tq; // FMA
    const Tw = KP951056516 * Tq - KP587785252 * Tp; // FNMS
    const Tt = KP250000000 * Tn + To; // FMA
    const Tu = Ts + Tt;
    const Tv = Tt - Ts;
    out[9] = -(Tr + Tu);
    out[6] = Tw + Tv;
    out[5] = Tr - Tu;
    out[8] = Tv - Tw;
  }

  return out;
}

module.exports = { r2cfII_10 };
