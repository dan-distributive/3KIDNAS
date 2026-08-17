'use strict';

// =============================================================================
// r2cfII_5.js -- faithful JS port of rdft/scalar/r2cf/r2cfII_5.c (non-FMA
// branch). R2HCII ("shifted") radix-5 direct codelet -- the "cldm"
// middle-column combine for an ODD outer radix (see r2cfII_7.js's header
// for the general odd-radix derivation: out[k]=Cr[k] for k=0..floor(r/2),
// out[r-1-k]=Ci[k] for k=0..floor(r/2)-1 -- for r=5 that's out[0..2]=Cr,
// out[3]=Ci[1], out[4]=Ci[0]).
// INPUT: ph[0..4], R0[k]=ph[2k] (k=0..2), R1[k]=ph[2k+1] (k=0..1).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.250000000000000000000000000000000000000000000;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP587785252 = 0.587785252292473129168705954639072768597652438;

function r2cfII_5(ph) {
  const T8 = ph[0];
  const T1 = ph[2];
  const T2 = ph[3];
  const T3 = T1 - T2;
  const T4 = ph[4];
  const T5 = ph[1];
  const T6 = T4 - T5;
  const T9 = T3 + T6;
  const Tc = T4 + T5;
  const Tb = T1 + T2;

  const out = new Float64Array(5);
  out[2] = T8 + T9;
  out[3] = KP587785252 * Tb - KP951056516 * Tc; // FNMS
  out[4] = -(KP951056516 * Tb + KP587785252 * Tc); // FMA
  const T7 = KP559016994 * (T3 - T6);
  const Ta = T8 - KP250000000 * T9; // FNMS
  out[0] = T7 + Ta;
  out[1] = Ta - T7;

  return out;
}

module.exports = { r2cfII_5 };
