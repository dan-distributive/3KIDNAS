'use strict';

// =============================================================================
// r2cfII_8.js -- faithful JS port of rdft/scalar/r2cf/r2cfII_8.c (non-FMA
// branch). R2HCII ("shifted") radix-8 direct codelet -- the "cldm" middle-
// column combine for an EVEN outer radix (see r2cfII_7.js's header for the
// general R2HCII derivation and the out[k]/out[r-1-k] convention).
//
// Radix parity changes the output SPLIT, not the convention: odd r has one
// self-paired half-integer frequency (an extra lone Cr entry with no Ci
// counterpart); even r has NONE (every half-integer frequency k+0.5 pairs
// with a DISTINCT (r-1-k)+0.5, since (r-1)/2 is not an integer when r is
// even) -- so for r=8 it's a clean 4-and-4 split: Cr[0..3]=out[0..3],
// Ci[0..3]=out[7..4] (out[r-1-k]=Ci[k], same formula, k now runs 0..r/2-1
// instead of 0..floor(r/2)-1 with a leftover Cr[floor(r/2)]). Confirmed
// directly against this codelet's own Cr/Ci index set (4+4, no overlap).
//
// INPUT: ph[0..7], R0[k]=ph[2k] (k=0..3), R1[k]=ph[2k+1] (k=0..3) -- same
// stride-2 convention as r2cf_r.js (see r2cf_6.js's header).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function r2cfII_8(ph) {
  const T1 = ph[0];
  const Tj = ph[4];
  const T2 = ph[2];
  const T3 = ph[6];
  const T4 = KP707106781 * (T2 - T3);
  const Ti = KP707106781 * (T2 + T3);

  const T6 = ph[1];
  const T7 = ph[5];
  const T8 = KP923879532 * T6 - KP382683432 * T7; // FNMS
  const Te = KP382683432 * T6 + KP923879532 * T7; // FMA
  const T9 = ph[3];
  const Ta = ph[7];
  const Tb = KP382683432 * T9 - KP923879532 * Ta; // FNMS
  const Tf = KP923879532 * T9 + KP382683432 * Ta; // FMA

  const out = new Float64Array(8);
  const T5 = T1 + T4;
  const Tc = T8 + Tb;
  out[3] = T5 - Tc;
  out[0] = T5 + Tc;
  const Th = Te + Tf;
  const Tk = Ti + Tj;
  out[7] = -(Th + Tk);
  out[4] = Tk - Th;

  const Td = T1 - T4;
  const Tg = Te - Tf;
  out[2] = Td - Tg;
  out[1] = Td + Tg;
  const Tl = Tb - T8;
  const Tm = Tj - Ti;
  out[5] = Tl - Tm;
  out[6] = Tl + Tm;

  return out;
}

module.exports = { r2cfII_8 };
