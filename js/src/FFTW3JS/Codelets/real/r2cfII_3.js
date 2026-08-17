'use strict';

// =============================================================================
// r2cfII_3.js -- faithful JS port of rdft/scalar/r2cf/r2cfII_3.c (non-FMA
// branch). R2HCII ("shifted") radix-3 direct codelet -- the "cldm"
// middle-column combine for an ODD outer radix (see r2cfII_7.js's header
// for the general odd-radix derivation and the out[k]=Cr[k] (k=0..
// floor(r/2)), out[r-1-k]=Ci[k] (k=0..floor(r/2)-1) convention -- for r=3
// that's out[0]=Cr[0], out[1]=Cr[1], out[2]=Ci[0]).
// INPUT: ph[0..2], R0[k]=ph[2k] (k=0), R1[k]=ph[2k+1] (k=0), plus
// R0[WS(rs,1)]=ph[2].
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function r2cfII_3(ph) {
  const T1 = ph[0];
  const T2 = ph[1];
  const T3 = ph[2];
  const T4 = T2 - T3;

  const out = new Float64Array(3);
  out[0] = KP500000000 * T4 + T1; // FMA
  out[1] = T1 - T4;
  out[2] = -(KP866025403 * (T2 + T3));
  return out;
}

module.exports = { r2cfII_3 };
