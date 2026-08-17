'use strict';

// =============================================================================
// r2cf_2.js -- faithful JS port of rdft/scalar/r2cf/r2cf_2.c (non-FMA branch),
// FFTW3's direct (base-case) radix-2 real->halfcomplex DFT codelet.
//
// Calling convention (derived from rdft/direct-r2c.c's apply_r2hc + mkplan:
// rs0=raw stride, R1=R0+rs0, codelet stride=2*raw -- so R0 walks EVEN input
// indices, R1 walks ODD indices; csi has stride -cs and ioffset=n*cs, so
// output is FFTW's packed "halfcomplex" format: O[0..floor(n/2)] = real
// parts forward from the start, O[n-1..ceil(n/2)] = imaginary parts
// BACKWARD from the end. Confirmed structurally via the stride/ioffset
// arithmetic AND numerically via a hand-derived DFT check on r2cf_3.):
//   x[0..n-1] (real input, length n) -> O[0..n-1] (packed halfcomplex).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

function r2cf_2(x) {
  const T1 = x[0], T2 = x[1];
  const O = new Float64Array(2);
  O[1] = T1 - T2;
  O[0] = T1 + T2;
  return O;
}

module.exports = { r2cf_2 };
