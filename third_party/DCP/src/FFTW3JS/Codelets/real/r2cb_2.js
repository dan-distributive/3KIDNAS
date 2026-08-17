'use strict';

// =============================================================================
// r2cb_2.js -- faithful JS port of rdft/scalar/r2cb/r2cb_2.c (non-FMA branch),
// FFTW3's direct (base-case) radix-2 halfcomplex->real inverse DFT codelet.
//
// Input O[0..n-1]: packed halfcomplex (same format r2cf_n/RealEngine1D
// produce -- see r2cf_2.js's header). Output x[0..n-1]: real, UNNORMALIZED
// (matches FFTW's own c2r convention -- caller divides by n).
// Convention confirmed via direct-r2c.c's apply_hc2r (same rs0/ioffset
// machinery as apply_r2hc, R0/R1 now the OUTPUT's even/odd-index split)
// AND numerically, by hand-deriving the inverse-DFT formula for r2cb_3 and
// matching its FMA/FNMS terms exactly.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

function r2cb_2(O) {
  const T1 = O[0], T2 = O[1];
  const x = new Float64Array(2);
  x[1] = T1 - T2;
  x[0] = T1 + T2;
  return x;
}

module.exports = { r2cb_2 };
