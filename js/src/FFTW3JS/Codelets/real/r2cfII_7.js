'use strict';

// =============================================================================
// r2cfII_7.js -- faithful JS port of rdft/scalar/r2cf/r2cfII_7.c (non-FMA
// branch). R2HCII ("shifted") radix-7 direct codelet: the "cldm" middle-
// column combine used by rdft/hc2hc-direct.c whenever the RECURSED
// sub-transform size m is EVEN (the standard hf_r/hf2_r twiddle codelets
// only handle columns 1..floor(m/2)-1; column m/2 itself -- the m-sized
// sub-transform's own Nyquist bin, a lone real value per phase, with no
// Im-pair to combine against -- needs this genuinely different transform).
//
// R2HCII(r)[k] = X[k+0.5] = sum_j x[j] * exp(-2*pi*i*j*(k+0.5)/r), the
// SAME r-point DFT but evaluated at HALF-INTEGER frequencies -- verified by
// hand against this exact codelet's arithmetic (k=0: matches Re/Im of
// X[0.5] exactly; k=3, r odd's own self-paired "1.5-analog" frequency
// (k+0.5=3.5, self-conjugate mod 7): matches the lone real Cr[3] with no
// Ci[3] counterpart, as expected -- exactly one self-paired frequency for
// ANY odd r, none for even r).
//
// INPUT/OUTPUT INDEX MAPPING (derived from rdft/direct-r2c.c's mkplan,
// NOT the same as r2cf_r.js's plain R2HC convention -- confirmed by
// tracing ego->k(I, I+rs0, O, O+ioffset, rs, csr, csi, ...)'s exact
// pointer arithmetic): input ph[0..r-1] (phase p's CURRENT value at the
// m/2 "middle" column) splits as R0[k]=ph[2k], R1[k]=ph[2k+1] (SAME
// even/odd split as r2cf_r.js). Output is WHERE THIS DIFFERS: Cr[k]
// writes to out[k] (k=0..floor(r/2)), Ci[k] writes to out[r-1-k] (k=0..
// floor(r/2)-1) -- NOT out[r-k] as r2cf_r.js's packed-halfcomplex uses.
// (ioffset(kind,n,cs) = cs*(n-1) for R2HCII, vs cs*n for plain R2HC --
// this one-less is exactly the out[r-1-k] vs out[r-k] difference.)
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP900968867 = 0.900968867902419126236102319507445051165919162;
const KP222520933 = 0.222520933956314404288902564496794759466355569;
const KP623489801 = 0.623489801858733530525004884004239810632274731;
const KP433883739 = 0.433883739117558120475768332848358754609990728;
const KP974927912 = 0.974927912181823607018131682993931217232785801;
const KP781831482 = 0.781831482468029808708444526674057750232334519;

function r2cfII_7(ph) {
  const T1 = ph[0], T8 = ph[1], T9 = ph[6];
  const Ta = T8 - T9, Td = T8 + T9;
  const T2 = ph[2], T3 = ph[5];
  const T4 = T2 - T3, Tb = T2 + T3;
  const T5 = ph[3], T6 = ph[4];
  const T7 = T5 - T6, Tc = T5 + T6;

  const M1 = KP781831482 * Tb + KP974927912 * Tc;
  const Ci0 = -(M1 + KP433883739 * Td);

  const M2 = KP781831482 * Tc - KP974927912 * Td;
  const Ci1 = M2 - KP433883739 * Tb;

  const M3 = KP623489801 * T4 + T1;
  const M4 = KP222520933 * T7 + KP900968867 * Ta;
  const Cr0 = M3 + M4;

  const M5 = KP974927912 * Tb - KP781831482 * Td;
  const Ci2 = M5 - KP433883739 * Tc;

  const M6 = KP900968867 * T7 + T1;
  const M7 = -(KP623489801 * Ta + KP222520933 * T4);
  const Cr2 = M6 + M7;

  const M8 = KP222520933 * Ta + T1;
  const M9 = -(KP623489801 * T7 + KP900968867 * T4);
  const Cr1 = M8 + M9;

  const N1 = T1 + T4, N2 = T7 + Ta;
  const Cr3 = N1 - N2;

  const out = new Float64Array(7);
  out[0] = Cr0; out[1] = Cr1; out[2] = Cr2; out[3] = Cr3;
  out[4] = Ci2; out[5] = Ci1; out[6] = Ci0;
  return out;
}

module.exports = { r2cfII_7 };
