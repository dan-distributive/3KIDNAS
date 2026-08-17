'use strict';

// =============================================================================
// r2cbIII_8.js -- faithful JS port of rdft/scalar/r2cb/r2cbIII_8.c (non-FMA
// branch). HC2RIII ("shifted") radix-8 direct codelet -- the BACKWARD
// counterpart of r2cfII_8.js (see that file's header for the R2HCII/
// HC2RIII math and the even-radix 4-and-4 Cr/Ci split).
//
// Calling convention swap vs forward (matches r2cbIII_7.js): INPUT
// in[0..7] uses r2cfII_8.js's OUTPUT convention (in[k]=Cr[k] for k=0..3,
// in[7-k]=Ci[k] for k=0..3); OUTPUT out[0..7] (phase p) uses r2cfII_8.js's
// INPUT convention (out[2k]=R0[k], out[2k+1]=R1[k]).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP1_414213562 = 1.414213562373095048801688724209698078569671875;
const KP765366864 = 0.765366864730179543456919968060797733522689125;
const KP1_847759065 = 1.847759065022573512256366378793576573644833252;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;

function r2cbIII_8(inArr) {
  const T1 = inArr[0];
  const T2 = inArr[3];
  const T3 = T1 + T2;
  const T7 = T1 - T2;
  const Td = inArr[7];
  const Te = inArr[4];
  const Tf = Td + Te;
  const Tl = Te - Td;

  const T4 = inArr[2];
  const T5 = inArr[1];
  const T6 = T4 + T5;
  const Tc = T4 - T5;
  const T8 = inArr[5];
  const T9 = inArr[6];
  const Ta = T8 + T9;
  const Tk = T8 - T9;

  const out = new Float64Array(8);
  out[0] = KP2_000000000 * (T3 + T6);
  out[4] = KP2_000000000 * (Tl - Tk);
  const Tb = T7 - Ta;
  const Tg = Tc + Tf;
  out[1] = KP1_847759065 * Tb - KP765366864 * Tg; // FNMS
  out[5] = -(KP765366864 * Tb + KP1_847759065 * Tg); // -FMA

  const Th = T7 + Ta;
  const Ti = Tc - Tf;
  out[3] = KP765366864 * Th + KP1_847759065 * Ti; // FMA
  out[7] = KP765366864 * Ti - KP1_847759065 * Th; // FNMS
  const Tj = T3 - T6;
  const Tm = Tk + Tl;
  out[2] = KP1_414213562 * (Tj + Tm);
  out[6] = KP1_414213562 * (Tm - Tj);

  return out;
}

module.exports = { r2cbIII_8 };
