'use strict';

// =============================================================================
// r2cbIII_7.js -- faithful JS port of rdft/scalar/r2cb/r2cbIII_7.c (non-FMA
// branch). HC2RIII ("shifted") radix-7 direct codelet: the BACKWARD
// counterpart of r2cfII_7.js -- see that file's header for the full
// derivation (R2HCII/HC2RIII math, the "cldm middle-column" role, and the
// out[r-1-k] vs out[r-k] physical-index distinction from plain r2cb_r.js).
//
// Calling convention swap vs the forward direction (matches how R2HC/HC2R
// themselves are related): rdft/direct-r2c.c's apply_hc2r calls
// ego->k(O, O+rs0, I, I+ioffset, ...) -- R0/R1 are now the OUTPUT (real
// phases), Cr/Ci are the INPUT (the shifted-frequency pairs), i.e. exactly
// reversed from r2cfII_7.js. So this function's INPUT convention (in[0..6])
// is identical to r2cfII_7.js's OUTPUT convention (in[k]=Cr[k] for k=0..3,
// in[r-1-k]=Ci[k] for k=0..2), and its OUTPUT convention (out[0..6], phase
// p) is identical to r2cfII_7.js's INPUT convention (out[2k]=R0[k],
// out[2k+1]=R1[k]).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP2_000000000 = 2.000000000000000000000000000000000000000000000;
const KP1_246979603 = 1.246979603717467061050009768008479621264549462;
const KP1_801937735 = 1.801937735804838252472204639014890102331838324;
const KP445041867 = 0.445041867912628808577805128993589518932711138;
const KP867767478 = 0.867767478235116240951536665696717509219981456;
const KP1_949855824 = 1.949855824363647214036263365987862434465571601;
const KP1_563662964 = 1.563662964936059617416889053348115500464669037;

function r2cbIII_7(inArr) {
  const T6 = inArr[4], T8 = inArr[6], T7 = inArr[5];
  const M1 = KP1_563662964 * T6 + KP1_949855824 * T7;
  const T9 = M1 + KP867767478 * T8;
  const M2 = KP1_563662964 * T7 - KP1_949855824 * T8;
  const Td = M2 - KP867767478 * T6;
  const M3 = KP1_949855824 * T6 - KP1_563662964 * T8;
  const Tb = M3 - KP867767478 * T7;

  const T1 = inArr[3], T4 = inArr[0], T2 = inArr[2], T3 = inArr[1];
  const M4 = KP445041867 * T3 + KP1_801937735 * T4;
  const M5 = -(KP1_246979603 * T2 + T1);
  const T5 = M4 + M5;
  const M6 = KP1_801937735 * T2 + KP445041867 * T4;
  const M7 = -(KP1_246979603 * T3 + T1);
  const Tc = M6 + M7;
  const M8 = KP1_246979603 * T4 + T1;
  const M9 = -(KP1_801937735 * T3 + KP445041867 * T2);
  const Ta = M8 + M9;

  const out = new Float64Array(7);
  out[1] = T5 - T9;
  out[6] = -(T5 + T9);
  out[4] = Td - Tc;
  out[3] = Tc + Td;
  out[5] = Tb - Ta;
  out[2] = Ta + Tb;
  out[0] = KP2_000000000 * (T2 + T3 + T4) + T1;
  return out;
}

module.exports = { r2cbIII_7 };
