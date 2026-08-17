'use strict';

// =============================================================================
// r2cb_7.js -- faithful JS port of rdft/scalar/r2cb/r2cb_7.c (non-FMA branch).
// O[0..6] packed halfcomplex (O[0..3]=Re0..Re3, O[4..6]=Im3..Im1, same
// convention as r2cf_7.js) -> x[0..6] real.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP2_000000000 = 2.000000000000000000000000000000000000000000000;
const KP1_801937735 = 1.801937735804838252472204639014890102331838324;
const KP445041867 = 0.445041867912628808577805128993589518932711138;
const KP1_246979603 = 1.246979603717467061050009768008479621264549462;
const KP867767478 = 0.867767478235116240951536665696717509219981456;
const KP1_949855824 = 1.949855824363647214036263365987862434465571601;
const KP1_563662964 = 1.563662964936059617416889053348115500464669037;

function r2cb_7(O) {
  const T6 = O[5], T8 = O[6], T7 = O[4];
  const M9 = KP1_563662964 * T6 - KP1_949855824 * T7;
  const T9 = M9 - KP867767478 * T8;
  const MD = KP867767478 * T6 + KP1_563662964 * T7;
  const Td = MD - KP1_949855824 * T8;
  const MB = KP1_563662964 * T8 + KP1_949855824 * T6;
  const Tb = MB + KP867767478 * T7;

  const T1 = O[0], T4 = O[3], T2 = O[1], T3 = O[2];
  const M5a = KP1_246979603 * T3 + T1, M5b = KP445041867 * T4 + KP1_801937735 * T2;
  const T5 = M5a - M5b;
  const MCa = KP1_246979603 * T4 + T1, MCb = KP1_801937735 * T3 + KP445041867 * T2;
  const Tc = MCa - MCb;
  const MAa = KP1_246979603 * T2 + T1, MAb = KP1_801937735 * T4 + KP445041867 * T3;
  const Ta = MAa - MAb;

  const x = new Float64Array(7);
  x[4] = T5 - T9;
  x[3] = T5 + T9;
  x[2] = Tc + Td;
  x[5] = Tc - Td;
  x[6] = Ta + Tb;
  x[1] = Ta - Tb;
  x[0] = KP2_000000000 * (T2 + T3 + T4) + T1;
  return x;
}

module.exports = { r2cb_7 };
