'use strict';

// =============================================================================
// r2cf_11.js -- faithful JS port of rdft/scalar/r2cf/r2cf_11.c (non-FMA
// branch). x[0..10] (real) -> O[0..10] packed halfcomplex: O[0..5]=Re0..Re5,
// O[6..10]=Im5..Im1 (O[n-k]=Im_k, same convention as r2cf_5/7/9.js).
//
// Several outputs here are 3-term FMA/FNMS/FNMA sums (matching n1_11.c's
// shape on the complex side) -- each macro's result is computed as its own
// named intermediate (M1/M2/M3) BEFORE combining via simple binary +/- in
// the SAME left-to-right order as the C source, never flattened into one
// expression -- see n1_11.js's header for why: floating-point addition
// isn't associative, so flattening silently changes rounding.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP654860733 = 0.654860733945285064056925072466293553183791199;
const KP142314838 = 0.142314838273285140443792668616369668791051361;
const KP959492973 = 0.959492973614497389890368057066327699062454848;
const KP415415013 = 0.415415013001886425529274149229623203524004910;
const KP841253532 = 0.841253532831181168861811648919367717513292498;
const KP989821441 = 0.989821441880932732376092037776718787376519372;
const KP909631995 = 0.909631995354518371411715383079028460060241051;
const KP281732556 = 0.281732556841429697711417915346616899035777899;
const KP540640817 = 0.540640817455597582107635954318691695431770608;
const KP755749574 = 0.755749574354258283774035843972344420179717445;

function r2cf_11(x) {
  const T1 = x[0];
  const T2 = x[2], T3 = x[9];
  const T4 = T2 + T3, Tl = T3 - T2;
  const Te = x[1], Tf = x[10];
  const Tg = Te + Tf, Th = Tf - Te;

  const Tb = x[3], Tc = x[8];
  const Td = Tb + Tc, Ti = Tc - Tb;

  const T8 = x[5], T9 = x[6];
  const Ta = T8 + T9, Tk = T9 - T8;
  const T5 = x[4], T6 = x[7];
  const T7 = T5 + T6, Tj = T6 - T5;

  const O = new Float64Array(11);

  const N7a = KP755749574 * Th + KP540640817 * Ti;
  const N7b = KP281732556 * Tj - KP909631995 * Tk;
  O[7] = N7a + N7b - KP989821441 * Tl;

  const N4a = KP841253532 * Td + T1;
  const N4b = KP415415013 * Ta - KP959492973 * T7;
  const N4c = -(KP142314838 * T4 + KP654860733 * Tg);
  O[4] = N4a + N4b + N4c;

  const N9a = KP909631995 * Th + KP755749574 * Tl;
  const N9b = -(KP540640817 * Tk + KP989821441 * Tj);
  O[9] = N9a + N9b - KP281732556 * Ti;

  const N6a = KP281732556 * Th + KP755749574 * Ti;
  const N6b = KP989821441 * Tk - KP909631995 * Tj;
  O[6] = N6a + N6b - KP540640817 * Tl;

  const N10a = KP540640817 * Th + KP909631995 * Tl;
  const N10b = KP989821441 * Ti + KP755749574 * Tj;
  O[10] = N10a + N10b + KP281732556 * Tk;

  const N8a = KP989821441 * Th + KP540640817 * Tj;
  const N8b = KP755749574 * Tk - KP909631995 * Ti;
  O[8] = N8a + N8b - KP281732556 * Tl;

  const N3a = KP415415013 * Td + T1;
  const N3b = KP841253532 * T7 - KP654860733 * Ta;
  const N3c = -(KP959492973 * T4 + KP142314838 * Tg);
  O[3] = N3a + N3b + N3c;

  const N1a = KP841253532 * Tg + T1;
  const N1b = KP415415013 * T4 - KP959492973 * Ta;
  const N1c = -(KP654860733 * T7 + KP142314838 * Td);
  O[1] = N1a + N1b + N1c;

  O[0] = T1 + Tg + T4 + Td + T7 + Ta;

  const N2a = KP415415013 * Tg + T1;
  const N2b = KP841253532 * Ta - KP142314838 * T7;
  const N2c = -(KP959492973 * Td + KP654860733 * T4);
  O[2] = N2a + N2b + N2c;

  const N5a = KP841253532 * T4 + T1;
  const N5b = KP415415013 * T7 - KP142314838 * Ta;
  const N5c = -(KP654860733 * Td + KP959492973 * Tg);
  O[5] = N5a + N5b + N5c;

  return O;
}

module.exports = { r2cf_11 };
