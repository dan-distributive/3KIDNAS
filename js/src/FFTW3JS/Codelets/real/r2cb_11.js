'use strict';

// =============================================================================
// r2cb_11.js -- faithful JS port of rdft/scalar/r2cb/r2cb_11.c (non-FMA
// branch). O[0..10] packed halfcomplex (O[0..5]=Re0..Re5, O[6..10]=Im5..
// Im1, same convention as r2cf_11.js) -> x[0..10] real.
//
// Several 3-term FMA/FNMS/FNMA sums here -- same grouping discipline as
// r2cf_11.js: each macro's result is its own named intermediate, combined
// via simple binary +/- in the C source's own left-to-right order.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP2_000000000 = 2.0;
const KP1_918985947 = 1.918985947228994779780736114132655398124909697;
const KP1_309721467 = 1.309721467890570128113850144932587106367582399;
const KP284629676 = 0.284629676546570280887585337232739337582102722;
const KP830830026 = 0.830830026003772851058548298459246407048009821;
const KP1_682507065 = 1.682507065662362337723623297838735435026584997;
const KP563465113 = 0.563465113682859395422835830693233798071555798;
const KP1_511499148 = 1.511499148708516567548071687944688840359434890;
const KP1_979642883 = 1.979642883761865464752184075553437574753038744;
const KP1_819263990 = 1.819263990709036742823430766158056920120482102;
const KP1_081281634 = 1.081281634911195164215271908637383390863541216;

function r2cb_11(O) {
  const T8 = O[9], Tc = O[10], T9 = O[7], Ta = O[6], Tb = O[8];

  const Da1 = KP1_081281634 * T8 + KP1_819263990 * T9;
  const Da2 = -(KP1_979642883 * Ta + KP1_511499148 * Tb);
  const Td = Da1 + Da2 - KP563465113 * Tc;

  const Dl1 = KP1_979642883 * T8 + KP1_819263990 * Ta;
  const Dl2 = -(KP563465113 * T9 + KP1_081281634 * Tb);
  const Tl = Dl1 + Dl2 - KP1_511499148 * Tc;

  const Df1 = KP563465113 * T8 + KP1_819263990 * Tb;
  const Df2 = -(KP1_511499148 * Ta + KP1_081281634 * T9);
  const Tf = Df1 + Df2 - KP1_979642883 * Tc;

  const Dh1 = KP1_081281634 * Tc + KP1_819263990 * T8;
  const Dh2 = KP1_979642883 * Tb + KP1_511499148 * T9;
  const Th = Dh1 + Dh2 + KP563465113 * Ta;

  const Dj1 = KP563465113 * Tb + KP1_979642883 * T9;
  const Dj2 = KP1_081281634 * Ta - KP1_511499148 * T8;
  const Tj = Dj1 + Dj2 - KP1_819263990 * Tc;

  const T1 = O[0], T2 = O[1], T6 = O[5], T5 = O[4], T4 = O[3], T3 = O[2];

  const D7a = KP1_682507065 * T3 + T1;
  const D7b = KP830830026 * T5 - KP284629676 * T6;
  const D7c = -(KP1_309721467 * T4 + KP1_918985947 * T2);
  const T7 = D7a + D7b + D7c;

  const Dka = KP1_682507065 * T4 + T1;
  const Dkb = KP830830026 * T6 - KP1_918985947 * T5;
  const Dkc = -(KP284629676 * T3 + KP1_309721467 * T2);
  const Tk = Dka + Dkb + Dkc;

  const Dea = KP830830026 * T4 + T1;
  const Deb = KP1_682507065 * T5 - KP1_309721467 * T6;
  const Dec = -(KP1_918985947 * T3 + KP284629676 * T2);
  const Te = Dea + Deb + Dec;

  const Dga = KP1_682507065 * T2 + T1;
  const Dgb = KP830830026 * T3 - KP1_918985947 * T6;
  const Dgc = -(KP1_309721467 * T5 + KP284629676 * T4);
  const Tg = Dga + Dgb + Dgc;

  const Dia = KP830830026 * T2 + T1;
  const Dib = KP1_682507065 * T6 - KP284629676 * T5;
  const Dic = -(KP1_918985947 * T4 + KP1_309721467 * T3);
  const Ti = Dia + Dib + Dic;

  const x = new Float64Array(11);
  x[6] = T7 - Td;
  x[8] = Te - Tf;
  x[4] = Tk + Tl;
  x[5] = T7 + Td;
  x[7] = Tk - Tl;
  x[2] = Ti + Tj;
  x[3] = Te + Tf;
  x[10] = Tg + Th;
  x[1] = Tg - Th;
  x[9] = Ti - Tj;
  x[0] = KP2_000000000 * (T2 + T3 + T4 + T5 + T6) + T1;

  return x;
}

module.exports = { r2cb_11 };
