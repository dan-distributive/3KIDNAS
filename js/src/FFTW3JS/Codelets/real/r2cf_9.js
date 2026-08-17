'use strict';

// =============================================================================
// r2cf_9.js -- faithful JS port of rdft/scalar/r2cf/r2cf_9.c (non-FMA branch).
// x[0..8] (real) -> O[0..8] packed halfcomplex: O[0..4]=Re0..Re4,
// O[5..8]=Im4..Im1 (O[n-k]=Im_k, same convention as r2cf_5/7.js).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP939692620 = 0.939692620785908384054109277324731469936208134;
const KP296198132 = 0.296198132726023843175338011893050938967728390;
const KP342020143 = 0.342020143325668733044099614682259580763083368;
const KP813797681 = 0.813797681349373692844693217248393223289101568;
const KP984807753 = 0.984807753012208059366743024589523013670643252;
const KP150383733 = 0.150383733180435296639271897612501926072238258;
const KP642787609 = 0.642787609686539326322643409907263432907559884;
const KP663413948 = 0.663413948168938396205421319635891297216863310;
const KP852868531 = 0.852868531952443209628250963940074071936020296;
const KP173648177 = 0.173648177666930348851716626769314796000375677;
const KP556670399 = 0.556670399226419366452912952047023132968291906;
const KP766044443 = 0.766044443118978035202392650555416673935832457;
const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP500000000 = 0.5;

function r2cf_9(x) {
  const T1 = x[0];
  const T2 = x[3], T3 = x[6];
  const T4 = T2 + T3, Tr = T3 - T2;

  const T6 = x[1], T7 = x[4], T8 = x[7];
  const T9 = T7 + T8;
  const Ta = T6 + T9;
  const Tl = T8 - T7;
  const Ti = T6 - KP500000000 * T9;

  const Tb = x[2], Tc = x[5], Td = x[8];
  const Te = Tc + Td;
  const Tf = Tb + Te;
  const Tk = Tb - KP500000000 * Te;
  const Tj = Td - Tc;

  const O = new Float64Array(9);
  O[6] = KP866025403 * (Tf - Ta);
  const T5 = T1 + T4;
  const Tg = Ta + Tf;
  O[3] = T5 - KP500000000 * Tg;
  O[0] = T5 + Tg;

  const Tt = KP866025403 * Tr;
  const Th = T1 - KP500000000 * T4;
  const Tm = KP766044443 * Ti + KP556670399 * Tl;
  const Tn = KP173648177 * Tk + KP852868531 * Tj;
  const To = Tm + Tn;
  const Tp = KP663413948 * Tl - KP642787609 * Ti;
  const Tq = KP150383733 * Tj - KP984807753 * Tk;
  const Ts = Tp + Tq;

  O[1] = Th + To;
  O[8] = Tt + Ts;
  const M4 = KP866025403 * (Tp - Tq) + Th;
  O[4] = M4 - KP500000000 * To;
  const M5 = KP866025403 * (Tr + (Tn - Tm));
  O[5] = M5 - KP500000000 * Ts;
  const M7a = KP813797681 * Tj - KP342020143 * Tk;
  const M7b = KP150383733 * Tl + KP984807753 * Ti;
  O[7] = M7a - M7b - Tt;
  const M2a = KP173648177 * Ti + Th;
  const M2b = KP296198132 * Tj + KP939692620 * Tk;
  O[2] = M2a - M2b - KP852868531 * Tl;

  return O;
}

module.exports = { r2cf_9 };
