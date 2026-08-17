'use strict';

// =============================================================================
// r2cf_10.js -- faithful JS port of rdft/scalar/r2cf/r2cf_10.c (non-FMA
// branch). x[0..9] (real) -> O[0..9] packed halfcomplex: O[0..5]=Re0..Re5,
// O[6..9]=Im4..Im1 (O[n-k]=Im_k). R0/R1 (stride-2) input convention -- see
// r2cf_6.js's header.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.250000000000000000000000000000000000000000000;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP587785252 = 0.587785252292473129168705954639072768597652438;

function r2cf_10(x) {
  const Tg = x[0];
  const Th = x[5];
  const Ti = Tg - Th;
  const Tt = Tg + Th;
  const T8 = x[4];
  const T9 = x[9];
  const Ta = T8 - T9;
  const Tn = T8 + T9;
  const Tb = x[6];
  const Tc = x[1];
  const Td = Tb - Tc;
  const To = Tb + Tc;
  const Te = Ta + Td;
  const Tv = Tn + To;
  const T1 = x[2];
  const T2 = x[7];
  const T3 = T1 - T2;
  const Tq = T1 + T2;
  const T4 = x[8];
  const T5 = x[3];
  const T6 = T4 - T5;
  const Tr = T4 + T5;
  const T7 = T3 + T6;
  const Tu = Tq + Tr;

  const O = new Float64Array(10);
  const Tl = Td - Ta;
  const Tm = T3 - T6;
  O[9] = KP587785252 * Tl - KP951056516 * Tm;
  O[7] = KP951056516 * Tl + KP587785252 * Tm;
  const Tf = KP559016994 * (T7 - Te);
  const Tj = T7 + Te;
  const Tk = Ti - KP250000000 * Tj;
  O[1] = Tf + Tk;
  O[5] = Ti + Tj;
  O[3] = Tk - Tf;

  const Tp = Tn - To;
  const Ts = Tq - Tr;
  O[8] = KP951056516 * Tp - KP587785252 * Ts;
  O[6] = KP587785252 * Tp + KP951056516 * Ts;
  const Ty = KP559016994 * (Tu - Tv);
  const Tw = Tu + Tv;
  const Tx = Tt - KP250000000 * Tw;
  O[2] = Tx - Ty;
  O[0] = Tt + Tw;
  O[4] = Ty + Tx;

  return O;
}

module.exports = { r2cf_10 };
