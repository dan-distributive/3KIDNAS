'use strict';

// =============================================================================
// r2cfII_12.js -- faithful JS port of rdft/scalar/r2cf/r2cfII_12.c
// (non-FMA branch). R2HCII ("shifted") radix-12 direct codelet -- the
// "cldm" middle-column combine for an EVEN outer radix (see r2cfII_8.js's
// header for the general even-radix derivation: clean 6-and-6 Cr/Ci
// split, out[11-k]=Ci[k]).
// INPUT: ph[0..11], R0[k]=ph[2k] (k=0..5), R1[k]=ph[2k+1] (k=0..5).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP353553390 = 0.353553390593273762200422181052424519642417969;
const KP707106781 = 0.707106781186547524400844362104849039284835938;
const KP612372435 = 0.612372435695794524549321018676472847991486870;
const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function r2cfII_12(ph) {
  const T1 = ph[0];
  const T3 = ph[4];
  const T2 = ph[8];
  const Tx = KP866025403 * (T2 + T3);
  const Tg = KP500000000 * (T3 - T2) + T1; // FMA
  const T4 = T1 + T2 - T3;
  const Tz = ph[6];
  const Th = ph[10];
  const Ti = ph[2];
  const Ty = Th + Ti;
  const Tj = KP866025403 * (Th - Ti);
  const TA = KP500000000 * Ty + Tz; // FMA

  const T5 = ph[3];
  const T6 = ph[11];
  const T7 = ph[7];
  const T8 = T6 - T7;
  const T9 = T5 + T8;
  const Tm = KP612372435 * (T6 + T7);
  const Tl = KP707106781 * T5 - KP353553390 * T8; // FNMS

  const Td = ph[9];
  const Ta = ph[5];
  const Tb = ph[1];
  const Tc = Ta - Tb;
  const Te = Tc - Td;
  const Tp = KP353553390 * Tc + KP707106781 * Td; // FMA
  const To = KP612372435 * (Ta + Tb);

  const out = new Float64Array(12);

  const Tf = KP707106781 * (T9 + Te);
  out[1] = T4 - Tf;
  out[4] = T4 + Tf;
  const TE = KP707106781 * (Te - T9);
  const TF = Tz - Ty;
  out[7] = TE - TF;
  out[10] = TE + TF;

  {
    const Tk = Tg - Tj;
    const TB = Tx - TA;
    const Tn = Tl - Tm;
    const Tq = To - Tp;
    const Tr = Tn + Tq;
    const Tw = Tn - Tq;
    out[5] = Tk - Tr;
    out[9] = Tw + TB;
    out[0] = Tk + Tr;
    out[8] = Tw - TB;
  }
  {
    const Ts = Tg + Tj;
    const TD = Tx + TA;
    const Tt = To + Tp;
    const Tu = Tm + Tl;
    const Tv = Tt - Tu;
    const TC = Tu + Tt;
    out[3] = Ts - Tv;
    out[6] = TD - TC;
    out[2] = Ts + Tv;
    out[11] = -(TC + TD);
  }

  return out;
}

module.exports = { r2cfII_12 };
