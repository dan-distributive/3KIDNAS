'use strict';

// =============================================================================
// r2cbIII_9.js -- faithful JS port of rdft/scalar/r2cb/r2cbIII_9.c
// (non-FMA branch). HC2RIII ("shifted") radix-9 direct codelet -- the
// BACKWARD counterpart of r2cfII_9.js. INPUT in[0..8] uses r2cfII_9.js's
// OUTPUT convention (in[0..4]=Cr, in[8-k]=Ci[k] for k=0..3); OUTPUT
// out[0..8] uses r2cfII_9.js's INPUT convention (out[2k]=R0[k],
// out[2k+1]=R1[k]).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP642787609 = 0.642787609686539326322643409907263432907559884;
const KP766044443 = 0.766044443118978035202392650555416673935832457;
const KP1_326827896 = 1.326827896337876792410842639271782594433726619;
const KP1_113340798 = 1.113340798452838732905825904094046265936583811;
const KP984807753 = 0.984807753012208059366743024589523013670643252;
const KP173648177 = 0.173648177666930348851716626769314796000375677;
const KP1_705737063 = 1.705737063904886419256501927880148143872040591;
const KP300767466 = 0.300767466360870593278543795225003852144476517;
const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;
const KP1_732050807 = 1.732050807568877293527446341505872366942805254;

function r2cbIII_9(inArr) {
  const Tg = inArr[7];
  const Th = KP1_732050807 * Tg;
  const T1 = inArr[4];
  const T2 = inArr[1];
  const Tf = T2 - T1;
  const T3 = KP2_000000000 * T2 + T1; // FMA
  const Ts = Tf - Th;
  const Ti = Tf + Th;

  const T4 = inArr[3];
  const Td = inArr[5];
  const T5 = inArr[0];
  const T6 = inArr[2];
  const T7 = T5 + T6;
  const Tm = KP866025403 * (T6 - T5);
  const Ta = inArr[6];
  const Tb = inArr[8];
  const Tc = Ta - Tb;
  const Tk = KP866025403 * (Tb + Ta);

  const T8 = T4 + T7;
  const Tn = KP500000000 * Tc + Td; // FMA
  const To = Tm - Tn;
  const Tu = Tm + Tn;
  const Tj = KP500000000 * T7 - T4; // FMS
  const Tl = Tj + Tk;
  const Tt = Tj - Tk;

  const out = new Float64Array(9);
  out[0] = KP2_000000000 * T8 + T3; // FMA
  const T9 = T8 - T3;
  const Te = KP1_732050807 * (Tc - Td);
  out[3] = T9 + Te;
  out[6] = Te - T9;

  {
    const Tr = KP300767466 * To - KP1_705737063 * Tl; // FNMS
    const Tp = KP173648177 * Tl + KP984807753 * To; // FMA
    const Tq = Ti - Tp;
    out[2] = -(KP2_000000000 * Tp + Ti); // FMA
    out[8] = Tr - Tq;
    out[5] = Tq + Tr;
    const Tx = KP1_113340798 * Tt + KP1_326827896 * Tu; // FMA
    const Tv = KP766044443 * Tt - KP642787609 * Tu; // FNMS
    const Tw = Tv - Ts;
    out[1] = KP2_000000000 * Tv + Ts; // FMA
    out[7] = Tx - Tw;
    out[4] = Tw + Tx;
  }

  return out;
}

module.exports = { r2cbIII_9 };
