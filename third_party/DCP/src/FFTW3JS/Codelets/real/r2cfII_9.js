'use strict';

// =============================================================================
// r2cfII_9.js -- faithful JS port of rdft/scalar/r2cf/r2cfII_9.c (non-FMA
// branch). R2HCII ("shifted") radix-9 direct codelet -- the "cldm"
// middle-column combine for an ODD outer radix (see r2cfII_7.js's header
// for the general odd-radix derivation: out[k]=Cr[k] for k=0..4,
// out[8-k]=Ci[k] for k=0..3).
// INPUT: ph[0..8], R0[k]=ph[2k] (k=0..4), R1[k]=ph[2k+1] (k=0..3).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP663413948 = 0.663413948168938396205421319635891297216863310;
const KP642787609 = 0.642787609686539326322643409907263432907559884;
const KP556670399 = 0.556670399226419366452912952047023132968291906;
const KP766044443 = 0.766044443118978035202392650555416673935832457;
const KP852868531 = 0.852868531952443209628250963940074071936020296;
const KP173648177 = 0.173648177666930348851716626769314796000375677;
const KP984807753 = 0.984807753012208059366743024589523013670643252;
const KP150383733 = 0.150383733180435296639271897612501926072238258;
const KP813797681 = 0.813797681349373692844693217248393223289101568;
const KP342020143 = 0.342020143325668733044099614682259580763083368;
const KP939692620 = 0.939692620785908384054109277324731469936208134;
const KP296198132 = 0.296198132726023843175338011893050938967728390;
const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP500000000 = 0.500000000000000000000000000000000000000000000;

function r2cfII_9(ph) {
  const T1 = ph[0];
  const T2 = ph[3];
  const T3 = ph[6];
  const T4 = T2 - T3;
  const To = T2 + T3;

  const T6 = ph[2];
  const T7 = ph[5];
  const T8 = ph[8];
  const T9 = T7 - T8;
  const Ta = T6 - T9;
  const Tl = T7 + T8;
  const Tk = KP500000000 * T9 + T6; // FMA

  const Tb = ph[4];
  const Tc = ph[1];
  const Td = ph[7];
  const Te = Tc + Td;
  const Tf = Tb - Te;
  const Ti = KP500000000 * Te + Tb; // FMA
  const Th = Tc - Td;

  const out = new Float64Array(9);
  out[7] = KP866025403 * (Tf - Ta);
  const T5 = T1 - T4;
  const Tg = Ta + Tf;
  out[1] = T5 - KP500000000 * Tg; // FNMS
  out[4] = T5 + Tg;

  const Tr = KP500000000 * T4 + T1; // FMA
  const Tt = KP296198132 * Th + KP939692620 * Ti; // FMA
  const Tw = KP342020143 * Ti - KP813797681 * Th; // FNMS
  const Tv = KP150383733 * Tl - KP984807753 * Tk; // FNMS
  const Tu = KP173648177 * Tk + KP852868531 * Tl; // FMA
  const Tp = KP766044443 * Tk - KP556670399 * Tl; // FNMS
  const Tq = KP852868531 * Th + KP173648177 * Ti; // FMA
  const Ts = Tp + Tq;
  const Tj = KP150383733 * Th - KP984807753 * Ti; // FNMS
  const Tm = KP642787609 * Tk + KP663413948 * Tl; // FMA
  const Tn = Tj - Tm;

  out[8] = Tn - KP866025403 * To; // FNMS
  out[0] = Tr + Ts;
  out[5] = KP866025403 * ((Tp - Tq) - To) - KP500000000 * Tn;
  out[3] = KP866025403 * (Tm + Tj) + Tr - KP500000000 * Ts;
  out[6] = KP866025403 * (To - (Tu + Tt)) + KP500000000 * (Tw - Tv);
  out[2] = KP500000000 * (Tt - Tu) + Tr + KP866025403 * (Tv + Tw);

  return out;
}

module.exports = { r2cfII_9 };
