'use strict';

// =============================================================================
// r2cfII_16.js -- faithful JS port of rdft/scalar/r2cf/r2cfII_16.c (non-FMA
// branch). R2HCII ("shifted") radix-16 direct codelet -- the "cldm"
// middle-column combine for an EVEN outer radix (see r2cfII_8.js's header
// for the general even-radix derivation: clean 8-and-8 Cr/Ci split, no
// self-paired frequency, out[15-k]=Ci[k]).
// INPUT: ph[0..15], R0[k]=ph[2k] (k=0..7), R1[k]=ph[2k+1] (k=0..7).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP555570233 = 0.555570233019602224742830813948532874374937191;
const KP831469612 = 0.831469612302545237078788377617905756738560812;
const KP980785280 = 0.980785280403230449126182236134239036973933731;
const KP195090322 = 0.195090322016128267848284868477022240927691618;
const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function r2cfII_16(ph) {
  const T1 = ph[0];
  const TU = ph[8];
  const T2 = ph[4];
  const T3 = ph[12];
  const T4 = KP707106781 * (T2 - T3);
  const TT = KP707106781 * (T2 + T3);
  const T5 = T1 + T4;
  const T11 = TU - TT;
  const TB = T1 - T4;
  const TV = TT + TU;

  const Tq = ph[15];
  const Tt = ph[7];
  const Tn = ph[3];
  const To = ph[11];
  const Tp = KP707106781 * (Tn - To);
  const Ts = KP707106781 * (Tn + To);
  const Tr = Tp - Tq;
  const TK = Tt - Ts;
  const Tu = Ts + Tt;
  const TJ = Tp + Tq;

  const Te = ph[1];
  const Tk = ph[9];
  const Tf = ph[5];
  const Tg = ph[13];
  const Th = KP707106781 * (Tf - Tg);
  const Tj = KP707106781 * (Tf + Tg);
  const Ti = Te + Th;
  const TH = Tk - Tj;
  const Tl = Tj + Tk;
  const TG = Te - Th;

  const T6 = ph[2];
  const T7 = ph[10];
  const T8 = KP923879532 * T6 - KP382683432 * T7; // FNMS
  const TC = KP382683432 * T6 + KP923879532 * T7; // FMA
  const T9 = ph[6];
  const Ta = ph[14];
  const Tb = KP382683432 * T9 - KP923879532 * Ta; // FNMS
  const TD = KP923879532 * T9 + KP382683432 * Ta; // FMA

  const Tc = T8 + Tb;
  const T10 = Tb - T8;
  const TE = TC - TD;
  const TS = TC + TD;

  const out = new Float64Array(16);

  {
    const Td = T5 - Tc;
    const TW = TS + TV;
    const Tm = KP195090322 * Ti + KP980785280 * Tl; // FMA
    const Tv = KP195090322 * Tr - KP980785280 * Tu; // FNMS
    const Tw = Tm + Tv;
    const TR = Tv - Tm;
    out[4] = Td - Tw;
    out[8] = TR + TW;
    out[3] = Td + Tw;
    out[15] = TR - TW;
  }
  {
    const Tx = T5 + Tc;
    const TY = TV - TS;
    const Ty = KP980785280 * Ti - KP195090322 * Tl; // FNMS
    const Tz = KP980785280 * Tr + KP195090322 * Tu; // FMA
    const TA = Ty + Tz;
    const TX = Tz - Ty;
    out[7] = Tx - TA;
    out[12] = TX + TY;
    out[0] = Tx + TA;
    out[11] = TX - TY;
  }
  {
    const TF = TB + TE;
    const T12 = T10 - T11;
    const TI = KP831469612 * TG + KP555570233 * TH; // FMA
    const TL = KP831469612 * TJ + KP555570233 * TK; // FMA
    const TM = TI - TL;
    const TZ = TI + TL;
    out[6] = TF - TM;
    out[13] = T12 - TZ;
    out[1] = TF + TM;
    out[10] = -(TZ + T12);
  }
  {
    const TN = TB - TE;
    const T14 = T10 + T11;
    const TO = KP831469612 * TK - KP555570233 * TJ; // FNMS
    const TP = KP831469612 * TH - KP555570233 * TG; // FNMS
    const TQ = TO - TP;
    const T13 = TP + TO;
    out[5] = TN - TQ;
    out[14] = T13 + T14;
    out[2] = TN + TQ;
    out[9] = T13 - T14;
  }

  return out;
}

module.exports = { r2cfII_16 };
