'use strict';

// =============================================================================
// r2cbIII_12.js -- faithful JS port of rdft/scalar/r2cb/r2cbIII_12.c
// (non-FMA branch). HC2RIII ("shifted") radix-12 direct codelet -- the
// BACKWARD counterpart of r2cfII_12.js. INPUT in[0..11] uses
// r2cfII_12.js's OUTPUT convention (in[k]=Cr[k] for k=0..5, in[11-k]=Ci[k]
// for k=0..5); OUTPUT out[0..11] uses r2cfII_12.js's INPUT convention
// (out[2k]=R0[k], out[2k+1]=R1[k]).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP1_414213562 = 1.414213562373095048801688724209698078569671875;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;
const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function r2cbIII_12(inArr) {
  const T1 = inArr[1];
  const T2 = inArr[5];
  const T3 = inArr[2];
  const T4 = T2 + T3;
  const T5 = T1 + T4;
  const Tw = KP866025403 * (T2 - T3);
  const Tb = T1 - KP500000000 * T4; // FNMS

  const Tq = inArr[10];
  const Tc = inArr[6];
  const Td = inArr[9];
  const Tr = Td - Tc;
  const Te = KP866025403 * (Tc + Td);
  const Tx = KP500000000 * Tr + Tq; // FMA
  const Ts = Tq - Tr;

  const T6 = inArr[4];
  const T7 = inArr[0];
  const T8 = inArr[3];
  const T9 = T7 + T8;
  const Ta = T6 + T9;
  const TA = KP866025403 * (T7 - T8);
  const Tg = T6 - KP500000000 * T9; // FNMS

  const To = inArr[7];
  const Th = inArr[11];
  const Ti = inArr[8];
  const Tn = Ti - Th;
  const Tj = KP866025403 * (Th + Ti);
  const Tz = KP500000000 * Tn + To; // FMA
  const Tp = Tn - To;

  const out = new Float64Array(12);
  out[0] = KP2_000000000 * (T5 + Ta);
  out[6] = KP2_000000000 * (Ts + Tp);
  const Tt = Tp - Ts;
  const Tu = T5 - Ta;
  out[3] = KP1_414213562 * (Tt - Tu);
  out[9] = KP1_414213562 * (Tu + Tt);

  {
    const Tf = Tb - Te;
    const Tk = Tg + Tj;
    const Tv = Tf - Tk;
    const Ty = Tw + Tx;
    const TB = Tz - TA;
    const TC = Ty + TB;
    out[4] = -(KP2_000000000 * (Tf + Tk));
    out[10] = KP2_000000000 * (TB - Ty);
    out[1] = KP1_414213562 * (Tv - TC);
    out[7] = KP1_414213562 * (Tv + TC);
  }
  {
    const Tl = Tb + Te;
    const Tm = Tg - Tj;
    const TF = Tm - Tl;
    const TD = TA + Tz;
    const TE = Tx - Tw;
    const TG = TE + TD;
    out[8] = KP2_000000000 * (Tl + Tm);
    out[5] = KP1_414213562 * (TF + TG);
    out[2] = KP2_000000000 * (TD - TE);
    out[11] = KP1_414213562 * (TF - TG);
  }

  return out;
}

module.exports = { r2cbIII_12 };
