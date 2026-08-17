'use strict';

// =============================================================================
// r2cbIII_20.js -- faithful JS port of rdft/scalar/r2cb/r2cbIII_20.c
// (non-FMA branch). HC2RIII ("shifted") radix-20 direct codelet -- the
// BACKWARD counterpart of r2cfII_20.js (see that file's header). INPUT
// in[0..19] uses r2cfII_20.js's OUTPUT convention (in[k]=Cr[k] for k=0..9,
// in[19-k]=Ci[k] for k=0..9); OUTPUT out[0..19] (phase p) uses
// r2cfII_20.js's INPUT convention (out[2k]=R0[k], out[2k+1]=R1[k]).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP1_414213562 = 1.414213562373095048801688724209698078569671875;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;
const KP250000000 = 0.250000000000000000000000000000000000000000000;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP559016994 = 0.559016994374947424102293417182819058860154590;

function r2cbIII_20(inArr) {
  const T1 = inArr[2];
  const T5 = inArr[9];
  const T6 = inArr[5];
  const T7 = T5 + T6;
  const T12 = T5 - T6;
  const T2 = inArr[6];
  const T3 = inArr[1];
  const T4 = T2 + T3;
  const T11 = T2 - T3;
  const Tj = KP559016994 * (T4 - T7);
  const T1k = KP587785252 * T11 - KP951056516 * T12; // FNMS
  const T13 = KP951056516 * T11 + KP587785252 * T12; // FMA
  const T8 = T4 + T7;
  const Tk = T1 - KP250000000 * T8; // FNMS

  const T17 = inArr[17];
  const Tp = inArr[14];
  const Tq = inArr[10];
  const Tr = Tp - Tq;
  const T15 = Tp + Tq;
  const Tm = inArr[13];
  const Tn = inArr[18];
  const To = Tm + Tn;
  const T14 = Tm - Tn;
  const Ts = KP951056516 * To + KP587785252 * Tr; // FMA
  const T16 = KP559016994 * (T14 + T15);
  const TI = KP587785252 * To - KP951056516 * Tr; // FNMS
  const T18 = T14 - T15;
  const T19 = T17 - KP250000000 * T18; // FNMS

  const Ta = inArr[7];
  const Te = inArr[0];
  const Tf = inArr[4];
  const Tg = Te + Tf;
  const TR = Te - Tf;
  const Tb = inArr[3];
  const Tc = inArr[8];
  const Td = Tb + Tc;
  const TQ = Tb - Tc;
  const Tu = KP559016994 * (Td - Tg);
  const T1i = KP587785252 * TQ - KP951056516 * TR; // FNMS
  const TS = KP951056516 * TQ + KP587785252 * TR; // FMA
  const Th = Td + Tg;
  const Tv = Ta - KP250000000 * Th; // FNMS

  const TX = inArr[12];
  const TA = inArr[15];
  const TB = inArr[19];
  const TC = TA - TB;
  const TU = TB + TA;
  const Tx = inArr[16];
  const Ty = inArr[11];
  const Tz = Tx + Ty;
  const TT = Ty - Tx;
  const TD = KP951056516 * Tz + KP587785252 * TC; // FMA
  const TV = KP559016994 * (TT - TU);
  const TL = KP951056516 * TC - KP587785252 * Tz; // FNMS
  const TW = TT + TU;
  const TY = KP250000000 * TW + TX; // FMA

  const out = new Float64Array(20);

  {
    const T9 = T1 + T8;
    const Ti = Ta + Th;
    const T1w = T9 - Ti;
    const T1t = T18 + T17;
    const T1u = TX - TW;
    const T1v = T1t + T1u;
    out[0] = KP2_000000000 * (T9 + Ti);
    out[10] = KP2_000000000 * (T1u - T1t);
    out[5] = KP1_414213562 * (T1v - T1w);
    out[15] = KP1_414213562 * (T1w + T1v);
  }
  {
    const TH = Tk - Tj;
    const TJ = TH + TI;
    const TO = TH - TI;
    const T1l = T19 - T16;
    const T1m = T1k + T1l;
    const T1q = T1l - T1k;
    const TK = Tv - Tu;
    const TM = TK + TL;
    const TN = TL - TK;
    const T1h = TV + TY;
    const T1j = T1h - T1i;
    const T1r = T1i + T1h;

    out[8] = KP2_000000000 * (TJ + TM);
    out[12] = KP2_000000000 * (TN - TO);
    out[18] = KP2_000000000 * (T1r - T1q);
    out[2] = KP2_000000000 * (T1j - T1m);

    const T1p = TM - TJ;
    const T1s = T1q + T1r;
    out[3] = KP1_414213562 * (T1p - T1s);
    out[13] = KP1_414213562 * (T1p + T1s);
    const T1n = TO + TN;
    const T1o = T1m + T1j;
    out[17] = KP1_414213562 * (T1n - T1o);
    out[7] = KP1_414213562 * (T1n + T1o);
  }
  {
    const Tl = Tj + Tk;
    const Tt = Tl - Ts;
    const TG = Tl + Ts;
    const T1a = T16 + T19;
    const T1b = T13 + T1a;
    const T1f = T1a - T13;
    const Tw = Tu + Tv;
    const TE = Tw + TD;
    const TF = TD - Tw;
    const TZ = TV - TY;
    const T10 = TS + TZ;
    const T1e = TZ - TS;

    out[16] = KP2_000000000 * (Tt + TE);
    out[4] = KP2_000000000 * (TF - TG);
    out[14] = KP2_000000000 * (T1f + T1e);
    out[6] = KP2_000000000 * (T1b + T10);

    const T1d = TG + TF;
    const T1g = T1e - T1f;
    out[9] = KP1_414213562 * (T1d + T1g);
    out[19] = KP1_414213562 * (T1g - T1d);
    const TP = Tt - TE;
    const T1c = T10 - T1b;
    out[1] = KP1_414213562 * (TP + T1c);
    out[11] = KP1_414213562 * (T1c - TP);
  }

  return out;
}

module.exports = { r2cbIII_20 };
