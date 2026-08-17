'use strict';

// =============================================================================
// hf_12.js -- faithful JS port of rdft/scalar/r2cf/hf_12.c (non-FMA
// branch). Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hf_6.js -- the
// PLAIN (non-alt-codegen) family, TW_FULL twinstr. Radix 12 has no hf2_12
// in real FFTW's own registered codelet set, so this plain family is the
// only way to cover it -- and real FFTW's actual choice for several N
// (e.g. n=24, 48, 132, 192) this port previously had to blacklist away
// via a different, wrong, unported-outer-radix workaround.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function hf_12(cr, ci, Wc, Ws) {
  const T1 = cr[0];
  const T1W = ci[0];
  const T3 = cr[4], T5 = ci[4];
  const T6 = Wc[4] * T3 + Ws[4] * T5; // FMA
  const T16 = Wc[4] * T5 - Ws[4] * T3; // FNMS
  const T8 = cr[8], Ta = ci[8];
  const Tb = Wc[8] * T8 + Ws[8] * Ta; // FMA
  const T17 = Wc[8] * Ta - Ws[8] * T8; // FNMS
  const T18 = KP866025403 * (T16 - T17);
  const T23 = KP866025403 * (Tb - T6);
  const Tc = T6 + Tb;
  const T15 = T1 - KP500000000 * Tc; // FNMS
  const T1V = T16 + T17;
  const T22 = T1W - KP500000000 * T1V; // FNMS

  const TO = cr[9], TQ = ci[9];
  const TR = Wc[9] * TO + Ws[9] * TQ; // FMA
  const T1E = Wc[9] * TQ - Ws[9] * TO; // FNMS
  const TY = cr[5], T10 = ci[5];
  const T11 = Wc[5] * TY + Ws[5] * T10; // FMA
  const T1n = Wc[5] * T10 - Ws[5] * TY; // FNMS
  const TT = cr[1], TV = ci[1];
  const TW = Wc[1] * TT + Ws[1] * TV; // FMA
  const T1m = Wc[1] * TV - Ws[1] * TT; // FNMS
  const T1o = KP866025403 * (T1m - T1n);
  const T1D = KP866025403 * (T11 - TW);
  const T12 = TW + T11;
  const T1l = TR - KP500000000 * T12; // FNMS
  const T1F = T1m + T1n;
  const T1G = T1E - KP500000000 * T1F; // FNMS

  const Tf = cr[6], Th = ci[6];
  const Ti = Wc[6] * Tf + Ws[6] * Th; // FMA
  const T1S = Wc[6] * Th - Ws[6] * Tf; // FNMS
  const Tp = cr[2], Tr = ci[2];
  const Ts = Wc[2] * Tp + Ws[2] * Tr; // FMA
  const T1c = Wc[2] * Tr - Ws[2] * Tp; // FNMS
  const Tk = cr[10], Tm = ci[10];
  const Tn = Wc[10] * Tk + Ws[10] * Tm; // FMA
  const T1b = Wc[10] * Tm - Ws[10] * Tk; // FNMS
  const T1d = KP866025403 * (T1b - T1c);
  const T26 = KP866025403 * (Ts - Tn);
  const Tt = Tn + Ts;
  const T1a = Ti - KP500000000 * Tt; // FNMS
  const T1T = T1b + T1c;
  const T25 = T1S - KP500000000 * T1T; // FNMS

  const Tx = cr[3], Tz = ci[3];
  const TA = Wc[3] * Tx + Ws[3] * Tz; // FMA
  const T1y = Wc[3] * Tz - Ws[3] * Tx; // FNMS
  const TH = cr[11], TJ = ci[11];
  const TK = Wc[11] * TH + Ws[11] * TJ; // FMA
  const T1i = Wc[11] * TJ - Ws[11] * TH; // FNMS
  const TC = cr[7], TE = ci[7];
  const TF = Wc[7] * TC + Ws[7] * TE; // FMA
  const T1h = Wc[7] * TE - Ws[7] * TC; // FNMS
  const T1j = KP866025403 * (T1h - T1i);
  const T1B = KP866025403 * (TK - TF);
  const TL = TF + TK;
  const T1g = TA - KP500000000 * TL; // FNMS
  const T1z = T1h + T1i;
  const T1A = T1y - KP500000000 * T1z; // FNMS

  const outCr = new Float64Array(12), outCi = new Float64Array(12);

  {
    const Td = T1 + Tc;
    const Tu = Ti + Tt;
    const Tv = Td + Tu;
    const T1N = Td - Tu;
    const T1U = T1S + T1T;
    const T1X = T1V + T1W;
    const T1Y = T1U + T1X;
    const T20 = T1X - T1U;

    const TM = TA + TL;
    const T13 = TR + T12;
    const T14 = TM + T13;
    const T1Z = TM - T13;
    const T1O = T1y + T1z;
    const T1P = T1E + T1F;
    const T1Q = T1O - T1P;
    const T1R = T1O + T1P;

    outCi[5] = Tv - T14;
    outCr[9] = T1Z - T20;
    outCi[8] = T1Z + T20;
    outCr[0] = Tv + T14;
    outCr[3] = T1N - T1Q;
    outCr[6] = T1R - T1Y;
    outCi[11] = T1R + T1Y;
    outCi[2] = T1N + T1Q;
  }
  {
    const T19 = T15 - T18;
    const T1e = T1a - T1d;
    const T1f = T19 + T1e;
    const T1x = T19 - T1e;
    const T24 = T22 - T23;
    const T27 = T25 - T26;
    const T28 = T24 - T27;
    const T2a = T27 + T24;

    const T1k = T1g - T1j;
    const T1p = T1l - T1o;
    const T1q = T1k + T1p;
    const T21 = T1p - T1k;
    const T1C = T1A - T1B;
    const T1H = T1D - T1G;
    const T1I = T1C + T1H;
    const T29 = T1H - T1C;

    outCr[2] = T1f - T1q;
    outCr[8] = T29 - T2a;
    outCi[9] = T29 + T2a;
    outCi[3] = T1f + T1q;
    outCi[0] = T1x - T1I;
    outCr[11] = T21 - T28;
    outCi[6] = T21 + T28;
    outCr[5] = T1x + T1I;
  }
  {
    const T1r = T15 + T18;
    const T1s = T1a + T1d;
    const T1t = T1r + T1s;
    const T1J = T1r - T1s;
    const T2c = T23 + T22;
    const T2d = T26 + T25;
    const T2e = T2c - T2d;
    const T2g = T2d + T2c;

    const T1u = T1g + T1j;
    const T1v = T1l + T1o;
    const T1w = T1u + T1v;
    const T2b = T1v - T1u;
    const T1K = T1B + T1A;
    const T1L = T1D + T1G;
    const T1M = T1K - T1L;
    const T2f = T1K + T1L;

    outCi[1] = T1t - T1w;
    outCr[1] = T1J + T1M;
    outCr[4] = T1t + T1w;
    outCi[4] = T1J - T1M;
    outCr[7] = T2b - T2e;
    outCi[7] = T2f + T2g;
    outCi[10] = T2b + T2e;
    outCr[10] = T2f - T2g;
  }

  return [outCr, outCi];
}

module.exports = { hf_12 };
