'use strict';

// =============================================================================
// hf_10.js -- faithful JS port of rdft/scalar/r2cf/hf_10.c (non-FMA
// branch). Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hf_6.js/
// hf_12.js -- the PLAIN (non-alt-codegen) family, TW_FULL twinstr. Radix
// 10 has no hf2_10 in real FFTW's own registered codelet set, so this
// plain family is the only way to cover it.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP250000000 = 0.250000000000000000000000000000000000000000000;
const KP559016994 = 0.559016994374947424102293417182819058860154590;

function hf_10(cr, ci, Wc, Ws) {
  const T1 = cr[0];
  const T1A = ci[0];
  const T3 = cr[5], T5 = ci[5];
  const T6 = Wc[5] * T3 + Ws[5] * T5; // FMA
  const T1B = Wc[5] * T5 - Ws[5] * T3; // FNMS
  const T7 = T1 - T6;
  const T1R = T1B + T1A;
  const TT = T1 + T6;
  const T1C = T1A - T1B;

  const Tw = cr[4], Ty = ci[4];
  const Tz = Wc[4] * Tw + Ws[4] * Ty; // FMA
  const T1b = Wc[4] * Ty - Ws[4] * Tw; // FNMS
  const TM = cr[1], TO = ci[1];
  const TP = Wc[1] * TM + Ws[1] * TO; // FMA
  const T1e = Wc[1] * TO - Ws[1] * TM; // FNMS
  const TB = cr[9], TD = ci[9];
  const TE = Wc[9] * TB + Ws[9] * TD; // FMA
  const T1c = Wc[9] * TD - Ws[9] * TB; // FNMS
  const TH = cr[6], TJ = ci[6];
  const TK = Wc[6] * TH + Ws[6] * TJ; // FMA
  const T1f = Wc[6] * TJ - Ws[6] * TH; // FNMS

  const TF = Tz - TE;
  const TQ = TK - TP;
  const TR = TF + TQ;
  const T1o = T1b + T1c;
  const T1p = T1f + T1e;
  const T1P = T1o + T1p;
  const TX = Tz + TE;
  const TY = TK + TP;
  const TZ = TX + TY;
  const T1d = T1b - T1c;
  const T1g = T1e - T1f;
  const T1x = T1g - T1d;

  const T9 = cr[2], Tb = ci[2];
  const Tc = Wc[2] * T9 + Ws[2] * Tb; // FMA
  const T14 = Wc[2] * Tb - Ws[2] * T9; // FNMS
  const Tp = cr[3], Tr = ci[3];
  const Ts = Wc[3] * Tp + Ws[3] * Tr; // FMA
  const T18 = Wc[3] * Tr - Ws[3] * Tp; // FNMS
  const Te = cr[7], Tg = ci[7];
  const Th = Wc[7] * Te + Ws[7] * Tg; // FMA
  const T15 = Wc[7] * Tg - Ws[7] * Te; // FNMS
  const Tk = cr[8], Tm = ci[8];
  const Tn = Wc[8] * Tk + Ws[8] * Tm; // FMA
  const T17 = Wc[8] * Tm - Ws[8] * Tk; // FNMS

  const Ti = Tc - Th;
  const Tt = Tn - Ts;
  const Tu = Ti + Tt;
  const T1r = T14 + T15;
  const T1s = T17 + T18;
  const T1O = T1r + T1s;
  const TU = Tc + Th;
  const TV = Tn + Ts;
  const TW = TU + TV;
  const T16 = T14 - T15;
  const T19 = T17 - T18;
  const T1y = T16 + T19;

  const outCr = new Float64Array(10), outCi = new Float64Array(10);

  {
    const T11 = KP559016994 * (Tu - TR);
    const TS = Tu + TR;
    const T12 = T7 - KP250000000 * TS; // FNMS
    const T1a = T16 - T19;
    const T1h = T1d + T1g;
    const T1i = KP951056516 * T1a + KP587785252 * T1h; // FMA
    const T1k = KP951056516 * T1h - KP587785252 * T1a; // FNMS
    outCi[4] = T7 + TS;
    const T1j = T12 - T11;
    outCi[2] = T1j - T1k;
    outCr[3] = T1j + T1k;
    const T13 = T11 + T12;
    outCi[0] = T13 - T1i;
    outCr[1] = T13 + T1i;
  }
  {
    const T1m = KP559016994 * (TW - TZ);
    const T10 = TW + TZ;
    const T1l = TT - KP250000000 * T10; // FNMS
    const T1q = T1o - T1p;
    const T1t = T1r - T1s;
    const T1u = KP951056516 * T1q - KP587785252 * T1t; // FNMS
    const T1w = KP951056516 * T1t + KP587785252 * T1q; // FMA
    outCr[0] = TT + T10;
    const T1v = T1m + T1l;
    outCr[4] = T1v - T1w;
    outCi[3] = T1v + T1w;
    const T1n = T1l - T1m;
    outCr[2] = T1n - T1u;
    outCi[1] = T1n + T1u;
  }
  {
    const T1H = KP559016994 * (T1y + T1x);
    const T1z = T1x - T1y;
    const T1G = KP250000000 * T1z + T1C; // FMA
    const T1D = Ti - Tt;
    const T1E = TQ - TF;
    const T1F = KP587785252 * T1D + KP951056516 * T1E; // FMA
    const T1J = KP587785252 * T1E - KP951056516 * T1D; // FNMS
    outCr[5] = T1z - T1C;
    const T1K = T1H + T1G;
    outCr[9] = T1J - T1K;
    outCi[8] = T1J + T1K;
    const T1I = T1G - T1H;
    outCr[7] = T1F - T1I;
    outCi[6] = T1F + T1I;
  }
  {
    const T1Q = KP559016994 * (T1O - T1P);
    const T1S = T1O + T1P;
    const T1T = T1R - KP250000000 * T1S; // FNMS
    const T1L = TU - TV;
    const T1M = TX - TY;
    const T1N = KP951056516 * T1L + KP587785252 * T1M; // FMA
    const T1V = KP951056516 * T1M - KP587785252 * T1L; // FNMS
    outCi[9] = T1S + T1R;
    const T1W = T1T - T1Q;
    outCr[8] = T1V - T1W;
    outCi[7] = T1V + T1W;
    const T1U = T1Q + T1T;
    outCr[6] = T1N - T1U;
    outCi[5] = T1N + T1U;
  }

  return [outCr, outCi];
}

module.exports = { hf_10 };
