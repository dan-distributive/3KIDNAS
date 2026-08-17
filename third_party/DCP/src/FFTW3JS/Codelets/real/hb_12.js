'use strict';

// =============================================================================
// hb_12.js -- faithful JS port of rdft/scalar/r2cb/hb_12.c (non-FMA
// branch). Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hb_6.js --
// the PLAIN (non-alt-codegen) backward family, BACKWARD counterpart of
// hf_12.js. Input cr[]/ci[] physical-index pairing transcribed literally
// from the C source's WS(rs,k) reads, not "cleaned up" to a uniform
// pattern (see hb_4.js's header for why that matters).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function hb_12(cr, ci, Wc, Ws) {
  const T1 = cr[0];
  const TD = ci[11];
  const T2 = cr[4], T3 = ci[3];
  const T4 = T2 + T3;
  const T1g = KP866025403 * (T2 - T3);
  const TE = ci[7], TF = cr[8];
  const TG = TE - TF;
  const T11 = KP866025403 * (TE + TF);
  const T5 = T1 + T4;
  const TH = TD + TG;
  const T10 = T1 - KP500000000 * T4; // FNMS
  const T12 = T10 - T11;
  const T1M = T10 + T11;
  const T1h = TD - KP500000000 * TG; // FNMS
  const T1i = T1g + T1h;
  const T1U = T1h - T1g;

  const Tc = cr[3];
  const Tp = ci[8];
  const Td = ci[4], Te = ci[0];
  const Tf = Td + Te;
  const T17 = KP866025403 * (Td - Te);
  const Tq = cr[7], Tr = cr[11];
  const Ts = Tq + Tr;
  const T1o = KP866025403 * (Tq - Tr);
  const Tg = Tc + Tf;
  const Tt = Tp - Ts;
  const T18 = KP500000000 * Ts + Tp; // FMA
  const T19 = T17 + T18;
  const T1X = T18 - T17;
  const T1n = Tc - KP500000000 * Tf; // FNMS
  const T1p = T1n + T1o;
  const T1P = T1n - T1o;

  const T6 = ci[5];
  const TL = cr[6];
  const T7 = ci[1], T8 = cr[2];
  const T9 = T7 + T8;
  const T1j = KP866025403 * (T7 - T8);
  const TI = ci[9], TJ = cr[10];
  const TK = TI - TJ;
  const T14 = KP866025403 * (TI + TJ);
  const Ta = T6 + T9;
  const TM = TK - TL;
  const T13 = T6 - KP500000000 * T9; // FNMS
  const T15 = T13 + T14;
  const T1N = T13 - T14;
  const T1k = KP500000000 * TK + TL; // FMA
  const T1l = T1j - T1k;
  const T1V = T1j + T1k;

  const Th = ci[2];
  const Tx = cr[9];
  const Ti = cr[1], Tj = cr[5];
  const Tk = Ti + Tj;
  const T1a = KP866025403 * (Ti - Tj);
  const Tu = ci[10], Tv = ci[6];
  const Tw = Tu + Tv;
  const T1r = KP866025403 * (Tv - Tu);
  const Tl = Th + Tk;
  const Ty = Tw - Tx;
  const T1b = KP500000000 * Tw + Tx; // FMA
  const T1c = T1a - T1b;
  const T1Y = T1a + T1b;
  const T1q = Th - KP500000000 * Tk; // FNMS
  const T1s = T1q + T1r;
  const T1Q = T1q - T1r;

  const outCr = new Float64Array(12), outCi = new Float64Array(12);

  {
    const Tb = T5 + Ta;
    const Tm = Tg + Tl;
    const TU = Tb - Tm;
    const TW = TH + TM;
    const TX = Tt + Ty;
    const TY = TW - TX;
    outCr[0] = Tb + Tm;
    outCi[0] = TW + TX;
    const TT = Wc[6], TV = Ws[6];
    outCr[6] = TT * TU - TV * TY; // FNMS
    outCi[6] = TV * TU + TT * TY; // FMA
  }
  {
    const To = T5 - Ta;
    const Tz = Tt - Ty;
    const TA = To - Tz;
    const TQ = To + Tz;
    const TC = Tg - Tl;
    const TN = TH - TM;
    const TO = TC + TN;
    const TS = TN - TC;

    const Tn = Wc[9], TB = Ws[9];
    outCr[9] = Tn * TA - TB * TO; // FNMS
    outCi[9] = Tn * TO + TB * TA; // FMA
    const TP = Wc[3], TR = Ws[3];
    outCr[3] = TP * TQ - TR * TS; // FNMS
    outCi[3] = TP * TS + TR * TQ; // FMA
  }
  {
    const T26 = T1M - T1N;
    const T27 = T1X + T1Y;
    const T28 = T26 - T27;
    const T2e = T26 + T27;
    const T2a = T1U + T1V;
    const T2b = T1P - T1Q;
    const T2c = T2a + T2b;
    const T2g = T2a - T2b;

    const T25 = Wc[5], T29 = Ws[5];
    outCr[5] = T25 * T28 - T29 * T2c; // FNMS
    outCi[5] = T25 * T2c + T29 * T28; // FMA
    const T2d = Wc[11], T2f = Ws[11];
    outCr[11] = T2d * T2e - T2f * T2g; // FNMS
    outCi[11] = T2d * T2g + T2f * T2e; // FMA
  }
  {
    const T1O = T1M + T1N;
    const T1R = T1P + T1Q;
    const T1S = T1O - T1R;
    const T22 = T1O + T1R;
    const T1W = T1U - T1V;
    const T1Z = T1X - T1Y;
    const T20 = T1W - T1Z;
    const T24 = T1W + T1Z;

    const T1L = Wc[2], T1T = Ws[2];
    outCr[2] = T1L * T1S - T1T * T20; // FNMS
    outCi[2] = T1T * T1S + T1L * T20; // FMA
    const T21 = Wc[8], T23 = Ws[8];
    outCr[8] = T21 * T22 - T23 * T24; // FNMS
    outCi[8] = T23 * T22 + T21 * T24; // FMA
  }
  {
    const T1A = T12 + T15;
    const T1B = T1p + T1s;
    const T1C = T1A - T1B;
    const T1I = T1A + T1B;
    const T1E = T1i + T1l;
    const T1F = T19 + T1c;
    const T1G = T1E - T1F;
    const T1K = T1E + T1F;

    const T1z = Wc[10], T1D = Ws[10];
    outCr[10] = T1z * T1C - T1D * T1G; // FNMS
    outCi[10] = T1D * T1C + T1z * T1G; // FMA
    const T1H = Wc[4], T1J = Ws[4];
    outCr[4] = T1H * T1I - T1J * T1K; // FNMS
    outCi[4] = T1J * T1I + T1H * T1K; // FMA
  }
  {
    const T16 = T12 - T15;
    const T1d = T19 - T1c;
    const T1e = T16 - T1d;
    const T1w = T16 + T1d;
    const T1m = T1i - T1l;
    const T1t = T1p - T1s;
    const T1u = T1m + T1t;
    const T1y = T1m - T1t;

    const TZ = Wc[1], T1f = Ws[1];
    outCr[1] = TZ * T1e - T1f * T1u; // FNMS
    outCi[1] = TZ * T1u + T1f * T1e; // FMA
    const T1v = Wc[7], T1x = Ws[7];
    outCr[7] = T1v * T1w - T1x * T1y; // FNMS
    outCi[7] = T1v * T1y + T1x * T1w; // FMA
  }

  return [outCr, outCi];
}

module.exports = { hb_12 };
