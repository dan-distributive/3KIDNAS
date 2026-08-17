'use strict';

// =============================================================================
// n1_16.js -- faithful JS port of dft/scalar/codelets/n1_16.c (non-FMA
// branch), FFTW3's direct (base-case) radix-16 complex DFT codelet.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function n1_16(ri, ii) {
  const T1 = ri[0], T2 = ri[8];
  const T3 = T1 + T2, TL = T1 - T2;
  const Tw = ii[0], Tx = ii[8];
  const Ty = Tw + Tx, T1k = Tw - Tx;

  const T4 = ri[4], T5 = ri[12];
  const T6 = T4 + T5, T1j = T4 - T5;
  const Tz = ii[4], TA = ii[12];
  const TB = Tz + TA, TM = Tz - TA;

  const T7 = T3 + T6, T1R = T3 - T6, T25 = Ty - TB, TC = Ty + TB;
  const TN = TL - TM, T1x = TL + TM, T1H = T1k - T1j, T1l = T1j + T1k;

  const Tn = ri[15], To = ri[7];
  const Tp = Tn + To, T17 = Tn - To;
  const T1d = ii[15], T1e = ii[7];
  const T1f = T1d - T1e, T20 = T1d + T1e;

  const Tq = ri[3], Tr = ri[11];
  const Ts = Tq + Tr, T1c = Tq - Tr;
  const T18 = ii[3], T19 = ii[11];
  const T1a = T18 - T19, T21 = T18 + T19;

  const Tt = Tp + Ts, T22 = T20 - T21, T2h = T20 + T21;
  const T1b = T17 - T1a, T1g = T1c + T1f, T1E = T1f - T1c, T1Z = Tp - Ts, T1D = T17 + T1a;

  const T8 = ri[2], T9 = ri[10];
  const Ta = T8 + T9, TP = T8 - T9;
  const TD = ii[2], TE = ii[10];
  const TF = TD + TE, TO = TD - TE;

  const Tb = ri[14], Tc = ri[6];
  const Td = Tb + Tc, TR = Tb - Tc;
  const TG = ii[14], TH = ii[6];
  const TI = TG + TH, TS = TG - TH;

  const Te = Ta + Td, T1S = TF - TI, T26 = Td - Ta, TJ = TF + TI;
  const TQ = TO - TP, T1m = TR - TS, T1n = TP + TO, TT = TR + TS;

  const Tg = ri[1], Th = ri[9];
  const Ti = Tg + Th, T11 = Tg - Th;
  const TX = ii[1], TY = ii[9];
  const TZ = TX - TY, T1V = TX + TY;

  const Tj = ri[5], Tk = ri[13];
  const Tl = Tj + Tk, TW = Tj - Tk;
  const T12 = ii[5], T13 = ii[13];
  const T14 = T12 - T13, T1W = T12 + T13;

  const Tm = Ti + Tl, T1X = T1V - T1W, T2g = T1V + T1W;
  const T10 = TW + TZ, T15 = T11 - T14, T1B = T11 + T14, T1U = Ti - Tl, T1A = TZ - TW;

  const ro = new Float64Array(16), io = new Float64Array(16);

  {
    const Tf = T7 + Te, Tu = Tm + Tt;
    ro[8] = Tf - Tu;
    ro[0] = Tf + Tu;
    const T2j = TC + TJ, T2k = T2g + T2h;
    io[8] = T2j - T2k;
    io[0] = T2j + T2k;
  }

  {
    const Tv = Tt - Tm, TK = TC - TJ;
    io[4] = Tv + TK;
    io[12] = TK - Tv;
    const T2f = T7 - Te, T2i = T2g - T2h;
    ro[12] = T2f - T2i;
    ro[4] = T2f + T2i;
  }

  {
    const T1T = T1R + T1S, T27 = T25 - T26;
    const T1Y = T1U + T1X, T23 = T1Z - T22;
    const T24 = KP707106781 * (T1Y + T23);
    const T28 = KP707106781 * (T23 - T1Y);
    ro[10] = T1T - T24;
    io[6] = T27 + T28;
    ro[2] = T1T + T24;
    io[14] = T27 - T28;
  }

  {
    const T29 = T1R - T1S, T2d = T26 + T25;
    const T2a = T1X - T1U, T2b = T1Z + T22;
    const T2c = KP707106781 * (T2a - T2b);
    const T2e = KP707106781 * (T2a + T2b);
    ro[14] = T29 - T2c;
    io[2] = T2d + T2e;
    ro[6] = T29 + T2c;
    io[10] = T2d - T2e;
  }

  {
    const TU = KP707106781 * (TQ - TT);
    const TV = TN + TU;
    const T1r = TN - TU;
    const T1o = KP707106781 * (T1m - T1n);
    const T1p = T1l - T1o;
    const T1v = T1l + T1o;

    const T16 = KP923879532 * T10 + KP382683432 * T15;
    const T1h = KP382683432 * T1b - KP923879532 * T1g;
    const T1i = T16 + T1h;
    const T1q = T1h - T16;
    const T1s = KP382683432 * T10 - KP923879532 * T15;
    const T1t = KP382683432 * T1g + KP923879532 * T1b;
    const T1u = T1s - T1t;
    const T1w = T1s + T1t;

    ro[11] = TV - T1i;
    io[11] = T1v - T1w;
    ro[3] = TV + T1i;
    io[3] = T1v + T1w;
    io[15] = T1p - T1q;
    ro[15] = T1r - T1u;
    io[7] = T1p + T1q;
    ro[7] = T1r + T1u;
  }

  {
    const T1y = KP707106781 * (T1n + T1m);
    const T1z = T1x + T1y;
    const T1L = T1x - T1y;
    const T1I = KP707106781 * (TQ + TT);
    const T1J = T1H - T1I;
    const T1P = T1H + T1I;

    const T1C = KP923879532 * T1B + KP382683432 * T1A;
    const T1F = KP923879532 * T1D - KP382683432 * T1E;
    const T1G = T1C + T1F;
    const T1K = T1F - T1C;
    const T1M = KP923879532 * T1A - KP382683432 * T1B;
    const T1N = KP923879532 * T1E + KP382683432 * T1D;
    const T1O = T1M - T1N;
    const T1Q = T1M + T1N;

    ro[9] = T1z - T1G;
    io[9] = T1P - T1Q;
    ro[1] = T1z + T1G;
    io[1] = T1P + T1Q;
    io[13] = T1J - T1K;
    ro[13] = T1L - T1O;
    io[5] = T1J + T1K;
    ro[5] = T1L + T1O;
  }

  return [ro, io];
}

module.exports = { n1_16 };
