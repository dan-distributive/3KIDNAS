'use strict';

// =============================================================================
// hf2_16.js -- faithful JS port of rdft/scalar/r2cf/hf2_16.c (non-FMA
// branch). Alternate-codegen EVEN-radix (r=16) twiddle codelet -- same role
// as hf2_8.js, larger radix. twinstr trig-generates W^1, W^3, W^9, W^15
// (four raw pairs); the rest are derived via complex products, exactly
// matching the C source's operation order. Same (cr,ci,Wc,Ws) ->
// (outCr,outCi) convention as hf2_8.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function hf2_16(cr, ci, Wc, Ws) {
  const T2 = Wc[1], T5 = Ws[1], Tg = Wc[3], Ti = Ws[3];
  const Th = T2 * Tg, Tn = T5 * Tg, Tj = T5 * Ti, Tm = T2 * Ti;
  const Tk = Th - Tj, To = Tm + Tn, TE = Tm - Tn, TC = Th + Tj;
  const T3 = Wc[9], T6 = Ws[9];
  const T7 = T5 * T6, Tv = Tg * T6, Ta = T2 * T6, Ts = Ti * T6;
  const T4 = T2 * T3, Tw = Ti * T3, Tb = T5 * T3, Tr = Tg * T3;

  const T8 = T4 + T7;
  const TW = Tv - Tw;
  const TJ = Ta + Tb;
  const Tt = Tr - Ts;
  const TU = Tr + Ts;
  const Tc = Ta - Tb;
  const Tx = Tv + Tw;
  const TH = T4 - T7;
  const TN = Wc[15], TO = Ws[15];
  const TP = T2 * TN + T5 * TO; // FMA
  const TR = T2 * TO - T5 * TN; // FNMS

  const T1d = Tk * T6, T1e = To * T3;
  const T1f = T1d - T1e, T1k = T1d + T1e;
  const T19 = Tk * T3, T1a = To * T6;
  const T1b = T19 + T1a, T1i = T19 - T1a;

  const T1w = TC * T6, T1x = TE * T3;
  const T1y = T1w - T1x, T1H = T1w + T1x;
  const T1s = TC * T3, T1t = TE * T6;
  const T1u = T1s + T1t, T1F = T1s - T1t;

  const T1 = cr[0], T3d = ci[0];
  const T9 = cr[8], Td = ci[8];
  const Te = T8 * T9 + Tc * Td; // FMA
  const T3c = T8 * Td - Tc * T9; // FNMS
  const Tf = T1 + Te;
  const T3s = T3d - T3c;
  const T1N = T1 - Te;
  const T3e = T3c + T3d;

  const Tl = cr[4], Tp = ci[4];
  const Tq = Tk * Tl + To * Tp; // FMA
  const T1O = Tk * Tp - To * Tl; // FNMS
  const Tu = cr[12], Ty = ci[12];
  const Tz = Tt * Tu + Tx * Ty; // FMA
  const T1P = Tt * Ty - Tx * Tu; // FNMS

  const TA = Tq + Tz;
  const T3r = Tq - Tz;
  const T1Q = T1O - T1P;
  const T3b = T1O + T1P;

  const TD = cr[2], TF = ci[2];
  const TG = TC * TD + TE * TF; // FMA
  const T1T = TC * TF - TE * TD; // FNMS
  const TI = cr[10], TK = ci[10];
  const TL = TH * TI + TJ * TK; // FMA
  const T1U = TH * TK - TJ * TI; // FNMS

  const TM = TG + TL;
  const T2N = T1T + T1U;
  const T1S = TG - TL;
  const T1V = T1T - T1U;
  const T1W = T1S - T1V;
  const T2w = T1S + T1V;

  const TQ = cr[14], TS = ci[14];
  const TT = TP * TQ + TR * TS; // FMA
  const T1Y = TP * TS - TR * TQ; // FNMS
  const TV = cr[6], TX = ci[6];
  const TY = TU * TV + TW * TX; // FMA
  const T1Z = TU * TX - TW * TV; // FNMS

  const TZ = TT + TY;
  const T2M = T1Y + T1Z;
  const T1X = TT - TY;
  const T20 = T1Y - T1Z;
  const T21 = T1X + T20;
  const T2x = T1X - T20;

  const T1p = cr[15], T1q = ci[15];
  const T1r = TN * T1p + TO * T1q; // FMA
  const T2f = TN * T1q - TO * T1p; // FNMS
  const T1G = cr[11], T1I = ci[11];
  const T1J = T1F * T1G + T1H * T1I; // FMA
  const T2m = T1F * T1I - T1H * T1G; // FNMS

  const T1v = cr[7], T1z = ci[7];
  const T1A = T1u * T1v + T1y * T1z; // FMA
  const T2g = T1u * T1z - T1y * T1v; // FNMS
  const T1C = cr[3], T1D = ci[3];
  const T1E = Tg * T1C + Ti * T1D; // FMA
  const T2l = Tg * T1D - Ti * T1C; // FNMS

  const T1B = T1r + T1A;
  const T1K = T1E + T1J;
  const T2V = T1B - T1K;
  const T2W = T2f + T2g;
  const T2X = T2l + T2m;
  const T2Y = T2W - T2X;

  const T2h = T2f - T2g;
  const T2i = T1E - T1J;
  const T2j = T2h + T2i;
  const T2E = T2h - T2i;
  const T2k = T1r - T1A;
  const T2n = T2l - T2m;
  const T2o = T2k - T2n;
  const T2D = T2k + T2n;

  const T12 = cr[1], T13 = ci[1];
  const T14 = T2 * T12 + T5 * T13; // FMA
  const T29 = T2 * T13 - T5 * T12; // FNMS
  const T1j = cr[13], T1l = ci[13];
  const T1m = T1i * T1j + T1k * T1l; // FMA
  const T26 = T1i * T1l - T1k * T1j; // FNMS

  const T15 = cr[9], T16 = ci[9];
  const T17 = T3 * T15 + T6 * T16; // FMA
  const T2a = T3 * T16 - T6 * T15; // FNMS
  const T1c = cr[5], T1g = ci[5];
  const T1h = T1b * T1c + T1f * T1g; // FMA
  const T25 = T1b * T1g - T1f * T1c; // FNMS

  const T18 = T14 + T17;
  const T1n = T1h + T1m;
  const T2Q = T18 - T1n;
  const T2R = T29 + T2a;
  const T2S = T25 + T26;
  const T2T = T2R - T2S;

  const T24 = T14 - T17;
  const T27 = T25 - T26;
  const T28 = T24 - T27;
  const T2B = T24 + T27;
  const T2b = T29 - T2a;
  const T2c = T1h - T1m;
  const T2d = T2b + T2c;
  const T2A = T2b - T2c;

  const outCr = new Float64Array(16), outCi = new Float64Array(16);

  {
    const T1R = T1N - T1Q;
    const T22 = KP707106781 * (T1W + T21);
    const T23 = T1R + T22;
    const T2r = T1R - T22;
    const T3q = KP707106781 * (T2w - T2x);
    const T3t = T3r + T3s;
    const T3u = T3q + T3t;
    const T3w = T3t - T3q;

    const T2e = KP923879532 * T28 - KP382683432 * T2d; // FNMS
    const T2p = KP382683432 * T2j + KP923879532 * T2o; // FMA
    const T2q = T2e + T2p;
    const T3v = T2p - T2e;
    const T2s = KP923879532 * T2d + KP382683432 * T28; // FMA
    const T2t = KP382683432 * T2o - KP923879532 * T2j; // FNMS
    const T2u = T2s + T2t;
    const T3p = T2t - T2s;

    outCr[7] = T23 - T2q;
    outCr[11] = T3v - T3w;
    outCi[12] = T3v + T3w;
    outCi[0] = T23 + T2q;
    outCi[4] = T2r - T2u;
    outCr[15] = T3p - T3u;
    outCi[8] = T3p + T3u;
    outCr[3] = T2r + T2u;
  }
  {
    const TB = Tf + TA;
    const T10 = TM + TZ;
    const T11 = TB + T10;
    const T35 = TB - T10;
    const T3a = T2N + T2M;
    const T3f = T3b + T3e;
    const T3g = T3a + T3f;
    const T3i = T3f - T3a;

    const T1o = T18 + T1n;
    const T1L = T1B + T1K;
    const T1M = T1o + T1L;
    const T3h = T1L - T1o;
    const T36 = T2W + T2X;
    const T37 = T2R + T2S;
    const T38 = T36 - T37;
    const T39 = T37 + T36;

    outCi[7] = T11 - T1M;
    outCr[12] = T3h - T3i;
    outCi[11] = T3h + T3i;
    outCr[0] = T11 + T1M;
    outCr[4] = T35 - T38;
    outCr[8] = T39 - T3g;
    outCi[15] = T39 + T3g;
    outCi[3] = T35 + T38;
  }
  {
    const T2v = T1N + T1Q;
    const T2y = KP707106781 * (T2w + T2x);
    const T2z = T2v + T2y;
    const T2H = T2v - T2y;
    const T3y = KP707106781 * (T21 - T1W);
    const T3z = T3s - T3r;
    const T3A = T3y + T3z;
    const T3C = T3z - T3y;

    const T2C = KP382683432 * T2A + KP923879532 * T2B; // FMA
    const T2F = KP923879532 * T2D - KP382683432 * T2E; // FNMS
    const T2G = T2C + T2F;
    const T3B = T2F - T2C;
    const T2I = KP382683432 * T2B - KP923879532 * T2A; // FNMS
    const T2J = KP923879532 * T2E + KP382683432 * T2D; // FMA
    const T2K = T2I + T2J;
    const T3x = T2J - T2I;

    outCi[6] = T2z - T2G;
    outCr[13] = T3B - T3C;
    outCi[10] = T3B + T3C;
    outCr[1] = T2z + T2G;
    outCr[5] = T2H - T2K;
    outCr[9] = T3x - T3A;
    outCi[14] = T3x + T3A;
    outCi[2] = T2H + T2K;
  }
  {
    const T2L = Tf - TA;
    const T2O = T2M - T2N;
    const T2P = T2L - T2O;
    const T31 = T2L + T2O;
    const T3k = TM - TZ;
    const T3l = T3e - T3b;
    const T3m = T3k + T3l;
    const T3o = T3l - T3k;

    const T2U = T2Q + T2T;
    const T2Z = T2V - T2Y;
    const T30 = KP707106781 * (T2U + T2Z);
    const T3j = KP707106781 * (T2Z - T2U);
    const T32 = T2Q - T2T;
    const T33 = T2V + T2Y;
    const T34 = KP707106781 * (T32 + T33);
    const T3n = KP707106781 * (T33 - T32);

    outCi[5] = T2P - T30;
    outCr[10] = T3n - T3o;
    outCi[13] = T3n + T3o;
    outCr[2] = T2P + T30;
    outCr[6] = T31 - T34;
    outCr[14] = T3j - T3m;
    outCi[9] = T3j + T3m;
    outCi[1] = T31 + T34;
  }

  return [outCr, outCi];
}

module.exports = { hf2_16 };
