'use strict';

// =============================================================================
// t2_16.js -- faithful JS port of dft/scalar/codelets/t2_16.c (non-FMA
// branch), FFTW3's "twiddle-log3 / precompute-twiddles" radix-16 twiddle
// codelet. twinstr only trig-generates W^1, W^3, W^9, W^15 (four raw
// pairs); the other 11 needed multiples (W^2,4,5,6,7,8,10,11,12,13,14) are
// DERIVED via one complex multiply each (same trick as t2_5/t2_8/t2_10.js):
//   (TC,TE)=W^2  (Tk,To)=W^4   (T1b,T1f)=W^5  (TU,TW)=W^6
//   (T1u,T1y)=W^7 (T8,Tc)=W^8  (TH,TJ)=W^10   (T1F,T1H)=W^11
//   (Tt,Tx)=W^12 (T1i,T1k)=W^13 (TP,TR)=W^14
// Same calling convention as every other twiddle codelet here (br/bi/Wc/Ws
// with Composite1D.js's full r-1 = 15 pair table already built) -- only
// indices 1, 3, 9, 15 are actually read.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function t2_16(br, bi, Wc, Ws) {
  const T2 = Wc[1], T5 = Ws[1];
  const Tg = Wc[3], Ti = Ws[3];
  const Th = T2 * Tg, Tn = T5 * Tg, Tj = T5 * Ti, Tm = T2 * Ti;
  const Tk = Th - Tj, To = Tm + Tn, TE = Tm - Tn, TC = Th + Tj;
  const T6 = Ws[9];
  const T7 = T5 * T6, Tv = Tg * T6, Ta = T2 * T6, Ts = Ti * T6;
  const T3 = Wc[9];
  const T4 = T2 * T3, Tw = Ti * T3, Tb = T5 * T3, Tr = Tg * T3;

  const T8 = T4 + T7, TW = Tv - Tw, TJ = Ta + Tb, Tt = Tr - Ts;
  const TU = Tr + Ts, Tc = Ta - Tb, Tx = Tv + Tw, TH = T4 - T7;
  const TN = Wc[15], TO = Ws[15];
  const TP = T2 * TN + T5 * TO;
  const TR = T2 * TO - T5 * TN;

  const T1d = Tk * T6, T1e = To * T3, T1f = T1d - T1e, T1k = T1d + T1e;
  const T19 = Tk * T3, T1a = To * T6, T1b = T19 + T1a, T1i = T19 - T1a;

  const T1w = TC * T6, T1x = TE * T3, T1y = T1w - T1x, T1H = T1w + T1x;
  const T1s = TC * T3, T1t = TE * T6, T1u = T1s + T1t, T1F = T1s - T1t;

  const T1 = br[0], T3d = bi[0];
  const T9 = br[8], Td = bi[8];
  const Te = T8 * T9 + Tc * Td;
  const T3c = T8 * Td - Tc * T9;
  const Tf = T1 + Te;
  const T3r = T3d - T3c;
  const T1N = T1 - Te;
  const T3e = T3c + T3d;

  const Tl = br[4], Tp = bi[4];
  const Tq = Tk * Tl + To * Tp;
  const T1O = Tk * Tp - To * Tl;
  const Tu = br[12], Ty = bi[12];
  const Tz = Tt * Tu + Tx * Ty;
  const T1P = Tt * Ty - Tx * Tu;

  const TA = Tq + Tz;
  const T3s = Tq - Tz;
  const T1Q = T1O - T1P;
  const T3b = T1O + T1P;

  const TD = br[2], TF = bi[2];
  const TG = TC * TD + TE * TF;
  const T1S = TC * TF - TE * TD;
  const TI = br[10], TK = bi[10];
  const TL = TH * TI + TJ * TK;
  const T1T = TH * TK - TJ * TI;

  const TM = TG + TL;
  const T2M = T1S + T1T;
  const T1U = T1S - T1T;
  const T1V = TG - TL;
  const T1W = T1U - T1V;
  const T2w = T1V + T1U;

  const TQ = br[14], TS = bi[14];
  const TT = TP * TQ + TR * TS;
  const T1Y = TP * TS - TR * TQ;
  const TV = br[6], TX = bi[6];
  const TY = TU * TV + TW * TX;
  const T1Z = TU * TX - TW * TV;

  const TZ = TT + TY;
  const T2N = T1Y + T1Z;
  const T1X = TT - TY;
  const T20 = T1Y - T1Z;
  const T21 = T1X + T20;
  const T2x = T1X - T20;

  const T1p = br[15], T1q = bi[15];
  const T1r = TN * T1p + TO * T1q;
  const T2k = TN * T1q - TO * T1p;
  const T1G = br[11], T1I = bi[11];
  const T1J = T1F * T1G + T1H * T1I;
  const T2h = T1F * T1I - T1H * T1G;

  const T1v = br[7], T1z = bi[7];
  const T1A = T1u * T1v + T1y * T1z;
  const T2l = T1u * T1z - T1y * T1v;
  const T1C = br[3], T1D = bi[3];
  const T1E = Tg * T1C + Ti * T1D;
  const T2g = Tg * T1D - Ti * T1C;

  const T1B = T1r + T1A;
  const T1K = T1E + T1J;
  const T2V = T1B - T1K;
  const T2W = T2k + T2l;
  const T2X = T2g + T2h;
  const T2Y = T2W - T2X;

  const T2f = T1r - T1A;
  const T2i = T2g - T2h;
  const T2j = T2f - T2i;
  const T2D = T2f + T2i;
  const T2m = T2k - T2l;
  const T2n = T1E - T1J;
  const T2o = T2m + T2n;
  const T2E = T2m - T2n;

  const T12 = br[1], T13 = bi[1];
  const T14 = T2 * T12 + T5 * T13;
  const T24 = T2 * T13 - T5 * T12;
  const T1j = br[13], T1l = bi[13];
  const T1m = T1i * T1j + T1k * T1l;
  const T2b = T1i * T1l - T1k * T1j;

  const T15 = br[9], T16 = bi[9];
  const T17 = T3 * T15 + T6 * T16;
  const T25 = T3 * T16 - T6 * T15;
  const T1c = br[5], T1g = bi[5];
  const T1h = T1b * T1c + T1f * T1g;
  const T2a = T1b * T1g - T1f * T1c;

  const T18 = T14 + T17;
  const T1n = T1h + T1m;
  const T2Q = T18 - T1n;
  const T2R = T24 + T25;
  const T2S = T2a + T2b;
  const T2T = T2R - T2S;

  const T26 = T24 - T25;
  const T27 = T1h - T1m;
  const T28 = T26 + T27;
  const T2A = T26 - T27;
  const T29 = T14 - T17;
  const T2c = T2a - T2b;
  const T2d = T29 - T2c;
  const T2B = T29 + T2c;

  const outR = new Float64Array(16), outI = new Float64Array(16);

  {
    const T1R = T1N - T1Q;
    const T22 = KP707106781 * (T1W - T21);
    const T23 = T1R + T22;
    const T2r = T1R - T22;
    const T3y = KP707106781 * (T2x - T2w);
    const T3z = T3s + T3r;
    const T3A = T3y + T3z;
    const T3C = T3z - T3y;

    const T2e = KP923879532 * T28 + KP382683432 * T2d;
    const T2p = KP382683432 * T2j - KP923879532 * T2o;
    const T2q = T2e + T2p;
    const T3B = T2p - T2e;
    const T2s = KP382683432 * T28 - KP923879532 * T2d;
    const T2t = KP382683432 * T2o + KP923879532 * T2j;
    const T2u = T2s - T2t;
    const T3x = T2s + T2t;

    outR[11] = T23 - T2q;
    outI[11] = T3A - T3x;
    outR[3] = T23 + T2q;
    outI[3] = T3x + T3A;
    outR[15] = T2r - T2u;
    outI[15] = T3C - T3B;
    outR[7] = T2r + T2u;
    outI[7] = T3B + T3C;
  }

  {
    const T2L = Tf - TA;
    const T2O = T2M - T2N;
    const T2P = T2L + T2O;
    const T31 = T2L - T2O;
    const T3k = TZ - TM;
    const T3l = T3e - T3b;
    const T3m = T3k + T3l;
    const T3o = T3l - T3k;

    const T2U = T2Q + T2T;
    const T2Z = T2V - T2Y;
    const T30 = KP707106781 * (T2U + T2Z);
    const T3n = KP707106781 * (T2Z - T2U);
    const T32 = T2T - T2Q;
    const T33 = T2V + T2Y;
    const T34 = KP707106781 * (T32 - T33);
    const T3j = KP707106781 * (T32 + T33);

    outR[10] = T2P - T30;
    outI[10] = T3m - T3j;
    outR[2] = T2P + T30;
    outI[2] = T3j + T3m;
    outR[14] = T31 - T34;
    outI[14] = T3o - T3n;
    outR[6] = T31 + T34;
    outI[6] = T3n + T3o;
  }

  {
    const T2v = T1N + T1Q;
    const T2y = KP707106781 * (T2w + T2x);
    const T2z = T2v + T2y;
    const T2H = T2v - T2y;
    const T3q = KP707106781 * (T1W + T21);
    const T3t = T3r - T3s;
    const T3u = T3q + T3t;
    const T3w = T3t - T3q;

    const T2C = KP382683432 * T2A + KP923879532 * T2B;
    const T2F = KP923879532 * T2D - KP382683432 * T2E;
    const T2G = T2C + T2F;
    const T3v = T2F - T2C;
    const T2I = KP923879532 * T2A - KP382683432 * T2B;
    const T2J = KP923879532 * T2E + KP382683432 * T2D;
    const T2K = T2I - T2J;
    const T3p = T2I + T2J;

    outR[9] = T2z - T2G;
    outI[9] = T3u - T3p;
    outR[1] = T2z + T2G;
    outI[1] = T3p + T3u;
    outR[13] = T2H - T2K;
    outI[13] = T3w - T3v;
    outR[5] = T2H + T2K;
    outI[5] = T3v + T3w;
  }

  {
    const TB = Tf + TA;
    const T10 = TM + TZ;
    const T11 = TB + T10;
    const T35 = TB - T10;
    const T3a = T2M + T2N;
    const T3f = T3b + T3e;
    const T3g = T3a + T3f;
    const T3i = T3f - T3a;

    const T1o = T18 + T1n;
    const T1L = T1B + T1K;
    const T1M = T1o + T1L;
    const T3h = T1L - T1o;
    const T36 = T2R + T2S;
    const T37 = T2W + T2X;
    const T38 = T36 - T37;
    const T39 = T36 + T37;

    outR[8] = T11 - T1M;
    outI[8] = T3g - T39;
    outR[0] = T11 + T1M;
    outI[0] = T39 + T3g;
    outR[12] = T35 - T38;
    outI[12] = T3i - T3h;
    outR[4] = T35 + T38;
    outI[4] = T3h + T3i;
  }

  return [outR, outI];
}

module.exports = { t2_16 };
