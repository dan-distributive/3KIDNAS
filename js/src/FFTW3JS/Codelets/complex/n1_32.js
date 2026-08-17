'use strict';

// =============================================================================
// n1_32.js -- faithful JS port of dft/scalar/codelets/n1_32.c (non-FMA
// branch), FFTW3's direct (base-case) radix-32 (split-radix-style) complex
// DFT codelet.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP831469612 = 0.831469612302545237078788377617905756738560812;
const KP555570233 = 0.555570233019602224742830813948532874374937191;
const KP195090322 = 0.195090322016128267848284868477022240927691618;
const KP980785280 = 0.980785280403230449126182236134239036973933731;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function n1_32(ri, ii) {
  const T1 = ri[0], T2 = ri[16];
  const T3 = T1 + T2, T1x = T1 - T2;
  const T12 = ii[0], T13 = ii[16];
  const T14 = T12 + T13, T2S = T12 - T13;
  const T4 = ri[8], T5 = ri[24];
  const T6 = T4 + T5, T2R = T4 - T5;
  const T15 = ii[8], T16 = ii[24];
  const T17 = T15 + T16, T1y = T15 - T16;
  const T7 = T3 + T6, T4r = T3 - T6, T4Z = T14 - T17, T18 = T14 + T17;
  const T1z = T1x - T1y, T3t = T1x + T1y, T3T = T2S - T2R, T2T = T2R + T2S;

  const T8 = ri[4], T9 = ri[20];
  const Ta = T8 + T9, T1B = T8 - T9;
  const T19 = ii[4], T1a = ii[20];
  const T1b = T19 + T1a, T1A = T19 - T1a;
  const Tb = ri[28], Tc = ri[12];
  const Td = Tb + Tc, T1D = Tb - Tc;
  const T1c = ii[28], T1d = ii[12];
  const T1e = T1c + T1d, T1E = T1c - T1d;
  const Te = Ta + Td, T1f = T1b + T1e, T50 = Td - Ta, T4s = T1b - T1e;
  const T2U = T1D - T1E, T2V = T1B + T1A;
  const T2W = KP707106781 * (T2U - T2V), T3u = KP707106781 * (T2V + T2U);
  const T1C = T1A - T1B, T1F = T1D + T1E;
  const T1G = KP707106781 * (T1C - T1F), T3U = KP707106781 * (T1C + T1F);

  const Tg = ri[2], Th = ri[18];
  const Ti = Tg + Th, T1L = Tg - Th;
  const T1h = ii[2], T1i = ii[18];
  const T1j = T1h + T1i, T1J = T1h - T1i;
  const Tj = ri[10], Tk = ri[26];
  const Tl = Tj + Tk, T1I = Tj - Tk;
  const T1k = ii[10], T1l = ii[26];
  const T1m = T1k + T1l, T1M = T1k - T1l;
  const Tm = Ti + Tl, T1n = T1j + T1m;
  const T1K = T1I + T1J, T1N = T1L - T1M;
  const T1O = KP382683432 * T1K - KP923879532 * T1N;
  const T2Z = KP923879532 * T1K + KP382683432 * T1N;
  const T3w = T1J - T1I, T3x = T1L + T1M;
  const T3y = KP923879532 * T3w - KP382683432 * T3x;
  const T3X = KP382683432 * T3w + KP923879532 * T3x;
  const T4u = T1j - T1m, T4v = Ti - Tl, T4w = T4u - T4v, T53 = T4v + T4u;

  const Tn = ri[30], To = ri[14];
  const Tp = Tn + To, T1S = Tn - To;
  const T1o = ii[30], T1p = ii[14];
  const T1q = T1o + T1p, T1Q = T1o - T1p;
  const Tq = ri[6], Tr = ri[22];
  const Ts = Tq + Tr, T1P = Tq - Tr;
  const T1r = ii[6], T1s = ii[22];
  const T1t = T1r + T1s, T1T = T1r - T1s;
  const Tt = Tp + Ts, T1u = T1q + T1t;
  const T1R = T1P + T1Q, T1U = T1S - T1T;
  const T1V = KP382683432 * T1R + KP923879532 * T1U;
  const T2Y = KP382683432 * T1U - KP923879532 * T1R;
  const T3z = T1Q - T1P, T3A = T1S + T1T;
  const T3B = KP923879532 * T3z + KP382683432 * T3A;
  const T3W = KP923879532 * T3A - KP382683432 * T3z;
  const T4x = Tp - Ts, T4y = T1q - T1t, T4z = T4x + T4y, T52 = T4x - T4y;

  const TL = ri[31], TM = ri[15];
  const TN = TL + TM, T2p = TL - TM;
  const T2H = ii[31], T2I = ii[15];
  const T2J = T2H - T2I, T4S = T2H + T2I;
  const TO = ri[7], TP = ri[23];
  const TQ = TO + TP, T2G = TO - TP;
  const T2q = ii[7], T2r = ii[23];
  const T2s = T2q - T2r, T4T = T2q + T2r;
  const TS = ri[3], TT = ri[19];
  const TU = TS + TT, T2x = TS - TT;
  const T2u = ii[3], T2v = ii[19];
  const T2w = T2u - T2v, T4O = T2u + T2v;
  const TV = ri[27], TW = ri[11];
  const TX = TV + TW, T2z = TV - TW;
  const T2A = ii[27], T2B = ii[11];
  const T2C = T2A - T2B, T4P = T2A + T2B;
  const T2t = T2p - T2s, T3L = T2p + T2s, T3O = T2J - T2G, T2K = T2G + T2J;
  const TR = TN + TQ, TY = TU + TX, T5F = TR - TY;
  const T5G = T4S + T4T, T5H = T4O + T4P, T5I = T5G - T5H;
  const T4N = TN - TQ, T4Q = T4O - T4P, T4R = T4N - T4Q, T5j = T4N + T4Q;
  const T2y = T2w - T2x, T2D = T2z + T2C;
  const T2E = KP707106781 * (T2y - T2D), T3P = KP707106781 * (T2y + T2D);
  const T4U = T4S - T4T, T4V = TX - TU, T4W = T4U - T4V, T5k = T4V + T4U;
  const T2L = T2z - T2C, T2M = T2x + T2w;
  const T2N = KP707106781 * (T2L - T2M), T3M = KP707106781 * (T2M + T2L);

  const Tw = ri[1], Tx = ri[17];
  const Ty = Tw + Tx, T2f = Tw - Tx;
  const T1Z = ii[1], T20 = ii[17];
  const T21 = T1Z - T20, T4C = T1Z + T20;
  const Tz = ri[9], TA = ri[25];
  const TB = Tz + TA, T1Y = Tz - TA;
  const T2g = ii[9], T2h = ii[25];
  const T2i = T2g - T2h, T4D = T2g + T2h;
  const TD = ri[5], TE = ri[21];
  const TF = TD + TE, T28 = TD - TE;
  const T29 = ii[5], T2a = ii[21];
  const T2b = T29 - T2a, T4I = T29 + T2a;
  const TG = ri[29], TH = ri[13];
  const TI = TG + TH, T23 = TG - TH;
  const T24 = ii[29], T25 = ii[13];
  const T26 = T24 - T25, T4J = T24 + T25;
  const T22 = T1Y + T21, T3E = T2f + T2i, T3H = T21 - T1Y, T2j = T2f - T2i;
  const TC = Ty + TB, TJ = TF + TI, T5A = TC - TJ;
  const T5B = T4C + T4D, T5C = T4I + T4J, T5D = T5B - T5C;
  const T4E = T4C - T4D, T4F = TI - TF, T4G = T4E - T4F, T5g = T4F + T4E;
  const T27 = T23 - T26, T2c = T28 + T2b;
  const T2d = KP707106781 * (T27 - T2c), T3F = KP707106781 * (T2c + T27);
  const T4H = Ty - TB, T4K = T4I - T4J, T4L = T4H - T4K, T5h = T4H + T4K;
  const T2k = T2b - T28, T2l = T23 + T26;
  const T2m = KP707106781 * (T2k - T2l), T3I = KP707106781 * (T2k + T2l);

  const outR = new Float64Array(32), outI = new Float64Array(32);

  {
    const T4t = T4r - T4s, T4A = KP707106781 * (T4w - T4z);
    const T4B = T4t + T4A, T57 = T4t - T4A;
    const T58 = KP382683432 * T4G - KP923879532 * T4L;
    const T59 = KP382683432 * T4W + KP923879532 * T4R;
    const T5a = T58 - T59, T5c = T58 + T59;
    const T4M = KP923879532 * T4G + KP382683432 * T4L;
    const T4X = KP382683432 * T4R - KP923879532 * T4W;
    const T4Y = T4M + T4X, T56 = T4X - T4M;
    const T51 = T4Z - T50, T54 = KP707106781 * (T52 - T53);
    const T55 = T51 - T54, T5b = T51 + T54;
    outR[22] = T4B - T4Y; outI[22] = T5b - T5c;
    outR[6] = T4B + T4Y; outI[6] = T5b + T5c;
    outI[30] = T55 - T56; outR[30] = T57 - T5a;
    outI[14] = T55 + T56; outR[14] = T57 + T5a;
  }

  {
    const T5d = T4r + T4s, T5e = KP707106781 * (T53 + T52);
    const T5f = T5d + T5e, T5r = T5d - T5e;
    const T5s = KP923879532 * T5g - KP382683432 * T5h;
    const T5t = KP923879532 * T5k + KP382683432 * T5j;
    const T5u = T5s - T5t, T5w = T5s + T5t;
    const T5i = KP382683432 * T5g + KP923879532 * T5h;
    const T5l = KP923879532 * T5j - KP382683432 * T5k;
    const T5m = T5i + T5l, T5q = T5l - T5i;
    const T5n = T50 + T4Z, T5o = KP707106781 * (T4w + T4z);
    const T5p = T5n - T5o, T5v = T5n + T5o;
    outR[18] = T5f - T5m; outI[18] = T5v - T5w;
    outR[2] = T5f + T5m; outI[2] = T5v + T5w;
    outI[26] = T5p - T5q; outR[26] = T5r - T5u;
    outI[10] = T5p + T5q; outR[10] = T5r + T5u;
  }

  {
    const T5x = T7 - Te, T5y = T1n - T1u;
    const T5z = T5x + T5y, T5P = T5x - T5y;
    const T5Q = T5D - T5A, T5R = T5F + T5I;
    const T5S = KP707106781 * (T5Q - T5R), T5U = KP707106781 * (T5Q + T5R);
    const T5E = T5A + T5D, T5J = T5F - T5I;
    const T5K = KP707106781 * (T5E + T5J), T5O = KP707106781 * (T5J - T5E);
    const T5L = T18 - T1f, T5M = Tt - Tm;
    const T5N = T5L - T5M, T5T = T5M + T5L;
    outR[20] = T5z - T5K; outI[20] = T5T - T5U;
    outR[4] = T5z + T5K; outI[4] = T5T + T5U;
    outI[28] = T5N - T5O; outR[28] = T5P - T5S;
    outI[12] = T5N + T5O; outR[12] = T5P + T5S;
  }

  {
    const Tf = T7 + Te, Tu = Tm + Tt;
    const Tv = Tf + Tu, T5V = Tf - Tu;
    const T5W = T5B + T5C, T5X = T5G + T5H;
    const T5Y = T5W - T5X, T60 = T5W + T5X;
    const TK = TC + TJ, TZ = TR + TY;
    const T10 = TK + TZ, T11 = TZ - TK;
    const T1g = T18 + T1f, T1v = T1n + T1u;
    const T1w = T1g - T1v, T5Z = T1g + T1v;
    outR[16] = Tv - T10; outI[16] = T5Z - T60;
    outR[0] = Tv + T10; outI[0] = T5Z + T60;
    outI[8] = T11 + T1w; outR[8] = T5V + T5Y;
    outI[24] = T1w - T11; outR[24] = T5V - T5Y;
  }

  {
    const T1H = T1z - T1G, T1W = T1O - T1V;
    const T1X = T1H + T1W, T33 = T1H - T1W;
    const T2X = T2T - T2W, T30 = T2Y - T2Z;
    const T31 = T2X - T30, T37 = T2X + T30;
    const T2e = T22 - T2d, T2n = T2j - T2m;
    const T2o = KP980785280 * T2e + KP195090322 * T2n;
    const T34 = KP195090322 * T2e - KP980785280 * T2n;
    const T2F = T2t - T2E, T2O = T2K - T2N;
    const T2P = KP195090322 * T2F - KP980785280 * T2O;
    const T35 = KP195090322 * T2O + KP980785280 * T2F;
    const T2Q = T2o + T2P;
    outR[23] = T1X - T2Q; outR[7] = T1X + T2Q;
    const T38 = T34 + T35;
    outI[23] = T37 - T38; outI[7] = T37 + T38;
    const T32 = T2P - T2o;
    outI[31] = T31 - T32; outI[15] = T31 + T32;
    const T36 = T34 - T35;
    outR[31] = T33 - T36; outR[15] = T33 + T36;
  }

  {
    const T3v = T3t - T3u, T3C = T3y - T3B;
    const T3D = T3v + T3C, T41 = T3v - T3C;
    const T3V = T3T - T3U, T3Y = T3W - T3X;
    const T3Z = T3V - T3Y, T45 = T3V + T3Y;
    const T3G = T3E - T3F, T3J = T3H - T3I;
    const T3K = KP555570233 * T3G + KP831469612 * T3J;
    const T42 = KP555570233 * T3J - KP831469612 * T3G;
    const T3N = T3L - T3M, T3Q = T3O - T3P;
    const T3R = KP555570233 * T3N - KP831469612 * T3Q;
    const T43 = KP831469612 * T3N + KP555570233 * T3Q;
    const T3S = T3K + T3R;
    outR[21] = T3D - T3S; outR[5] = T3D + T3S;
    const T46 = T42 + T43;
    outI[21] = T45 - T46; outI[5] = T45 + T46;
    const T40 = T3R - T3K;
    outI[29] = T3Z - T40; outI[13] = T3Z + T40;
    const T44 = T42 - T43;
    outR[29] = T41 - T44; outR[13] = T41 + T44;
  }

  {
    const T47 = T3t + T3u, T48 = T3X + T3W;
    const T49 = T47 + T48, T4l = T47 - T48;
    const T4h = T3T + T3U, T4i = T3y + T3B;
    const T4j = T4h - T4i, T4p = T4h + T4i;
    const T4a = T3E + T3F, T4b = T3H + T3I;
    const T4c = KP980785280 * T4a + KP195090322 * T4b;
    const T4m = KP980785280 * T4b - KP195090322 * T4a;
    const T4d = T3L + T3M, T4e = T3O + T3P;
    const T4f = KP980785280 * T4d - KP195090322 * T4e;
    const T4n = KP195090322 * T4d + KP980785280 * T4e;
    const T4g = T4c + T4f;
    outR[17] = T49 - T4g; outR[1] = T49 + T4g;
    const T4q = T4m + T4n;
    outI[17] = T4p - T4q; outI[1] = T4p + T4q;
    const T4k = T4f - T4c;
    outI[25] = T4j - T4k; outI[9] = T4j + T4k;
    const T4o = T4m - T4n;
    outR[25] = T4l - T4o; outR[9] = T4l + T4o;
  }

  {
    const T39 = T1z + T1G, T3a = T2Z + T2Y;
    const T3b = T39 + T3a, T3n = T39 - T3a;
    const T3j = T2T + T2W, T3k = T1O + T1V;
    const T3l = T3j - T3k, T3r = T3j + T3k;
    const T3c = T22 + T2d, T3d = T2j + T2m;
    const T3e = KP555570233 * T3c + KP831469612 * T3d;
    const T3o = KP831469612 * T3c - KP555570233 * T3d;
    const T3f = T2t + T2E, T3g = T2K + T2N;
    const T3h = KP831469612 * T3f - KP555570233 * T3g;
    const T3p = KP831469612 * T3g + KP555570233 * T3f;
    const T3i = T3e + T3h;
    outR[19] = T3b - T3i; outR[3] = T3b + T3i;
    const T3s = T3o + T3p;
    outI[19] = T3r - T3s; outI[3] = T3r + T3s;
    const T3m = T3h - T3e;
    outI[27] = T3l - T3m; outI[11] = T3l + T3m;
    const T3q = T3o - T3p;
    outR[27] = T3n - T3q; outR[11] = T3n + T3q;
  }

  return [outR, outI];
}

module.exports = { n1_32 };
