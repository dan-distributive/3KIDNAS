'use strict';

// =============================================================================
// t2_32.js -- faithful JS port of dft/scalar/codelets/t2_32.c (non-FMA
// branch), FFTW3's "twiddle-log3 / precompute-twiddles" radix-32 twiddle
// codelet. twinstr only trig-generates W^1, W^3, W^9, W^27 (four raw
// pairs); the other 27 needed multiples (2, 4, 5, ..., 31 excluding
// 1/3/9/27) are DERIVED via complex products of those four, exactly
// mirroring the C source's own derivation chain bit-for-bit (same trick
// as t2_5/t2_8/t2_10/t2_16/t2_20/t2_25.js). Same calling convention as
// every other twiddle codelet here (br/bi/Wc/Ws with Composite1D.js's
// full r-1 = 31 pair table already built) -- only indices 1, 3, 9, 27 are
// actually read; every other Wc[k]/Ws[k] is deliberately ignored in favor
// of recomputing the equivalent product chain, to match FFTW's exact
// rounding (real FFTW never has a full table either -- it only has the 4
// raw generator pairs from twinstr).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP195090322 = 0.195090322016128267848284868477022240927691618;
const KP980785280 = 0.980785280403230449126182236134239036973933731;
const KP555570233 = 0.555570233019602224742830813948532874374937191;
const KP831469612 = 0.831469612302545237078788377617905756738560812;
const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function t2_32(br, bi, Wc, Ws) {
  const outR = new Float64Array(32), outI = new Float64Array(32);

  // -- precompute all needed twiddle products from the 4 raw generators --
  const T2 = Wc[1], T5 = Ws[1], T3 = Wc[3], T6 = Ws[3];
  const T4 = T2 * T3, Tc = T5 * T3, T7 = T5 * T6, Tb = T2 * T6;
  const T8 = T4 + T7, TM = T4 - T7, TO = Tb + Tc, Td = Tb - Tc;
  const T9 = Wc[9];
  const Ts = T2 * T9, T1d = T6 * T9, Tx = T5 * T9, T18 = T3 * T9;
  const Te = Ws[9];
  const Tt = T5 * Te, T1c = T3 * Te, Tw = T2 * Te, T19 = T6 * Te;
  const Th = Wc[27];
  const TB = T3 * Th, T14 = T5 * Th, TG = T6 * Th, TZ = T2 * Th;
  const Tl = Ws[27];
  const TC = T6 * Tl, T13 = T2 * Tl, TF = T3 * Tl, T10 = T5 * Tl;

  const TD = TB + TC, TH = TF - TG, T1y = TZ + T10, T1H = TF + TG;
  const T15 = T13 + T14, T1A = T13 - T14, T11 = TZ - T10, T1F = TB - TC;
  const T1n = T9 * Th + Te * Tl, T1p = T9 * Tl - Te * Th;

  const T2o = T8 * Th, T2p = Td * Tl;
  const T2q = T2o + T2p, T2I = T2o - T2p;
  const T2s = T8 * Tl, T2t = Td * Th;
  const T2u = T2s - T2t, T2K = T2s + T2t;

  const T2T = TM * Th, T2U = TO * Tl;
  const T2V = T2T - T2U, T3b = T2T + T2U;
  const T2X = TM * Tl, T2Y = TO * Th;
  const T2Z = T2X + T2Y, T3d = T2X - T2Y;
  const Tu = Ts + Tt, Ty = Tw - Tx;
  const T3l = Tu * Th + Ty * Tl, T3n = Tu * Tl - Ty * Th;

  const T1t = Ts - Tt, T1v = Tw + Tx;
  const T2f = T1t * Th + T1v * Tl, T2h = T1t * Tl - T1v * Th;
  const T1a = T18 - T19, T1e = T1c + T1d;
  const T32 = T1a * Th + T1e * Tl, T34 = T1a * Tl - T1e * Th;
  const T1W = T18 + T19, T1Y = T1c - T1d;
  const T2C = T1W * Th + T1Y * Tl, T2E = T1W * Tl - T1Y * Th;

  const Ta = T8 * T9, Tf = Td * Te;
  const Tg = Ta - Tf, TR = Ta + Tf;
  const Ti = T8 * Te, Tj = Td * T9;
  const Tk = Ti + Tj, TS = Ti - Tj;

  const Tm = Tg * Th + Tk * Tl, TV = TR * Tl - TS * Th;
  const To = Tg * Tl - Tk * Th, TT = TR * Th + TS * Tl;

  const T1K = TM * T9, T1L = TO * Te;
  const T1M = T1K - T1L, T21 = T1K + T1L;
  const T1N = TM * Te, T1O = TO * T9;
  const T1P = T1N + T1O, T22 = T1N - T1O;

  const T1Q = T1M * Th + T1P * Tl, T25 = T21 * Tl - T22 * Th;
  const T1S = T1M * Tl - T1P * Th, T23 = T21 * Th + T22 * Tl;

  // -- butterfly 0: k = 0, 16, 8, 24 --
  const T1 = br[0], T7G = bi[0];
  const Tn = br[16], Tp = bi[16];
  const Tq = Tm * Tn + To * Tp, T7F = Tm * Tp - To * Tn;
  const Tv = br[8], Tz = bi[8];
  const TA = Tu * Tv + Ty * Tz, T3C = Tu * Tz - Ty * Tv;
  const TE = br[24], TI = bi[24];
  const TJ = TD * TE + TH * TI, T3D = TD * TI - TH * TE;

  const Tr = T1 + Tq, TK = TA + TJ;
  const TL = Tr + TK, T6f = Tr - TK;
  const T8a = T7G - T7F, T8b = TA - TJ;
  const T8c = T8a - T8b, T8q = T8b + T8a;

  const T3B = T1 - Tq, T3E = T3C - T3D;
  const T3F = T3B - T3E, T5t = T3B + T3E;
  const T7E = T3C + T3D, T7H = T7F + T7G;
  const T7I = T7E + T7H, T7W = T7H - T7E;

  // -- butterfly 1: k = 1, 25, 17, 9 --
  const T2c = br[1], T2d = bi[1];
  const T2e = T2 * T2c + T5 * T2d, T4g = T2 * T2d - T5 * T2c;
  const T2r = br[25], T2v = bi[25];
  const T2w = T2q * T2r + T2u * T2v, T4z = T2q * T2v - T2u * T2r;

  const T2g = br[17], T2i = bi[17];
  const T2j = T2f * T2g + T2h * T2i, T4h = T2f * T2i - T2h * T2g;
  const T2l = br[9], T2m = bi[9];
  const T2n = T9 * T2l + Te * T2m, T4y = T9 * T2m - Te * T2l;

  const T2k = T2e + T2j, T2x = T2n + T2w;
  const T2y = T2k + T2x, T6B = T2k - T2x;
  const T6w = T4g + T4h, T6x = T4y + T4z;
  const T6y = T6w - T6x, T7j = T6w + T6x;

  const T4i = T4g - T4h, T4j = T2n - T2w;
  const T4k = T4i + T4j, T5J = T4i - T4j;
  const T4x = T2e - T2j, T4A = T4y - T4z;
  const T4B = T4x - T4A, T5G = T4x + T4A;

  // -- butterfly 2: k = 31, 23, 15, 7 --
  const T2W = br[31], T30 = bi[31];
  const T31 = T2V * T2W + T2Z * T30, T4Y = T2V * T30 - T2Z * T2W;
  const T3c = br[23], T3e = bi[23];
  const T3f = T3b * T3c + T3d * T3e, T4J = T3b * T3e - T3d * T3c;

  const T33 = br[15], T35 = bi[15];
  const T36 = T32 * T33 + T34 * T35, T4Z = T32 * T35 - T34 * T33;
  const T38 = br[7], T39 = bi[7];
  const T3a = TR * T38 + TS * T39, T4I = TR * T39 - TS * T38;

  const T37 = T31 + T36, T3g = T3a + T3f;
  const T3h = T37 + T3g, T6H = T37 - T3g;
  const T6M = T4Y + T4Z, T6N = T4I + T4J;
  const T6O = T6M - T6N, T7o = T6M + T6N;

  const T4H = T31 - T36, T4K = T4I - T4J;
  const T4L = T4H - T4K, T5N = T4H + T4K;
  const T50 = T4Y - T4Z, T51 = T3a - T3f;
  const T52 = T50 + T51, T5Q = T50 - T51;

  // -- butterfly 3: k = 4, 12, 20, 28 --
  const TN = br[4], TP = bi[4];
  const TQ = TM * TN + TO * TP, T3G = TM * TP - TO * TN;
  const T1b = br[12], T1f = bi[12];
  const T1g = T1a * T1b + T1e * T1f, T3N = T1a * T1f - T1e * T1b;

  const TU = br[20], TW = bi[20];
  const TX = TT * TU + TV * TW, T3H = TT * TW - TV * TU;
  const T12 = br[28], T16 = bi[28];
  const T17 = T11 * T12 + T15 * T16, T3M = T11 * T16 - T15 * T12;

  const TY = TQ + TX, T1h = T17 + T1g;
  const T1i = TY + T1h, T7V = T1h - TY;
  const T6g = T3G + T3H, T6h = T3M + T3N;
  const T6i = T6g - T6h, T7D = T6g + T6h;

  const T3I = T3G - T3H, T3J = TQ - TX;
  const T3K = T3I - T3J, T5u = T3J + T3I;
  const T3L = T17 - T1g, T3O = T3M - T3N;
  const T3P = T3L + T3O, T5v = T3L - T3O;

  // -- butterfly 4: k = 2, 26, 18, 10 --
  const T1k = br[2], T1l = bi[2];
  const T1m = T8 * T1k + Td * T1l, T3S = T8 * T1l - Td * T1k;
  const T1z = br[26], T1B = bi[26];
  const T1C = T1y * T1z + T1A * T1B, T3Z = T1y * T1B - T1A * T1z;

  const T1o = br[18], T1q = bi[18];
  const T1r = T1n * T1o + T1p * T1q, T3T = T1n * T1q - T1p * T1o;
  const T1u = br[10], T1w = bi[10];
  const T1x = T1t * T1u + T1v * T1w, T3Y = T1t * T1w - T1v * T1u;

  const T1s = T1m + T1r, T1D = T1x + T1C;
  const T1E = T1s + T1D, T6n = T1s - T1D;
  const T6k = T3S + T3T, T6l = T3Y + T3Z;
  const T6m = T6k - T6l, T7e = T6k + T6l;

  const T3U = T3S - T3T, T3V = T1x - T1C;
  const T3W = T3U + T3V, T5y = T3U - T3V;
  const T3X = T1m - T1r, T40 = T3Y - T3Z;
  const T41 = T3X - T40, T5z = T3X + T40;

  // -- butterfly 5: k = 30, 22, 14, 6 --
  const T1G = br[30], T1I = bi[30];
  const T1J = T1F * T1G + T1H * T1I, T43 = T1F * T1I - T1H * T1G;
  const T24 = br[22], T26 = bi[22];
  const T27 = T23 * T24 + T25 * T26, T4a = T23 * T26 - T25 * T24;

  const T1R = br[14], T1T = bi[14];
  const T1U = T1Q * T1R + T1S * T1T, T44 = T1Q * T1T - T1S * T1R;
  const T1X = br[6], T1Z = bi[6];
  const T20 = T1W * T1X + T1Y * T1Z, T49 = T1W * T1Z - T1Y * T1X;

  const T1V = T1J + T1U, T28 = T20 + T27;
  const T29 = T1V + T28, T6p = T1V - T28;
  const T6q = T43 + T44, T6r = T49 + T4a;
  const T6s = T6q - T6r, T7f = T6q + T6r;

  const T45 = T43 - T44, T46 = T20 - T27;
  const T47 = T45 + T46, T5B = T45 - T46;
  const T48 = T1J - T1U, T4b = T49 - T4a;
  const T4c = T48 - T4b, T5C = T48 + T4b;

  // -- butterfly 6: k = 5, 21, 29, 13 --
  const T2z = br[5], T2A = bi[5];
  const T2B = T21 * T2z + T22 * T2A, T4r = T21 * T2A - T22 * T2z;
  const T2D = br[21], T2F = bi[21];
  const T2G = T2C * T2D + T2E * T2F, T4s = T2C * T2F - T2E * T2D;

  const T4q = T2B - T2G, T4t = T4r - T4s;

  const T2J = br[29], T2L = bi[29];
  const T2M = T2I * T2J + T2K * T2L, T4m = T2I * T2L - T2K * T2J;
  const T2N = br[13], T2O = bi[13];
  const T2P = T1M * T2N + T1P * T2O, T4n = T1M * T2O - T1P * T2N;

  const T4l = T2M - T2P, T4o = T4m - T4n;

  const T2H = T2B + T2G, T2Q = T2M + T2P;
  const T2R = T2H + T2Q, T6z = T2Q - T2H;
  const T6C = T4r + T4s, T6D = T4m + T4n;
  const T6E = T6C - T6D, T7k = T6C + T6D;

  const T4p = T4l - T4o, T4u = T4q + T4t;
  const T4v = KP707106781 * (T4p - T4u), T5H = KP707106781 * (T4u + T4p);
  const T4C = T4t - T4q, T4D = T4l + T4o;
  const T4E = KP707106781 * (T4C - T4D), T5K = KP707106781 * (T4C + T4D);

  // -- butterfly 7: k = 3, 19, 27, 11 --
  const T3i = br[3], T3j = bi[3];
  const T3k = T3 * T3i + T6 * T3j, T4M = T3 * T3j - T6 * T3i;
  const T3m = br[19], T3o = bi[19];
  const T3p = T3l * T3m + T3n * T3o, T4N = T3l * T3o - T3n * T3m;

  const T4O = T4M - T4N, T4P = T3k - T3p;

  const T3r = br[27], T3s = bi[27];
  const T3t = Th * T3r + Tl * T3s, T4S = Th * T3s - Tl * T3r;
  const T3u = br[11], T3v = bi[11];
  const T3w = Tg * T3u + Tk * T3v, T4T = Tg * T3v - Tk * T3u;

  const T4R = T3t - T3w, T4U = T4S - T4T;

  const T3q = T3k + T3p, T3x = T3t + T3w;
  const T3y = T3q + T3x, T6P = T3x - T3q;
  const T6I = T4M + T4N, T6J = T4S + T4T;
  const T6K = T6I - T6J, T7p = T6I + T6J;

  const T4Q = T4O - T4P, T4V = T4R + T4U;
  const T4W = KP707106781 * (T4Q - T4V), T5R = KP707106781 * (T4Q + T4V);
  const T53 = T4R - T4U, T54 = T4P + T4O;
  const T55 = KP707106781 * (T53 - T54), T5O = KP707106781 * (T54 + T53);

  // -- combine group A: outputs 16, 0, 24, 8 --
  const T1j = TL + T1i, T2a = T1E + T29;
  const T2b = T1j + T2a, T7x = T1j - T2a;
  const T7C = T7e + T7f, T7J = T7D + T7I;
  const T7K = T7C + T7J, T7M = T7J - T7C;

  const T2S = T2y + T2R, T3z = T3h + T3y;
  const T3A = T2S + T3z, T7L = T3z - T2S;
  const T7y = T7j + T7k, T7z = T7o + T7p;
  const T7A = T7y - T7z, T7B = T7y + T7z;

  outR[16] = T2b - T3A;
  outI[16] = T7K - T7B;
  outR[0] = T2b + T3A;
  outI[0] = T7B + T7K;
  outR[24] = T7x - T7A;
  outI[24] = T7M - T7L;
  outR[8] = T7x + T7A;
  outI[8] = T7L + T7M;

  // -- combine group B: outputs 20, 4, 28, 12 --
  const T7d = TL - T1i, T7g = T7e - T7f;
  const T7h = T7d + T7g, T7t = T7d - T7g;
  const T7O = T29 - T1E, T7P = T7I - T7D;
  const T7Q = T7O + T7P, T7S = T7P - T7O;

  const T7i = T2y - T2R, T7l = T7j - T7k;
  const T7m = T7i + T7l, T7u = T7l - T7i;
  const T7n = T3h - T3y, T7q = T7o - T7p;
  const T7r = T7n - T7q, T7v = T7n + T7q;

  const T7s = KP707106781 * (T7m + T7r);
  outR[20] = T7h - T7s;
  outR[4] = T7h + T7s;
  const T7N = KP707106781 * (T7u + T7v);
  outI[4] = T7N + T7Q;
  outI[20] = T7Q - T7N;
  const T7w = KP707106781 * (T7u - T7v);
  outR[28] = T7t - T7w;
  outR[12] = T7t + T7w;
  const T7R = KP707106781 * (T7r - T7m);
  outI[12] = T7R + T7S;
  outI[28] = T7S - T7R;

  // -- combine group C: outputs 22, 6, 30, 14, 18, 2, 26, 10 --
  const T6o = T6m - T6n, T6t = T6p + T6s;
  const T6u = KP707106781 * (T6o - T6t), T7U = KP707106781 * (T6o + T6t);

  const T6j = T6f - T6i, T7X = T7V + T7W, T83 = T7W - T7V, T6X = T6f + T6i;

  const T75 = T6H + T6K, T76 = T6O + T6P;
  const T77 = KP923879532 * T75 - KP382683432 * T76;
  const T7b = KP382683432 * T75 + KP923879532 * T76;
  const T6Y = T6n + T6m, T6Z = T6p - T6s;
  const T70 = KP707106781 * (T6Y + T6Z), T82 = KP707106781 * (T6Z - T6Y);

  const T6A = T6y - T6z, T6F = T6B - T6E;
  const T6G = KP923879532 * T6A + KP382683432 * T6F;
  const T6U = KP382683432 * T6A - KP923879532 * T6F;

  const T72 = T6y + T6z, T73 = T6B + T6E;
  const T74 = KP382683432 * T72 + KP923879532 * T73;
  const T7a = KP923879532 * T72 - KP382683432 * T73;
  const T6L = T6H - T6K, T6Q = T6O - T6P;
  const T6R = KP382683432 * T6L - KP923879532 * T6Q;
  const T6V = KP923879532 * T6L + KP382683432 * T6Q;

  const T6v = T6j + T6u, T6S = T6G + T6R;
  outR[22] = T6v - T6S;
  outR[6] = T6v + T6S;
  const T81 = T6U + T6V, T84 = T82 + T83;
  outI[6] = T81 + T84;
  outI[22] = T84 - T81;

  const T6T = T6j - T6u, T6W = T6U - T6V;
  outR[30] = T6T - T6W;
  outR[14] = T6T + T6W;
  const T85 = T6R - T6G, T86 = T83 - T82;
  outI[14] = T85 + T86;
  outI[30] = T86 - T85;

  const T71 = T6X + T70, T78 = T74 + T77;
  outR[18] = T71 - T78;
  outR[2] = T71 + T78;
  const T7T = T7a + T7b, T7Y = T7U + T7X;
  outI[2] = T7T + T7Y;
  outI[18] = T7Y - T7T;

  const T79 = T6X - T70, T7c = T7a - T7b;
  outR[26] = T79 - T7c;
  outR[10] = T79 + T7c;
  const T7Z = T77 - T74, T80 = T7X - T7U;
  outI[10] = T7Z + T80;
  outI[26] = T80 - T7Z;

  // -- combine group D: outputs 23, 7, 31, 15, 19, 3, 27, 11 --
  const T3Q = KP707106781 * (T3K - T3P);
  const T3R = T3F - T3Q, T5d = T3F + T3Q;
  const T8p = KP707106781 * (T5v - T5u);
  const T8r = T8p + T8q, T8x = T8q - T8p;

  const T42 = KP382683432 * T3W - KP923879532 * T41;
  const T4d = KP382683432 * T47 + KP923879532 * T4c;
  const T4e = T42 - T4d, T8o = T42 + T4d;
  const T5l = T4L + T4W, T5m = T52 + T55;
  const T5n = KP831469612 * T5l - KP555570233 * T5m;
  const T5r = KP555570233 * T5l + KP831469612 * T5m;

  const T4w = T4k - T4v, T4F = T4B - T4E;
  const T4G = KP980785280 * T4w + KP195090322 * T4F;
  const T5a = KP195090322 * T4w - KP980785280 * T4F;
  const T5e = KP923879532 * T3W + KP382683432 * T41;
  const T5f = KP382683432 * T4c - KP923879532 * T47;
  const T5g = T5e + T5f, T8w = T5f - T5e;

  const T5i = T4k + T4v, T5j = T4B + T4E;
  const T5k = KP555570233 * T5i + KP831469612 * T5j;
  const T5q = KP831469612 * T5i - KP555570233 * T5j;
  const T4X = T4L - T4W, T56 = T52 - T55;
  const T57 = KP195090322 * T4X - KP980785280 * T56;
  const T5b = KP980785280 * T4X + KP195090322 * T56;

  const T4f = T3R + T4e, T58 = T4G + T57;
  outR[23] = T4f - T58;
  outR[7] = T4f + T58;
  const T8v = T5a + T5b, T8y = T8w + T8x;
  outI[7] = T8v + T8y;
  outI[23] = T8y - T8v;

  const T59 = T3R - T4e, T5c = T5a - T5b;
  outR[31] = T59 - T5c;
  outR[15] = T59 + T5c;
  const T8z = T57 - T4G, T8A = T8x - T8w;
  outI[15] = T8z + T8A;
  outI[31] = T8A - T8z;

  const T5h = T5d + T5g, T5o = T5k + T5n;
  outR[19] = T5h - T5o;
  outR[3] = T5h + T5o;
  const T8n = T5q + T5r, T8s = T8o + T8r;
  outI[3] = T8n + T8s;
  outI[19] = T8s - T8n;

  const T5p = T5d - T5g, T5s = T5q - T5r;
  outR[27] = T5p - T5s;
  outR[11] = T5p + T5s;
  const T8t = T5n - T5k, T8u = T8r - T8o;
  outI[11] = T8t + T8u;
  outI[27] = T8u - T8t;

  // -- combine group E: outputs 21, 5, 29, 13, 17, 1, 25, 9 --
  const T5w = KP707106781 * (T5u + T5v);
  const T5x = T5t - T5w, T5Z = T5t + T5w;
  const T89 = KP707106781 * (T3K + T3P);
  const T8d = T89 + T8c, T8j = T8c - T89;

  const T5A = KP923879532 * T5y - KP382683432 * T5z;
  const T5D = KP382683432 * T5C + KP923879532 * T5B;
  const T5E = T5A - T5D, T88 = T5A + T5D;
  const T67 = T5N + T5O, T68 = T5Q + T5R;
  const T69 = KP980785280 * T67 - KP195090322 * T68;
  const T6d = KP195090322 * T67 + KP980785280 * T68;

  const T5I = T5G - T5H, T5L = T5J - T5K;
  const T5M = KP555570233 * T5I + KP831469612 * T5L;
  const T5W = KP555570233 * T5L - KP831469612 * T5I;
  const T60 = KP382683432 * T5y + KP923879532 * T5z;
  const T61 = KP923879532 * T5C - KP382683432 * T5B;
  const T62 = T60 + T61, T8i = T61 - T60;

  const T64 = T5G + T5H, T65 = T5J + T5K;
  const T66 = KP980785280 * T64 + KP195090322 * T65;
  const T6c = KP980785280 * T65 - KP195090322 * T64;
  const T5P = T5N - T5O, T5S = T5Q - T5R;
  const T5T = KP555570233 * T5P - KP831469612 * T5S;
  const T5X = KP831469612 * T5P + KP555570233 * T5S;

  const T5F = T5x + T5E, T5U = T5M + T5T;
  outR[21] = T5F - T5U;
  outR[5] = T5F + T5U;
  const T8h = T5W + T5X, T8k = T8i + T8j;
  outI[5] = T8h + T8k;
  outI[21] = T8k - T8h;

  const T5V = T5x - T5E, T5Y = T5W - T5X;
  outR[29] = T5V - T5Y;
  outR[13] = T5V + T5Y;
  const T8l = T5T - T5M, T8m = T8j - T8i;
  outI[13] = T8l + T8m;
  outI[29] = T8m - T8l;

  const T63 = T5Z + T62, T6a = T66 + T69;
  outR[17] = T63 - T6a;
  outR[1] = T63 + T6a;
  const T87 = T6c + T6d, T8e = T88 + T8d;
  outI[1] = T87 + T8e;
  outI[17] = T8e - T87;

  const T6b = T5Z - T62, T6e = T6c - T6d;
  outR[25] = T6b - T6e;
  outR[9] = T6b + T6e;
  const T8f = T69 - T66, T8g = T8d - T88;
  outI[9] = T8f + T8g;
  outI[25] = T8g - T8f;

  return [outR, outI];
}

module.exports = { t2_32 };
