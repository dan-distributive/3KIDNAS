'use strict';

// =============================================================================
// hf2_32.js -- faithful JS port of rdft/scalar/r2cf/hf2_32.c (non-FMA
// branch). Alternate-codegen EVEN-radix (r=32) twiddle codelet -- same role
// as hf2_16.js, larger radix. twinstr trig-generates W^1, W^3, W^9, W^27
// (four raw pairs); the rest are derived via complex products, exactly
// matching the C source's operation order. Same (cr,ci,Wc,Ws) ->
// (outCr,outCi) convention as hf2_16.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP555570233 = 0.555570233019602224742830813948532874374937191;
const KP831469612 = 0.831469612302545237078788377617905756738560812;
const KP980785280 = 0.980785280403230449126182236134239036973933731;
const KP195090322 = 0.195090322016128267848284868477022240927691618;
const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function hf2_32(cr, ci, Wc, Ws) {
  const T2 = Wc[1], T5 = Ws[1], T3 = Wc[3], T6 = Ws[3];
  const T4 = T2 * T3, Tc = T5 * T3, T7 = T5 * T6, Tb = T2 * T6;
  const T8 = T4 + T7, TM = T4 - T7, TO = Tb + Tc, Td = Tb - Tc;
  const T9 = Wc[9], Te = Ws[9];
  const Ts = T2 * T9, T1d = T6 * T9, Tx = T5 * T9, T18 = T3 * T9;
  const Tt = T5 * Te, T1c = T3 * Te, Tw = T2 * Te, T19 = T6 * Te;
  const Th = Wc[27], Tl = Ws[27];
  const TB = T3 * Th, T14 = T5 * Th, TG = T6 * Th, TZ = T2 * Th;
  const TC = T6 * Tl, T13 = T2 * Tl, TF = T3 * Tl, T10 = T5 * Tl;

  const TD = TB + TC;
  const TH = TF - TG;
  const T1y = TZ + T10;
  const T1H = TF + TG;
  const T15 = T13 + T14;
  const T1A = T13 - T14;
  const T11 = TZ - T10;
  const T1F = TB - TC;
  const T1n = T9 * Th + Te * Tl; // FMA
  const T1p = T9 * Tl - Te * Th; // FNMS

  const T2o = T8 * Th, T2p = Td * Tl;
  const T2q = T2o + T2p, T2I = T2o - T2p;
  const T2s = T8 * Tl, T2t = Td * Th;
  const T2u = T2s - T2t, T2K = T2s + T2t;

  const T2T = TM * Th, T2U = TO * Tl;
  const T2V = T2T - T2U, T3b = T2T + T2U;
  const T2X = TM * Tl, T2Y = TO * Th;
  const T2Z = T2X + T2Y, T3d = T2X - T2Y;
  const Tu = Ts + Tt;
  const Ty = Tw - Tx;
  const T3l = Tu * Th + Ty * Tl; // FMA
  const T3n = Tu * Tl - Ty * Th; // FNMS

  const T1t = Ts - Tt;
  const T1v = Tw + Tx;
  const T2f = T1t * Th + T1v * Tl; // FMA
  const T2h = T1t * Tl - T1v * Th; // FNMS
  const T1a = T18 - T19;
  const T1e = T1c + T1d;
  const T32 = T1a * Th + T1e * Tl; // FMA
  const T34 = T1a * Tl - T1e * Th; // FNMS
  const T1W = T18 + T19;
  const T1Y = T1c - T1d;
  const T2C = T1W * Th + T1Y * Tl; // FMA
  const T2E = T1W * Tl - T1Y * Th; // FNMS

  const Ta = T8 * T9, Tf = Td * Te;
  const Tg = Ta - Tf, TR = Ta + Tf;
  const Ti = T8 * Te, Tj = Td * T9;
  const Tk = Ti + Tj, TS = Ti - Tj;

  const Tm = Tg * Th + Tk * Tl; // FMA
  const TV = TR * Tl - TS * Th; // FNMS
  const To = Tg * Tl - Tk * Th; // FNMS
  const TT = TR * Th + TS * Tl; // FMA

  const T1K = TM * T9, T1L = TO * Te;
  const T1M = T1K - T1L, T21 = T1K + T1L;
  const T1N = TM * Te, T1O = TO * T9;
  const T1P = T1N + T1O, T22 = T1N - T1O;

  const T1Q = T1M * Th + T1P * Tl; // FMA
  const T25 = T21 * Tl - T22 * Th; // FNMS
  const T1S = T1M * Tl - T1P * Th; // FNMS
  const T23 = T21 * Th + T22 * Tl; // FMA

  const T1 = cr[0];
  const T7G = ci[0];
  const Tn = cr[16], Tp = ci[16];
  const Tq = Tm * Tn + To * Tp; // FMA
  const T7F = Tm * Tp - To * Tn; // FNMS

  const Tv = cr[8], Tz = ci[8];
  const TA = Tu * Tv + Ty * Tz; // FMA
  const T3C = Tu * Tz - Ty * Tv; // FNMS
  const TE = cr[24], TI = ci[24];
  const TJ = TD * TE + TH * TI; // FMA
  const T3D = TD * TI - TH * TE; // FNMS

  const Tr = T1 + Tq;
  const TK = TA + TJ;
  const TL = Tr + TK;
  const T6f = Tr - TK;
  const T8a = TA - TJ;
  const T8b = T7G - T7F;
  const T8c = T8a + T8b;
  const T8q = T8b - T8a;

  const T3B = T1 - Tq;
  const T3E = T3C - T3D;
  const T3F = T3B + T3E;
  const T5t = T3B - T3E;
  const T7E = T3C + T3D;
  const T7H = T7F + T7G;
  const T7I = T7E + T7H;
  const T7W = T7H - T7E;

  const T2c = cr[1], T2d = ci[1];
  const T2e = T2 * T2c + T5 * T2d; // FMA
  const T4x = T2 * T2d - T5 * T2c; // FNMS
  const T2r = cr[25], T2v = ci[25];
  const T2w = T2q * T2r + T2u * T2v; // FMA
  const T4i = T2q * T2v - T2u * T2r; // FNMS

  const T2g = cr[17], T2i = ci[17];
  const T2j = T2f * T2g + T2h * T2i; // FMA
  const T4y = T2f * T2i - T2h * T2g; // FNMS
  const T2l = cr[9], T2m = ci[9];
  const T2n = T9 * T2l + Te * T2m; // FMA
  const T4h = T9 * T2m - Te * T2l; // FNMS

  const T2k = T2e + T2j;
  const T2x = T2n + T2w;
  const T2y = T2k + T2x;
  const T6B = T2k - T2x;
  const T6w = T4x + T4y;
  const T6x = T4h + T4i;
  const T6y = T6w - T6x;
  const T7j = T6w + T6x;

  const T4g = T2e - T2j;
  const T4j = T4h - T4i;
  const T4k = T4g + T4j;
  const T5G = T4g - T4j;
  const T4z = T4x - T4y;
  const T4A = T2n - T2w;
  const T4B = T4z - T4A;
  const T5J = T4z + T4A;

  const T2W = cr[31], T30 = ci[31];
  const T31 = T2V * T2W + T2Z * T30; // FMA
  const T4H = T2V * T30 - T2Z * T2W; // FNMS
  const T3c = cr[23], T3e = ci[23];
  const T3f = T3b * T3c + T3d * T3e; // FMA
  const T50 = T3b * T3e - T3d * T3c; // FNMS

  const T33 = cr[15], T35 = ci[15];
  const T36 = T32 * T33 + T34 * T35; // FMA
  const T4I = T32 * T35 - T34 * T33; // FNMS
  const T38 = cr[7], T39 = ci[7];
  const T3a = TR * T38 + TS * T39; // FMA
  const T4Z = TR * T39 - TS * T38; // FNMS

  const T37 = T31 + T36;
  const T3g = T3a + T3f;
  const T3h = T37 + T3g;
  const T6H = T37 - T3g;
  const T6M = T4H + T4I;
  const T6N = T4Z + T50;
  const T6O = T6M - T6N;
  const T7o = T6M + T6N;

  const T4J = T4H - T4I;
  const T4K = T3a - T3f;
  const T4L = T4J - T4K;
  const T5Q = T4J + T4K;
  const T4Y = T31 - T36;
  const T51 = T4Z - T50;
  const T52 = T4Y + T51;
  const T5N = T4Y - T51;

  const TN = cr[4], TP = ci[4];
  const TQ = TM * TN + TO * TP; // FMA
  const T3H = TM * TP - TO * TN; // FNMS
  const T1b = cr[12], T1f = ci[12];
  const T1g = T1a * T1b + T1e * T1f; // FMA
  const T3N = T1a * T1f - T1e * T1b; // FNMS

  const TU = cr[20], TW = ci[20];
  const TX = TT * TU + TV * TW; // FMA
  const T3I = TT * TW - TV * TU; // FNMS
  const T12 = cr[28], T16 = ci[28];
  const T17 = T11 * T12 + T15 * T16; // FMA
  const T3M = T11 * T16 - T15 * T12; // FNMS

  const TY = TQ + TX;
  const T1h = T17 + T1g;
  const T1i = TY + T1h;
  const T7V = TY - T1h;
  const T6g = T3M + T3N;
  const T6h = T3H + T3I;
  const T6i = T6g - T6h;
  const T7D = T6h + T6g;

  const T3G = TQ - TX;
  const T3J = T3H - T3I;
  const T3K = T3G + T3J;
  const T5u = T3G - T3J;
  const T3L = T17 - T1g;
  const T3O = T3M - T3N;
  const T3P = T3L - T3O;
  const T5v = T3L + T3O;

  const T1k = cr[2], T1l = ci[2];
  const T1m = T8 * T1k + Td * T1l; // FMA
  const T3X = T8 * T1l - Td * T1k; // FNMS
  const T1z = cr[26], T1B = ci[26];
  const T1C = T1y * T1z + T1A * T1B; // FMA
  const T3U = T1y * T1B - T1A * T1z; // FNMS

  const T1o = cr[18], T1q = ci[18];
  const T1r = T1n * T1o + T1p * T1q; // FMA
  const T3Y = T1n * T1q - T1p * T1o; // FNMS
  const T1u = cr[10], T1w = ci[10];
  const T1x = T1t * T1u + T1v * T1w; // FMA
  const T3T = T1t * T1w - T1v * T1u; // FNMS

  const T1s = T1m + T1r;
  const T1D = T1x + T1C;
  const T1E = T1s + T1D;
  const T6k = T1s - T1D;
  const T6l = T3X + T3Y;
  const T6m = T3T + T3U;
  const T6n = T6l - T6m;
  const T7f = T6l + T6m;

  const T3S = T1m - T1r;
  const T3V = T3T - T3U;
  const T3W = T3S + T3V;
  const T5z = T3S - T3V;
  const T3Z = T3X - T3Y;
  const T40 = T1x - T1C;
  const T41 = T3Z - T40;
  const T5y = T3Z + T40;

  const T1G = cr[30], T1I = ci[30];
  const T1J = T1F * T1G + T1H * T1I; // FMA
  const T43 = T1F * T1I - T1H * T1G; // FNMS
  const T24 = cr[22], T26 = ci[22];
  const T27 = T23 * T24 + T25 * T26; // FMA
  const T4a = T23 * T26 - T25 * T24; // FNMS

  const T1R = cr[14], T1T = ci[14];
  const T1U = T1Q * T1R + T1S * T1T; // FMA
  const T44 = T1Q * T1T - T1S * T1R; // FNMS
  const T1X = cr[6], T1Z = ci[6];
  const T20 = T1W * T1X + T1Y * T1Z; // FMA
  const T49 = T1W * T1Z - T1Y * T1X; // FNMS

  const T1V = T1J + T1U;
  const T28 = T20 + T27;
  const T29 = T1V + T28;
  const T6p = T1V - T28;
  const T6q = T43 + T44;
  const T6r = T49 + T4a;
  const T6s = T6q - T6r;
  const T7e = T6q + T6r;

  const T45 = T43 - T44;
  const T46 = T20 - T27;
  const T47 = T45 - T46;
  const T5C = T45 + T46;
  const T48 = T1J - T1U;
  const T4b = T49 - T4a;
  const T4c = T48 + T4b;
  const T5B = T48 - T4b;

  const T2z = cr[5], T2A = ci[5];
  const T2B = T21 * T2z + T22 * T2A; // FMA
  const T4m = T21 * T2A - T22 * T2z; // FNMS
  const T2D = cr[21], T2F = ci[21];
  const T2G = T2C * T2D + T2E * T2F; // FMA
  const T4n = T2C * T2F - T2E * T2D; // FNMS

  const T4l = T2B - T2G;
  const T4o = T4m - T4n;

  const T2J = cr[29], T2L = ci[29];
  const T2M = T2I * T2J + T2K * T2L; // FMA
  const T4q = T2I * T2L - T2K * T2J; // FNMS
  const T2N = cr[13], T2O = ci[13];
  const T2P = T1M * T2N + T1P * T2O; // FMA
  const T4r = T1M * T2O - T1P * T2N; // FNMS

  const T4s = T4q - T4r;
  const T4t = T2M - T2P;

  const T2H = T2B + T2G;
  const T2Q = T2M + T2P;
  const T2R = T2H + T2Q;
  const T6z = T2H - T2Q;
  const T6C = T4q + T4r;
  const T6D = T4m + T4n;
  const T6E = T6C - T6D;
  const T7k = T6D + T6C;

  const T4p = T4l + T4o;
  const T4u = T4s - T4t;
  const T4v = KP707106781 * (T4p - T4u);
  const T5K = KP707106781 * (T4p + T4u);
  const T4C = T4t + T4s;
  const T4D = T4l - T4o;
  const T4E = KP707106781 * (T4C - T4D);
  const T5H = KP707106781 * (T4D + T4C);

  const T3i = cr[3], T3j = ci[3];
  const T3k = T3 * T3i + T6 * T3j; // FMA
  const T4S = T3 * T3j - T6 * T3i; // FNMS
  const T3m = cr[19], T3o = ci[19];
  const T3p = T3l * T3m + T3n * T3o; // FMA
  const T4T = T3l * T3o - T3n * T3m; // FNMS

  const T4R = T3k - T3p;
  const T4U = T4S - T4T;

  const T3r = cr[27], T3s = ci[27];
  const T3t = Th * T3r + Tl * T3s; // FMA
  const T4N = Th * T3s - Tl * T3r; // FNMS
  const T3u = cr[11], T3v = ci[11];
  const T3w = Tg * T3u + Tk * T3v; // FMA
  const T4O = Tg * T3v - Tk * T3u; // FNMS

  const T4M = T3t - T3w;
  const T4P = T4N - T4O;

  const T3q = T3k + T3p;
  const T3x = T3t + T3w;
  const T3y = T3q + T3x;
  const T6P = T3q - T3x;
  const T6I = T4N + T4O;
  const T6J = T4S + T4T;
  const T6K = T6I - T6J;
  const T7p = T6J + T6I;

  const T4Q = T4M + T4P;
  const T4V = T4R - T4U;
  const T4W = KP707106781 * (T4Q - T4V);
  const T5O = KP707106781 * (T4V + T4Q);
  const T53 = T4R + T4U;
  const T54 = T4P - T4M;
  const T55 = KP707106781 * (T53 - T54);
  const T5R = KP707106781 * (T53 + T54);

  const outCr = new Float64Array(32), outCi = new Float64Array(32);

  {
    const T1j = TL + T1i;
    const T2a = T1E + T29;
    const T2b = T1j + T2a;
    const T7x = T1j - T2a;
    const T7C = T7f + T7e;
    const T7J = T7D + T7I;
    const T7K = T7C + T7J;
    const T7M = T7J - T7C;

    const T2S = T2y + T2R;
    const T3z = T3h + T3y;
    const T3A = T2S + T3z;
    const T7L = T3z - T2S;
    const T7y = T7o + T7p;
    const T7z = T7j + T7k;
    const T7A = T7y - T7z;
    const T7B = T7z + T7y;

    outCi[15] = T2b - T3A;
    outCr[24] = T7L - T7M;
    outCi[23] = T7L + T7M;
    outCr[0] = T2b + T3A;
    outCr[8] = T7x - T7A;
    outCr[16] = T7B - T7K;
    outCi[31] = T7B + T7K;
    outCi[7] = T7x + T7A;
  }
  {
    const T5w = KP707106781 * (T5u + T5v);
    const T5x = T5t - T5w;
    const T5Z = T5t + T5w;
    const T89 = KP707106781 * (T3K - T3P);
    const T8d = T89 + T8c;
    const T8j = T8c - T89;

    const T5A = KP923879532 * T5y + KP382683432 * T5z; // FMA
    const T5D = KP382683432 * T5B - KP923879532 * T5C; // FNMS
    const T5E = T5A + T5D;
    const T88 = T5A - T5D;
    const T67 = T5N + T5O;
    const T68 = T5Q + T5R;
    const T69 = KP195090322 * T67 - KP980785280 * T68; // FNMS
    const T6d = KP980785280 * T67 + KP195090322 * T68; // FMA

    const T5I = T5G - T5H;
    const T5L = T5J - T5K;
    const T5M = KP831469612 * T5I + KP555570233 * T5L; // FMA
    const T5W = KP555570233 * T5I - KP831469612 * T5L; // FNMS
    const T60 = KP923879532 * T5z - KP382683432 * T5y; // FNMS
    const T61 = KP382683432 * T5C + KP923879532 * T5B; // FMA
    const T62 = T60 + T61;
    const T8i = T61 - T60;

    const T64 = T5G + T5H;
    const T65 = T5J + T5K;
    const T66 = KP195090322 * T64 + KP980785280 * T65; // FMA
    const T6c = KP980785280 * T64 - KP195090322 * T65; // FNMS
    const T5P = T5N - T5O;
    const T5S = T5Q - T5R;
    const T5T = KP831469612 * T5P - KP555570233 * T5S; // FNMS
    const T5X = KP555570233 * T5P + KP831469612 * T5S; // FMA

    const T5F = T5x + T5E;
    const T5U = T5M + T5T;
    outCi[12] = T5F - T5U;
    outCr[3] = T5F + T5U;
    const T8h = T5X - T5W;
    const T8k = T8i + T8j;
    outCr[19] = T8h - T8k;
    outCi[28] = T8h + T8k;

    const T8l = T5T - T5M;
    const T8m = T8j - T8i;
    outCr[27] = T8l - T8m;
    outCi[20] = T8l + T8m;
    const T5V = T5x - T5E;
    const T5Y = T5W + T5X;
    outCr[11] = T5V - T5Y;
    outCi[4] = T5V + T5Y;

    const T63 = T5Z - T62;
    const T6a = T66 + T69;
    outCi[8] = T63 - T6a;
    outCr[7] = T63 + T6a;
    const T87 = T69 - T66;
    const T8e = T88 + T8d;
    outCr[31] = T87 - T8e;
    outCi[16] = T87 + T8e;

    const T8f = T6d - T6c;
    const T8g = T8d - T88;
    outCr[23] = T8f - T8g;
    outCi[24] = T8f + T8g;
    const T6b = T5Z + T62;
    const T6e = T6c + T6d;
    outCr[15] = T6b - T6e;
    outCi[0] = T6b + T6e;
  }
  {
    const T7d = TL - T1i;
    const T7g = T7e - T7f;
    const T7h = T7d - T7g;
    const T7t = T7d + T7g;
    const T7O = T1E - T29;
    const T7P = T7I - T7D;
    const T7Q = T7O + T7P;
    const T7S = T7P - T7O;

    const T7i = T2y - T2R;
    const T7l = T7j - T7k;
    const T7m = T7i + T7l;
    const T7u = T7i - T7l;
    const T7n = T3h - T3y;
    const T7q = T7o - T7p;
    const T7r = T7n - T7q;
    const T7v = T7n + T7q;

    const T7s = KP707106781 * (T7m + T7r);
    outCi[11] = T7h - T7s;
    outCr[4] = T7h + T7s;
    const T7R = KP707106781 * (T7v - T7u);
    outCr[20] = T7R - T7S;
    outCi[27] = T7R + T7S;
    const T7w = KP707106781 * (T7u + T7v);
    outCr[12] = T7t - T7w;
    outCi[3] = T7t + T7w;
    const T7N = KP707106781 * (T7r - T7m);
    outCr[28] = T7N - T7Q;
    outCi[19] = T7N + T7Q;
  }
  {
    const T6j = T6f - T6i;
    const T7X = T7V + T7W;
    const T83 = T7W - T7V;
    const T6X = T6f + T6i;
    const T6o = T6k + T6n;
    const T6t = T6p - T6s;
    const T6u = KP707106781 * (T6o + T6t);
    const T7U = KP707106781 * (T6o - T6t);

    const T75 = T6O + T6P;
    const T76 = T6H + T6K;
    const T77 = KP382683432 * T75 + KP923879532 * T76; // FMA
    const T7b = KP382683432 * T76 - KP923879532 * T75; // FNMS
    const T6Y = T6k - T6n;
    const T6Z = T6p + T6s;
    const T70 = KP707106781 * (T6Y + T6Z);
    const T82 = KP707106781 * (T6Z - T6Y);

    const T6A = T6y - T6z;
    const T6F = T6B - T6E;
    const T6G = KP382683432 * T6A + KP923879532 * T6F; // FMA
    const T6U = KP382683432 * T6F - KP923879532 * T6A; // FNMS

    const T72 = T6B + T6E;
    const T73 = T6y + T6z;
    const T74 = KP923879532 * T72 - KP382683432 * T73; // FNMS
    const T7a = KP923879532 * T73 + KP382683432 * T72; // FMA
    const T6L = T6H - T6K;
    const T6Q = T6O - T6P;
    const T6R = KP923879532 * T6L - KP382683432 * T6Q; // FNMS
    const T6V = KP923879532 * T6Q + KP382683432 * T6L; // FMA

    const T6v = T6j + T6u;
    const T6S = T6G + T6R;
    outCi[13] = T6v - T6S;
    outCr[2] = T6v + T6S;
    const T81 = T6V - T6U;
    const T84 = T82 + T83;
    outCr[18] = T81 - T84;
    outCi[29] = T81 + T84;

    const T85 = T6R - T6G;
    const T86 = T83 - T82;
    outCr[26] = T85 - T86;
    outCi[21] = T85 + T86;
    const T6T = T6j - T6u;
    const T6W = T6U + T6V;
    outCr[10] = T6T - T6W;
    outCi[5] = T6T + T6W;

    const T71 = T6X + T70;
    const T78 = T74 + T77;
    outCr[14] = T71 - T78;
    outCi[1] = T71 + T78;
    const T7T = T7b - T7a;
    const T7Y = T7U + T7X;
    outCr[30] = T7T - T7Y;
    outCi[17] = T7T + T7Y;

    const T7Z = T77 - T74;
    const T80 = T7X - T7U;
    outCr[22] = T7Z - T80;
    outCi[25] = T7Z + T80;
    const T79 = T6X - T70;
    const T7c = T7a + T7b;
    outCi[9] = T79 - T7c;
    outCr[6] = T79 + T7c;
  }
  {
    const T3Q = KP707106781 * (T3K + T3P);
    const T3R = T3F - T3Q;
    const T5d = T3F + T3Q;
    const T8p = KP707106781 * (T5v - T5u);
    const T8r = T8p + T8q;
    const T8x = T8q - T8p;

    const T42 = KP382683432 * T3W - KP923879532 * T41; // FNMS
    const T4d = KP923879532 * T47 + KP382683432 * T4c; // FMA
    const T4e = T42 + T4d;
    const T8o = T4d - T42;
    const T5l = T52 + T55;
    const T5m = T4L + T4W;
    const T5n = KP980785280 * T5l - KP195090322 * T5m; // FNMS
    const T5r = KP195090322 * T5l + KP980785280 * T5m; // FMA

    const T4w = T4k - T4v;
    const T4F = T4B - T4E;
    const T4G = KP831469612 * T4w - KP555570233 * T4F; // FNMS
    const T5a = KP555570233 * T4w + KP831469612 * T4F; // FMA
    const T5e = KP382683432 * T41 + KP923879532 * T3W; // FMA
    const T5f = KP923879532 * T4c - KP382683432 * T47; // FNMS
    const T5g = T5e + T5f;
    const T8w = T5e - T5f;

    const T5i = T4B + T4E;
    const T5j = T4k + T4v;
    const T5k = KP195090322 * T5i + KP980785280 * T5j; // FMA
    const T5q = KP195090322 * T5j - KP980785280 * T5i; // FNMS
    const T4X = T4L - T4W;
    const T56 = T52 - T55;
    const T57 = KP555570233 * T4X + KP831469612 * T56; // FMA
    const T5b = KP555570233 * T56 - KP831469612 * T4X; // FNMS

    const T4f = T3R + T4e;
    const T58 = T4G + T57;
    outCr[13] = T4f - T58;
    outCi[2] = T4f + T58;
    const T8v = T5b - T5a;
    const T8y = T8w + T8x;
    outCr[29] = T8v - T8y;
    outCi[18] = T8v + T8y;

    const T8z = T57 - T4G;
    const T8A = T8x - T8w;
    outCr[21] = T8z - T8A;
    outCi[26] = T8z + T8A;
    const T59 = T3R - T4e;
    const T5c = T5a + T5b;
    outCi[10] = T59 - T5c;
    outCr[5] = T59 + T5c;

    const T5h = T5d + T5g;
    const T5o = T5k + T5n;
    outCi[14] = T5h - T5o;
    outCr[1] = T5h + T5o;
    const T8n = T5r - T5q;
    const T8s = T8o + T8r;
    outCr[17] = T8n - T8s;
    outCi[30] = T8n + T8s;

    const T8t = T5n - T5k;
    const T8u = T8r - T8o;
    outCr[25] = T8t - T8u;
    outCi[22] = T8t + T8u;
    const T5p = T5d - T5g;
    const T5s = T5q + T5r;
    outCr[9] = T5p - T5s;
    outCi[6] = T5p + T5s;
  }

  return [outCr, outCi];
}

module.exports = { hf2_32 };
