'use strict';

// =============================================================================
// hb2_32.js -- faithful JS port of rdft/scalar/r2cb/hb2_32.c (non-FMA
// branch). Alternate-codegen EVEN-radix (r=32) backward twiddle codelet --
// same role as hb2_16.js, larger radix (see hf2_32.js's header for the
// forward-direction analogue). twinstr trig-generates W^1, W^3, W^9, W^27.
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hb2_16.js -- note the
// cross-indexed cr[k]/ci[k'] pairing pattern, transcribed literally.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP555570233 = 0.555570233019602224742830813948532874374937191;
const KP831469612 = 0.831469612302545237078788377617905756738560812;
const KP980785280 = 0.980785280403230449126182236134239036973933731;
const KP195090322 = 0.195090322016128267848284868477022240927691618;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function hb2_32(cr, ci, Wc, Ws) {
  const T11 = Wc[1], T14 = Ws[1], T12 = Wc[3], T15 = Ws[3];
  const T13 = T11 * T12, T1b = T14 * T12, T16 = T14 * T15, T1a = T11 * T15;
  const T17 = T13 + T16, T2z = T13 - T16, T2B = T1a + T1b, T1c = T1a - T1b;
  const T18 = Wc[9], T1d = Ws[9];
  const T2P = T12 * T18, T3q = T14 * T18, T2U = T15 * T18, T3l = T11 * T18;
  const T2Q = T15 * T1d, T3p = T11 * T1d, T2T = T12 * T1d, T3m = T14 * T1d;
  const T1g = Wc[27], T1k = Ws[27];
  const T2D = T11 * T1g, T3g = T15 * T1g, T2K = T14 * T1g, T39 = T12 * T1g;
  const T2E = T14 * T1k, T3f = T12 * T1k, T2J = T11 * T1k, T3a = T15 * T1k;

  const T2F = T2D - T2E;
  const T2L = T2J + T2K;
  const T3t = T39 - T3a;
  const T4H = T2J - T2K;
  const T3h = T3f - T3g;
  const T3V = T3f + T3g;
  const T3b = T39 + T3a;
  const T4v = T2D + T2E;
  const T4T = T18 * T1g + T1d * T1k; // FMA
  const T4X = T18 * T1k - T1d * T1g; // FNMS

  const T6r = T17 * T1g, T6s = T1c * T1k;
  const T6t = T6r - T6s, T71 = T6r + T6s;
  const T6x = T17 * T1k, T6y = T1c * T1g;
  const T6z = T6x + T6y, T75 = T6x - T6y;

  const T7Z = T2z * T1g, T80 = T2B * T1k;
  const T81 = T7Z + T80, T8x = T7Z - T80;
  const T8d = T2z * T1k, T8e = T2B * T1g;
  const T8f = T8d - T8e, T8z = T8d + T8e;
  const T2R = T2P - T2Q;
  const T2V = T2T + T2U;
  const T8p = T2R * T1g + T2V * T1k; // FMA
  const T8t = T2R * T1k - T2V * T1g; // FNMS

  const T4r = T2P + T2Q;
  const T4t = T2T - T2U;
  const T53 = T4r * T1g + T4t * T1k; // FMA
  const T69 = T4r * T1k - T4t * T1g; // FNMS
  const T3n = T3l + T3m;
  const T3r = T3p - T3q;
  const T7P = T3n * T1g + T3r * T1k; // FMA
  const T7T = T3n * T1k - T3r * T1g; // FNMS
  const T4P = T3l - T3m;
  const T4R = T3p + T3q;
  const T6F = T4P * T1g + T4R * T1k; // FMA
  const T6R = T4P * T1k - T4R * T1g; // FNMS

  const T19 = T17 * T18, T1e = T1c * T1d;
  const T1f = T19 + T1e, T2X = T19 - T1e;
  const T1h = T17 * T1d, T1i = T1c * T18;
  const T1j = T1h - T1i, T2Y = T1h + T1i;

  const T1l = T1f * T1g + T1j * T1k; // FMA
  const T31 = T2X * T1k - T2Y * T1g; // FNMS
  const T2d = T1f * T1k - T1j * T1g; // FNMS
  const T2Z = T2X * T1g + T2Y * T1k; // FMA

  const T47 = T2z * T18, T48 = T2B * T1d;
  const T49 = T47 - T48, T4h = T47 + T48;
  const T4a = T2z * T1d, T4b = T2B * T18;
  const T4c = T4a + T4b, T4i = T4a - T4b;

  const T4d = T49 * T1g + T4c * T1k; // FMA
  const T4n = T4h * T1k - T4i * T1g; // FNMS
  const T4f = T49 * T1k - T4c * T1g; // FNMS
  const T4j = T4h * T1g + T4i * T1k; // FMA

  const T1 = cr[0];
  const T2 = ci[15];
  const T3 = T1 + T2;
  const T54 = T1 - T2;
  const T2m = ci[27];
  const T2n = cr[20];
  const T2o = T2m - T2n;
  const T58 = T2m + T2n;

  const T2p = ci[19];
  const T2q = cr[28];
  const T2r = T2p - T2q;
  const T5b = T2p + T2q;
  const T4 = cr[8];
  const T5 = ci[7];
  const T6 = T4 + T5;
  const T6a = T4 - T5;

  const T8 = cr[4];
  const T9 = ci[11];
  const Ta = T8 + T9;
  const T57 = T8 - T9;
  const T2f = ci[31];
  const T2g = cr[16];
  const T2h = T2f - T2g;
  const T6b = T2f + T2g;

  const T2i = ci[23];
  const T2j = cr[24];
  const T2k = T2i - T2j;
  const T55 = T2i + T2j;
  const Tb = ci[3];
  const Tc = cr[12];
  const Td = Tb + Tc;
  const T5a = Tb - Tc;

  const T56 = T54 - T55;
  const T7b = T54 + T55;
  const T7C = T6b - T6a;
  const T6c = T6a + T6b;
  const T7 = T3 + T6;
  const Te = Ta + Td;
  const Tf = T7 + Te;
  const T1m = T7 - Te;

  const T6d = T57 + T58;
  const T6e = T5a + T5b;
  const T6f = KP707106781 * (T6d - T6e);
  const T7c = KP707106781 * (T6d + T6e);
  const T3W = T2h - T2k;
  const T3X = Ta - Td;
  const T3Y = T3W - T3X;
  const T4I = T3X + T3W;

  const T2l = T2h + T2k;
  const T2s = T2o + T2r;
  const T2t = T2l - T2s;
  const T32 = T2l + T2s;

  const T59 = T57 - T58;
  const T5c = T5a - T5b;
  const T5d = KP707106781 * (T59 + T5c);
  const T7D = KP707106781 * (T59 - T5c);
  const T3u = T3 - T6;
  const T3v = T2r - T2o;
  const T3w = T3u - T3v;
  const T4w = T3u + T3v;

  const Tg = cr[2];
  const Th = ci[13];
  const Ti = Tg + Th;
  const T5p = Tg - Th;
  const T1u = ci[29];
  const T1v = cr[18];
  const T1w = T1u - T1v;
  const T5n = T1u + T1v;

  const T1x = ci[21];
  const T1y = cr[26];
  const T1z = T1x - T1y;
  const T5q = T1x + T1y;
  const Tj = cr[10];
  const Tk = ci[5];
  const Tl = Tj + Tk;
  const T5m = Tj - Tk;

  const Tn = ci[1];
  const To = cr[14];
  const Tp = Tn + To;
  const T5i = Tn - To;
  const T1n = ci[17];
  const T1o = cr[30];
  const T1p = T1n - T1o;
  const T5g = T1n + T1o;

  const T1q = ci[25];
  const T1r = cr[22];
  const T1s = T1q - T1r;
  const T5j = T1q + T1r;
  const Tq = cr[6];
  const Tr = ci[9];
  const Ts = Tq + Tr;
  const T5f = Tq - Tr;

  const Tm = Ti + Tl;
  const Tt = Tp + Ts;
  const Tu = Tm + Tt;
  const T2e = Tm - Tt;
  const T7e = T5p + T5q;
  const T7f = T5n - T5m;
  const T7g = KP382683432 * T7e - KP923879532 * T7f; // FNMS
  const T7F = KP923879532 * T7e + KP382683432 * T7f; // FMA

  const T7h = T5i + T5j;
  const T7i = T5f + T5g;
  const T7j = KP382683432 * T7h - KP923879532 * T7i; // FNMS
  const T7G = KP923879532 * T7h + KP382683432 * T7i; // FMA
  const T1t = T1p + T1s;
  const T1A = T1w + T1z;
  const T1B = T1t - T1A;
  const T33 = T1A + T1t;

  const T3x = T1p - T1s;
  const T3y = Tp - Ts;
  const T3z = T3x - T3y;
  const T40 = T3y + T3x;
  const T5h = T5f - T5g;
  const T5k = T5i - T5j;
  const T5l = KP923879532 * T5h - KP382683432 * T5k; // FNMS
  const T6i = KP923879532 * T5k + KP382683432 * T5h; // FMA

  const T5o = T5m + T5n;
  const T5r = T5p - T5q;
  const T5s = KP923879532 * T5o + KP382683432 * T5r; // FMA
  const T6h = KP923879532 * T5r - KP382683432 * T5o; // FNMS
  const T3A = Ti - Tl;
  const T3B = T1w - T1z;
  const T3C = T3A + T3B;
  const T3Z = T3A - T3B;

  const Tw = cr[1];
  const Tx = ci[14];
  const Ty = Tw + Tx;
  const T5v = Tw - Tx;
  const Tz = cr[9];
  const TA = ci[6];
  const TB = Tz + TA;
  const T5G = Tz - TA;
  const T1H = ci[22];
  const T1I = cr[25];
  const T1J = T1H - T1I;
  const T5w = T1H + T1I;

  const T1E = ci[30];
  const T1F = cr[17];
  const T1G = T1E - T1F;
  const T5H = T1E + T1F;

  const TG = ci[2];
  const TH = cr[13];
  const T5B = TG - TH;
  const T1O = ci[18];
  const T1P = cr[29];
  const T5C = T1O + T1P;
  const TI = TG + TH;
  const T5K = T5B + T5C;
  const T1Q = T1O - T1P;
  const T5D = T5B - T5C;

  const TD = cr[5];
  const TE = ci[10];
  const T5y = TD - TE;
  const T1L = ci[26];
  const T1M = cr[21];
  const T5z = T1L + T1M;
  const TF = TD + TE;
  const T5J = T5y + T5z;
  const T1N = T1L - T1M;
  const T5A = T5y - T5z;

  const TC = Ty + TB;
  const TJ = TF + TI;
  const TK = TC + TJ;
  const T1D = TC - TJ;
  const T7t = T5H - T5G;
  const T7u = KP707106781 * (T5A - T5D);
  const T7v = T7t + T7u;
  const T86 = T7t - T7u;

  const T7w = T5v + T5w;
  const T7x = KP707106781 * (T5J + T5K);
  const T7y = T7w - T7x;
  const T85 = T7w + T7x;
  const T1K = T1G + T1J;
  const T1R = T1N + T1Q;
  const T1S = T1K - T1R;
  const T35 = T1K + T1R;

  const T3M = T1G - T1J;
  const T3N = TF - TI;
  const T3O = T3M - T3N;
  const T4C = T3N + T3M;
  const T5x = T5v - T5w;
  const T5E = KP707106781 * (T5A + T5D);
  const T5F = T5x - T5E;
  const T6J = T5x + T5E;

  const T5I = T5G + T5H;
  const T5L = KP707106781 * (T5J - T5K);
  const T5M = T5I - T5L;
  const T6K = T5I + T5L;
  const T3P = Ty - TB;
  const T3Q = T1Q - T1N;
  const T3R = T3P - T3Q;
  const T4D = T3P + T3Q;

  const TL = ci[0];
  const TM = cr[15];
  const TN = TL + TM;
  const T5O = TL - TM;

  const TO = cr[7];
  const TP = ci[8];
  const TQ = TO + TP;
  const T5Z = TO - TP;
  const T1Y = ci[24];
  const T1Z = cr[23];
  const T20 = T1Y - T1Z;
  const T5P = T1Y + T1Z;

  const T1V = ci[16];
  const T1W = cr[31];
  const T1X = T1V - T1W;
  const T60 = T1V + T1W;

  const TV = ci[4];
  const TW = cr[11];
  const T5U = TV - TW;
  const T25 = ci[20];
  const T26 = cr[27];
  const T5V = T25 + T26;
  const TX = TV + TW;
  const T63 = T5U + T5V;
  const T27 = T25 - T26;
  const T5W = T5U - T5V;

  const TS = cr[3];
  const TT = ci[12];
  const T5R = TS - TT;
  const T22 = ci[28];
  const T23 = cr[19];
  const T5S = T22 + T23;
  const TU = TS + TT;
  const T62 = T5R + T5S;
  const T24 = T22 - T23;
  const T5T = T5R - T5S;

  const TR = TN + TQ;
  const TY = TU + TX;
  const TZ = TR + TY;
  const T1U = TR - TY;
  const T7m = KP707106781 * (T5T - T5W);
  const T7n = T5Z + T60;
  const T7o = T7m - T7n;
  const T89 = T7n + T7m;

  const T7p = T5O + T5P;
  const T7q = KP707106781 * (T62 + T63);
  const T7r = T7p - T7q;
  const T88 = T7p + T7q;
  const T21 = T1X + T20;
  const T28 = T24 + T27;
  const T29 = T21 - T28;
  const T36 = T21 + T28;

  const T3F = T1X - T20;
  const T3G = TU - TX;
  const T3H = T3F - T3G;
  const T4z = T3G + T3F;
  const T5Q = T5O - T5P;
  const T5X = KP707106781 * (T5T + T5W);
  const T5Y = T5Q - T5X;
  const T6M = T5Q + T5X;

  const T61 = T5Z - T60;
  const T64 = KP707106781 * (T62 - T63);
  const T65 = T61 - T64;
  const T6N = T61 + T64;
  const T3I = TN - TQ;
  const T3J = T27 - T24;
  const T3K = T3I - T3J;
  const T4A = T3I + T3J;

  const outCr = new Float64Array(32), outCi = new Float64Array(32);

  {
    const Tv = Tf + Tu;
    const T10 = TK + TZ;
    const T30 = Tv - T10;
    const T34 = T32 + T33;
    const T37 = T35 + T36;
    const T38 = T34 - T37;
    outCr[0] = Tv + T10;
    outCi[0] = T34 + T37;
    outCr[16] = T2Z * T30 - T31 * T38; // FNMS
    outCi[16] = T31 * T30 + T2Z * T38; // FMA
  }
  {
    const T3c = Tf - Tu;
    const T3d = T36 - T35;
    const T3e = T3c - T3d;
    const T3o = T3c + T3d;
    const T3i = T32 - T33;
    const T3j = TK - TZ;
    const T3k = T3i - T3j;
    const T3s = T3j + T3i;

    outCr[24] = T3b * T3e - T3h * T3k; // FNMS
    outCi[24] = T3b * T3k + T3h * T3e; // FMA
    outCr[8] = T3n * T3o - T3r * T3s; // FNMS
    outCi[8] = T3n * T3s + T3r * T3o; // FMA
  }
  {
    const T1C = T1m + T1B;
    const T2u = T2e + T2t;
    const T2M = T2t - T2e;
    const T2G = T1m - T1B;

    const T2v = T1D + T1S;
    const T2w = T29 - T1U;
    const T2x = KP707106781 * (T2v + T2w);
    const T2H = KP707106781 * (T2w - T2v);
    const T1T = T1D - T1S;
    const T2a = T1U + T29;
    const T2b = KP707106781 * (T1T + T2a);
    const T2N = KP707106781 * (T1T - T2a);

    const T2c = T1C - T2b;
    const T2y = T2u - T2x;
    outCr[20] = T1l * T2c - T2d * T2y; // FNMS
    outCi[20] = T2d * T2c + T1l * T2y; // FMA
    const T2S = T2G + T2H;
    const T2W = T2M + T2N;
    outCr[12] = T2R * T2S - T2V * T2W; // FNMS
    outCi[12] = T2R * T2W + T2V * T2S; // FMA

    const T2A = T1C + T2b;
    const T2C = T2u + T2x;
    outCr[4] = T2z * T2A - T2B * T2C; // FNMS
    outCi[4] = T2B * T2A + T2z * T2C; // FMA
    const T2I = T2G - T2H;
    const T2O = T2M - T2N;
    outCr[28] = T2F * T2I - T2L * T2O; // FNMS
    outCi[28] = T2F * T2O + T2L * T2I; // FMA
  }
  {
    const T4x = KP707106781 * (T3Z + T40);
    const T4y = T4w - T4x;
    const T4U = T4w + T4x;
    const T4J = KP707106781 * (T3C + T3z);
    const T4K = T4I - T4J;
    const T4Y = T4I + T4J;

    const T4B = KP923879532 * T4z - KP382683432 * T4A; // FNMS
    const T4E = KP923879532 * T4C + KP382683432 * T4D; // FMA
    const T4F = T4B - T4E;
    const T4Z = T4E + T4B;
    const T4L = KP923879532 * T4D - KP382683432 * T4C; // FNMS
    const T4M = KP382683432 * T4z + KP923879532 * T4A; // FMA
    const T4N = T4L - T4M;
    const T4V = T4L + T4M;

    const T4G = T4y - T4F;
    const T4O = T4K - T4N;
    outCr[26] = T4v * T4G - T4H * T4O; // FNMS
    outCi[26] = T4H * T4G + T4v * T4O; // FMA
    const T51 = T4U + T4V;
    const T52 = T4Y + T4Z;
    outCr[2] = T17 * T51 - T1c * T52; // FNMS
    outCi[2] = T17 * T52 + T1c * T51; // FMA

    const T4Q = T4y + T4F;
    const T4S = T4K + T4N;
    outCr[10] = T4P * T4Q - T4R * T4S; // FNMS
    outCi[10] = T4R * T4Q + T4P * T4S; // FMA
    const T4W = T4U - T4V;
    const T50 = T4Y - T4Z;
    outCr[18] = T4T * T4W - T4X * T50; // FNMS
    outCi[18] = T4T * T50 + T4X * T4W; // FMA
  }
  {
    const T3D = KP707106781 * (T3z - T3C);
    const T3E = T3w - T3D;
    const T4k = T3w + T3D;
    const T41 = KP707106781 * (T3Z - T40);
    const T42 = T3Y - T41;
    const T4o = T3Y + T41;

    const T3L = KP382683432 * T3H - KP923879532 * T3K; // FNMS
    const T3S = KP382683432 * T3O + KP923879532 * T3R; // FMA
    const T3T = T3L - T3S;
    const T4p = T3S + T3L;
    const T43 = KP382683432 * T3R - KP923879532 * T3O; // FNMS
    const T44 = KP923879532 * T3H + KP382683432 * T3K; // FMA
    const T45 = T43 - T44;
    const T4l = T43 + T44;

    const T3U = T3E - T3T;
    const T46 = T42 - T45;
    outCr[30] = T3t * T3U - T3V * T46; // FNMS
    outCi[30] = T3V * T3U + T3t * T46; // FMA
    const T4s = T4k + T4l;
    const T4u = T4o + T4p;
    outCr[6] = T4r * T4s - T4t * T4u; // FNMS
    outCi[6] = T4r * T4u + T4t * T4s; // FMA

    const T4e = T3E + T3T;
    const T4g = T42 + T45;
    outCr[14] = T4d * T4e - T4f * T4g; // FNMS
    outCi[14] = T4f * T4e + T4d * T4g; // FMA
    const T4m = T4k - T4l;
    const T4q = T4o - T4p;
    outCr[22] = T4j * T4m - T4n * T4q; // FNMS
    outCi[22] = T4j * T4q + T4n * T4m; // FMA
  }
  {
    const T6G = T56 + T5d;
    const T6H = T6h + T6i;
    const T6I = T6G + T6H;
    const T72 = T6G - T6H;
    const T6V = KP195090322 * T6J + KP980785280 * T6K; // FMA
    const T6W = KP980785280 * T6N - KP195090322 * T6M; // FNMS
    const T6X = T6V + T6W;
    const T73 = T6W - T6V;

    const T6L = KP980785280 * T6J - KP195090322 * T6K; // FNMS
    const T6O = KP980785280 * T6M + KP195090322 * T6N; // FMA
    const T6P = T6L + T6O;
    const T77 = T6L - T6O;
    const T6S = T6c + T6f;
    const T6T = T5s + T5l;
    const T6U = T6S + T6T;
    const T76 = T6S - T6T;

    const T6Q = T6I - T6P;
    const T6Y = T6U - T6X;
    outCr[17] = T6F * T6Q - T6R * T6Y; // FNMS
    outCi[17] = T6R * T6Q + T6F * T6Y; // FMA
    const T79 = T72 + T73;
    const T7a = T76 + T77;
    outCr[9] = T18 * T79 - T1d * T7a; // FNMS
    outCi[9] = T18 * T7a + T1d * T79; // FMA

    const T6Z = T6I + T6P;
    const T70 = T6U + T6X;
    outCr[1] = T11 * T6Z - T14 * T70; // FNMS
    outCi[1] = T14 * T6Z + T11 * T70; // FMA
    const T74 = T72 - T73;
    const T78 = T76 - T77;
    outCr[25] = T71 * T74 - T75 * T78; // FNMS
    outCi[25] = T71 * T78 + T75 * T74; // FMA
  }
  {
    const T82 = T7b + T7c;
    const T83 = T7F + T7G;
    const T84 = T82 - T83;
    const T8q = T82 + T83;
    const T8j = KP195090322 * T86 + KP980785280 * T85; // FMA
    const T8k = KP195090322 * T89 + KP980785280 * T88; // FMA
    const T8l = T8j - T8k;
    const T8r = T8j + T8k;

    const T87 = KP195090322 * T85 - KP980785280 * T86; // FNMS
    const T8a = KP195090322 * T88 - KP980785280 * T89; // FNMS
    const T8b = T87 + T8a;
    const T8v = T87 - T8a;
    const T8g = T7C - T7D;
    const T8h = T7g - T7j;
    const T8i = T8g + T8h;
    const T8u = T8g - T8h;

    const T8c = T84 - T8b;
    const T8m = T8i - T8l;
    outCr[23] = T81 * T8c - T8f * T8m; // FNMS
    outCi[23] = T8f * T8c + T81 * T8m; // FMA
    const T8y = T8q + T8r;
    const T8A = T8u - T8v;
    outCr[31] = T8x * T8y - T8z * T8A; // FNMS
    outCi[31] = T8x * T8A + T8z * T8y; // FMA

    const T8n = T84 + T8b;
    const T8o = T8i + T8l;
    outCr[7] = T1f * T8n - T1j * T8o; // FNMS
    outCi[7] = T1j * T8n + T1f * T8o; // FMA
    const T8s = T8q - T8r;
    const T8w = T8u + T8v;
    outCr[15] = T8p * T8s - T8t * T8w; // FNMS
    outCi[15] = T8p * T8w + T8t * T8s; // FMA
  }
  {
    const T5e = T56 - T5d;
    const T5t = T5l - T5s;
    const T5u = T5e + T5t;
    const T6u = T5e - T5t;
    const T6l = KP831469612 * T5F + KP555570233 * T5M; // FMA
    const T6m = KP555570233 * T65 - KP831469612 * T5Y; // FNMS
    const T6n = T6l + T6m;
    const T6v = T6m - T6l;

    const T5N = KP555570233 * T5F - KP831469612 * T5M; // FNMS
    const T66 = KP555570233 * T5Y + KP831469612 * T65; // FMA
    const T67 = T5N + T66;
    const T6B = T5N - T66;
    const T6g = T6c - T6f;
    const T6j = T6h - T6i;
    const T6k = T6g + T6j;
    const T6A = T6g - T6j;

    const T68 = T5u - T67;
    const T6o = T6k - T6n;
    outCr[21] = T53 * T68 - T69 * T6o; // FNMS
    outCi[21] = T69 * T68 + T53 * T6o; // FMA
    const T6D = T6u + T6v;
    const T6E = T6A + T6B;
    outCr[13] = T49 * T6D - T4c * T6E; // FNMS
    outCi[13] = T49 * T6E + T4c * T6D; // FMA

    const T6p = T5u + T67;
    const T6q = T6k + T6n;
    outCr[5] = T4h * T6p - T4i * T6q; // FNMS
    outCi[5] = T4i * T6p + T4h * T6q; // FMA
    const T6w = T6u - T6v;
    const T6C = T6A - T6B;
    outCr[29] = T6t * T6w - T6z * T6C; // FNMS
    outCi[29] = T6t * T6C + T6z * T6w; // FMA
  }
  {
    const T7d = T7b - T7c;
    const T7k = T7g + T7j;
    const T7l = T7d - T7k;
    const T7Q = T7d + T7k;
    const T7J = KP831469612 * T7y - KP555570233 * T7v; // FNMS
    const T7K = KP555570233 * T7o + KP831469612 * T7r; // FMA
    const T7L = T7J - T7K;
    const T7R = T7J + T7K;

    const T7s = KP831469612 * T7o - KP555570233 * T7r; // FNMS
    const T7z = KP831469612 * T7v + KP555570233 * T7y; // FMA
    const T7A = T7s - T7z;
    const T7V = T7z + T7s;
    const T7E = T7C + T7D;
    const T7H = T7F - T7G;
    const T7I = T7E - T7H;
    const T7U = T7E + T7H;

    const T7B = T7l - T7A;
    const T7M = T7I - T7L;
    outCr[27] = T1g * T7B - T1k * T7M; // FNMS
    outCi[27] = T1k * T7B + T1g * T7M; // FMA
    const T7X = T7Q + T7R;
    const T7Y = T7U + T7V;
    outCr[3] = T12 * T7X - T15 * T7Y; // FNMS
    outCi[3] = T12 * T7Y + T15 * T7X; // FMA

    const T7N = T7l + T7A;
    const T7O = T7I + T7L;
    outCr[11] = T2X * T7N - T2Y * T7O; // FNMS
    outCi[11] = T2Y * T7N + T2X * T7O; // FMA
    const T7S = T7Q - T7R;
    const T7W = T7U - T7V;
    outCr[19] = T7P * T7S - T7T * T7W; // FNMS
    outCi[19] = T7P * T7W + T7T * T7S; // FMA
  }

  return [outCr, outCi];
}

module.exports = { hb2_32 };
