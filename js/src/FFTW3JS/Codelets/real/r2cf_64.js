'use strict';

// =============================================================================
// r2cf_64.js -- faithful JS port of rdft/scalar/r2cf/r2cf_64.c (non-FMA
// branch). x[0..63] (real) -> O[0..63] packed halfcomplex: O[0..32]=Re0..
// Re32, O[33..63]=Im31..Im1 (O[n-k]=Im_k). R0/R1 (stride-2) input
// convention -- see r2cf_6.js's header.
//
// NOT WIRED into RealEngine1D.js's r2cf registry -- deliberately. Checked
// fftw_fprint_plan for n=64 and every multiple of 64 tried (192, 320, 448,
// 576, 704, 1600): real FFTW's FFTW_ESTIMATE planner NEVER selects the bare
// direct r2cf_64 codelet -- it always prefers an even-radix-outer CT split
// (hf2_8/hf2_16/hf2_20/hf2_25/hf2_32, an alternate-codegen EVEN-radix
// twiddle-codelet family this port hasn't built at all) over treating 64
// as an opaque direct block, since 64=2^6 is too composite for the direct
// codelet to ever win FFTW's own cost comparison. This means r2cf_64 can
// NEVER be validated bit-exact against real compiled FFTW3 through the
// normal planning path (there is no N where ground truth actually invokes
// it), and wiring it into chooseReal's candidate pool would make
// chooseReal('R2HC',64) wrongly predict "r2cf_64 direct" (confirmed: our
// current cost formula, which doesn't know about the hf2 family, ranks it
// as cheapest) -- an unverifiable, likely-wrong structural claim. Kept in
// the tree unwired (same "park it, document it" treatment as hf2_25) in
// case a future session's hf2/even-radix-outer work turns up a genuine
// use as a recursion leaf; do not wire this in without first confirming a
// specific N where fftw_fprint_plan shows real FFTW actually choosing it.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP773010453 = 0.773010453362736960810906609758469800971041293;
const KP634393284 = 0.634393284163645498215171613225493370675687095;
const KP098017140 = 0.098017140329560601994195563888641845861136673;
const KP995184726 = 0.995184726672196886244836953109479921575474869;
const KP290284677 = 0.290284677254462367636192375817395274691476278;
const KP956940335 = 0.956940335732208864935797886980269969482849206;
const KP471396736 = 0.471396736825997648556387625905254377657460319;
const KP881921264 = 0.881921264348355029712756863660388349508442621;
const KP195090322 = 0.195090322016128267848284868477022240927691618;
const KP980785280 = 0.980785280403230449126182236134239036973933731;
const KP555570233 = 0.555570233019602224742830813948532874374937191;
const KP831469612 = 0.831469612302545237078788377617905756738560812;
const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function r2cf_64(x) {
  const T1 = x[0];
  const T2 = x[32];
  const T3 = T1 + T2;
  const T11 = T1 - T2;
  const Tb = x[56];
  const Tc = x[24];
  const Td = Tb + Tc;
  const T13 = Tb - Tc;

  const T4 = x[16];
  const T5 = x[48];
  const T6 = T4 + T5;
  const T2S = T4 - T5;
  const T8 = x[8];
  const T9 = x[40];
  const Ta = T8 + T9;
  const T12 = T8 - T9;

  const T4l = T3 - T6;
  const T5a = Td - Ta;
  const T14 = KP707106781 * (T12 + T13);
  const T15 = T11 + T14;
  const T3n = T11 - T14;
  const T2R = KP707106781 * (T13 - T12);
  const T2T = T2R - T2S;
  const T3Q = T2S + T2R;
  const T7 = T3 + T6;
  const Te = Ta + Td;
  const Tf = T7 + Te;

  const T1N = x[57];
  const T1O = x[25];
  const T1P = T1N - T1O;
  const T4J = T1N + T1O;
  const T1Z = x[1];
  const T20 = x[33];
  const T21 = T1Z - T20;
  const T4y = T1Z + T20;

  const T1Q = x[9];
  const T1R = x[41];
  const T1S = T1Q - T1R;
  const T4K = T1Q + T1R;
  const T1U = x[17];
  const T1V = x[49];
  const T1W = T1U - T1V;
  const T4z = T1U + T1V;

  const T4A = T4y - T4z;
  const T4L = T4J - T4K;

  const T1T = KP707106781 * (T1P - T1S);
  const T1X = T1T - T1W;
  const T3B = T1W + T1T;
  const T22 = KP707106781 * (T1S + T1P);
  const T23 = T21 + T22;
  const T3y = T21 - T22;
  const T5G = T4y + T4z;
  const T5H = T4K + T4J;
  const T5I = T5G + T5H;
  const T66 = T5G - T5H;

  const T29 = x[63];
  const T2a = x[31];
  const T2b = T29 - T2a;
  const T4P = T29 + T2a;
  const T2E = x[15];
  const T2F = x[47];
  const T2G = T2E - T2F;
  const T4Q = T2E + T2F;

  const T2c = x[7];
  const T2d = x[39];
  const T2e = T2c - T2d;
  const T51 = T2c + T2d;
  const T2f = x[55];
  const T2g = x[23];
  const T2h = T2f - T2g;
  const T50 = T2f + T2g;

  const T4R = T4P - T4Q;
  const T52 = T50 - T51;

  const T2i = KP707106781 * (T2e + T2h);
  const T2j = T2b + T2i;
  const T3F = T2b - T2i;
  const T2D = KP707106781 * (T2h - T2e);
  const T2H = T2D - T2G;
  const T3I = T2G + T2D;
  const T5N = T4P + T4Q;
  const T5O = T51 + T50;
  const T5P = T5N + T5O;
  const T69 = T5N - T5O;

  const TL = x[62];
  const TM = x[30];
  const TN = TL + TM;
  const T1e = TL - TM;
  const TV = x[54];
  const TW = x[22];
  const TX = TV + TW;
  const T1g = TV - TW;

  const TO = x[14];
  const TP = x[46];
  const TQ = TO + TP;
  const T1k = TO - TP;
  const TS = x[6];
  const TT = x[38];
  const TU = TS + TT;
  const T1f = TS - TT;

  const T1h = KP707106781 * (T1f + T1g);
  const T1i = T1e + T1h;
  const T3t = T1e - T1h;
  const T1j = KP707106781 * (T1g - T1f);
  const T1l = T1j - T1k;
  const T3u = T1k + T1j;

  const TR = TN + TQ;
  const TY = TU + TX;
  const TZ = TR + TY;
  const T63 = TR - TY;
  const T4t = TN - TQ;
  const T4u = TX - TU;
  const T4v = KP923879532 * T4t - KP382683432 * T4u; // FNMS
  const T58 = KP382683432 * T4t + KP923879532 * T4u; // FMA

  const Tw = x[2];
  const Tx = x[34];
  const Ty = Tw + Tx;
  const T1s = Tw - Tx;
  const TG = x[58];
  const TH = x[26];
  const TI = TG + TH;
  const T1n = TG - TH;

  const Tz = x[18];
  const TA = x[50];
  const TB = Tz + TA;
  const T1q = Tz - TA;
  const TD = x[10];
  const TE = x[42];
  const TF = TD + TE;
  const T1o = TD - TE;

  const T1p = KP707106781 * (T1n - T1o);
  const T1r = T1p - T1q;
  const T3r = T1q + T1p;
  const T1t = KP707106781 * (T1o + T1n);
  const T1u = T1s + T1t;
  const T3q = T1s - T1t;

  const TC = Ty + TB;
  const TJ = TF + TI;
  const TK = TC + TJ;
  const T62 = TC - TJ;
  const T4q = Ty - TB;
  const T4r = TI - TF;
  const T4s = KP923879532 * T4q + KP382683432 * T4r; // FMA
  const T57 = KP923879532 * T4r - KP382683432 * T4q; // FNMS

  const Tg = x[4];
  const Th = x[36];
  const Ti = Tg + Th;
  const T16 = Tg - Th;
  const Tq = x[12];
  const Tr = x[44];
  const Ts = Tq + Tr;
  const T1a = Tq - Tr;

  const Tj = x[20];
  const Tk = x[52];
  const Tl = Tj + Tk;
  const T17 = Tj - Tk;
  const Tn = x[60];
  const To = x[28];
  const Tp = Tn + To;
  const T19 = Tn - To;

  const Tm = Ti + Tl;
  const Tt = Tp + Ts;
  const Tu = Tm + Tt;
  const T4m = Ti - Tl;
  const T4n = Tp - Ts;
  const T4o = KP707106781 * (T4m + T4n);
  const T5b = KP707106781 * (T4n - T4m);

  const T18 = KP923879532 * T16 - KP382683432 * T17; // FNMS
  const T1b = KP923879532 * T19 + KP382683432 * T1a; // FMA
  const T1c = T18 + T1b;
  const T3R = T1b - T18;
  const T2O = KP382683432 * T19 - KP923879532 * T1a; // FNMS
  const T2P = KP382683432 * T16 + KP923879532 * T17; // FMA
  const T2Q = T2O - T2P;
  const T3o = T2P + T2O;

  const T1y = x[61];
  const T1z = x[29];
  const T1A = T1y - T1z;
  const T4E = T1y + T1z;
  const T1I = x[21];
  const T1J = x[53];
  const T1K = T1I - T1J;
  const T4C = T1I + T1J;

  const T1B = x[13];
  const T1C = x[45];
  const T1D = T1B - T1C;
  const T4F = T1B + T1C;
  const T1F = x[5];
  const T1G = x[37];
  const T1H = T1F - T1G;
  const T4B = T1F + T1G;

  const T1E = KP382683432 * T1A - KP923879532 * T1D; // FNMS
  const T1L = KP382683432 * T1H + KP923879532 * T1K; // FMA
  const T1M = T1E - T1L;
  const T3z = T1L + T1E;
  const T5J = T4B + T4C;
  const T5K = T4E + T4F;
  const T5L = T5J + T5K;
  const T67 = T5K - T5J;

  const T24 = KP923879532 * T1H - KP382683432 * T1K; // FNMS
  const T25 = KP923879532 * T1A + KP382683432 * T1D; // FMA
  const T26 = T24 + T25;
  const T3C = T25 - T24;
  const T4D = T4B - T4C;
  const T4G = T4E - T4F;
  const T4H = KP707106781 * (T4D + T4G);
  const T4M = KP707106781 * (T4G - T4D);

  const T2k = x[3];
  const T2l = x[35];
  const T2m = T2k - T2l;
  const T4S = T2k + T2l;
  const T2u = x[11];
  const T2v = x[43];
  const T2w = T2u - T2v;
  const T4W = T2u + T2v;

  const T2n = x[19];
  const T2o = x[51];
  const T2p = T2n - T2o;
  const T4T = T2n + T2o;
  const T2r = x[59];
  const T2s = x[27];
  const T2t = T2r - T2s;
  const T4V = T2r + T2s;

  const T2q = KP923879532 * T2m - KP382683432 * T2p; // FNMS
  const T2x = KP923879532 * T2t + KP382683432 * T2w; // FMA
  const T2y = T2q + T2x;
  const T3J = T2x - T2q;
  const T5Q = T4S + T4T;
  const T5R = T4V + T4W;
  const T5S = T5Q + T5R;
  const T6a = T5R - T5Q;

  const T2A = KP382683432 * T2t - KP923879532 * T2w; // FNMS
  const T2B = KP382683432 * T2m + KP923879532 * T2p; // FMA
  const T2C = T2A - T2B;
  const T3G = T2B + T2A;
  const T4U = T4S - T4T;
  const T4X = T4V - T4W;
  const T4Y = KP707106781 * (T4U + T4X);
  const T53 = KP707106781 * (T4X - T4U);

  const O = new Float64Array(64);

  {
    const Tv = Tf + Tu;
    const T10 = TK + TZ;
    const T5X = Tv + T10;
    const T5Y = T5I + T5L;
    const T5Z = T5P + T5S;
    const T60 = T5Y + T5Z;
    O[16] = Tv - T10;
    O[48] = T5Z - T5Y;
    O[32] = T5X - T60;
    O[0] = T5X + T60;
  }
  {
    const T5F = Tf - Tu;
    const T5V = TZ - TK;
    const T5M = T5I - T5L;
    const T5T = T5P - T5S;
    const T5U = KP707106781 * (T5M + T5T);
    const T5W = KP707106781 * (T5T - T5M);
    O[24] = T5F - T5U;
    O[40] = T5W - T5V;
    O[8] = T5F + T5U;
    O[56] = T5V + T5W;
  }
  {
    const T61 = T7 - Te;
    const T64 = KP707106781 * (T62 + T63);
    const T65 = T61 + T64;
    const T6l = T61 - T64;
    const T6i = KP923879532 * T67 - KP382683432 * T66; // FNMS
    const T6j = KP382683432 * T69 + KP923879532 * T6a; // FMA
    const T6k = T6i + T6j;
    const T6m = T6j - T6i;

    const T68 = KP923879532 * T66 + KP382683432 * T67; // FMA
    const T6b = KP923879532 * T69 - KP382683432 * T6a; // FNMS
    const T6c = T68 + T6b;
    const T6g = T6b - T68;
    const T6d = KP707106781 * (T63 - T62);
    const T6e = Tt - Tm;
    const T6f = T6d - T6e;
    const T6h = T6e + T6d;

    O[28] = T65 - T6c;
    O[36] = T6k - T6h;
    O[4] = T65 + T6c;
    O[60] = T6h + T6k;
    O[52] = T6f + T6g;
    O[12] = T6l + T6m;
    O[44] = T6g - T6f;
    O[20] = T6l - T6m;
  }
  {
    const T5l = T4l - T4o;
    const T5m = T58 - T57;
    const T5n = T5l + T5m;
    const T5D = T5l - T5m;
    const T5v = T4v - T4s;
    const T5w = T5b - T5a;
    const T5x = T5v - T5w;
    const T5z = T5w + T5v;

    const T5o = T4A - T4H;
    const T5p = T4M - T4L;
    const T5q = KP831469612 * T5o + KP555570233 * T5p; // FMA
    const T5A = KP831469612 * T5p - KP555570233 * T5o; // FNMS
    const T5r = T4R - T4Y;
    const T5s = T53 - T52;
    const T5t = KP831469612 * T5r - KP555570233 * T5s; // FNMS
    const T5B = KP555570233 * T5r + KP831469612 * T5s; // FMA

    const T5u = T5q + T5t;
    O[26] = T5n - T5u;
    O[6] = T5n + T5u;
    const T5C = T5A + T5B;
    O[58] = T5z + T5C;
    O[38] = T5C - T5z;
    const T5y = T5t - T5q;
    O[54] = T5x + T5y;
    O[42] = T5y - T5x;
    const T5E = T5B - T5A;
    O[22] = T5D - T5E;
    O[10] = T5D + T5E;
  }
  {
    const T4p = T4l + T4o;
    const T4w = T4s + T4v;
    const T4x = T4p + T4w;
    const T5j = T4p - T4w;
    const T59 = T57 + T58;
    const T5c = T5a + T5b;
    const T5d = T59 - T5c;
    const T5f = T5c + T59;

    const T4I = T4A + T4H;
    const T4N = T4L + T4M;
    const T4O = KP980785280 * T4I + KP195090322 * T4N; // FMA
    const T5g = KP980785280 * T4N - KP195090322 * T4I; // FNMS
    const T4Z = T4R + T4Y;
    const T54 = T52 + T53;
    const T55 = KP980785280 * T4Z - KP195090322 * T54; // FNMS
    const T5h = KP195090322 * T4Z + KP980785280 * T54; // FMA

    const T56 = T4O + T55;
    O[30] = T4x - T56;
    O[2] = T4x + T56;
    const T5i = T5g + T5h;
    O[62] = T5f + T5i;
    O[34] = T5i - T5f;
    const T5e = T55 - T4O;
    O[50] = T5d + T5e;
    O[46] = T5e - T5d;
    const T5k = T5h - T5g;
    O[18] = T5j - T5k;
    O[14] = T5j + T5k;
  }
  {
    const T3p = T3n + T3o;
    const T41 = T3n - T3o;
    const T4c = T3R - T3Q;
    const T3S = T3Q + T3R;
    const T3s = KP831469612 * T3q + KP555570233 * T3r; // FMA
    const T3v = KP831469612 * T3t - KP555570233 * T3u; // FNMS
    const T3w = T3s + T3v;
    const T4b = T3v - T3s;

    const T47 = T3F - T3G;
    const T48 = T3J - T3I;
    const T49 = KP881921264 * T47 - KP471396736 * T48; // FNMS
    const T4h = KP471396736 * T47 + KP881921264 * T48; // FMA
    const T3N = KP831469612 * T3r - KP555570233 * T3q; // FNMS
    const T3O = KP555570233 * T3t + KP831469612 * T3u; // FMA
    const T3P = T3N + T3O;
    const T42 = T3O - T3N;

    const T3A = T3y + T3z;
    const T3D = T3B + T3C;
    const T3E = KP956940335 * T3A + KP290284677 * T3D; // FMA
    const T3W = KP956940335 * T3D - KP290284677 * T3A; // FNMS

    const T44 = T3y - T3z;
    const T45 = T3C - T3B;
    const T46 = KP881921264 * T44 + KP471396736 * T45; // FMA
    const T4g = KP881921264 * T45 - KP471396736 * T44; // FNMS
    const T3H = T3F + T3G;
    const T3K = T3I + T3J;
    const T3L = KP956940335 * T3H - KP290284677 * T3K; // FNMS
    const T3X = KP290284677 * T3H + KP956940335 * T3K; // FMA

    const T3x = T3p + T3w;
    const T3M = T3E + T3L;
    O[29] = T3x - T3M;
    O[3] = T3x + T3M;
    const T3V = T3S + T3P;
    const T3Y = T3W + T3X;
    O[61] = T3V + T3Y;
    O[35] = T3Y - T3V;

    const T3T = T3P - T3S;
    const T3U = T3L - T3E;
    O[51] = T3T + T3U;
    O[45] = T3U - T3T;
    const T3Z = T3p - T3w;
    const T40 = T3X - T3W;
    O[19] = T3Z - T40;
    O[13] = T3Z + T40;

    const T43 = T41 + T42;
    const T4a = T46 + T49;
    O[27] = T43 - T4a;
    O[5] = T43 + T4a;
    const T4f = T4c + T4b;
    const T4i = T4g + T4h;
    O[59] = T4f + T4i;
    O[37] = T4i - T4f;

    const T4d = T4b - T4c;
    const T4e = T49 - T46;
    O[53] = T4d + T4e;
    O[43] = T4e - T4d;
    const T4j = T41 - T42;
    const T4k = T4h - T4g;
    O[21] = T4j - T4k;
    O[11] = T4j + T4k;
  }
  {
    const T1d = T15 - T1c;
    const T33 = T15 + T1c;
    const T3e = T2T + T2Q;
    const T2U = T2Q - T2T;
    const T1m = KP195090322 * T1i + KP980785280 * T1l; // FMA
    const T1v = KP980785280 * T1r - KP195090322 * T1u; // FNMS
    const T1w = T1m - T1v;
    const T3d = T1v + T1m;

    const T39 = T2j + T2y;
    const T3a = T2H + T2C;
    const T3b = KP995184726 * T39 - KP098017140 * T3a; // FNMS
    const T3j = KP098017140 * T39 + KP995184726 * T3a; // FMA
    const T2L = KP980785280 * T1i - KP195090322 * T1l; // FNMS
    const T2M = KP195090322 * T1r + KP980785280 * T1u; // FMA
    const T2N = T2L - T2M;
    const T34 = T2M + T2L;

    const T1Y = T1M - T1X;
    const T27 = T23 - T26;
    const T28 = KP634393284 * T1Y + KP773010453 * T27; // FMA
    const T2Y = KP773010453 * T1Y - KP634393284 * T27; // FNMS

    const T36 = T1X + T1M;
    const T37 = T23 + T26;
    const T38 = KP098017140 * T36 + KP995184726 * T37; // FMA
    const T3i = KP995184726 * T36 - KP098017140 * T37; // FNMS
    const T2z = T2j - T2y;
    const T2I = T2C - T2H;
    const T2J = KP773010453 * T2z - KP634393284 * T2I; // FNMS
    const T2Z = KP634393284 * T2z + KP773010453 * T2I; // FMA

    const T1x = T1d + T1w;
    const T2K = T28 + T2J;
    O[25] = T1x - T2K;
    O[7] = T1x + T2K;
    const T2X = T2U + T2N;
    const T30 = T2Y + T2Z;
    O[57] = T2X + T30;
    O[39] = T30 - T2X;

    const T2V = T2N - T2U;
    const T2W = T2J - T28;
    O[55] = T2V + T2W;
    O[41] = T2W - T2V;
    const T31 = T1d - T1w;
    const T32 = T2Z - T2Y;
    O[23] = T31 - T32;
    O[9] = T31 + T32;

    const T35 = T33 + T34;
    const T3c = T38 + T3b;
    O[31] = T35 - T3c;
    O[1] = T35 + T3c;
    const T3h = T3e + T3d;
    const T3k = T3i + T3j;
    O[63] = T3h + T3k;
    O[33] = T3k - T3h;

    const T3f = T3d - T3e;
    const T3g = T3b - T38;
    O[49] = T3f + T3g;
    O[47] = T3g - T3f;
    const T3l = T33 - T34;
    const T3m = T3j - T3i;
    O[17] = T3l - T3m;
    O[15] = T3l + T3m;
  }

  return O;
}

module.exports = { r2cf_64 };
