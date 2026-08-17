'use strict';

// =============================================================================
// t2_64.js -- faithful JS port of dft/scalar/codelets/t2_64.c (non-FMA
// branch), FFTW3's "twiddle-log3 / precompute-twiddles" radix-64 twiddle
// codelet. twinstr only trig-generates W^1, W^3, W^9, W^27, W^63 (five raw
// pairs); every other needed multiple is DERIVED via complex products of
// those five, exactly mirroring the C source's own derivation chain
// bit-for-bit (same trick as t2_5/t2_8/t2_10/t2_16/t2_20/t2_25/t2_32.js).
// Same calling convention as every other twiddle codelet here (br/bi/Wc/Ws
// with Composite1D.js's full r-1 = 63 pair table already built) -- only
// indices 1, 3, 9, 27, 63 are actually read; every other Wc[k]/Ws[k] is
// deliberately ignored in favor of recomputing the equivalent product
// chain, to match FFTW's exact rounding (real FFTW never has a full table
// either -- it only has the 5 raw generator pairs from twinstr).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP471396736 = 0.471396736825997648556387625905254377657460319;
const KP881921264 = 0.881921264348355029712756863660388349508442621;
const KP290284677 = 0.290284677254462367636192375817395274691476278;
const KP956940335 = 0.956940335732208864935797886980269969482849206;
const KP634393284 = 0.634393284163645498215171613225493370675687095;
const KP773010453 = 0.773010453362736960810906609758469800971041293;
const KP098017140 = 0.098017140329560601994195563888641845861136673;
const KP995184726 = 0.995184726672196886244836953109479921575474869;
const KP555570233 = 0.555570233019602224742830813948532874374937191;
const KP831469612 = 0.831469612302545237078788377617905756738560812;
const KP980785280 = 0.980785280403230449126182236134239036973933731;
const KP195090322 = 0.195090322016128267848284868477022240927691618;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function t2_64(br, bi, Wc, Ws) {
  const outR = new Float64Array(64), outI = new Float64Array(64);

  // -- precompute all needed twiddle products from the 5 raw generators --
  const T2 = Wc[1], T5 = Ws[1], T3 = Wc[3], T6 = Ws[3];
  const Te = Ws[9], T9 = Wc[9];
  const T4 = T2 * T3, T1d = T5 * T9, T19 = T5 * Te, Tb = T2 * T6;
  const T1c = T2 * Te, T7 = T5 * T6, Tc = T5 * T3, T18 = T2 * T9;
  const TR = T3 * Te, TO = T6 * Te, TS = T6 * T9, TN = T3 * T9;
  const TP = TN - TO, T3e = TR - TS, T1e = T1c - T1d, T39 = T1c + T1d;
  const T3c = TN + TO, TT = TR + TS, T1a = T18 + T19, T37 = T18 - T19;
  const T8 = T4 - T7, Ta = T8 * T9, Tj = T8 * Te;
  const Tw = T4 + T7, Tx = Tw * T9, TC = Tw * Te;
  const Td = Tb + Tc, Tf = Td * Te, Tk = Td * T9;
  const Ty = Tb - Tc, Tz = Ty * Te, TD = Ty * T9;
  const Tm = Ws[27];
  const T1B = T6 * Tm, T1E = T3 * Tm, T2o = T2 * Tm, T2l = T5 * Tm, T1T = T9 * Tm, T1Q = Te * Tm;
  const Th = Wc[27];
  const T1A = T3 * Th, T1F = T6 * Th, T2p = T5 * Th, T2k = T2 * Th, T1U = Te * Th, T1P = T9 * Th;

  const T1C = T1A + T1B, T3K = T1E + T1F, T1V = T1T + T1U, T3x = T2o - T2p;
  const T3I = T1A - T1B, T1G = T1E - T1F, T1R = T1P - T1Q;
  const T3v = T2k + T2l, T2m = T2k - T2l, T2q = T2o + T2p;
  const T5W = T8 * Th, T5X = Td * Tm, T5Y = T5W - T5X, T6u = T5W + T5X;
  const T51 = Tw * Th, T52 = Ty * Tm, T53 = T51 + T52, T5B = T51 - T52;
  const T60 = T8 * Tm, T61 = Td * Th, T62 = T60 + T61, T6w = T60 - T61;
  const T55 = Tw * Tm, T56 = Ty * Th, T57 = T55 - T56, T5D = T55 + T56;

  const T2V = T1P + T1Q, T2X = T1T - T1U;
  const Tg = Ta + Tf, Ti = Tg * Th, Tq = Tg * Tm;
  const TE = TC + TD, TF = TE * Tm, TJ = TE * Th;
  const T3W = T37 * Tm, T3X = T39 * Th, T3Y = T3W - T3X;
  const T3T = T37 * Th, T3U = T39 * Tm, T3V = T3T + T3U;
  const T3h = T3c * Tm, T3i = T3e * Th, T3j = T3h - T3i;
  const Tl = Tj - Tk, Tn = Tl * Tm, Tr = Tl * Th;
  const TA = Tx - Tz, TB = TA * Th, TI = TA * Tm;
  const T3d = T3c * Th, T3f = T3e * Tm, T3g = T3d + T3f;
  const T1j = Tj + Tk, T1k = T1j * Tm, T1o = T1j * Th;
  const T1t = Tx + Tz, T1Z = T1t * Th, T23 = T1t * Tm;
  const TQ = TP * Th, TU = TT * Tm, TV = TQ + TU;
  const T2A = T1a * Tm, T2B = T1e * Th, T2C = T2A - T2B;
  const T2x = T1a * Th, T2y = T1e * Tm, T2z = T2x + T2y;
  const T1u = TC - TD, T20 = T1u * Tm, T24 = T1u * Th;
  const TX = TP * Tm, TY = TT * Th, TZ = TX - TY;
  const T1h = Ta - Tf, T1i = T1h * Th, T1n = T1h * Tm;
  const To = Ti - Tn, T1p = T1n + T1o;
  const T6j = TQ - TU, T6H = T2A + T2B;
  const Ts = Tq + Tr, T1l = T1i - T1k;
  const T6l = TX + TY, T6F = T2x - T2y;
  const T2P = T1Z - T20, T4b = TI + TJ;
  const T4x = T3d - T3f, T5i = T3W + T3X;
  const T2R = T23 + T24, T49 = TB - TF;
  const T4z = T3h + T3i, T5g = T3T - T3U;
  const TG = TB + TF, T4k = Ti + Tn;
  const T4m = Tq - Tr, TK = TI - TJ;
  const T21 = T1Z + T20, T3O = T1i + T1k;
  const T3Q = T1n - T1o, T25 = T23 - T24;
  const TW = Wc[63], T10 = Ws[63];

  const T11 = TV * TW + TZ * T10, T79 = T21 * T10 - T25 * TW;
  const T6X = T8 * T10 - Td * TW, T5M = T2V * T10 - T2X * TW;
  const T6b = TG * T10 - TK * TW, T1v = T1t * TW + T1u * T10;
  const T30 = T1h * TW + T1j * T10, T69 = TG * TW + TK * T10;
  const T77 = T21 * TW + T25 * T10, T13 = TV * T10 - TZ * TW;
  const T2F = T2z * T10 - T2C * TW, T2D = T2z * TW + T2C * T10;
  const T6p = T1a * TW + T1e * T10, T6O = TP * TW + TT * T10;
  const T1x = T1t * T10 - T1u * TW, T2a = TA * T10 - TE * TW;
  const T2f = T3 * TW + T6 * T10, T6V = T8 * TW + Td * T10;
  const T28 = TA * TW + TE * T10, T6r = T1a * T10 - T1e * TW;
  const T2h = T3 * T10 - T6 * TW, T6Q = TP * T10 - TT * TW;
  const T32 = T1h * T10 - T1j * TW, T5K = T2V * TW + T2X * T10;
  const T5w = Tw * TW + Ty * T10, T4G = T3O * TW + T3Q * T10;
  const T4Q = T4k * TW + T4m * T10, T3m = T3g * T10 - T3j * TW;
  const T4h = T9 * T10 - Te * TW, T4I = T3O * T10 - T3Q * TW;
  const T5y = Tw * T10 - Ty * TW, T3k = T3g * TW + T3j * T10;
  const T4f = T9 * TW + Te * T10, T41 = T3V * T10 - T3Y * TW;
  const T4S = T4k * T10 - T4m * TW, T4Y = T3c * T10 - T3e * TW;
  const T3q = Tg * TW + Tl * T10, T3D = T2 * TW + T5 * T10;
  const T3F = T2 * T10 - T5 * TW, T5r = T37 * T10 - T39 * TW;
  const T3s = Tg * T10 - Tl * TW, T4W = T3c * TW + T3e * T10;
  const T3Z = T3V * TW + T3Y * T10, T5p = T37 * TW + T39 * T10;

  // -- butterfly group 1: k = 0, 32, 16, 48 --
  const T1 = br[0], Ti1 = bi[0];
  const Tp = br[32], Tt = bi[32];
  const Tu = To * Tp + Ts * Tt, Ti0 = To * Tt - Ts * Tp;
  const TH = br[16], TL = bi[16];
  const TM = TG * TH + TK * TL, T7i = TG * TL - TK * TH;
  const T12 = br[48], T14 = bi[48];
  const T15 = T11 * T12 + T13 * T14, T7j = T11 * T14 - T13 * T12;

  const Tv = T1 + Tu, T16 = TM + T15;
  const T17 = Tv + T16, TdV = Tv - T16;
  const Tj1 = Ti1 - Ti0, Tj2 = TM - T15;
  const Tj3 = Tj1 - Tj2, Tjx = Tj2 + Tj1;

  const T7h = T1 - Tu, T7k = T7i - T7j;
  const T7l = T7h - T7k, TbJ = T7h + T7k;
  const ThZ = T7i + T7j, Ti2 = Ti0 + Ti1;
  const Ti3 = ThZ + Ti2, Tix = Ti2 - ThZ;

  // -- butterfly group 2: k = 8, 40, 56, 24 --
  const T1b = br[8], T1f = bi[8];
  const T1g = T1a * T1b + T1e * T1f, T7m = T1a * T1f - T1e * T1b;
  const T1m = br[40], T1q = bi[40];
  const T1r = T1l * T1m + T1p * T1q, T7n = T1l * T1q - T1p * T1m;

  const T7o = T7m - T7n, T7p = T1g - T1r;

  const T1w = br[56], T1y = bi[56];
  const T1z = T1v * T1w + T1x * T1y, T7s = T1v * T1y - T1x * T1w;
  const T1D = br[24], T1H = bi[24];
  const T1I = T1C * T1D + T1G * T1H, T7t = T1C * T1H - T1G * T1D;

  const T7r = T1z - T1I, T7u = T7s - T7t;

  const T1s = T1g + T1r, T1J = T1z + T1I;
  const T1K = T1s + T1J, Tiw = T1J - T1s;
  const TdW = T7m + T7n, TdX = T7s + T7t;
  const TdY = TdW - TdX, ThY = TdW + TdX;

  const T7q = T7o - T7p, T7v = T7r + T7u;
  const T7w = KP707106781 * (T7q - T7v), Tj0 = KP707106781 * (T7q + T7v);
  const TbK = T7p + T7o, TbL = T7r - T7u;
  const TbM = KP707106781 * (TbK + TbL), Tjw = KP707106781 * (TbL - TbK);

  // -- butterfly group 3: k = 4, 36, 20, 52 --
  const T1M = br[4], T1N = bi[4];
  const T1O = T8 * T1M + Td * T1N, T7y = T8 * T1N - Td * T1M;
  const T1S = br[36], T1W = bi[36];
  const T1X = T1R * T1S + T1V * T1W, T7z = T1R * T1W - T1V * T1S;

  const T1Y = T1O + T1X, Te0 = T7y + T7z;
  const T7A = T7y - T7z, T7D = T1O - T1X;

  const T22 = br[20], T26 = bi[20];
  const T27 = T21 * T22 + T25 * T26, T7E = T21 * T26 - T25 * T22;
  const T29 = br[52], T2b = bi[52];
  const T2c = T28 * T29 + T2a * T2b, T7F = T28 * T2b - T2a * T29;

  const T2d = T27 + T2c, Te1 = T7E + T7F;
  const T7B = T27 - T2c, T7G = T7E - T7F;

  const T2e = T1Y + T2d, TgA = Te0 + Te1;
  const T7C = T7A + T7B, T7H = T7D - T7G;
  const T7I = KP382683432 * T7C - KP923879532 * T7H;
  const TaY = KP923879532 * T7C + KP382683432 * T7H;

  const TbO = T7A - T7B, TbP = T7D + T7G;
  const TbQ = KP923879532 * TbO - KP382683432 * TbP;
  const Tda = KP382683432 * TbO + KP923879532 * TbP;
  const Te2 = Te0 - Te1, Te3 = T1Y - T2d;
  const Te4 = Te2 - Te3, TfO = Te3 + Te2;

  // -- butterfly group 4: k = 60, 28, 12, 44 --
  const T2g = br[60], T2i = bi[60];
  const T2j = T2f * T2g + T2h * T2i, T7J = T2f * T2i - T2h * T2g;
  const T2n = br[28], T2r = bi[28];
  const T2s = T2m * T2n + T2q * T2r, T7K = T2m * T2r - T2q * T2n;

  const T2t = T2j + T2s, Te6 = T7J + T7K;
  const T7L = T7J - T7K, T7O = T2j - T2s;

  const T2u = br[12], T2v = bi[12];
  const T2w = TP * T2u + TT * T2v, T7P = TP * T2v - TT * T2u;
  const T2E = br[44], T2G = bi[44];
  const T2H = T2D * T2E + T2F * T2G, T7Q = T2D * T2G - T2F * T2E;

  const T2I = T2w + T2H, Te7 = T7P + T7Q;
  const T7M = T2w - T2H, T7R = T7P - T7Q;

  const T2J = T2t + T2I, TgB = Te6 + Te7;
  const T7N = T7L + T7M, T7S = T7O - T7R;
  const T7T = KP382683432 * T7N + KP923879532 * T7S;
  const TaZ = KP382683432 * T7S - KP923879532 * T7N;

  const TbR = T7L - T7M, TbS = T7O + T7R;
  const TbT = KP923879532 * TbR + KP382683432 * TbS;
  const Tdb = KP923879532 * TbS - KP382683432 * TbR;
  const Te5 = T2t - T2I, Te8 = Te6 - Te7;
  const Te9 = Te5 + Te8, TfP = Te5 - Te8;

  // -- butterfly group 5: k = 2, 34, 18, 50, 10, 42, 58, 26 --
  const T2M = br[2], T2N = bi[2];
  const T2O = Tw * T2M + Ty * T2N, T7W = Tw * T2N - Ty * T2M;
  const T2Q = br[34], T2S = bi[34];
  const T2T = T2P * T2Q + T2R * T2S, T7X = T2P * T2S - T2R * T2Q;

  const T2U = T2O + T2T, Tec = T7W + T7X;

  const T2W = br[18], T2Y = bi[18];
  const T2Z = T2V * T2W + T2X * T2Y, T8e = T2V * T2Y - T2X * T2W;
  const T31 = br[50], T33 = bi[50];
  const T34 = T30 * T31 + T32 * T33, T8f = T30 * T33 - T32 * T31;

  const T35 = T2Z + T34, Ted = T8e + T8f;

  const T38 = br[10], T3a = bi[10];
  const T3b = T37 * T38 + T39 * T3a, T87 = T37 * T3a - T39 * T38;
  const T3l = br[42], T3n = bi[42];
  const T3o = T3k * T3l + T3m * T3n, T88 = T3k * T3n - T3m * T3l;

  const T3p = T3b + T3o, Tei = T87 + T88;
  const T86 = T3b - T3o, T89 = T87 - T88;

  const T3r = br[58], T3t = bi[58];
  const T3u = T3q * T3r + T3s * T3t, T82 = T3q * T3t - T3s * T3r;
  const T3w = br[26], T3y = bi[26];
  const T3z = T3v * T3w + T3x * T3y, T83 = T3v * T3y - T3x * T3w;

  const T3A = T3u + T3z, Tej = T82 + T83;
  const T81 = T3u - T3z, T84 = T82 - T83;

  const T36 = T2U + T35, T3B = T3p + T3A;
  const TgH = T36 - T3B, TgE = Tec + Ted;
  const TgF = Tei + Tej, TgG = TgE - TgF;

  const T7Y = T7W - T7X, T7Z = T2Z - T34;
  const T80 = T7Y + T7Z, TbW = T7Y - T7Z;
  const Teh = T2U - T35, Tek = Tei - Tej;
  const Tel = Teh - Tek, TfT = Teh + Tek;

  const T85 = T81 - T84, T8a = T86 + T89;
  const T8b = KP707106781 * (T85 - T8a), Tc0 = KP707106781 * (T8a + T85);
  const T8i = T89 - T86, T8j = T81 + T84;
  const T8k = KP707106781 * (T8i - T8j), TbX = KP707106781 * (T8i + T8j);

  const Tee = Tec - Ted, Tef = T3A - T3p;
  const Teg = Tee - Tef, TfS = Tee + Tef;
  const T8d = T2O - T2T, T8g = T8e - T8f;
  const T8h = T8d - T8g, TbZ = T8d + T8g;

  // -- butterfly group 6: k = 62, 30, 14, 46, 6, 38, 54, 22 --
  const T3E = br[62], T3G = bi[62];
  const T3H = T3D * T3E + T3F * T3G, T8n = T3D * T3G - T3F * T3E;
  const T3J = br[30], T3L = bi[30];
  const T3M = T3I * T3J + T3K * T3L, T8o = T3I * T3L - T3K * T3J;

  const T3N = T3H + T3M, Ten = T8n + T8o;

  const T3P = br[14], T3R = bi[14];
  const T3S = T3O * T3P + T3Q * T3R, T8F = T3O * T3R - T3Q * T3P;
  const T40 = br[46], T42 = bi[46];
  const T43 = T3Z * T40 + T41 * T42, T8G = T3Z * T42 - T41 * T40;

  const T44 = T3S + T43, Teo = T8F + T8G;

  const T46 = br[6], T47 = bi[6];
  const T48 = T3c * T46 + T3e * T47, T8y = T3c * T47 - T3e * T46;
  const T4a = br[38], T4c = bi[38];
  const T4d = T49 * T4a + T4b * T4c, T8z = T49 * T4c - T4b * T4a;

  const T4e = T48 + T4d, Tet = T8y + T8z;
  const T8x = T48 - T4d, T8A = T8y - T8z;

  const T4g = br[54], T4i = bi[54];
  const T4j = T4f * T4g + T4h * T4i, T8t = T4f * T4i - T4h * T4g;
  const T4l = br[22], T4n = bi[22];
  const T4o = T4k * T4l + T4m * T4n, T8u = T4k * T4n - T4m * T4l;

  const T4p = T4j + T4o, Teu = T8t + T8u;
  const T8s = T4j - T4o, T8v = T8t - T8u;

  const T45 = T3N + T44, T4q = T4e + T4p;
  const TgJ = T45 - T4q, TgK = Ten + Teo;
  const TgL = Tet + Teu, TgM = TgK - TgL;

  const T8p = T8n - T8o, T8q = T3S - T43;
  const T8r = T8p + T8q, Tc6 = T8p - T8q;
  const Tes = T3N - T44, Tev = Tet - Teu;
  const Tew = Tes - Tev, TfW = Tes + Tev;

  const T8w = T8s - T8v, T8B = T8x + T8A;
  const T8C = KP707106781 * (T8w - T8B), Tc4 = KP707106781 * (T8B + T8w);
  const T8J = T8A - T8x, T8K = T8s + T8v;
  const T8L = KP707106781 * (T8J - T8K), Tc7 = KP707106781 * (T8J + T8K);

  const Tep = Ten - Teo, Teq = T4p - T4e;
  const Ter = Tep - Teq, TfV = Tep + Teq;
  const T8E = T3H - T3M, T8H = T8F - T8G;
  const T8I = T8E - T8H, Tc3 = T8E + T8H;

  // -- butterfly group 7: k = 63, 31, 15, 47, 7, 39, 55, 23 --
  const T5T = br[63], T5U = bi[63];
  const T5V = TW * T5T + T10 * T5U, Tao = TW * T5U - T10 * T5T;
  const T5Z = br[31], T63 = bi[31];
  const T64 = T5Y * T5Z + T62 * T63, Tap = T5Y * T63 - T62 * T5Z;

  const T65 = T5V + T64, Tfi = Tao + Tap;

  const T66 = br[15], T67 = bi[15];
  const T68 = TV * T66 + TZ * T67, T9K = TV * T67 - TZ * T66;
  const T6a = br[47], T6c = bi[47];
  const T6d = T69 * T6a + T6b * T6c, T9L = T69 * T6c - T6b * T6a;

  const T6e = T68 + T6d, Tfj = T9K + T9L;

  const T6g = br[7], T6h = bi[7];
  const T6i = T1t * T6g + T1u * T6h, T9O = T1t * T6h - T1u * T6g;
  const T6k = br[39], T6m = bi[39];
  const T6n = T6j * T6k + T6l * T6m, T9P = T6j * T6m - T6l * T6k;

  const T6o = T6i + T6n, Tf2 = T9O + T9P;
  const T9Q = T9O - T9P, T9R = T6i - T6n;

  const T6q = br[55], T6s = bi[55];
  const T6t = T6p * T6q + T6r * T6s, T9U = T6p * T6s - T6r * T6q;
  const T6v = br[23], T6x = bi[23];
  const T6y = T6u * T6v + T6w * T6x, T9V = T6u * T6x - T6w * T6v;

  const T6z = T6t + T6y, Tf3 = T9U + T9V;
  const T9T = T6t - T6y, T9W = T9U - T9V;

  const T6f = T65 + T6e, T6A = T6o + T6z;
  const T6B = T6f + T6A, Th1 = T6f - T6A;
  const Tfk = Tfi - Tfj, Tfl = T6z - T6o;
  const Tfm = Tfk - Tfl, Tga = Tfk + Tfl;

  const Th6 = Tfi + Tfj, Th7 = Tf2 + Tf3;
  const Th8 = Th6 - Th7, ThI = Th6 + Th7;
  const T9J = T5V - T64, T9M = T9K - T9L;
  const T9N = T9J - T9M, Tcv = T9J + T9M;

  const T9S = T9Q - T9R, T9X = T9T + T9W;
  const T9Y = KP707106781 * (T9S - T9X), TcH = KP707106781 * (T9S + T9X);
  const Tat = T9T - T9W, Tau = T9R + T9Q;
  const Tav = KP707106781 * (Tat - Tau), Tcw = KP707106781 * (Tau + Tat);

  const Tf1 = T65 - T6e, Tf4 = Tf2 - Tf3;
  const Tf5 = Tf1 - Tf4, Tg7 = Tf1 + Tf4;
  const Taq = Tao - Tap, Tar = T68 - T6d;
  const Tas = Taq + Tar, TcG = Taq - Tar;

  // -- butterfly group 8: k = 1, 33, 17, 49, 9, 41, 57, 25 --
  const T4u = br[1], T4v = bi[1];
  const T4w = T2 * T4u + T5 * T4v, T8Q = T2 * T4v - T5 * T4u;
  const T4y = br[33], T4A = bi[33];
  const T4B = T4x * T4y + T4z * T4A, T8R = T4x * T4A - T4z * T4y;

  const T4C = T4w + T4B, TeA = T8Q + T8R;

  const T4D = br[17], T4E = bi[17];
  const T4F = T3V * T4D + T3Y * T4E, T9w = T3V * T4E - T3Y * T4D;
  const T4H = br[49], T4J = bi[49];
  const T4K = T4G * T4H + T4I * T4J, T9x = T4G * T4J - T4I * T4H;

  const T4L = T4F + T4K, TeB = T9w + T9x;

  const T4N = br[9], T4O = bi[9];
  const T4P = T9 * T4N + Te * T4O, T91 = T9 * T4O - Te * T4N;
  const T4R = br[41], T4T = bi[41];
  const T4U = T4Q * T4R + T4S * T4T, T92 = T4Q * T4T - T4S * T4R;

  const T4V = T4P + T4U, TeS = T91 + T92;
  const T90 = T4P - T4U, T93 = T91 - T92;

  const T4X = br[57], T4Z = bi[57];
  const T50 = T4W * T4X + T4Y * T4Z, T8W = T4W * T4Z - T4Y * T4X;
  const T54 = br[25], T58 = bi[25];
  const T59 = T53 * T54 + T57 * T58, T8X = T53 * T58 - T57 * T54;

  const T5a = T50 + T59, TeT = T8W + T8X;
  const T8V = T50 - T59, T8Y = T8W - T8X;

  const T4M = T4C + T4L, T5b = T4V + T5a;
  const T5c = T4M + T5b, TgV = T4M - T5b;
  const TeR = T4C - T4L, TeU = TeS - TeT;
  const TeV = TeR - TeU, Tg0 = TeR + TeU;

  const TgQ = TeA + TeB, TgR = TeS + TeT;
  const TgS = TgQ - TgR, ThD = TgQ + TgR;
  const T8S = T8Q - T8R, T8T = T4F - T4K;
  const T8U = T8S + T8T, Tcc = T8S - T8T;

  const T8Z = T8V - T8Y, T94 = T90 + T93;
  const T95 = KP707106781 * (T8Z - T94), Tco = KP707106781 * (T94 + T8Z);
  const T9A = T93 - T90, T9B = T8V + T8Y;
  const T9C = KP707106781 * (T9A - T9B), Tcd = KP707106781 * (T9A + T9B);

  const TeC = TeA - TeB, TeD = T5a - T4V;
  const TeE = TeC - TeD, Tg3 = TeC + TeD;
  const T9v = T4w - T4B, T9y = T9w - T9x;
  const T9z = T9v - T9y, Tcn = T9v + T9y;

  // -- butterfly group 9: k = 5, 37, 13, 45, 21, 53, 61, 29 --
  const T5d = br[5], T5e = bi[5];
  const T5f = Tg * T5d + Tl * T5e, T9i = Tg * T5e - Tl * T5d;
  const T5h = br[37], T5j = bi[37];
  const T5k = T5g * T5h + T5i * T5j, T9j = T5g * T5j - T5i * T5h;

  const T5l = T5f + T5k, TeL = T9i + T9j;
  const T9k = T9i - T9j, T9n = T5f - T5k;

  const T5H = br[13], T5I = bi[13];
  const T5J = T1h * T5H + T1j * T5I, T98 = T1h * T5I - T1j * T5H;
  const T5L = br[45], T5N = bi[45];
  const T5O = T5K * T5L + T5M * T5N, T99 = T5K * T5N - T5M * T5L;

  const T5P = T5J + T5O, TeH = T98 + T99;
  const T9a = T98 - T99, T9f = T5J - T5O;

  const T5m = br[21], T5n = bi[21];
  const T5o = T3g * T5m + T3j * T5n, T9o = T3g * T5n - T3j * T5m;
  const T5q = br[53], T5s = bi[53];
  const T5t = T5p * T5q + T5r * T5s, T9p = T5p * T5s - T5r * T5q;

  const T5u = T5o + T5t, TeM = T9o + T9p;
  const T9l = T5o - T5t, T9q = T9o - T9p;

  const T5x = br[61], T5z = bi[61];
  const T5A = T5w * T5x + T5y * T5z, T9c = T5w * T5z - T5y * T5x;
  const T5C = br[29], T5E = bi[29];
  const T5F = T5B * T5C + T5D * T5E, T9d = T5B * T5E - T5D * T5C;

  const T5G = T5A + T5F, TeG = T9c + T9d;
  const T97 = T5A - T5F, T9e = T9c - T9d;

  const T5v = T5l + T5u, T5Q = T5G + T5P;
  const T5R = T5v + T5Q, TgT = T5Q - T5v;
  const TeK = T5l - T5u, TeN = TeL - TeM;
  const TeO = TeK + TeN, TeW = TeN - TeK;

  const TgW = TeL + TeM, TgX = TeG + TeH;
  const TgY = TgW - TgX, ThE = TgW + TgX;
  const T9b = T97 - T9a, T9g = T9e + T9f;
  const T9h = KP382683432 * T9b - KP923879532 * T9g;
  const T9F = KP382683432 * T9g + KP923879532 * T9b;

  const T9m = T9k + T9l, T9r = T9n - T9q;
  const T9s = KP923879532 * T9m + KP382683432 * T9r;
  const T9E = KP382683432 * T9m - KP923879532 * T9r;
  const Tci = T9k - T9l, Tcj = T9n + T9q;
  const Tck = KP382683432 * Tci + KP923879532 * Tcj;
  const Tcq = KP923879532 * Tci - KP382683432 * Tcj;

  const TeF = T5G - T5P, TeI = TeG - TeH;
  const TeJ = TeF - TeI, TeX = TeF + TeI;
  const Tcf = T97 + T9a, Tcg = T9e - T9f;
  const Tch = KP923879532 * Tcf - KP382683432 * Tcg;
  const Tcr = KP923879532 * Tcg + KP382683432 * Tcf;

  // -- butterfly group 10: k = 3, 35, 11, 43, 19, 51, 59, 27 --
  const T6C = br[3], T6D = bi[3];
  const T6E = T3 * T6C + T6 * T6D, Ta0 = T3 * T6D - T6 * T6C;
  const T6G = br[35], T6I = bi[35];
  const T6J = T6F * T6G + T6H * T6I, Ta1 = T6F * T6I - T6H * T6G;

  const T6K = T6E + T6J, Tf6 = Ta0 + Ta1;
  const Ta2 = Ta0 - Ta1, Ta5 = T6E - T6J;

  const T74 = br[11], T75 = bi[11];
  const T76 = TA * T74 + TE * T75, Tah = TA * T75 - TE * T74;
  const T78 = br[43], T7a = bi[43];
  const T7b = T77 * T78 + T79 * T7a, Tai = T77 * T7a - T79 * T78;

  const T7c = T76 + T7b, Tfd = Tah + Tai;
  const Tae = T76 - T7b, Taj = Tah - Tai;

  const T6L = br[19], T6M = bi[19];
  const T6N = T2z * T6L + T2C * T6M, Ta6 = T2z * T6M - T2C * T6L;
  const T6P = br[51], T6R = bi[51];
  const T6S = T6O * T6P + T6Q * T6R, Ta7 = T6O * T6R - T6Q * T6P;

  const T6T = T6N + T6S, Tf7 = Ta6 + Ta7;
  const Ta3 = T6N - T6S, Ta8 = Ta6 - Ta7;

  const T6W = br[59], T6Y = bi[59];
  const T6Z = T6V * T6W + T6X * T6Y, Tab = T6V * T6Y - T6X * T6W;
  const T70 = br[27], T71 = bi[27];
  const T72 = Th * T70 + Tm * T71, Tac = Th * T71 - Tm * T70;

  const T73 = T6Z + T72, Tfc = Tab + Tac;
  const Tad = Tab - Tac, Tag = T6Z - T72;

  const T6U = T6K + T6T, T7d = T73 + T7c;
  const T7e = T6U + T7d, Th9 = T7d - T6U;
  const Tfb = T73 - T7c, Tfe = Tfc - Tfd;
  const Tff = Tfb + Tfe, Tfn = Tfb - Tfe;

  const Th2 = Tf6 + Tf7, Th3 = Tfc + Tfd;
  const Th4 = Th2 - Th3, ThJ = Th2 + Th3;
  const Ta4 = Ta2 + Ta3, Ta9 = Ta5 - Ta8;
  const Taa = KP382683432 * Ta4 - KP923879532 * Ta9;
  const Tay = KP923879532 * Ta4 + KP382683432 * Ta9;

  const Taf = Tad + Tae, Tak = Tag - Taj;
  const Tal = KP382683432 * Taf + KP923879532 * Tak;
  const Tax = KP382683432 * Tak - KP923879532 * Taf;
  const TcB = Tad - Tae, TcC = Tag + Taj;
  const TcD = KP923879532 * TcB + KP382683432 * TcC;
  const TcJ = KP923879532 * TcC - KP382683432 * TcB;

  const Tf8 = Tf6 - Tf7, Tf9 = T6K - T6T;
  const Tfa = Tf8 - Tf9, Tfo = Tf9 + Tf8;
  const Tcy = Ta2 - Ta3, Tcz = Ta5 + Ta8;
  const TcA = KP923879532 * Tcy - KP382683432 * Tcz;
  const TcK = KP382683432 * Tcy + KP923879532 * Tcz;

  // -- combine group A: outputs 32, 0, 48, 16, 40, 8, 56, 24 --
  const T1L = T17 + T1K, T2K = T2e + T2J;
  const T2L = T1L + T2K, Thx = T1L - T2K;
  const ThS = ThD + ThE, ThT = ThI + ThJ;
  const ThU = ThS - ThT, ThV = ThS + ThT;

  const ThX = TgA + TgB, Ti4 = ThY + Ti3;
  const Ti5 = ThX + Ti4, Tib = Ti4 - ThX;
  const T3C = T36 + T3B, T4r = T45 + T4q;
  const T4s = T3C + T4r, Tia = T4r - T3C;

  const T5S = T5c + T5R, T7f = T6B + T7e;
  const T7g = T5S + T7f, Ti7 = T7f - T5S;
  const ThC = T5c - T5R, ThF = ThD - ThE;
  const ThG = ThC + ThF, ThO = ThF - ThC;

  const ThH = T6B - T7e, ThK = ThI - ThJ;
  const ThL = ThH - ThK, ThP = ThH + ThK;
  const Thy = TgE + TgF, Thz = TgK + TgL;
  const ThA = Thy - Thz, ThW = Thy + Thz;

  const T4t = T2L + T4s;
  outR[32] = T4t - T7g;
  outR[0] = T4t + T7g;
  const Ti6 = ThW + Ti5;
  outI[0] = ThV + Ti6;
  outI[32] = Ti6 - ThV;
  const ThR = T2L - T4s;
  outR[48] = ThR - ThU;
  outR[16] = ThR + ThU;
  const Ti8 = Ti5 - ThW;
  outI[16] = Ti7 + Ti8;
  outI[48] = Ti8 - Ti7;

  const ThB = Thx + ThA, ThM = KP707106781 * (ThG + ThL);
  outR[40] = ThB - ThM;
  outR[8] = ThB + ThM;
  const Ti9 = KP707106781 * (ThO + ThP), Tic = Tia + Tib;
  outI[8] = Ti9 + Tic;
  outI[40] = Tic - Ti9;

  const ThN = Thx - ThA, ThQ = KP707106781 * (ThO - ThP);
  outR[56] = ThN - ThQ;
  outR[24] = ThN + ThQ;
  const Tid = KP707106781 * (ThL - ThG), Tie = Tib - Tia;
  outI[24] = Tid + Tie;
  outI[56] = Tie - Tid;

  // -- combine group B: outputs 44, 12, 60, 28, 36, 4, 52, 20 --
  const Tgz = T17 - T1K, TgC = TgA - TgB;
  const TgD = Tgz - TgC, Thh = Tgz + TgC;
  const Thp = Th1 + Th4, Thq = Th8 + Th9;
  const Thr = KP923879532 * Thp - KP382683432 * Thq;
  const Thv = KP382683432 * Thp + KP923879532 * Thq;

  const Tih = T2J - T2e, Tii = Ti3 - ThY;
  const Tij = Tih + Tii, Tip = Tii - Tih;
  const TgI = TgG - TgH, TgN = TgJ + TgM;
  const TgO = KP707106781 * (TgI - TgN), Tig = KP707106781 * (TgI + TgN);

  const TgU = TgS - TgT, TgZ = TgV - TgY;
  const Th0 = KP923879532 * TgU + KP382683432 * TgZ;
  const The = KP382683432 * TgU - KP923879532 * TgZ;
  const Thi = TgH + TgG, Thj = TgJ - TgM;
  const Thk = KP707106781 * (Thi + Thj), Tio = KP707106781 * (Thj - Thi);

  const Thm = TgS + TgT, Thn = TgV + TgY;
  const Tho = KP382683432 * Thm + KP923879532 * Thn;
  const Thu = KP923879532 * Thm - KP382683432 * Thn;
  const Th5 = Th1 - Th4, Tha = Th8 - Th9;
  const Thb = KP382683432 * Th5 - KP923879532 * Tha;
  const Thf = KP382683432 * Tha + KP923879532 * Th5;

  const TgP = TgD + TgO, Thc = Th0 + Thb;
  outR[44] = TgP - Thc;
  outR[12] = TgP + Thc;
  const Tin = The + Thf, Tiq = Tio + Tip;
  outI[12] = Tin + Tiq;
  outI[44] = Tiq - Tin;

  const Thd = TgD - TgO, Thg = The - Thf;
  outR[60] = Thd - Thg;
  outR[28] = Thd + Thg;
  const Tir = Thb - Th0, Tis = Tip - Tio;
  outI[28] = Tir + Tis;
  outI[60] = Tis - Tir;

  const Thl = Thh + Thk, Ths = Tho + Thr;
  outR[36] = Thl - Ths;
  outR[4] = Thl + Ths;
  const Tif = Thu + Thv, Tik = Tig + Tij;
  outI[4] = Tif + Tik;
  outI[36] = Tik - Tif;

  const Tht = Thh - Thk, Thw = Thu - Thv;
  outR[52] = Tht - Thw;
  outR[20] = Tht + Thw;
  const Til = Thr - Tho, Tim = Tij - Tig;
  outI[20] = Til + Tim;
  outI[52] = Tim - Til;

  // -- combine group C: outputs 46, 14, 62, 30, 38, 6, 54, 22 --
  const TdZ = TdV - TdY, Tea = KP707106781 * (Te4 - Te9);
  const Teb = TdZ - Tea, Tfx = TdZ + Tea;

  const Tem = KP382683432 * Teg - KP923879532 * Tel;
  const Tex = KP382683432 * Ter + KP923879532 * Tew;
  const Tey = Tem - Tex, TiK = Tem + Tex;
  const TiL = KP707106781 * (TfP - TfO), TiM = Tix - Tiw;
  const TiN = TiL + TiM, TiT = TiM - TiL;

  const Tfy = KP923879532 * Teg + KP382683432 * Tel;
  const Tfz = KP382683432 * Tew - KP923879532 * Ter;
  const TfA = Tfy + Tfz, TiS = Tfz - Tfy;

  const Tfg = KP707106781 * (Tfa - Tff);
  const Tfh = Tf5 - Tfg, TfF = Tf5 + Tfg;
  const Tfp = KP707106781 * (Tfn - Tfo);
  const Tfq = Tfm - Tfp, TfG = Tfm + Tfp;
  const Tfr = KP195090322 * Tfh - KP980785280 * Tfq;
  const TfL = KP555570233 * TfF + KP831469612 * TfG;
  const Tfv = KP195090322 * Tfq + KP980785280 * Tfh;
  const TfH = KP831469612 * TfF - KP555570233 * TfG;

  const TeP = KP707106781 * (TeJ - TeO);
  const TeQ = TeE - TeP, TfC = TeE + TeP;
  const TeY = KP707106781 * (TeW - TeX);
  const TeZ = TeV - TeY, TfD = TeV + TeY;
  const Tf0 = KP980785280 * TeQ + KP195090322 * TeZ;
  const TfK = KP831469612 * TfC - KP555570233 * TfD;
  const Tfu = KP195090322 * TeQ - KP980785280 * TeZ;
  const TfE = KP555570233 * TfC + KP831469612 * TfD;

  const Tez = Teb + Tey, Tfs = Tf0 + Tfr;
  outR[46] = Tez - Tfs;
  outR[14] = Tez + Tfs;
  const TiR = Tfu + Tfv, TiU = TiS + TiT;
  outI[14] = TiR + TiU;
  outI[46] = TiU - TiR;

  const Tft = Teb - Tey, Tfw = Tfu - Tfv;
  outR[62] = Tft - Tfw;
  outR[30] = Tft + Tfw;
  const TiV = Tfr - Tf0, TiW = TiT - TiS;
  outI[30] = TiV + TiW;
  outI[62] = TiW - TiV;

  const TfB = Tfx + TfA, TfI = TfE + TfH;
  outR[38] = TfB - TfI;
  outR[6] = TfB + TfI;
  const TiJ = TfK + TfL, TiO = TiK + TiN;
  outI[6] = TiJ + TiO;
  outI[38] = TiO - TiJ;

  const TfJ = Tfx - TfA, TfM = TfK - TfL;
  outR[54] = TfJ - TfM;
  outR[22] = TfJ + TfM;
  const TiP = TfH - TfE, TiQ = TiN - TiK;
  outI[22] = TiP + TiQ;
  outI[54] = TiQ - TiP;

  // -- combine group D: outputs 42, 10, 58, 26, 34, 2, 50, 18 --
  const TfN = TdV + TdY, TfQ = KP707106781 * (TfO + TfP);
  const TfR = TfN - TfQ, Tgj = TfN + TfQ;

  const TfU = KP923879532 * TfS - KP382683432 * TfT;
  const TfX = KP923879532 * TfV + KP382683432 * TfW;
  const TfY = TfU - TfX, Tiu = TfU + TfX;
  const Tiv = KP707106781 * (Te4 + Te9), Tiy = Tiw + Tix;
  const Tiz = Tiv + Tiy, TiF = Tiy - Tiv;

  const Tgk = KP382683432 * TfS + KP923879532 * TfT;
  const Tgl = KP923879532 * TfW - KP382683432 * TfV;
  const Tgm = Tgk + Tgl, TiE = Tgl - Tgk;

  const Tg8 = KP707106781 * (Tfo + Tfn);
  const Tg9 = Tg7 - Tg8, Tgr = Tg7 + Tg8;
  const Tgb = KP707106781 * (Tfa + Tff);
  const Tgc = Tga - Tgb, Tgs = Tga + Tgb;
  const Tgd = KP555570233 * Tg9 - KP831469612 * Tgc;
  const Tgx = KP980785280 * Tgs + KP195090322 * Tgr;
  const Tgh = KP555570233 * Tgc + KP831469612 * Tg9;
  const Tgt = KP980785280 * Tgr - KP195090322 * Tgs;

  const Tg1 = KP707106781 * (TeO + TeJ);
  const Tg2 = Tg0 - Tg1, Tgo = Tg0 + Tg1;
  const Tg4 = KP707106781 * (TeW + TeX);
  const Tg5 = Tg3 - Tg4, Tgp = Tg3 + Tg4;
  const Tg6 = KP555570233 * Tg2 + KP831469612 * Tg5;
  const Tgw = KP980785280 * Tgp - KP195090322 * Tgo;
  const Tgg = KP555570233 * Tg5 - KP831469612 * Tg2;
  const Tgq = KP980785280 * Tgo + KP195090322 * Tgp;

  const TfZ = TfR + TfY, Tge = Tg6 + Tgd;
  outR[42] = TfZ - Tge;
  outR[10] = TfZ + Tge;
  const TiD = Tgg + Tgh, TiG = TiE + TiF;
  outI[10] = TiD + TiG;
  outI[42] = TiG - TiD;

  const Tgf = TfR - TfY, Tgi = Tgg - Tgh;
  outR[58] = Tgf - Tgi;
  outR[26] = Tgf + Tgi;
  const TiH = Tgd - Tg6, TiI = TiF - TiE;
  outI[26] = TiH + TiI;
  outI[58] = TiI - TiH;

  const Tgn = Tgj + Tgm, Tgu = Tgq + Tgt;
  outR[34] = Tgn - Tgu;
  outR[2] = Tgn + Tgu;
  const Tit = Tgw + Tgx, TiA = Tiu + Tiz;
  outI[2] = Tit + TiA;
  outI[34] = TiA - Tit;

  const Tgv = Tgj - Tgm, Tgy = Tgw - Tgx;
  outR[50] = Tgv - Tgy;
  outR[18] = Tgv + Tgy;
  const TiB = Tgt - Tgq, TiC = Tiz - Tiu;
  outI[18] = TiB + TiC;
  outI[50] = TiC - TiB;

  // -- combine group E: outputs 47, 15, 63, 31, 39, 7, 55, 23 --
  const T7x = T7l - T7w, T7U = T7I - T7T;
  const T7V = T7x - T7U, TaH = T7x + T7U;
  const TjL = TaZ - TaY, TjM = Tjx - Tjw;
  const TjN = TjL + TjM, TjT = TjM - TjL;

  const T8c = T80 - T8b, T8l = T8h - T8k;
  const T8m = KP195090322 * T8c - KP980785280 * T8l;
  const TaI = KP980785280 * T8c + KP195090322 * T8l;
  const T8D = T8r - T8C, T8M = T8I - T8L;
  const T8N = KP195090322 * T8D + KP980785280 * T8M;
  const TaJ = KP195090322 * T8M - KP980785280 * T8D;

  const T8O = T8m - T8N, TjS = TaJ - TaI;
  const TaK = TaI + TaJ, TjK = T8m + T8N;

  const T96 = T8U - T95, T9t = T9h - T9s;
  const T9u = T96 - T9t, TaM = T96 + T9t;
  const T9D = T9z - T9C, T9G = T9E - T9F;
  const T9H = T9D - T9G, TaN = T9D + T9G;

  const T9I = KP995184726 * T9u + KP098017140 * T9H;
  const TaU = KP773010453 * TaM - KP634393284 * TaN;
  const TaE = KP098017140 * T9u - KP995184726 * T9H;
  const TaO = KP634393284 * TaM + KP773010453 * TaN;

  const T9Z = T9N - T9Y, Tam = Taa - Tal;
  const Tan = T9Z - Tam, TaP = T9Z + Tam;
  const Taw = Tas - Tav, Taz = Tax - Tay;
  const TaA = Taw - Taz, TaQ = Taw + Taz;

  const TaB = KP098017140 * Tan - KP995184726 * TaA;
  const TaV = KP634393284 * TaP + KP773010453 * TaQ;
  const TaF = KP995184726 * Tan + KP098017140 * TaA;
  const TaR = KP773010453 * TaP - KP634393284 * TaQ;

  const T8P = T7V + T8O, TaC = T9I + TaB;
  outR[47] = T8P - TaC;
  outR[15] = T8P + TaC;
  const TjR = TaE + TaF, TjU = TjS + TjT;
  outI[15] = TjR + TjU;
  outI[47] = TjU - TjR;

  const TaD = T7V - T8O, TaG = TaE - TaF;
  outR[63] = TaD - TaG;
  outR[31] = TaD + TaG;
  const TjV = TaB - T9I, TjW = TjT - TjS;
  outI[31] = TjV + TjW;
  outI[63] = TjW - TjV;

  const TaL = TaH + TaK, TaS = TaO + TaR;
  outR[39] = TaL - TaS;
  outR[7] = TaL + TaS;
  const TjJ = TaU + TaV, TjO = TjK + TjN;
  outI[7] = TjJ + TjO;
  outI[39] = TjO - TjJ;

  const TaT = TaH - TaK, TaW = TaU - TaV;
  outR[55] = TaT - TaW;
  outR[23] = TaT + TaW;
  const TjP = TaR - TaO, TjQ = TjN - TjK;
  outI[23] = TjP + TjQ;
  outI[55] = TjQ - TjP;

  // -- combine group F: outputs 45, 13, 61, 29, 37, 5, 53, 21 --
  const TbN = TbJ - TbM, TbU = TbQ - TbT;
  const TbV = TbN - TbU, TcT = TbN + TbU;
  const Tjh = Tdb - Tda, Tji = Tj3 - Tj0;
  const Tjj = Tjh + Tji, Tjp = Tji - Tjh;

  const TbY = TbW - TbX, Tc1 = TbZ - Tc0;
  const Tc2 = KP555570233 * TbY - KP831469612 * Tc1;
  const TcU = KP555570233 * Tc1 + KP831469612 * TbY;
  const Tc5 = Tc3 - Tc4, Tc8 = Tc6 - Tc7;
  const Tc9 = KP831469612 * Tc5 + KP555570233 * Tc8;
  const TcV = KP555570233 * Tc5 - KP831469612 * Tc8;

  const Tca = Tc2 - Tc9, Tjo = TcV - TcU;
  const TcW = TcU + TcV, Tjg = Tc2 + Tc9;

  const Tce = Tcc - Tcd, Tcl = Tch - Tck;
  const Tcm = Tce - Tcl, TcY = Tce + Tcl;
  const Tcp = Tcn - Tco, Tcs = Tcq - Tcr;
  const Tct = Tcp - Tcs, TcZ = Tcp + Tcs;

  const Tcu = KP956940335 * Tcm + KP290284677 * Tct;
  const Td6 = KP881921264 * TcY - KP471396736 * TcZ;
  const TcQ = KP290284677 * Tcm - KP956940335 * Tct;
  const Td0 = KP471396736 * TcY + KP881921264 * TcZ;

  const Tcx = Tcv - Tcw, TcE = TcA - TcD;
  const TcF = Tcx - TcE, Td1 = Tcx + TcE;
  const TcI = TcG - TcH, TcL = TcJ - TcK;
  const TcM = TcI - TcL, Td2 = TcI + TcL;

  const TcN = KP290284677 * TcF - KP956940335 * TcM;
  const Td7 = KP471396736 * Td1 + KP881921264 * Td2;
  const TcR = KP956940335 * TcF + KP290284677 * TcM;
  const Td3 = KP881921264 * Td1 - KP471396736 * Td2;

  const Tcb = TbV + Tca, TcO = Tcu + TcN;
  outR[45] = Tcb - TcO;
  outR[13] = Tcb + TcO;
  const Tjn = TcQ + TcR, Tjq = Tjo + Tjp;
  outI[13] = Tjn + Tjq;
  outI[45] = Tjq - Tjn;

  const TcP = TbV - Tca, TcS = TcQ - TcR;
  outR[61] = TcP - TcS;
  outR[29] = TcP + TcS;
  const Tjr = TcN - Tcu, Tjs = Tjp - Tjo;
  outI[29] = Tjr + Tjs;
  outI[61] = Tjs - Tjr;

  const TcX = TcT + TcW, Td4 = Td0 + Td3;
  outR[37] = TcX - Td4;
  outR[5] = TcX + Td4;
  const Tjf = Td6 + Td7, Tjk = Tjg + Tjj;
  outI[5] = Tjf + Tjk;
  outI[37] = Tjk - Tjf;

  const Td5 = TcT - TcW, Td8 = Td6 - Td7;
  outR[53] = Td5 - Td8;
  outR[21] = Td5 + Td8;
  const Tjl = Td3 - Td0, Tjm = Tjj - Tjg;
  outI[21] = Tjl + Tjm;
  outI[53] = Tjm - Tjl;

  // -- combine group G: outputs 41, 9, 57, 25, 33, 1, 49, 17 --
  const Td9 = TbJ + TbM, Tdc = Tda + Tdb;
  const Tdd = Td9 - Tdc, TdF = Td9 + Tdc;
  const TiZ = TbQ + TbT, Tj4 = Tj0 + Tj3;
  const Tj5 = TiZ + Tj4, Tjb = Tj4 - TiZ;

  const Tde = TbW + TbX, Tdf = TbZ + Tc0;
  const Tdg = KP980785280 * Tde - KP195090322 * Tdf;
  const TdG = KP195090322 * Tde + KP980785280 * Tdf;
  const Tdh = Tc3 + Tc4, Tdi = Tc6 + Tc7;
  const Tdj = KP195090322 * Tdh + KP980785280 * Tdi;
  const TdH = KP980785280 * Tdh - KP195090322 * Tdi;

  const Tdk = Tdg - Tdj, Tja = TdH - TdG;
  const TdI = TdG + TdH, TiY = Tdg + Tdj;

  const Tdm = Tcn + Tco, Tdn = Tck + Tch;
  const Tdo = Tdm - Tdn, TdK = Tdm + Tdn;
  const Tdp = Tcc + Tcd, Tdq = Tcq + Tcr;
  const Tdr = Tdp - Tdq, TdL = Tdp + Tdq;

  const Tds = KP634393284 * Tdo + KP773010453 * Tdr;
  const TdS = KP995184726 * TdL - KP098017140 * TdK;
  const TdC = KP634393284 * Tdr - KP773010453 * Tdo;
  const TdM = KP995184726 * TdK + KP098017140 * TdL;

  const Tdt = Tcv + Tcw, Tdu = TcK + TcJ;
  const Tdv = Tdt - Tdu, TdN = Tdt + Tdu;
  const Tdw = TcG + TcH, Tdx = TcA + TcD;
  const Tdy = Tdw - Tdx, TdO = Tdw + Tdx;

  const Tdz = KP634393284 * Tdv - KP773010453 * Tdy;
  const TdT = KP995184726 * TdO + KP098017140 * TdN;
  const TdD = KP634393284 * Tdy + KP773010453 * Tdv;
  const TdP = KP995184726 * TdN - KP098017140 * TdO;

  const Tdl = Tdd + Tdk, TdA = Tds + Tdz;
  outR[41] = Tdl - TdA;
  outR[9] = Tdl + TdA;
  const Tj9 = TdC + TdD, Tjc = Tja + Tjb;
  outI[9] = Tj9 + Tjc;
  outI[41] = Tjc - Tj9;

  const TdB = Tdd - Tdk, TdE = TdC - TdD;
  outR[57] = TdB - TdE;
  outR[25] = TdB + TdE;
  const Tjd = Tdz - Tds, Tje = Tjb - Tja;
  outI[25] = Tjd + Tje;
  outI[57] = Tje - Tjd;

  const TdJ = TdF + TdI, TdQ = TdM + TdP;
  outR[33] = TdJ - TdQ;
  outR[1] = TdJ + TdQ;
  const TiX = TdS + TdT, Tj6 = TiY + Tj5;
  outI[1] = TiX + Tj6;
  outI[33] = Tj6 - TiX;

  const TdR = TdF - TdI, TdU = TdS - TdT;
  outR[49] = TdR - TdU;
  outR[17] = TdR + TdU;
  const Tj7 = TdP - TdM, Tj8 = Tj5 - TiY;
  outI[17] = Tj7 + Tj8;
  outI[49] = Tj8 - Tj7;

  // -- combine group H: outputs 43, 11, 59, 27, 35, 3, 51, 19 --
  const TaX = T7l + T7w, Tb0 = TaY + TaZ;
  const Tb1 = TaX - Tb0, Tbt = TaX + Tb0;
  const Tjv = T7I + T7T, Tjy = Tjw + Tjx;
  const Tjz = Tjv + Tjy, TjF = Tjy - Tjv;

  const Tb2 = T80 + T8b, Tb3 = T8h + T8k;
  const Tb4 = KP831469612 * Tb2 - KP555570233 * Tb3;
  const Tbu = KP555570233 * Tb2 + KP831469612 * Tb3;
  const Tb5 = T8r + T8C, Tb6 = T8I + T8L;
  const Tb7 = KP831469612 * Tb5 + KP555570233 * Tb6;
  const Tbv = KP831469612 * Tb6 - KP555570233 * Tb5;

  const Tb8 = Tb4 - Tb7, TjE = Tbv - Tbu;
  const Tbw = Tbu + Tbv, Tju = Tb4 + Tb7;

  const Tba = T9z + T9C, Tbb = T9s + T9h;
  const Tbc = Tba - Tbb, Tby = Tba + Tbb;
  const Tbd = T8U + T95, Tbe = T9E + T9F;
  const Tbf = Tbd - Tbe, Tbz = Tbd + Tbe;

  const Tbg = KP471396736 * Tbc + KP881921264 * Tbf;
  const TbG = KP956940335 * Tbz - KP290284677 * Tby;
  const Tbq = KP471396736 * Tbf - KP881921264 * Tbc;
  const TbA = KP956940335 * Tby + KP290284677 * Tbz;

  const Tbh = T9N + T9Y, Tbi = Tay + Tax;
  const Tbj = Tbh - Tbi, TbB = Tbh + Tbi;
  const Tbk = Tas + Tav, Tbl = Taa + Tal;
  const Tbm = Tbk - Tbl, TbC = Tbk + Tbl;

  const Tbn = KP471396736 * Tbj - KP881921264 * Tbm;
  const TbH = KP956940335 * TbC + KP290284677 * TbB;
  const Tbr = KP881921264 * Tbj + KP471396736 * Tbm;
  const TbD = KP956940335 * TbB - KP290284677 * TbC;

  const Tb9 = Tb1 + Tb8, Tbo = Tbg + Tbn;
  outR[43] = Tb9 - Tbo;
  outR[11] = Tb9 + Tbo;
  const TjD = Tbq + Tbr, TjG = TjE + TjF;
  outI[11] = TjD + TjG;
  outI[43] = TjG - TjD;

  const Tbp = Tb1 - Tb8, Tbs = Tbq - Tbr;
  outR[59] = Tbp - Tbs;
  outR[27] = Tbp + Tbs;
  const TjH = Tbn - Tbg, TjI = TjF - TjE;
  outI[27] = TjH + TjI;
  outI[59] = TjI - TjH;

  const Tbx = Tbt + Tbw, TbE = TbA + TbD;
  outR[35] = Tbx - TbE;
  outR[3] = Tbx + TbE;
  const Tjt = TbG + TbH, TjA = Tju + Tjz;
  outI[3] = Tjt + TjA;
  outI[35] = TjA - Tjt;

  const TbF = Tbt - Tbw, TbI = TbG - TbH;
  outR[51] = TbF - TbI;
  outR[19] = TbF + TbI;
  const TjB = TbD - TbA, TjC = Tjz - Tju;
  outI[19] = TjB + TjC;
  outI[51] = TjC - TjB;

  return [outR, outI];
}

module.exports = { t2_64 };
