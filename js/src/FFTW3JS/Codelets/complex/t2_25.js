'use strict';

// =============================================================================
// t2_25.js -- faithful JS port of dft/scalar/codelets/t2_25.c (non-FMA
// branch), FFTW3's "twiddle-log3 / precompute-twiddles" radix-25 twiddle
// codelet. twinstr only trig-generates W^1, W^3, W^9, W^24 (four raw
// pairs); the other 20 needed multiples are DERIVED via one complex
// multiply each (same trick as t2_5/t2_8/t2_10/t2_16/t2_20.js). Same
// calling convention as every other twiddle codelet here (br/bi/Wc/Ws with
// Composite1D.js's full r-1 = 24 pair table already built) -- only indices
// 1, 3, 9, 24 are actually read.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP998026728 = 0.998026728428271561952336806863450553336905220;
const KP062790519 = 0.062790519529313376076178224565631133122484832;
const KP425779291 = 0.425779291565072648862502445744251703979973042;
const KP904827052 = 0.904827052466019527713668647932697593970413911;
const KP992114701 = 0.992114701314477831049793042785778521453036709;
const KP125333233 = 0.125333233564304245373118759816508793942918247;
const KP637423989 = 0.637423989748689710176712811676016195434917298;
const KP770513242 = 0.770513242775789230803009636396177847271667672;
const KP684547105 = 0.684547105928688673732283357621209269889519233;
const KP728968627 = 0.728968627421411523146730319055259111372571664;
const KP481753674 = 0.481753674101715274987191502872129653528542010;
const KP876306680 = 0.876306680043863587308115903922062583399064238;
const KP844327925 = 0.844327925502015078548558063966681505381659241;
const KP535826794 = 0.535826794978996618271308767867639978063575346;
const KP248689887 = 0.248689887164854788242283746006447968417567406;
const KP968583161 = 0.968583161128631119490168375464735813836012403;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;

function t2_25(br, bi, Wc, Ws) {
  const T2 = Wc[1], T5 = Ws[1], T3 = Wc[3], T6 = Ws[3];
  const T4 = T2 * T3, Tc = T5 * T3, T7 = T5 * T6, Tb = T2 * T6;
  const T8 = T4 - T7, Td = Tb + Tc, T16 = Tb - Tc, T14 = T4 + T7;
  const Te = Ws[9];
  const Tw = T5 * Te, TT = T3 * Te, Tz = T2 * Te, TQ = T6 * Te;
  const T9 = Wc[9];
  const Tv = T2 * T9, TU = T6 * T9, TA = T5 * T9, TP = T3 * T9;

  const T21 = TP - TQ, T23 = TT + TU;

  const Tx = Tv - Tw, TR = TP + TQ, T1g = Tz - TA, TB = Tz + TA, T1f = Tv + Tw, TV = TT - TU;
  const T15 = T14 * T9, T17 = T16 * Te, T1Q = T15 + T17;
  const Ta = T8 * T9, Tf = Td * Te, Tg = Ta + Tf;
  const T1a = T14 * Te, T1b = T16 * T9, T1S = T1a - T1b;
  const Ti = T8 * Te, Tj = Td * T9, Tk = Ti - Tj;
  const T18 = T15 - T17, T2s = Ti + Tj, T1c = T1a + T1b, T2q = Ta - Tf;

  const Tn = Wc[24], To = Ws[24];
  const Tp = T8 * Tn + Td * To;
  const Tr = T8 * To - Td * Tn;
  const T28 = T1Q * To - T1S * Tn;
  const T2x = TR * To - TV * Tn;
  const TY = T3 * Tn + T6 * To;
  const T2k = T2 * Tn + T5 * To;
  const T2m = T2 * To - T5 * Tn;
  const T2v = TR * Tn + TV * To;
  const TG = T9 * To - Te * Tn;
  const TE = T9 * Tn + Te * To;
  const T10 = T3 * To - T6 * Tn;
  const T1h = T1f * Tn + T1g * To;
  const T1E = Tg * Tn + Tk * To;
  const T26 = T1Q * Tn + T1S * To;
  const T1B = Tx * To - TB * Tn;
  const T1G = Tg * To - Tk * Tn;
  const T1V = T14 * Tn + T16 * To;
  const T1X = T14 * To - T16 * Tn;
  const T1z = Tx * Tn + TB * To;
  const T1j = T1f * To - T1g * Tn;

  const T1 = br[0], T6v = bi[0];

  const Th = br[5], Tl = bi[5];
  const Tm = Tg * Th + Tk * Tl;
  const T2I = Tg * Tl - Tk * Th;
  const Tq = br[20], Ts = bi[20];
  const Tt = Tp * Tq + Tr * Ts;
  const T2J = Tp * Ts - Tr * Tq;

  const Tu = Tm + Tt, T6s = T2I + T2J;

  const Ty = br[10], TC = bi[10];
  const TD = Tx * Ty + TB * TC;
  const T2L = Tx * TC - TB * Ty;
  const TF = br[15], TH = bi[15];
  const TI = TE * TF + TG * TH;
  const T2M = TE * TH - TG * TF;

  const TJ = TD + TI, T6t = T2L + T2M;
  const T2F = KP559016994 * (Tu - TJ), T6I = KP559016994 * (T6s - T6t);
  const TK = Tu + TJ, T2G = T1 - KP250000000 * TK;
  const T6u = T6s + T6t, T6J = T6v - KP250000000 * T6u;

  const T6L = Tm - Tt, T6M = TD - TI;
  const T6N = KP951056516 * T6L + KP587785252 * T6M;
  const T7c = KP951056516 * T6M - KP587785252 * T6L;
  const T2K = T2I - T2J, T2N = T2L - T2M;
  const T2O = KP951056516 * T2K + KP587785252 * T2N;
  const T52 = KP951056516 * T2N - KP587785252 * T2K;

  const T2e = br[3], T2f = bi[3];
  const T2g = T3 * T2e + T6 * T2f;
  const T4c = T3 * T2f - T6 * T2e;

  const T2h = br[8], T2i = bi[8];
  const T2j = T1f * T2h + T1g * T2i;
  const T41 = T1f * T2i - T1g * T2h;
  const T2w = br[18], T2y = bi[18];
  const T2z = T2v * T2w + T2x * T2y;
  const T45 = T2v * T2y - T2x * T2w;

  const T2l = br[23], T2n = bi[23];
  const T2o = T2k * T2l + T2m * T2n;
  const T42 = T2k * T2n - T2m * T2l;
  const T2r = br[13], T2t = bi[13];
  const T2u = T2q * T2r + T2s * T2t;
  const T44 = T2q * T2t - T2s * T2r;

  const T43 = T41 - T42, T46 = T44 - T45, T4h = T2u - T2z, T4g = T2j - T2o;
  const T49 = T41 + T42, T4a = T44 + T45, T4d = T49 + T4a;
  const T2p = T2j + T2o, T2A = T2u + T2z, T2B = T2p + T2A;

  const T2C = T2g + T2B, T6k = T4c + T4d;

  const T47 = KP951056516 * T43 + KP587785252 * T46;
  const T5r = KP951056516 * T46 - KP587785252 * T43;
  const T3Y = KP559016994 * (T2p - T2A), T3Z = T2g - KP250000000 * T2B;
  const T40 = T3Y + T3Z, T5q = T3Z - T3Y;
  const T48 = T40 + T47, T5X = T5q + T5r;
  const T4L = T40 - T47, T5s = T5q - T5r;

  const T4i = KP951056516 * T4g + KP587785252 * T4h;
  const T5t = KP951056516 * T4h - KP587785252 * T4g;
  const T4b = KP559016994 * (T49 - T4a), T4e = T4c - KP250000000 * T4d;
  const T4f = T4b + T4e, T5u = T4e - T4b;
  const T4j = T4f - T4i, T5W = T5u - T5t;
  const T4K = T4i + T4f, T5v = T5t + T5u;

  const TM = br[1], TN = bi[1];
  const TO = T2 * TM + T5 * TN;
  const T34 = T2 * TN - T5 * TM;

  const TS = br[6], TW = bi[6];
  const TX = TR * TS + TV * TW;
  const T2T = TR * TW - TV * TS;
  const T1i = br[16], T1k = bi[16];
  const T1l = T1h * T1i + T1j * T1k;
  const T2X = T1h * T1k - T1j * T1i;

  const TZ = br[21], T11 = bi[21];
  const T12 = TY * TZ + T10 * T11;
  const T2U = TY * T11 - T10 * TZ;
  const T19 = br[11], T1d = bi[11];
  const T1e = T18 * T19 + T1c * T1d;
  const T2W = T18 * T1d - T1c * T19;

  const T2V = T2T - T2U, T2Y = T2W - T2X, T39 = T1e - T1l, T38 = TX - T12;
  const T31 = T2T + T2U, T32 = T2W + T2X, T35 = T31 + T32;
  const T13 = TX + T12, T1m = T1e + T1l, T1n = T13 + T1m;

  const T1o = TO + T1n, T6g = T34 + T35;

  const T2Z = KP951056516 * T2V + KP587785252 * T2Y;
  const T55 = KP951056516 * T2Y - KP587785252 * T2V;
  const T2Q = KP559016994 * (T13 - T1m), T2R = TO - KP250000000 * T1n;
  const T2S = T2Q + T2R, T54 = T2R - T2Q;
  const T30 = T2S + T2Z, T5M = T54 + T55;
  const T4A = T2S - T2Z, T56 = T54 - T55;

  const T3a = KP951056516 * T38 + KP587785252 * T39;
  const T57 = KP951056516 * T39 - KP587785252 * T38;
  const T33 = KP559016994 * (T31 - T32), T36 = T34 - KP250000000 * T35;
  const T37 = T33 + T36, T58 = T36 - T33;
  const T3b = T37 - T3a, T5N = T58 - T57;
  const T4B = T3a + T37, T59 = T57 + T58;

  const T1p = br[4], T1q = bi[4];
  const T1r = T8 * T1p + Td * T1q;
  const T3r = T8 * T1q - Td * T1p;

  const T1s = br[9], T1t = bi[9];
  const T1u = T9 * T1s + Te * T1t;
  const T3g = T9 * T1t - Te * T1s;
  const T1F = br[19], T1H = bi[19];
  const T1I = T1E * T1F + T1G * T1H;
  const T3k = T1E * T1H - T1G * T1F;

  const T1v = br[24], T1w = bi[24];
  const T1x = Tn * T1v + To * T1w;
  const T3h = Tn * T1w - To * T1v;
  const T1A = br[14], T1C = bi[14];
  const T1D = T1z * T1A + T1B * T1C;
  const T3j = T1z * T1C - T1B * T1A;

  const T3i = T3g - T3h, T3l = T3j - T3k, T3w = T1D - T1I, T3v = T1u - T1x;
  const T3o = T3g + T3h, T3p = T3j + T3k, T3s = T3o + T3p;
  const T1y = T1u + T1x, T1J = T1D + T1I, T1K = T1y + T1J;

  const T1L = T1r + T1K, T6h = T3r + T3s;

  const T3m = KP951056516 * T3i + KP587785252 * T3l;
  const T5f = KP951056516 * T3l - KP587785252 * T3i;
  const T3d = KP559016994 * (T1y - T1J), T3e = T1r - KP250000000 * T1K;
  const T3f = T3d + T3e, T5e = T3e - T3d;
  const T3n = T3f + T3m, T5Q = T5e + T5f;
  const T4D = T3f - T3m, T5g = T5e - T5f;

  const T3x = KP951056516 * T3v + KP587785252 * T3w;
  const T5b = KP951056516 * T3w - KP587785252 * T3v;
  const T3q = KP559016994 * (T3o - T3p), T3t = T3r - KP250000000 * T3s;
  const T3u = T3q + T3t, T5c = T3t - T3q;
  const T3y = T3u - T3x, T5P = T5c - T5b;
  const T4E = T3x + T3u, T5d = T5b + T5c;

  const T1N = br[2], T1O = bi[2];
  const T1P = T14 * T1N + T16 * T1O;
  const T3P = T14 * T1O - T16 * T1N;

  const T1R = br[7], T1T = bi[7];
  const T1U = T1Q * T1R + T1S * T1T;
  const T3E = T1Q * T1T - T1S * T1R;
  const T27 = br[17], T29 = bi[17];
  const T2a = T26 * T27 + T28 * T29;
  const T3I = T26 * T29 - T28 * T27;

  const T1W = br[22], T1Y = bi[22];
  const T1Z = T1V * T1W + T1X * T1Y;
  const T3F = T1V * T1Y - T1X * T1W;
  const T22 = br[12], T24 = bi[12];
  const T25 = T21 * T22 + T23 * T24;
  const T3H = T21 * T24 - T23 * T22;

  const T3G = T3E - T3F, T3J = T3H - T3I, T3U = T25 - T2a, T3T = T1U - T1Z;
  const T3M = T3E + T3F, T3N = T3H + T3I, T3Q = T3M + T3N;
  const T20 = T1U + T1Z, T2b = T25 + T2a, T2c = T20 + T2b;

  const T2d = T1P + T2c, T6j = T3P + T3Q;

  const T3K = KP951056516 * T3G + KP587785252 * T3J;
  const T5k = KP951056516 * T3J - KP587785252 * T3G;
  const T3B = KP559016994 * (T20 - T2b), T3C = T1P - KP250000000 * T2c;
  const T3D = T3B + T3C, T5j = T3C - T3B;
  const T3L = T3D + T3K, T5T = T5j + T5k;
  const T4I = T3D - T3K, T5l = T5j - T5k;

  const T3V = KP951056516 * T3T + KP587785252 * T3U;
  const T5m = KP951056516 * T3U - KP587785252 * T3T;
  const T3O = KP559016994 * (T3M - T3N), T3R = T3P - KP250000000 * T3Q;
  const T3S = T3O + T3R, T5n = T3R - T3O;
  const T3W = T3S - T3V, T5U = T5n - T5m;
  const T4H = T3V + T3S, T5o = T5m + T5n;

  const outR = new Float64Array(25), outI = new Float64Array(25);

  {
    const T6i = T6g - T6h, T6l = T6j - T6k;
    const T6m = KP951056516 * T6i + KP587785252 * T6l;
    const T6o = KP951056516 * T6l - KP587785252 * T6i;
    const TL = T1 + TK;
    const T1M = T1o + T1L;
    const T2D = T2d + T2C;
    const T2E = T1M + T2D;
    const T6d = KP559016994 * (T1M - T2D);
    const T6e = TL - KP250000000 * T2E;

    outR[0] = TL + T2E;
    const T6n = T6e - T6d;
    outR[10] = T6n - T6o;
    outR[15] = T6n + T6o;
    const T6f = T6d + T6e;
    outR[20] = T6f - T6m;
    outR[5] = T6f + T6m;
  }

  {
    const T6A = T1o - T1L, T6B = T2d - T2C;
    const T6C = KP951056516 * T6A + KP587785252 * T6B;
    const T6D = KP951056516 * T6B - KP587785252 * T6A;
    const T6w = T6u + T6v;
    const T6p = T6g + T6h;
    const T6q = T6j + T6k;
    const T6r = T6p + T6q;
    const T6x = KP559016994 * (T6p - T6q);
    const T6y = T6w - KP250000000 * T6r;

    outI[0] = T6r + T6w;
    const T6E = T6y - T6x;
    outI[10] = T6D + T6E;
    outI[15] = T6E - T6D;
    const T6z = T6x + T6y;
    outI[5] = T6z - T6C;
    outI[20] = T6C + T6z;
  }

  {
    const T2H = T2F + T2G;
    const T2P = T2H + T2O;
    const T4z = T2H - T2O;
    const T6K = T6I + T6J;
    const T6O = T6K - T6N;
    const T70 = T6N + T6K;

    const T3c = KP968583161 * T30 + KP248689887 * T3b;
    const T3z = KP535826794 * T3n + KP844327925 * T3y;
    const T3A = T3c + T3z;
    const T3X = KP876306680 * T3L + KP481753674 * T3W;
    const T4k = KP728968627 * T48 + KP684547105 * T4j;
    const T4l = T3X + T4k;
    const T4m = T3A + T4l;
    const T6T = T3X - T4k;
    const T4n = KP559016994 * (T3A - T4l);
    const T6S = T3c - T3z;

    const T4S = KP535826794 * T4B - KP844327925 * T4A;
    const T4T = KP770513242 * T4D - KP637423989 * T4E;
    const T6X = T4S + T4T;
    const T4V = KP125333233 * T4L + KP992114701 * T4K;
    const T4W = KP904827052 * T4I + KP425779291 * T4H;
    const T6Y = T4W + T4V;
    const T4U = T4S - T4T;
    const T71 = KP559016994 * (T6X + T6Y);
    const T4X = T4V - T4W;
    const T6Z = T6X - T6Y;

    const T4C = KP535826794 * T4A + KP844327925 * T4B;
    const T4F = KP637423989 * T4D + KP770513242 * T4E;
    const T4G = T4C - T4F;
    const T4J = KP904827052 * T4H - KP425779291 * T4I;
    const T4M = KP125333233 * T4K - KP992114701 * T4L;
    const T4N = T4J + T4M;
    const T4O = T4G + T4N;
    const T75 = T4J - T4M;
    const T4P = KP559016994 * (T4G - T4N);
    const T74 = T4C + T4F;

    const T4q = KP968583161 * T3b - KP248689887 * T30;
    const T4r = KP535826794 * T3y - KP844327925 * T3n;
    const T6F = T4q + T4r;
    const T4t = KP876306680 * T3W - KP481753674 * T3L;
    const T4u = KP728968627 * T4j - KP684547105 * T48;
    const T6G = T4t + T4u;
    const T4s = T4q - T4r;
    const T6P = KP559016994 * (T6F - T6G);
    const T4v = T4t - T4u;
    const T6H = T6F + T6G;

    outR[1] = T2P + T4m;
    outI[1] = T6H + T6O;
    outR[4] = T4z + T4O;
    outI[4] = T6Z + T70;

    const T4w = KP951056516 * T4s + KP587785252 * T4v;
    const T4y = KP951056516 * T4v - KP587785252 * T4s;
    const T4o = T2P - KP250000000 * T4m;
    const T4p = T4n + T4o;
    const T4x = T4o - T4n;
    outR[21] = T4p - T4w;
    outR[16] = T4x + T4y;
    outR[6] = T4p + T4w;
    outR[11] = T4x - T4y;

    const T6U = KP951056516 * T6S + KP587785252 * T6T;
    const T6V = KP951056516 * T6T - KP587785252 * T6S;
    const T6Q = T6O - KP250000000 * T6H;
    const T6R = T6P + T6Q;
    const T6W = T6Q - T6P;
    outI[6] = T6R - T6U;
    outI[16] = T6W - T6V;
    outI[21] = T6U + T6R;
    outI[11] = T6V + T6W;

    const T4Y = KP951056516 * T4U + KP587785252 * T4X;
    const T50 = KP951056516 * T4X - KP587785252 * T4U;
    const T4Q = T4z - KP250000000 * T4O;
    const T4R = T4P + T4Q;
    const T4Z = T4Q - T4P;
    outR[24] = T4R - T4Y;
    outR[19] = T4Z + T50;
    outR[9] = T4R + T4Y;
    outR[14] = T4Z - T50;

    const T76 = KP951056516 * T74 + KP587785252 * T75;
    const T77 = KP951056516 * T75 - KP587785252 * T74;
    const T72 = T70 - KP250000000 * T6Z;
    const T73 = T71 + T72;
    const T78 = T72 - T71;
    outI[9] = T73 - T76;
    outI[19] = T78 - T77;
    outI[24] = T76 + T73;
    outI[14] = T77 + T78;
  }

  {
    const T51 = T2G - T2F;
    const T53 = T51 - T52;
    const T5L = T51 + T52;
    const T7d = T6J - T6I;
    const T7e = T7c + T7d;
    const T7q = T7d - T7c;

    const T5a = KP876306680 * T56 + KP481753674 * T59;
    const T5h = KP904827052 * T5d - KP425779291 * T5g;
    const T5i = T5a + T5h;
    const T5p = KP535826794 * T5l + KP844327925 * T5o;
    const T5w = KP062790519 * T5s + KP998026728 * T5v;
    const T5x = T5p + T5w;
    const T5y = T5i + T5x;
    const T7j = T5p - T5w;
    const T5z = KP559016994 * (T5i - T5x);
    const T7i = T5a - T5h;

    const T64 = KP728968627 * T5N - KP684547105 * T5M;
    const T65 = KP125333233 * T5Q + KP992114701 * T5P;
    const T7n = T64 - T65;
    const T67 = KP062790519 * T5U - KP998026728 * T5T;
    const T68 = KP770513242 * T5X + KP637423989 * T5W;
    const T7o = T67 - T68;
    const T66 = T64 + T65;
    const T7r = KP559016994 * (T7n - T7o);
    const T69 = T67 + T68;
    const T7p = T7n + T7o;

    const T5O = KP728968627 * T5M + KP684547105 * T5N;
    const T5R = KP125333233 * T5P - KP992114701 * T5Q;
    const T5S = T5O + T5R;
    const T5V = KP062790519 * T5T + KP998026728 * T5U;
    const T5Y = KP770513242 * T5W - KP637423989 * T5X;
    const T5Z = T5V + T5Y;
    const T60 = T5S + T5Z;
    const T7v = T5V - T5Y;
    const T61 = KP559016994 * (T5S - T5Z);
    const T7u = T5O - T5R;

    const T5C = KP876306680 * T59 - KP481753674 * T56;
    const T5D = KP904827052 * T5g + KP425779291 * T5d;
    const T79 = T5C - T5D;
    const T5F = KP535826794 * T5o - KP844327925 * T5l;
    const T5G = KP062790519 * T5v - KP998026728 * T5s;
    const T7a = T5F + T5G;
    const T5E = T5C + T5D;
    const T7f = KP559016994 * (T79 - T7a);
    const T5H = T5F - T5G;
    const T7b = T79 + T7a;

    outR[2] = T53 + T5y;
    outI[2] = T7b + T7e;
    outR[3] = T5L + T60;
    outI[3] = T7p + T7q;

    const T5I = KP951056516 * T5E + KP587785252 * T5H;
    const T5K = KP951056516 * T5H - KP587785252 * T5E;
    const T5A = T53 - KP250000000 * T5y;
    const T5B = T5z + T5A;
    const T5J = T5A - T5z;
    outR[22] = T5B - T5I;
    outR[17] = T5J + T5K;
    outR[7] = T5B + T5I;
    outR[12] = T5J - T5K;

    const T7k = KP951056516 * T7i + KP587785252 * T7j;
    const T7l = KP951056516 * T7j - KP587785252 * T7i;
    const T7g = T7e - KP250000000 * T7b;
    const T7h = T7f + T7g;
    const T7m = T7g - T7f;
    outI[7] = T7h - T7k;
    outI[17] = T7m - T7l;
    outI[22] = T7k + T7h;
    outI[12] = T7l + T7m;

    const T6a = KP951056516 * T66 + KP587785252 * T69;
    const T6c = KP951056516 * T69 - KP587785252 * T66;
    const T62 = T5L - KP250000000 * T60;
    const T63 = T61 + T62;
    const T6b = T62 - T61;
    outR[23] = T63 - T6a;
    outR[18] = T6b + T6c;
    outR[8] = T63 + T6a;
    outR[13] = T6b - T6c;

    const T7w = KP951056516 * T7u + KP587785252 * T7v;
    const T7x = KP951056516 * T7v - KP587785252 * T7u;
    const T7s = T7q - KP250000000 * T7p;
    const T7t = T7r + T7s;
    const T7y = T7s - T7r;
    outI[8] = T7t - T7w;
    outI[18] = T7y - T7x;
    outI[23] = T7w + T7t;
    outI[13] = T7x + T7y;
  }

  return [outR, outI];
}

module.exports = { t2_25 };
