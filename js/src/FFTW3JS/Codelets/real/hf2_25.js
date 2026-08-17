'use strict';

// =============================================================================
// hf2_25.js -- faithful JS port of rdft/scalar/r2cf/hf2_25.c (non-FMA
// branch). Alternate-codegen ("twiddle-log3/precompute-twiddles") sibling
// of the (unported, dormant) hf_25.js -- same math, different rounding
// (see hf2_5.js's header for the pattern). twinstr only trig-generates
// W^1, W^3, W^9, W^24 (four raw pairs, same generator set as the complex
// side's t2_25.js); every other needed multiple is DERIVED via complex
// products of those four. Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention
// as hf_25/hf2_5.js -- only indices 1, 3, 9, 24 are actually read from
// Wc/Ws.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP998026728 = 0.998026728428271561952336806863450553336905220;
const KP062790519 = 0.062790519529313376076178224565631133122484832;
const KP684547105 = 0.684547105928688673732283357621209269889519233;
const KP728968627 = 0.728968627421411523146730319055259111372571664;
const KP481753674 = 0.481753674101715274987191502872129653528542010;
const KP876306680 = 0.876306680043863587308115903922062583399064238;
const KP248689887 = 0.248689887164854788242283746006447968417567406;
const KP968583161 = 0.968583161128631119490168375464735813836012403;
const KP992114701 = 0.992114701314477831049793042785778521453036709;
const KP125333233 = 0.125333233564304245373118759816508793942918247;
const KP425779291 = 0.425779291565072648862502445744251703979973042;
const KP904827052 = 0.904827052466019527713668647932697593970413911;
const KP637423989 = 0.637423989748689710176712811676016195434917298;
const KP770513242 = 0.770513242775789230803009636396177847271667672;
const KP844327925 = 0.844327925502015078548558063966681505381659241;
const KP535826794 = 0.535826794978996618271308767867639978063575346;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;

function hf2_25(cr, ci, Wc, Ws) {
  const T2 = Wc[1], T5 = Ws[1], T3 = Wc[3], T6 = Ws[3];
  const T4 = T2 * T3, Tc = T5 * T3, T7 = T5 * T6, Tb = T2 * T6;
  const T8 = T4 - T7, Td = Tb + Tc, T16 = Tb - Tc, T14 = T4 + T7;
  const Te = Ws[9];
  const Tw = T5 * Te, TT = T3 * Te, Tz = T2 * Te, TQ = T6 * Te;
  const T9 = Wc[9];
  const Tv = T2 * T9, TU = T6 * T9, TA = T5 * T9, TP = T3 * T9;

  const T21 = TP - TQ, T23 = TT + TU;

  const T15 = T14 * T9, T17 = T16 * Te, T1Q = T15 + T17;
  const Ta = T8 * T9, Tf = Td * Te, Tg = Ta + Tf;
  const T1a = T14 * Te, T1b = T16 * T9, T1S = T1a - T1b;
  const Ti = T8 * Te, Tj = Td * T9, Tk = Ti - Tj;
  const T18 = T15 - T17, T2s = Ti + Tj, T1c = T1a + T1b, T2q = Ta - Tf;

  const Tx = Tv - Tw, TR = TP + TQ, T1g = Tz - TA, TB = Tz + TA, T1f = Tv + Tw, TV = TT - TU;

  const Tn = Wc[24], To = Ws[24];

  const Tp = T8 * Tn + Td * To, Tr = T8 * To - Td * Tn;
  const T28 = T1Q * To - T1S * Tn, T2x = TR * To - TV * Tn;
  const TY = T3 * Tn + T6 * To, T2k = T2 * Tn + T5 * To;
  const T2m = T2 * To - T5 * Tn, T2v = TR * Tn + TV * To;
  const TG = T9 * To - Te * Tn, TE = T9 * Tn + Te * To;
  const T10 = T3 * To - T6 * Tn;
  const T1h = T1f * Tn + T1g * To, T1E = Tg * Tn + Tk * To;
  const T26 = T1Q * Tn + T1S * To, T1B = Tx * To - TB * Tn;
  const T1G = Tg * To - Tk * Tn, T1V = T14 * Tn + T16 * To;
  const T1X = T14 * To - T16 * Tn, T1z = Tx * Tn + TB * To;
  const T1j = T1f * To - T1g * Tn;

  const T1 = cr[0], T6v = ci[0];

  const Th = cr[5], Tl = ci[5];
  const Tm = Tg * Th + Tk * Tl;
  const T2I = Tg * Tl - Tk * Th;
  const Tq = cr[20], Ts = ci[20];
  const Tt = Tp * Tq + Tr * Ts;
  const T2J = Tp * Ts - Tr * Tq;

  const Tu = Tm + Tt, T6w = T2I + T2J;

  const Ty = cr[10], TC = ci[10];
  const TD = Tx * Ty + TB * TC;
  const T2L = Tx * TC - TB * Ty;
  const TF = cr[15], TH = ci[15];
  const TI = TE * TF + TG * TH;
  const T2M = TE * TH - TG * TF;

  const TJ = TD + TI, T6x = T2L + T2M;
  const T2F = KP559016994 * (Tu - TJ), T6A = KP559016994 * (T6w - T6x);
  const TK = Tu + TJ, T2G = T1 - KP250000000 * TK;
  const T6y = T6w + T6x, T6z = T6v - KP250000000 * T6y;

  const T6s = TD - TI, T6t = Tm - Tt;
  const T6u = KP951056516 * T6s - KP587785252 * T6t;
  const T71 = KP951056516 * T6t + KP587785252 * T6s;
  const T2K = T2I - T2J, T2N = T2L - T2M;
  const T2O = KP951056516 * T2K + KP587785252 * T2N;
  const T52 = KP951056516 * T2N - KP587785252 * T2K;

  const T2e = cr[3], T2f = ci[3];
  const T2g = T3 * T2e + T6 * T2f;
  const T48 = T3 * T2f - T6 * T2e;

  const T2h = cr[8], T2i = ci[8];
  const T2j = T1f * T2h + T1g * T2i;
  const T41 = T1f * T2i - T1g * T2h;
  const T2w = cr[18], T2y = ci[18];
  const T2z = T2v * T2w + T2x * T2y;
  const T45 = T2v * T2y - T2x * T2w;

  const T2l = cr[23], T2n = ci[23];
  const T2o = T2k * T2l + T2m * T2n;
  const T42 = T2k * T2n - T2m * T2l;
  const T2r = cr[13], T2t = ci[13];
  const T2u = T2q * T2r + T2s * T2t;
  const T44 = T2q * T2t - T2s * T2r;

  const T3Y = T2j - T2o, T3Z = T2u - T2z, T4h = T44 - T45, T4g = T41 - T42;
  const T43 = T41 + T42, T46 = T44 + T45, T49 = T43 + T46;
  const T2p = T2j + T2o, T2A = T2u + T2z, T2B = T2p + T2A;

  const T2C = T2g + T2B, T6k = T48 + T49;

  const T40 = KP951056516 * T3Y + KP587785252 * T3Z;
  const T5r = KP951056516 * T3Z - KP587785252 * T3Y;
  const T47 = KP559016994 * (T43 - T46), T4a = T48 - KP250000000 * T49;
  const T4b = T47 + T4a, T5q = T4a - T47;
  const T4c = T40 + T4b, T5X = T5r + T5q;
  const T4L = T4b - T40, T5s = T5q - T5r;

  const T4i = KP951056516 * T4g + KP587785252 * T4h;
  const T5u = KP951056516 * T4h - KP587785252 * T4g;
  const T4d = KP559016994 * (T2p - T2A), T4e = T2g - KP250000000 * T2B;
  const T4f = T4d + T4e, T5t = T4e - T4d;
  const T4j = T4f - T4i, T5W = T5t - T5u;
  const T4K = T4f + T4i, T5v = T5t + T5u;

  const TM = cr[1], TN = ci[1];
  const TO = T2 * TM + T5 * TN;
  const T37 = T2 * TN - T5 * TM;

  const TS = cr[6], TW = ci[6];
  const TX = TR * TS + TV * TW;
  const T2T = TR * TW - TV * TS;
  const T1i = cr[16], T1k = ci[16];
  const T1l = T1h * T1i + T1j * T1k;
  const T2X = T1h * T1k - T1j * T1i;

  const TZ = cr[21], T11 = ci[21];
  const T12 = TY * TZ + T10 * T11;
  const T2U = TY * T11 - T10 * TZ;
  const T19 = cr[11], T1d = ci[11];
  const T1e = T18 * T19 + T1c * T1d;
  const T2W = T18 * T1d - T1c * T19;

  const T2V = T2T - T2U, T2Y = T2W - T2X, T32 = T1e - T1l, T31 = TX - T12;
  const T34 = T2T + T2U, T35 = T2W + T2X, T38 = T34 + T35;
  const T13 = TX + T12, T1m = T1e + T1l, T1n = T13 + T1m;

  const T1o = TO + T1n, T6g = T37 + T38;

  const T2Z = KP951056516 * T2V + KP587785252 * T2Y;
  const T55 = KP951056516 * T2Y - KP587785252 * T2V;
  const T2Q = KP559016994 * (T13 - T1m), T2R = TO - KP250000000 * T1n;
  const T2S = T2Q + T2R, T54 = T2R - T2Q;
  const T30 = T2S - T2Z, T5M = T54 - T55;
  const T4A = T2S + T2Z, T56 = T54 + T55;

  const T33 = KP951056516 * T31 + KP587785252 * T32;
  const T58 = KP951056516 * T32 - KP587785252 * T31;
  const T36 = KP559016994 * (T34 - T35), T39 = T37 - KP250000000 * T38;
  const T3a = T36 + T39, T57 = T39 - T36;
  const T3b = T33 + T3a, T5N = T58 + T57;
  const T4B = T3a - T33, T59 = T57 - T58;

  const T1p = cr[4], T1q = ci[4];
  const T1r = T8 * T1p + Td * T1q;
  const T3n = T8 * T1q - Td * T1p;

  const T1s = cr[9], T1t = ci[9];
  const T1u = T9 * T1s + Te * T1t;
  const T3g = T9 * T1t - Te * T1s;
  const T1F = cr[19], T1H = ci[19];
  const T1I = T1E * T1F + T1G * T1H;
  const T3k = T1E * T1H - T1G * T1F;

  const T1v = cr[24], T1w = ci[24];
  const T1x = Tn * T1v + To * T1w;
  const T3h = Tn * T1w - To * T1v;
  const T1A = cr[14], T1C = ci[14];
  const T1D = T1z * T1A + T1B * T1C;
  const T3j = T1z * T1C - T1B * T1A;

  const T3d = T1x - T1u, T3e = T1D - T1I, T3w = T3j - T3k, T3v = T3g - T3h;
  const T3i = T3g + T3h, T3l = T3j + T3k, T3o = T3i + T3l;
  const T1y = T1u + T1x, T1J = T1D + T1I, T1K = T1y + T1J;

  const T1L = T1r + T1K, T6h = T3n + T3o;

  const T3f = KP951056516 * T3d - KP587785252 * T3e;
  const T5c = KP951056516 * T3e + KP587785252 * T3d;
  const T3m = KP559016994 * (T3i - T3l), T3p = T3n - KP250000000 * T3o;
  const T3q = T3m + T3p, T5b = T3p - T3m;
  const T3r = T3f - T3q, T5P = T5c + T5b;
  const T4E = T3f + T3q, T5d = T5b - T5c;

  const T3x = KP951056516 * T3v + KP587785252 * T3w;
  const T5f = KP951056516 * T3w - KP587785252 * T3v;
  const T3s = KP559016994 * (T1y - T1J), T3t = T1r - KP250000000 * T1K;
  const T3u = T3s + T3t, T5e = T3t - T3s;
  const T3y = T3u - T3x, T5Q = T5e - T5f;
  const T4D = T3u + T3x, T5g = T5e + T5f;

  const T1N = cr[2], T1O = ci[2];
  const T1P = T14 * T1N + T16 * T1O;
  const T3L = T14 * T1O - T16 * T1N;

  const T1R = cr[7], T1T = ci[7];
  const T1U = T1Q * T1R + T1S * T1T;
  const T3E = T1Q * T1T - T1S * T1R;
  const T27 = cr[17], T29 = ci[17];
  const T2a = T26 * T27 + T28 * T29;
  const T3I = T26 * T29 - T28 * T27;

  const T1W = cr[22], T1Y = ci[22];
  const T1Z = T1V * T1W + T1X * T1Y;
  const T3F = T1V * T1Y - T1X * T1W;
  const T22 = cr[12], T24 = ci[12];
  const T25 = T21 * T22 + T23 * T24;
  const T3H = T21 * T24 - T23 * T22;

  const T3B = T1U - T1Z, T3C = T25 - T2a, T3U = T3H - T3I, T3T = T3E - T3F;
  const T3G = T3E + T3F, T3J = T3H + T3I, T3M = T3G + T3J;
  const T20 = T1U + T1Z, T2b = T25 + T2a, T2c = T20 + T2b;

  const T2d = T1P + T2c, T6j = T3L + T3M;

  const T3D = KP951056516 * T3B + KP587785252 * T3C;
  const T5n = KP951056516 * T3C - KP587785252 * T3B;
  const T3K = KP559016994 * (T3G - T3J), T3N = T3L - KP250000000 * T3M;
  const T3O = T3K + T3N, T5m = T3N - T3K;
  const T3P = T3D + T3O, T5U = T5n + T5m;
  const T4I = T3O - T3D, T5o = T5m - T5n;

  const T3V = KP951056516 * T3T + KP587785252 * T3U;
  const T5k = KP951056516 * T3U - KP587785252 * T3T;
  const T3Q = KP559016994 * (T20 - T2b), T3R = T1P - KP250000000 * T2c;
  const T3S = T3Q + T3R, T5j = T3R - T3Q;
  const T3W = T3S - T3V, T5T = T5j - T5k;
  const T4H = T3S + T3V, T5l = T5j + T5k;

  const outCr = new Float64Array(25), outCi = new Float64Array(25);

  const T6i = T6g - T6h, T6l = T6j - T6k;
  const T6m = KP951056516 * T6i + KP587785252 * T6l;
  const T6o = KP951056516 * T6l - KP587785252 * T6i;
  const TL = T1 + TK;
  const T1M = T1o + T1L, T2D = T2d + T2C, T2E = T1M + T2D;
  const T6d = KP559016994 * (T1M - T2D), T6e = TL - KP250000000 * T2E;

  outCr[0] = TL + T2E;
  const T6n = T6e - T6d;
  outCr[10] = T6n - T6o;
  outCi[9] = T6n + T6o;
  const T6f = T6d + T6e;
  outCi[4] = T6f - T6m;
  outCr[5] = T6f + T6m;

  const T2H = T2F + T2G;
  const T2P = T2H - T2O, T4z = T2H + T2O;
  const T70 = T6A + T6z;
  const T72 = T70 - T71, T7e = T71 + T70;

  const T3c = KP535826794 * T30 + KP844327925 * T3b;
  const T3z = KP770513242 * T3r - KP637423989 * T3y;
  const T3A = T3c + T3z;
  const T3X = KP904827052 * T3P - KP425779291 * T3W;
  const T4k = KP125333233 * T4c - KP992114701 * T4j;
  const T4l = T3X + T4k;
  const T4m = T3A + T4l, T7j = T3X - T4k;
  const T4n = KP559016994 * (T3A - T4l), T7i = T3z - T3c;

  const T4S = KP968583161 * T4B - KP248689887 * T4A;
  const T4T = KP535826794 * T4E - KP844327925 * T4D;
  const T73 = T4S + T4T;
  const T4V = KP876306680 * T4I - KP481753674 * T4H;
  const T4W = KP728968627 * T4L - KP684547105 * T4K;
  const T74 = T4V + T4W;
  const T4U = T4S - T4T, T77 = KP559016994 * (T73 - T74);
  const T4X = T4V - T4W, T75 = T73 + T74;

  const T4C = KP968583161 * T4A + KP248689887 * T4B;
  const T4F = KP535826794 * T4D + KP844327925 * T4E;
  const T4G = T4C + T4F;
  const T4J = KP876306680 * T4H + KP481753674 * T4I;
  const T4M = KP728968627 * T4K + KP684547105 * T4L;
  const T4N = T4J + T4M;
  const T4O = T4G + T4N, T6Y = T4J - T4M;
  const T4P = KP559016994 * (T4G - T4N), T6X = T4F - T4C;

  const T4q = KP535826794 * T3b - KP844327925 * T30;
  const T4r = KP770513242 * T3y + KP637423989 * T3r;
  const T7b = T4q + T4r;
  const T4t = KP125333233 * T4j + KP992114701 * T4c;
  const T4u = KP904827052 * T3W + KP425779291 * T3P;
  const T7c = T4u + T4t;
  const T4s = T4q - T4r, T7f = T7b - T7c;
  const T4v = T4t - T4u, T7d = KP559016994 * (T7b + T7c);

  outCr[4] = T2P + T4m;
  outCi[23] = T75 + T72;
  outCi[20] = T7f + T7e;
  outCr[1] = T4z + T4O;

  const T4w = KP951056516 * T4s + KP587785252 * T4v;
  const T4y = KP951056516 * T4v - KP587785252 * T4s;
  const T4o = T2P - KP250000000 * T4m;
  const T4p = T4n + T4o, T4x = T4o - T4n;
  outCi[0] = T4p - T4w;
  outCi[5] = T4x + T4y;
  outCr[9] = T4p + T4w;
  outCi[10] = T4x - T4y;

  const T6Z = KP587785252 * T6X + KP951056516 * T6Y;
  const T79 = KP951056516 * T6X - KP587785252 * T6Y;
  const T76 = T72 - KP250000000 * T75;
  const T78 = T76 - T77, T7a = T77 + T76;
  outCr[16] = T6Z - T78;
  outCi[18] = T79 + T7a;
  outCi[13] = T6Z + T78;
  outCr[21] = T79 - T7a;

  const T7k = KP587785252 * T7i + KP951056516 * T7j;
  const T7l = KP951056516 * T7i - KP587785252 * T7j;
  const T7g = T7e - KP250000000 * T7f;
  const T7h = T7d - T7g, T7m = T7d + T7g;
  outCr[14] = T7h - T7k;
  outCi[15] = T7l + T7m;
  outCr[19] = T7k + T7h;
  outCr[24] = T7l - T7m;

  const T4Y = KP951056516 * T4U + KP587785252 * T4X;
  const T50 = KP951056516 * T4X - KP587785252 * T4U;
  const T4Q = T4z - KP250000000 * T4O;
  const T4R = T4P + T4Q, T4Z = T4Q - T4P;
  outCi[3] = T4R - T4Y;
  outCi[8] = T4Z + T50;
  outCr[6] = T4R + T4Y;
  outCr[11] = T4Z - T50;

  const T7n = T1L - T1o, T7o = T2d - T2C;
  const T7p = KP587785252 * T7n + KP951056516 * T7o;
  const T7x = KP951056516 * T7n - KP587785252 * T7o;
  const T7q = T6y + T6v;
  const T7r = T6g + T6h, T7s = T6j + T6k, T7t = T7r + T7s;
  const T7u = T7q - KP250000000 * T7t, T7v = KP559016994 * (T7r - T7s);

  outCi[24] = T7t + T7q;
  const T7y = T7v + T7u;
  outCr[20] = T7x - T7y;
  outCi[19] = T7x + T7y;
  const T7w = T7u - T7v;
  outCr[15] = T7p - T7w;
  outCi[14] = T7p + T7w;

  const T51 = T2G - T2F;
  const T53 = T51 + T52, T5L = T51 - T52;
  const T6B = T6z - T6A;
  const T6C = T6u + T6B, T6O = T6B - T6u;

  const T5a = KP728968627 * T56 + KP684547105 * T59;
  const T5h = KP125333233 * T5d - KP992114701 * T5g;
  const T5i = T5a + T5h;
  const T5p = KP062790519 * T5l + KP998026728 * T5o;
  const T5w = KP770513242 * T5s - KP637423989 * T5v;
  const T5x = T5p + T5w;
  const T5y = T5i + T5x, T6T = T5p - T5w;
  const T5z = KP559016994 * (T5i - T5x), T6S = T5h - T5a;

  const T64 = KP876306680 * T5N - KP481753674 * T5M;
  const T65 = KP904827052 * T5Q + KP425779291 * T5P;
  const T6D = T64 - T65;
  const T67 = KP535826794 * T5U - KP844327925 * T5T;
  const T68 = KP062790519 * T5X - KP998026728 * T5W;
  const T6E = T67 + T68;
  const T66 = T64 + T65, T6H = KP559016994 * (T6D - T6E);
  const T69 = T67 - T68, T6F = T6D + T6E;

  const T5O = KP876306680 * T5M + KP481753674 * T5N;
  const T5R = KP904827052 * T5P - KP425779291 * T5Q;
  const T5S = T5O + T5R;
  const T5V = KP535826794 * T5T + KP844327925 * T5U;
  const T5Y = KP062790519 * T5W + KP998026728 * T5X;
  const T5Z = T5V + T5Y;
  const T60 = T5S + T5Z, T6q = T5V - T5Y;
  const T61 = KP559016994 * (T5S - T5Z), T6p = T5R - T5O;

  const T5C = KP728968627 * T59 - KP684547105 * T56;
  const T5D = KP125333233 * T5g + KP992114701 * T5d;
  const T6L = T5C - T5D;
  const T5F = KP062790519 * T5o - KP998026728 * T5l;
  const T5G = KP770513242 * T5v + KP637423989 * T5s;
  const T6M = T5F - T5G;
  const T5E = T5C + T5D, T6P = T6L + T6M;
  const T5H = T5F + T5G, T6N = KP559016994 * (T6L - T6M);

  outCr[3] = T53 + T5y;
  outCi[22] = T6F + T6C;
  outCi[21] = T6P + T6O;
  outCr[2] = T5L + T60;

  const T6r = KP587785252 * T6p + KP951056516 * T6q;
  const T6J = KP951056516 * T6p - KP587785252 * T6q;
  const T6G = T6C - KP250000000 * T6F;
  const T6I = T6G - T6H, T6K = T6H + T6G;
  outCr[17] = T6r - T6I;
  outCi[17] = T6J + T6K;
  outCi[12] = T6r + T6I;
  outCr[22] = T6J - T6K;

  const T6a = KP951056516 * T66 + KP587785252 * T69;
  const T6c = KP951056516 * T69 - KP587785252 * T66;
  const T62 = T5L - KP250000000 * T60;
  const T63 = T61 + T62, T6b = T62 - T61;
  outCi[2] = T63 - T6a;
  outCi[7] = T6b + T6c;
  outCr[7] = T63 + T6a;
  outCr[12] = T6b - T6c;

  const T5I = KP951056516 * T5E + KP587785252 * T5H;
  const T5K = KP951056516 * T5H - KP587785252 * T5E;
  const T5A = T53 - KP250000000 * T5y;
  const T5B = T5z + T5A, T5J = T5A - T5z;
  outCi[1] = T5B - T5I;
  outCi[6] = T5J + T5K;
  outCr[8] = T5B + T5I;
  outCi[11] = T5J - T5K;

  const T6U = KP587785252 * T6S + KP951056516 * T6T;
  const T6V = KP951056516 * T6S - KP587785252 * T6T;
  const T6Q = T6O - KP250000000 * T6P;
  const T6R = T6N - T6Q, T6W = T6N + T6Q;
  outCr[13] = T6R - T6U;
  outCi[16] = T6V + T6W;
  outCr[18] = T6U + T6R;
  outCr[23] = T6V - T6W;

  return [outCr, outCi];
}

module.exports = { hf2_25 };
