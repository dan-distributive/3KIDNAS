'use strict';

// =============================================================================
// hb2_25.js -- JS port of rdft/scalar/r2cb/hb2_25.c (non-FMA branch).
// Alternate-codegen ("twiddle-log3/precompute-twiddles") sibling of the
// (unported, dormant) hb_25.js -- same math, different rounding (see
// hb2_5.js's header for the pattern, and hf2_25.js's header for the
// forward-direction analogue). twinstr only trig-generates W^1, W^3, W^9,
// W^24 (four raw pairs, same generator set as hf2_25.js/t2_25.js); every
// other needed multiple is DERIVED via complex products of those four.
// Same (cr,ci,Wc,Ws) -> [outCr,outCi] convention as hb2_5.js -- only
// indices 1, 3, 9, 24 are actually read from Wc/Ws.
//
// **NOT VERIFIED CORRECT -- KNOWN BUGGY, NOT REGISTERED.** An automated
// semantic diff against hb2_25.c (every assignment evaluated with random
// inputs, C vs JS, in source order) shows the dependency graph below is a
// byte-for-byte faithful transcription of the C source -- ruling out a
// stray macro/sign/index typo. Despite that, composing this codelet with
// the independently-proven-correct hf2_25.js (hb2_25(hf2_25(x)) should
// equal exactly 25*x, the same exact property confirmed for the sibling
// hf2_5/hb2_5 pair) fails at exactly the phases p where p%5 is 0 or 3
// (excluding p=0) -- a large (20-300%), reproducible, NOT rounding-scale
// error. Ground-truth-confirmed wrong too (real N=625 HC2R differs from
// this engine's output by ~1300 when this codelet is used). Root cause
// NOT isolated this session -- see RealEngine1D.js's hb2 registry comment
// (this codelet is deliberately left OUT of that registry, unlike the
// registered-but-excluded hf2_25) for the full investigation writeup and
// candidate next steps. Kept in the tree for a future debugging pass, not
// deleted -- do not register or call this function until root-caused.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP998026728 = 0.998026728428271561952336806863450553336905220;
const KP062790519 = 0.062790519529313376076178224565631133122484832;
const KP992114701 = 0.992114701314477831049793042785778521453036709;
const KP125333233 = 0.125333233564304245373118759816508793942918247;
const KP425779291 = 0.425779291565072648862502445744251703979973042;
const KP904827052 = 0.904827052466019527713668647932697593970413911;
const KP248689887 = 0.248689887164854788242283746006447968417567406;
const KP968583161 = 0.968583161128631119490168375464735813836012403;
const KP770513242 = 0.770513242775789230803009636396177847271667672;
const KP637423989 = 0.637423989748689710176712811676016195434917298;
const KP844327925 = 0.844327925502015078548558063966681505381659241;
const KP535826794 = 0.535826794978996618271308767867639978063575346;
const KP684547105 = 0.684547105928688673732283357621209269889519233;
const KP728968627 = 0.728968627421411523146730319055259111372571664;
const KP481753674 = 0.481753674101715274987191502872129653528542010;
const KP876306680 = 0.876306680043863587308115903922062583399064238;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP250000000 = 0.25;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;

function hb2_25(cr, ci, Wc, Ws) {
  const TN = Wc[1], TQ = Ws[1], TO = Wc[3], TR = Ws[3];
  const TP = TN * TO, TX = TQ * TO, TS = TQ * TR, TW = TN * TR;
  const TT = TP - TS, TY = TW + TX, T2t = TW - TX, T2r = TP + TS;
  const TZ = Ws[9];
  const T2c = TQ * TZ, T4j = TO * TZ, T2h = TN * TZ, T4e = TR * TZ;
  const TU = Wc[9];
  const T2b = TN * TU, T4k = TR * TU, T2i = TQ * TU, T4d = TO * TU;

  const T4f = T4d - T4e, T4l = T4j + T4k;

  const T2d = T2b - T2c, T4v = T2b + T2c, T5m = T4j - T4k;
  const T2j = T2h + T2i, T5l = T4d + T4e, T4X = T2h - T2i;
  const T2s = T2r * TU, T2u = T2t * TZ, T2v = T2s + T2u;
  const TV = TT * TU, T10 = TY * TZ, T11 = TV + T10;
  const T3P = T2r * TZ, T3Q = T2t * TU, T3R = T3P - T3Q;
  const T1J = TT * TZ, T1K = TY * TU, T1L = T1J - T1K;
  const T5d = TV - T10, T6x = T3P + T3Q, T5h = T1J + T1K, T6t = T2s - T2u;

  const T25 = Wc[24], T26 = Ws[24];
  const T27 = TT * T25 + TY * T26;
  const T29 = TT * T26 - TY * T25;
  const T6D = T4v * T26 - T4X * T25;
  const T7v = T11 * T26 - T1L * T25;
  const T49 = T2r * T25 + T2t * T26;
  const T7l = T2d * T25 + T2j * T26;
  const T7p = T2d * T26 - T2j * T25;
  const T7t = T11 * T25 + T1L * T26;
  const T2p = TU * T26 - TZ * T25;
  const T2n = TU * T25 + TZ * T26;
  const T4b = T2r * T26 - T2t * T25;
  const T4p = T2v * T25 + T3R * T26;
  const T5n = T5l * T25 + T5m * T26;
  const T6B = T4v * T25 + T4X * T26;
  const T5b = TN * T26 - TQ * T25;
  const T5p = T5l * T26 - T5m * T25;
  const T6p = TO * T25 + TR * T26;
  const T6r = TO * T26 - TR * T25;
  const T59 = TN * T25 + TQ * T26;
  const T4r = T2v * T26 - T3R * T25;

  const T1 = cr[0];
  const T2 = cr[5], T3 = ci[4];
  const T4 = T2 + T3;
  const T5 = cr[10], T6 = ci[9];
  const T7 = T5 + T6;
  const T8 = T4 + T7;
  const T3Z = T5 - T6;
  const T3Y = T2 - T3;

  const T9 = T1 + T8;
  const T6i = KP951056516 * T3Y + KP587785252 * T3Z;
  const T40 = KP587785252 * T3Y - KP951056516 * T3Z;
  const T3x = T1 - KP250000000 * T8;
  const T3y = KP559016994 * (T4 - T7);
  const T3z = T3x - T3y;
  const T5Y = T3y + T3x;

  const Ta = cr[1];
  const Te = cr[11], Tf = ci[8];
  const Tg = Te + Tf;
  const T2E = Te - Tf;
  const Tb = cr[6], Tc = ci[3];
  const Td = Tb + Tc;
  const T2D = Tb - Tc;

  const T2x = KP559016994 * (Td - Tg);
  const T5w = KP951056516 * T2D + KP587785252 * T2E;
  const T2F = KP587785252 * T2D - KP951056516 * T2E;
  const Th = Td + Tg;
  const T2w = Ta - KP250000000 * Th;

  const T1e = ci[20];
  const T1i = cr[14], T1j = cr[19];
  const T1k = T1i + T1j;
  const T2W = T1j - T1i;
  const T1f = ci[15], T1g = cr[24];
  const T1h = T1f - T1g;
  const T2V = T1f + T1g;

  const T2P = KP559016994 * (T1h + T1k);
  const T5B = KP951056516 * T2V + KP587785252 * T2W;
  const T2X = KP587785252 * T2V - KP951056516 * T2W;
  const T1l = T1h - T1k;
  const T2O = T1e - KP250000000 * T1l;

  const Tj = cr[4];
  const Tn = ci[10], To = ci[5];
  const Tp = Tn + To;
  const T2M = Tn - To;
  const Tk = cr[9], Tl = ci[0];
  const Tm = Tk + Tl;
  const T2L = Tk - Tl;

  const T2N = KP587785252 * T2L - KP951056516 * T2M;
  const T5D = KP951056516 * T2L + KP587785252 * T2M;
  const T2T = KP559016994 * (Tm - Tp);
  const Tq = Tm + Tp;
  const T2S = Tj - KP250000000 * Tq;

  const T15 = ci[23];
  const T19 = ci[13], T1a = cr[16];
  const T1b = T19 - T1a;
  const T2A = T19 + T1a;
  const T16 = ci[18], T17 = cr[21];
  const T18 = T16 - T17;
  const T2z = T16 + T17;

  const T2B = KP587785252 * T2z - KP951056516 * T2A;
  const T5u = KP951056516 * T2z + KP587785252 * T2A;
  const T2H = KP559016994 * (T18 - T1b);
  const T1c = T18 + T1b;
  const T2G = T15 - KP250000000 * T1c;

  const Ti = Ta + Th;
  const Tr = Tj + Tq;
  const Ts = Ti + Tr;
  const T1d = T15 + T1c;
  const T1m = T1e + T1l;
  const T1P = T1d + T1m;

  const T2y = T2w - T2x;
  const T2C = T2y - T2B;
  const T4w = T2y + T2B;
  const T2I = T2G - T2H;
  const T2J = T2F + T2I;
  const T4x = T2I - T2F;
  const T2K = KP876306680 * T2C - KP481753674 * T2J;
  const T4P = KP728968627 * T4w + KP684547105 * T4x;
  const T3H = KP876306680 * T2J + KP481753674 * T2C;
  const T4y = KP728968627 * T4w - KP684547105 * T4x;

  const T5A = T2T + T2S;
  const T5C = T5A - T5B;
  const T6M = T5A + T5B;
  const T5E = T2O + T2P;
  const T5F = T5D + T5E;
  const T6L = T5E - T5D;
  const T5G = KP535826794 * T5C - KP844327925 * T5F;
  const T71 = KP770513242 * T6M + KP637423989 * T6L;
  const T65 = KP535826794 * T5F + KP844327925 * T5C;
  const T6N = KP770513242 * T6L - KP637423989 * T6M;

  const T5t = T2x + T2w;
  const T5v = T5t - T5u;
  const T6I = T5t + T5u;
  const T5x = T2H + T2G;
  const T5y = T5w + T5x;
  const T6J = T5x - T5w;
  const T5z = KP968583161 * T5v - KP248689887 * T5y;
  const T70 = KP844327925 * T6I + KP535826794 * T6J;
  const T64 = KP968583161 * T5y + KP248689887 * T5v;
  const T6K = KP535826794 * T6I - KP844327925 * T6J;

  const T2Q = T2O - T2P;
  const T2R = T2N + T2Q;
  const T4z = T2Q - T2N;
  const T2U = T2S - T2T;
  const T2Y = T2U - T2X;
  const T4A = T2U + T2X;
  const T2Z = KP904827052 * T2R + KP425779291 * T2Y;
  const T4Q = KP125333233 * T4z - KP992114701 * T4A;
  const T3I = KP904827052 * T2Y - KP425779291 * T2R;
  const T4B = KP125333233 * T4A + KP992114701 * T4z;

  const T1S = ci[24];
  const T1T = ci[19], T1U = cr[20];
  const T1V = T1T - T1U;
  const T1W = ci[14], T1X = cr[15];
  const T1Y = T1W - T1X;
  const T1Z = T1V + T1Y;
  const T3B = T1W + T1X;
  const T3A = T1T + T1U;

  const T20 = T1S + T1Z;
  const T5Z = KP951056516 * T3A + KP587785252 * T3B;
  const T3C = KP587785252 * T3A - KP951056516 * T3B;
  const T41 = T1S - KP250000000 * T1Z;
  const T42 = KP559016994 * (T1V - T1Y);
  const T43 = T41 - T42;
  const T6j = T42 + T41;

  const Tt = cr[2];
  const Tu = cr[7], Tv = ci[2];
  const Tw = Tu + Tv;
  const T38 = Tu - Tv;
  const Tx = cr[12], Ty = ci[7];
  const Tz = Tx + Ty;
  const T39 = Tx - Ty;

  const T32 = KP559016994 * (Tw - Tz);
  const T5L = KP951056516 * T38 + KP587785252 * T39;
  const T3a = KP587785252 * T38 - KP951056516 * T39;
  const TA = Tw + Tz;
  const T31 = Tt - KP250000000 * TA;

  const T1o = ci[22];
  const T1p = ci[17], T1q = cr[22];
  const T1r = T1p - T1q;
  const T34 = T1p + T1q;
  const T1s = ci[12], T1t = cr[17];
  const T1u = T1s - T1t;
  const T35 = T1s + T1t;

  const T36 = KP587785252 * T34 - KP951056516 * T35;
  const T5J = KP951056516 * T34 + KP587785252 * T35;
  const T3c = KP559016994 * (T1r - T1u);
  const T1v = T1r + T1u;
  const T3b = T1o - KP250000000 * T1v;

  const TC = cr[3];
  const TG = ci[11], TH = ci[6];
  const TI = TG + TH;
  const T3o = TG - TH;
  const TD = cr[8], TE = ci[1];
  const TF = TD + TE;
  const T3n = TD - TE;

  const T3h = KP559016994 * (TF - TI);
  const T5S = KP951056516 * T3n + KP587785252 * T3o;
  const T3p = KP587785252 * T3n - KP951056516 * T3o;
  const TJ = TF + TI;
  const T3g = TC - KP250000000 * TJ;

  const T1x = ci[21];
  const T1B = cr[13], T1C = cr[18];
  const T1D = T1B + T1C;
  const T3k = T1C - T1B;
  const T1y = ci[16], T1z = cr[23];
  const T1A = T1y - T1z;
  const T3j = T1y + T1z;

  const T3l = KP587785252 * T3j - KP951056516 * T3k;
  const T5Q = KP951056516 * T3j + KP587785252 * T3k;
  const T3r = KP559016994 * (T1A + T1D);
  const T1E = T1A - T1D;
  const T3q = T1x - KP250000000 * T1E;

  const TB = Tt + TA;
  const TK = TC + TJ;
  const TL = TB + TK;
  const T1w = T1o + T1v;
  const T1F = T1x + T1E;
  const T1Q = T1w + T1F;

  const T33 = T31 - T32;
  const T37 = T33 - T36;
  const T4D = T33 + T36;
  const T3d = T3b - T3c;
  const T3e = T3a + T3d;
  const T4E = T3d - T3a;
  const T3f = KP535826794 * T37 - KP844327925 * T3e;
  const T4S = KP998026728 * T4D + KP062790519 * T4E;
  const T3K = KP535826794 * T3e + KP844327925 * T37;
  const T4F = KP062790519 * T4D - KP998026728 * T4E;

  const T5P = T3h + T3g;
  const T5R = T5P - T5Q;
  const T6T = T5P + T5Q;
  const T5T = T3q + T3r;
  const T5U = T5S + T5T;
  const T6S = T5T - T5S;
  const T5V = KP728968627 * T5R - KP684547105 * T5U;
  const T74 = KP125333233 * T6T - KP992114701 * T6S;
  const T68 = KP728968627 * T5U + KP684547105 * T5R;
  const T6U = KP992114701 * T6T + KP125333233 * T6S;

  const T5I = T32 + T31;
  const T5K = T5I - T5J;
  const T6Q = T5I + T5J;
  const T5M = T3c + T3b;
  const T5N = T5L + T5M;
  const T6P = T5M - T5L;
  const T5O = KP876306680 * T5K - KP481753674 * T5N;
  const T73 = KP904827052 * T6Q - KP425779291 * T6P;
  const T67 = KP876306680 * T5N + KP481753674 * T5K;
  const T6R = KP425779291 * T6Q + KP904827052 * T6P;

  const T3i = T3g - T3h;
  const T3m = T3i - T3l;
  const T4H = T3i + T3l;
  const T3s = T3q - T3r;
  const T3t = T3p + T3s;
  const T4G = T3s - T3p;
  const T3u = KP062790519 * T3m - KP998026728 * T3t;
  const T4T = KP770513242 * T4H - KP637423989 * T4G;
  const T3L = KP998026728 * T3m + KP062790519 * T3t;
  const T4I = KP637423989 * T4H + KP770513242 * T4G;

  const T12 = KP559016994 * (Ts - TL);
  const TM = Ts + TL;
  const T13 = TM - KP250000000 * T9;
  const T14 = T12 + T13;
  const T2e = T13 - T12;
  const T1R = KP559016994 * (T1P - T1Q);
  const T21 = T1P + T1Q;
  const T22 = T20 - KP250000000 * T21;
  const T23 = T1R + T22;
  const T2l = T22 - T1R;

  const T1n = T1d - T1m;
  const T1G = T1w - T1F;
  const T1H = KP951056516 * T1n + KP587785252 * T1G;
  const T2f = KP587785252 * T1n - KP951056516 * T1G;
  const T1M = Ti - Tr;
  const T1N = TB - TK;
  const T1O = KP951056516 * T1M + KP587785252 * T1N;
  const T2k = KP587785252 * T1M - KP951056516 * T1N;

  const outCr = new Float64Array(25), outCi = new Float64Array(25);

  outCr[0] = T9 + TM;
  outCi[0] = T20 + T21;
  const T1I = T14 - T1H;
  const T24 = T1O + T23;
  outCr[5] = T11 * T1I - T1L * T24;
  outCi[5] = T1L * T1I + T11 * T24;
  const T2o = T2e + T2f;
  const T2q = T2l - T2k;
  outCr[15] = T2n * T2o - T2p * T2q;
  outCi[15] = T2p * T2o + T2n * T2q;

  const T2g = T2e - T2f;
  const T2m = T2k + T2l;
  outCr[10] = T2d * T2g - T2j * T2m;
  outCi[10] = T2j * T2g + T2d * T2m;
  const T28 = T14 + T1H;
  const T2a = T23 - T1O;
  outCr[20] = T27 * T28 - T29 * T2a;
  outCi[20] = T29 * T28 + T27 * T2a;

  const T72 = T70 + T71;
  const T75 = T73 - T74;
  const T76 = KP951056516 * T72 + KP587785252 * T75;
  const T7n = KP587785252 * T72 - KP951056516 * T75;
  const T78 = T6K - T6N;
  const T79 = T6U - T6R;
  const T7a = KP951056516 * T78 + KP587785252 * T79;
  const T7q = KP587785252 * T78 - KP951056516 * T79;

  const T6H = T5Y + T5Z;
  const T6O = T6K + T6N;
  const T6V = T6R + T6U;
  const T6W = T6O - T6V;
  const T6X = T6H - KP250000000 * T6W;
  const T6Y = KP559016994 * (T6O + T6V);
  const T7e = T6j - T6i;
  const T7b = T70 - T71;
  const T7c = T73 + T74;
  const T7f = T7b + T7c;
  const T7d = KP559016994 * (T7b - T7c);
  const T7g = T7e - KP250000000 * T7f;

  const T7x = T6H + T6W;
  const T7y = T7e + T7f;
  outCr[4] = TT * T7x - TY * T7y;
  outCi[4] = TY * T7x + TT * T7y;

  const T7m = T6X - T6Y;
  const T7o = T7m - T7n;
  const T7u = T7m + T7n;
  const T7r = T7g - T7d;
  const T7s = T7q + T7r;
  const T7w = T7r - T7q;
  outCr[14] = T7l * T7o - T7p * T7s;
  outCi[14] = T7p * T7o + T7l * T7s;
  outCr[19] = T7t * T7u - T7v * T7w;
  outCi[19] = T7v * T7u + T7t * T7w;

  const T6Z = T6X + T6Y;
  const T77 = T6Z - T76;
  const T7j = T6Z + T76;
  const T7h = T7d + T7g;
  const T7i = T7a + T7h;
  const T7k = T7h - T7a;
  outCr[9] = TU * T77 - TZ * T7i;
  outCi[9] = TZ * T77 + TU * T7i;
  outCr[24] = T25 * T7j - T26 * T7k;
  outCi[24] = T26 * T7j + T25 * T7k;

  const T3J = T3H - T3I;
  const T3M = T3K - T3L;
  const T3N = KP951056516 * T3J + KP587785252 * T3M;
  const T4h = KP587785252 * T3J - KP951056516 * T3M;
  const T3S = T2K + T2Z;
  const T3T = T3f - T3u;
  const T3U = KP951056516 * T3S + KP587785252 * T3T;
  const T4m = KP587785252 * T3S - KP951056516 * T3T;

  const T3D = T3z - T3C;
  const T30 = T2K - T2Z;
  const T3v = T3f + T3u;
  const T3E = T30 + T3v;
  const T3w = KP559016994 * (T30 - T3v);
  const T3F = T3D - KP250000000 * T3E;
  const T44 = T40 + T43;
  const T3V = T3H + T3I;
  const T3W = T3K + T3L;
  const T45 = T3V + T3W;
  const T3X = KP559016994 * (T3V - T3W);
  const T46 = T44 - KP250000000 * T45;

  const T4t = T3D + T3E;
  const T4u = T44 + T45;
  outCr[2] = T2r * T4t - T2t * T4u;
  outCi[2] = T2t * T4t + T2r * T4u;

  const T4g = T3F - T3w;
  const T4i = T4g - T4h;
  const T4q = T4g + T4h;
  const T4n = T46 - T3X;
  const T4o = T4m + T4n;
  const T4s = T4n - T4m;
  outCr[12] = T4f * T4i - T4l * T4o;
  outCi[12] = T4l * T4i + T4f * T4o;
  outCr[17] = T4p * T4q - T4r * T4s;
  outCi[17] = T4r * T4q + T4p * T4s;

  const T3G = T3w + T3F;
  const T3O = T3G - T3N;
  const T4a = T3G + T3N;
  const T47 = T3X + T46;
  const T48 = T3U + T47;
  const T4c = T47 - T3U;
  outCr[7] = T2v * T3O - T3R * T48;
  outCi[7] = T3R * T3O + T2v * T48;
  outCr[22] = T49 * T4a - T4b * T4c;
  outCi[22] = T4b * T4a + T49 * T4c;

  const T4R = T4P - T4Q;
  const T4U = T4S - T4T;
  const T4V = KP951056516 * T4R + KP587785252 * T4U;
  const T5f = KP587785252 * T4R - KP951056516 * T4U;
  const T4Y = T4y + T4B;
  const T4Z = T4F + T4I;
  const T50 = KP951056516 * T4Y + KP587785252 * T4Z;
  const T5i = KP587785252 * T4Y - KP951056516 * T4Z;

  const T4L = T3z + T3C;
  const T4C = T4y - T4B;
  const T4J = T4F - T4I;
  const T4M = T4C + T4J;
  const T4K = KP559016994 * (T4C - T4J);
  const T4N = T4L - KP250000000 * T4M;
  const T54 = T43 - T40;
  const T51 = T4P + T4Q;
  const T52 = T4S + T4T;
  const T55 = T51 + T52;
  const T53 = KP559016994 * (T51 - T52);
  const T56 = T54 - KP250000000 * T55;

  const T5r = T4L + T4M;
  const T5s = T54 + T55;
  outCr[3] = TO * T5r - TR * T5s;
  outCi[3] = TR * T5r + TO * T5s;

  const T5e = T4N - T4K;
  const T5g = T5e - T5f;
  const T5o = T5e + T5f;
  const T5j = T56 - T53;
  const T5k = T5i + T5j;
  const T5q = T5j - T5i;
  outCr[13] = T5d * T5g - T5h * T5k;
  outCi[13] = T5h * T5g + T5d * T5k;
  outCr[18] = T5n * T5o - T5p * T5q;
  outCi[18] = T5p * T5o + T5n * T5q;

  const T4O = T4K + T4N;
  const T4W = T4O - T4V;
  const T5a = T4O + T4V;
  const T57 = T53 + T56;
  const T58 = T50 + T57;
  const T5c = T57 - T50;
  outCr[8] = T4v * T4W - T4X * T58;
  outCi[8] = T4X * T4W + T4v * T58;
  outCr[23] = T59 * T5a - T5b * T5c;
  outCi[23] = T5b * T5a + T59 * T5c;

  const T66 = T64 - T65;
  const T69 = T67 - T68;
  const T6a = KP951056516 * T66 + KP587785252 * T69;
  const T6v = KP587785252 * T66 - KP951056516 * T69;
  const T6c = T5z - T5G;
  const T6d = T5O - T5V;
  const T6e = KP951056516 * T6c + KP587785252 * T6d;
  const T6y = KP587785252 * T6c - KP951056516 * T6d;

  const T60 = T5Y - T5Z;
  const T5H = T5z + T5G;
  const T5W = T5O + T5V;
  const T61 = T5H + T5W;
  const T5X = KP559016994 * (T5H - T5W);
  const T62 = T60 - KP250000000 * T61;
  const T6k = T6i + T6j;
  const T6f = T64 + T65;
  const T6g = T67 + T68;
  const T6l = T6f + T6g;
  const T6h = KP559016994 * (T6f - T6g);
  const T6m = T6k - KP250000000 * T6l;

  const T6F = T60 + T61;
  const T6G = T6k + T6l;
  outCr[1] = TN * T6F - TQ * T6G;
  outCi[1] = TQ * T6F + TN * T6G;

  const T6u = T62 - T5X;
  const T6w = T6u - T6v;
  const T6C = T6u + T6v;
  const T6z = T6m - T6h;
  const T6A = T6y + T6z;
  const T6E = T6z - T6y;
  outCr[11] = T6t * T6w - T6x * T6A;
  outCi[11] = T6x * T6w + T6t * T6A;
  outCr[16] = T6B * T6C - T6D * T6E;
  outCi[16] = T6D * T6C + T6B * T6E;

  const T63 = T5X + T62;
  const T6b = T63 - T6a;
  const T6q = T63 + T6a;
  const T6n = T6h + T6m;
  const T6o = T6e + T6n;
  const T6s = T6n - T6e;
  outCr[6] = T5l * T6b - T5m * T6o;
  outCi[6] = T5m * T6b + T5l * T6o;
  outCr[21] = T6p * T6q - T6r * T6s;
  outCi[21] = T6r * T6q + T6p * T6s;

  return [outCr, outCi];
}

module.exports = { hb2_25 };
