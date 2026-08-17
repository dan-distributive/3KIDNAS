'use strict';

// =============================================================================
// r2cb_25.js -- faithful JS port of rdft/scalar/r2cb/r2cb_25.c (non-FMA
// branch). O[0..24] packed halfcomplex (O[0..12]=Re0..Re12, O[13..24]=
// Im12..Im1, same convention as r2cf_25.js) -> x[0..24] real.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP425779291 = 0.425779291565072648862502445744251703979973042;
const KP904827052 = 0.904827052466019527713668647932697593970413911;
const KP535826794 = 0.535826794978996618271308767867639978063575346;
const KP844327925 = 0.844327925502015078548558063966681505381659241;
const KP876306680 = 0.876306680043863587308115903922062583399064238;
const KP481753674 = 0.481753674101715274987191502872129653528542010;
const KP968583161 = 0.968583161128631119490168375464735813836012403;
const KP248689887 = 0.248689887164854788242283746006447968417567406;
const KP062790519 = 0.062790519529313376076178224565631133122484832;
const KP998026728 = 0.998026728428271561952336806863450553336905220;
const KP728968627 = 0.728968627421411523146730319055259111372571664;
const KP684547105 = 0.684547105928688673732283357621209269889519233;
const KP250000000 = 0.25;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP500000000 = 0.5;
const KP2_000000000 = 2.0;
const KP1_118033988 = 1.118033988749894848204586834365638117720309180;
const KP1_175570504 = 1.175570504584946258337411909278145537195304875;
const KP1_902113032 = 1.902113032590307144232878666758764286811397268;

function r2cb_25(O) {
  const Ts = O[20], Tt = O[15];
  const Tu = KP1_902113032 * Ts + KP1_175570504 * Tt;
  const T1G = KP1_175570504 * Ts - KP1_902113032 * Tt;

  const T1 = O[0], T2 = O[5], T3 = O[10];
  const T4 = T2 + T3;
  const Tp = KP1_118033988 * (T2 - T3);
  const T5 = KP2_000000000 * T4 + T1;
  const Tq = T1 - KP500000000 * T4;
  const Tr = Tp + Tq, T1F = Tq - Tp;

  const T6 = O[1], TN = O[24];

  const T7 = O[6], T8 = O[4];
  const T9 = T7 + T8;
  const Ta = O[11], Tb = O[9];
  const Tc = Ta + Tb;
  const Td = T9 + Tc;
  const TI = Ta - Tb;
  const Tw = KP559016994 * (T9 - Tc);
  const TH = T7 - T8;

  const Tz = O[19], TA = O[21];
  const TK = Tz - TA;
  const TC = O[14], TD = O[16];
  const TL = TC - TD;
  const TB = Tz + TA;
  const TO = TK + TL;
  const TE = TC + TD;
  const TM = KP559016994 * (TK - TL);

  const Te = T6 + Td;

  const TJ = KP951056516 * TH + KP587785252 * TI;
  const T1L = KP587785252 * TH - KP951056516 * TI;
  const TP = TN - KP250000000 * TO;
  const TQ = TM + TP;
  const T1M = TP - TM;
  const TR = TJ + TQ;
  const T27 = T1M - T1L;
  const T1r = TQ - TJ;
  const T1N = T1L + T1M;

  const TF = KP951056516 * TB + KP587785252 * TE;
  const T1J = KP587785252 * TB - KP951056516 * TE;
  const Tx = T6 - KP250000000 * Td;
  const Ty = Tw + Tx;
  const T1I = Tx - Tw;
  const TG = Ty - TF;
  const T26 = T1I + T1J;
  const T1q = Ty + TF;
  const T1K = T1I - T1J;

  const Tf = O[2], T1a = O[23];

  const Tg = O[7], Th = O[3];
  const Ti = Tg + Th;
  const Tj = O[12], Tk = O[8];
  const Tl = Tj + Tk;
  const Tm = Ti + Tl;
  const T15 = Tj - Tk;
  const TT = KP559016994 * (Ti - Tl);
  const T14 = Tg - Th;

  const TW = O[18], TX = O[22];
  const T17 = TW - TX;
  const TZ = O[13], T10 = O[17];
  const T18 = TZ - T10;
  const TY = TW + TX;
  const T1b = T17 + T18;
  const T11 = TZ + T10;
  const T19 = KP559016994 * (T17 - T18);

  const Tn = Tf + Tm;

  const T16 = KP951056516 * T14 + KP587785252 * T15;
  const T1S = KP587785252 * T14 - KP951056516 * T15;
  const T1c = T1a - KP250000000 * T1b;
  const T1d = T19 + T1c;
  const T1T = T1c - T19;
  const T1e = T16 + T1d;
  const T2a = T1T - T1S;
  const T1u = T1d - T16;
  const T1U = T1S + T1T;

  const T12 = KP951056516 * TY + KP587785252 * T11;
  const T1Q = KP587785252 * TY - KP951056516 * T11;
  const TU = Tf - KP250000000 * Tm;
  const TV = TT + TU;
  const T1P = TU - TT;
  const T13 = TV - T12;
  const T29 = T1P + T1Q;
  const T1t = TV + T12;
  const T1R = T1P - T1Q;

  const x = new Float64Array(25);

  const T2m = KP1_118033988 * (Te - Tn);
  const To = Te + Tn;
  const T2l = T5 - KP500000000 * To;
  const T2o = TO + TN;
  const T2p = T1b + T1a;
  const T2q = KP1_175570504 * T2o - KP1_902113032 * T2p;
  const T2s = KP1_902113032 * T2o + KP1_175570504 * T2p;
  x[0] = KP2_000000000 * To + T5;
  const T2r = T2m + T2l;
  x[5] = T2r - T2s;
  x[20] = T2r + T2s;
  const T2n = T2l - T2m;
  x[10] = T2n - T2q;
  x[15] = T2n + T2q;

  const T2g = KP684547105 * T26 + KP728968627 * T27;
  const T2h = KP998026728 * T29 + KP062790519 * T2a;
  const T2i = KP1_175570504 * T2g - KP1_902113032 * T2h;
  const T2k = KP1_902113032 * T2g + KP1_175570504 * T2h;
  const T25 = T1F + T1G;
  const T28 = KP728968627 * T26 - KP684547105 * T27;
  const T2b = KP062790519 * T29 - KP998026728 * T2a;
  const T2c = T28 + T2b;
  const T2d = T25 - KP500000000 * T2c;
  const T2e = KP1_118033988 * (T28 - T2b);

  x[3] = KP2_000000000 * T2c + T25;
  const T2j = T2e + T2d;
  x[8] = T2j - T2k;
  x[23] = T2j + T2k;
  const T2f = T2d - T2e;
  x[13] = T2f - T2i;
  x[18] = T2f + T2i;

  const T1k = KP248689887 * TG + KP968583161 * TR;
  const T1l = KP481753674 * T13 + KP876306680 * T1e;
  const T1m = KP1_175570504 * T1k - KP1_902113032 * T1l;
  const T1o = KP1_902113032 * T1k + KP1_175570504 * T1l;
  const Tv = Tr - Tu;
  const TS = KP968583161 * TG - KP248689887 * TR;
  const T1f = KP876306680 * T13 - KP481753674 * T1e;
  const T1g = TS + T1f;
  const T1h = Tv - KP500000000 * T1g;
  const T1i = KP1_118033988 * (TS - T1f);

  x[1] = KP2_000000000 * T1g + Tv;
  const T1n = T1i + T1h;
  x[6] = T1n - T1o;
  x[21] = T1n + T1o;
  const T1j = T1h - T1i;
  x[11] = T1j - T1m;
  x[16] = T1j + T1m;

  const T1A = KP844327925 * T1q + KP535826794 * T1r;
  const T1B = KP904827052 * T1t - KP425779291 * T1u;
  const T1C = KP1_175570504 * T1A - KP1_902113032 * T1B;
  const T1E = KP1_902113032 * T1A + KP1_175570504 * T1B;
  const T1p = Tr + Tu;
  const T1s = KP535826794 * T1q - KP844327925 * T1r;
  const T1v = KP425779291 * T1t + KP904827052 * T1u;
  const T1w = T1s - T1v;
  const T1x = T1p - KP500000000 * T1w;
  const T1y = KP1_118033988 * (T1s + T1v);

  x[4] = KP2_000000000 * T1w + T1p;
  const T1D = T1x + T1y;
  x[9] = T1D - T1E;
  x[24] = T1E + T1D;
  const T1z = T1x - T1y;
  x[14] = T1z - T1C;
  x[19] = T1C + T1z;

  const T20 = KP481753674 * T1K + KP876306680 * T1N;
  const T21 = KP844327925 * T1R + KP535826794 * T1U;
  const T22 = KP1_175570504 * T20 - KP1_902113032 * T21;
  const T24 = KP1_902113032 * T20 + KP1_175570504 * T21;
  const T1H = T1F - T1G;
  const T1O = KP876306680 * T1K - KP481753674 * T1N;
  const T1V = KP535826794 * T1R - KP844327925 * T1U;
  const T1W = T1O + T1V;
  const T1X = T1H - KP500000000 * T1W;
  const T1Y = KP1_118033988 * (T1O - T1V);

  x[2] = KP2_000000000 * T1W + T1H;
  const T23 = T1Y + T1X;
  x[7] = T23 - T24;
  x[22] = T23 + T24;
  const T1Z = T1X - T1Y;
  x[12] = T1Z - T22;
  x[17] = T1Z + T22;

  return x;
}

module.exports = { r2cb_25 };
