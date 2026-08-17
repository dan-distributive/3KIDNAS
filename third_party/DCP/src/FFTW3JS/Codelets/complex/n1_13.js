'use strict';

// =============================================================================
// n1_13.js -- faithful JS port of dft/scalar/codelets/n1_13.c (non-FMA
// branch), FFTW3's direct (base-case) radix-13 complex DFT codelet.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP2_000000000 = 2.0;
const KP083333333 = 0.083333333333333333333333333333333333333333333;
const KP251768516 = 0.251768516431883313623436926934233488546674281;
const KP075902986 = 0.075902986037193865983102897245103540356428373;
const KP132983124 = 0.132983124607418643793760531921092974399165133;
const KP258260390 = 0.258260390311744861420450644284508567852516811;
const KP1_732050807 = 1.732050807568877293527446341505872366942805254;
const KP300238635 = 0.300238635966332641462884626667381504676006424;
const KP011599105 = 0.011599105605768290721655456654083252189827041;
const KP156891391 = 0.156891391051584611046832726756003269660212636;
const KP256247671 = 0.256247671582936600958684654061725059144125175;
const KP174138601 = 0.174138601152135905005660794929264742616964676;
const KP575140729 = 0.575140729474003121368385547455453388461001608;
const KP503537032 = 0.503537032863766627246873853868466977093348562;
const KP113854479 = 0.113854479055790798974654345867655310534642560;
const KP265966249 = 0.265966249214837287587521063842185948798330267;
const KP387390585 = 0.387390585467617292130675966426762851778775217;
const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP300462606 = 0.300462606288665774426601772289207995520941381;
const KP500000000 = 0.5;

function n1_13(ri, ii) {
  const T1 = ri[0], T1q = ii[0];

  const Td = ri[8], Te = ri[5];
  const Tf = Td + Te, Tp = Td - Te;

  const T7 = ri[12], T8 = ri[10], T9 = ri[4];
  const Ta = T8 + T9;
  const Tb = T7 + Ta, TC = T8 - T9;
  const Tx = T7 - KP500000000 * Ta;

  const T2 = ri[1], T3 = ri[3], T4 = ri[9];
  const T5 = T3 + T4;
  const T6 = T2 + T5, TB = T3 - T4;
  const Tw = T2 - KP500000000 * T5;

  const Tg = ri[11], Th = ri[6];
  const Ti = Tg + Th, Tq = Tg - Th;
  const Tj = ri[7], Tk = ri[2];
  const Tl = Tj + Tk, Tr = Tj - Tk;

  const Tm = Ti + Tl;
  const Ts = Tq + Tr;
  const Tt = Tp + Ts;
  const Tu = T6 - Tb;
  const Tc = T6 + Tb;
  const Tn = Tf + Tm;
  const To = Tc + Tn;
  const T22 = KP300462606 * (Tc - Tn);

  const T1Y = TB + TC;
  const T1Z = Tq - Tr;
  const T20 = T1Y - T1Z;
  const T24 = T1Y + T1Z;
  const TD = KP866025403 * (TB - TC);
  const TE = Tp - KP500000000 * Ts;
  const TF = TD - TE;
  const TH = TD + TE;

  const Ty = Tw - Tx;
  const Tz = KP866025403 * (Ti - Tl);
  const TA = Ty + Tz;
  const TI = Ty - Tz;
  const T1V = Tw + Tx;
  const T1W = Tf - KP500000000 * Tm;
  const T1X = T1V - T1W;
  const T25 = T1V + T1W;

  const TX = ii[8], TY = ii[5];
  const TZ = TX + TY, T2b = TX - TY;

  const TR = ii[12], TS = ii[10], TT = ii[4];
  const TU = TS + TT;
  const TV = TR - KP500000000 * TU;
  const T1i = TR + TU;
  const T1a = TS - TT;

  const TM = ii[1], TN = ii[3], TO = ii[9];
  const TP = TN + TO;
  const TQ = TM - KP500000000 * TP;
  const T1h = TM + TP;
  const T19 = TN - TO;

  const T10 = ii[11], T11 = ii[6];
  const T12 = T10 + T11, T1d = T10 - T11;
  const T13 = ii[7], T14 = ii[2];
  const T15 = T13 + T14, T1c = T13 - T14;

  const T16 = T12 + T15;
  const T2c = T1d + T1c;
  const T2a = T1h - T1i;
  const T2d = T2b + T2c;
  const TW = TQ + TV;
  const T17 = TZ - KP500000000 * T16;
  const T18 = TW - T17;
  const T1n = TW + T17;

  const T2i = TQ - TV;
  const T2j = KP866025403 * (T15 - T12);
  const T2k = T2i + T2j;
  const T2n = T2i - T2j;
  const T1j = T1h + T1i;
  const T1k = TZ + T16;
  const T1l = KP300462606 * (T1j - T1k);
  const T1r = T1j + T1k;

  const T1b = T19 + T1a;
  const T1e = T1c - T1d;
  const T1f = T1b + T1e;
  const T1o = T1e - T1b;
  const T2f = T2b - KP500000000 * T2c;
  const T2g = KP866025403 * (T1a - T19);
  const T2h = T2f - T2g;
  const T2m = T2g + T2f;

  const ro = new Float64Array(13), io = new Float64Array(13);

  ro[0] = T1 + To;
  io[0] = T1q + T1r;

  const T1B = KP387390585 * T1f + KP265966249 * T18;
  const T1C = KP113854479 * T1o + KP503537032 * T1n;
  const T1D = T1B + T1C;
  const T1N = T1C - T1B;
  const T1y = KP575140729 * Tu + KP174138601 * Tt;
  const T1v = KP256247671 * TI - KP156891391 * TH;
  const T1w = KP011599105 * TF + KP300238635 * TA;
  const T1x = T1v - T1w;
  const T1E = T1y + T1x;
  const T1O = KP1_732050807 * (T1v + T1w);

  const Tv = KP575140729 * Tt - KP174138601 * Tu;
  const TG = KP011599105 * TA - KP300238635 * TF;
  const TJ = KP256247671 * TH + KP156891391 * TI;
  const TK = TG - TJ;
  const T1J = KP1_732050807 * (TJ + TG);
  const T1Q = Tv - TK;

  const T1g = KP258260390 * T18 - KP132983124 * T1f;
  const T1H = T1l - T1g;
  const T1p = KP075902986 * T1n - KP251768516 * T1o;
  const T1s = T1q - KP083333333 * T1r;
  const T1G = T1s - T1p;
  const T1m = KP2_000000000 * T1g + T1l;
  const T1R = T1H + T1G;
  const T1t = KP2_000000000 * T1p + T1s;
  const T1I = T1G - T1H;

  const TL = KP2_000000000 * TK + Tv;
  const T1u = T1m + T1t;
  io[1] = TL + T1u;
  io[12] = T1u - TL;

  const T1z = KP2_000000000 * T1x - T1y;
  const T1A = T1t - T1m;
  io[5] = T1z + T1A;
  io[8] = T1A - T1z;
  const T1T = T1R - T1Q;
  const T1U = T1O + T1N;
  io[4] = T1T - T1U;
  io[10] = T1U + T1T;

  const T1P = T1N - T1O;
  const T1S = T1Q + T1R;
  io[3] = T1P + T1S;
  io[9] = T1S - T1P;

  const T1L = T1J + T1I;
  const T1M = T1E + T1D;
  io[6] = T1L - T1M;
  io[11] = T1M + T1L;
  const T1F = T1D - T1E;
  const T1K = T1I - T1J;
  io[2] = T1F + T1K;
  io[7] = T1K - T1F;

  const T2w = KP387390585 * T20 + KP265966249 * T1X;
  const T2x = KP113854479 * T24 - KP503537032 * T25;
  const T2y = T2w + T2x;
  const T2I = T2w - T2x;
  const T2J = KP575140729 * T2a + KP174138601 * T2d;
  const T2z = KP011599105 * T2m - KP300238635 * T2n;
  const T2A = KP256247671 * T2k - KP156891391 * T2h;
  const T2K = T2z + T2A;
  const T2B = KP1_732050807 * (T2z - T2A);
  const T2L = T2J + T2K;

  const T2e = KP174138601 * T2a - KP575140729 * T2d;
  const T2l = KP256247671 * T2h + KP156891391 * T2k;
  const T2o = KP300238635 * T2m + KP011599105 * T2n;
  const T2p = T2l - T2o;
  const T2u = T2e - T2p;
  const T2G = KP1_732050807 * (T2o + T2l);

  const T21 = KP258260390 * T1X - KP132983124 * T20;
  const T2r = T22 - T21;
  const T26 = KP251768516 * T24 + KP075902986 * T25;
  const T27 = T1 - KP083333333 * To;
  const T2s = T27 - T26;
  const T23 = KP2_000000000 * T21 + T22;
  const T2F = T2s - T2r;
  const T28 = KP2_000000000 * T26 + T27;
  const T2t = T2r + T2s;

  const T29 = T23 + T28;
  const T2q = KP2_000000000 * T2p + T2e;
  ro[12] = T29 - T2q;
  ro[1] = T29 + T2q;

  const T2v = T2t - T2u;
  const T2C = T2y - T2B;
  ro[10] = T2v - T2C;
  ro[4] = T2v + T2C;
  const T2P = T28 - T23;
  const T2Q = KP2_000000000 * T2K - T2J;
  ro[5] = T2P - T2Q;
  ro[8] = T2P + T2Q;

  const T2N = T2F - T2G;
  const T2O = T2L - T2I;
  ro[11] = T2N - T2O;
  ro[6] = T2N + T2O;

  const T2H = T2F + T2G;
  const T2M = T2I + T2L;
  ro[7] = T2H - T2M;
  ro[2] = T2H + T2M;
  const T2D = T2t + T2u;
  const T2E = T2y + T2B;
  ro[3] = T2D - T2E;
  ro[9] = T2D + T2E;

  return [ro, io];
}

module.exports = { n1_13 };
