'use strict';

// =============================================================================
// n1_15.js -- faithful JS port of dft/scalar/codelets/n1_15.c (non-FMA
// branch), FFTW3's direct (base-case) radix-15 (3x5) complex DFT codelet.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP500000000 = 0.5;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function n1_15(ri, ii) {
  const T1 = ri[0], T1z = ii[0];
  const T2 = ri[5], T3 = ri[10];
  const T4 = T2 + T3, T1y = KP866025403 * (T3 - T2);
  const Tu = ii[5], Tv = ii[10];
  const Tw = KP866025403 * (Tu - Tv), T1A = Tu + Tv;
  const T5 = T1 + T4, T2l = T1z + T1A;
  const Tt = T1 - KP500000000 * T4;
  const Tx = Tt - Tw, TV = Tt + Tw;
  const T1B = T1z - KP500000000 * T1A;
  const T1C = T1y + T1B, T20 = T1B - T1y;

  const Th = ri[6], Ti = ri[11], Tj = ri[1];
  const Tk = Ti + Tj;
  const TJ = Th - KP500000000 * Tk;
  const T1h = KP866025403 * (Tj - Ti);
  const T1i = ii[6], TK = ii[11], TL = ii[1];
  const T1j = TK + TL;
  const TM = KP866025403 * (TK - TL);
  const T1k = T1i - KP500000000 * T1j;

  const Tm = ri[9], Tn = ri[14], To = ri[4];
  const Tp = Tn + To;
  const TO = Tm - KP500000000 * Tp;
  const T1m = KP866025403 * (To - Tn);
  const T1n = ii[9], TP = ii[14], TQ = ii[4];
  const T1o = TP + TQ;
  const TR = KP866025403 * (TP - TQ);
  const T1p = T1n - KP500000000 * T1o;

  const Tl = Th + Tk, Tq = Tm + Tp, Tr = Tl + Tq;
  const TN = TJ - TM, TS = TO - TR, TT = TN + TS;
  const T2c = T1i + T1j, T2d = T1n + T1o, T2n = T2c + T2d;
  const T1O = T1k - T1h, T1P = T1p - T1m, T22 = T1O + T1P;
  const T1l = T1h + T1k, T1q = T1m + T1p, T1w = T1l + T1q;
  const TZ = TJ + TM, T10 = TO + TR, T11 = TZ + T10;

  const T6 = ri[3], T7 = ri[8], T8 = ri[13];
  const T9 = T7 + T8;
  const Ty = T6 - KP500000000 * T9;
  const T16 = KP866025403 * (T8 - T7);
  const T17 = ii[3], Tz = ii[8], TA = ii[13];
  const T18 = Tz + TA;
  const TB = KP866025403 * (Tz - TA);
  const T19 = T17 - KP500000000 * T18;

  const Tb = ri[12], Tc = ri[2], Td = ri[7];
  const Te = Tc + Td;
  const TD = Tb - KP500000000 * Te;
  const T1b = KP866025403 * (Td - Tc);
  const T1c = ii[12], TE = ii[2], TF = ii[7];
  const T1d = TE + TF;
  const TG = KP866025403 * (TE - TF);
  const T1e = T1c - KP500000000 * T1d;

  const Ta = T6 + T9, Tf = Tb + Te, Tg = Ta + Tf;
  const TC = Ty - TB, TH = TD - TG, TI = TC + TH;
  const T2f = T17 + T18, T2g = T1c + T1d, T2m = T2f + T2g;
  const T1R = T19 - T16, T1S = T1e - T1b, T21 = T1R + T1S;
  const T1a = T16 + T19, T1f = T1b + T1e, T1v = T1a + T1f;
  const TW = Ty + TB, TX = TD + TG, TY = TW + TX;

  const ro = new Float64Array(15), io = new Float64Array(15);

  {
    const T2a = KP559016994 * (Tg - Tr);
    const Ts = Tg + Tr;
    const T29 = T5 - KP250000000 * Ts;
    const T2e = T2c - T2d;
    const T2h = T2f - T2g;
    const T2i = KP951056516 * T2e - KP587785252 * T2h;
    const T2k = KP951056516 * T2h + KP587785252 * T2e;
    ro[0] = T5 + Ts;
    const T2j = T2a + T29;
    ro[9] = T2j - T2k;
    ro[6] = T2j + T2k;
    const T2b = T29 - T2a;
    ro[12] = T2b - T2i;
    ro[3] = T2b + T2i;
  }

  {
    const T2q = KP559016994 * (T2m - T2n);
    const T2o = T2m + T2n;
    const T2p = T2l - KP250000000 * T2o;
    const T2s = Tl - Tq;
    const T2t = Ta - Tf;
    const T2u = KP951056516 * T2s - KP587785252 * T2t;
    const T2w = KP951056516 * T2t + KP587785252 * T2s;
    io[0] = T2l + T2o;
    const T2v = T2q + T2p;
    io[6] = T2v - T2w;
    io[9] = T2w + T2v;
    const T2r = T2p - T2q;
    io[3] = T2r - T2u;
    io[12] = T2u + T2r;
  }

  {
    const T1M = KP559016994 * (TI - TT);
    const TU = TI + TT;
    const T1L = Tx - KP250000000 * TU;
    const T1Q = T1O - T1P;
    const T1T = T1R - T1S;
    const T1U = KP951056516 * T1Q - KP587785252 * T1T;
    const T1W = KP951056516 * T1T + KP587785252 * T1Q;
    ro[5] = Tx + TU;
    const T1V = T1M + T1L;
    ro[14] = T1V - T1W;
    ro[11] = T1V + T1W;
    const T1N = T1L - T1M;
    ro[2] = T1N - T1U;
    ro[8] = T1N + T1U;
  }

  {
    const T25 = KP559016994 * (T21 - T22);
    const T23 = T21 + T22;
    const T24 = T20 - KP250000000 * T23;
    const T1X = TN - TS;
    const T1Y = TC - TH;
    const T1Z = KP951056516 * T1X - KP587785252 * T1Y;
    const T28 = KP951056516 * T1Y + KP587785252 * T1X;
    io[5] = T20 + T23;
    const T27 = T25 + T24;
    io[11] = T27 - T28;
    io[14] = T28 + T27;
    const T26 = T24 - T25;
    io[2] = T1Z + T26;
    io[8] = T26 - T1Z;
  }

  {
    const T1x = KP559016994 * (T1v - T1w);
    const T1D = T1v + T1w;
    const T1E = T1C - KP250000000 * T1D;
    const T1G = TW - TX;
    const T1H = TZ - T10;
    const T1I = KP951056516 * T1G + KP587785252 * T1H;
    const T1J = KP951056516 * T1H - KP587785252 * T1G;
    io[10] = T1C + T1D;
    const T1K = T1E - T1x;
    io[7] = T1J + T1K;
    io[13] = T1K - T1J;
    const T1F = T1x + T1E;
    io[1] = T1F - T1I;
    io[4] = T1I + T1F;
  }

  {
    const T13 = KP559016994 * (TY - T11);
    const T12 = TY + T11;
    const T14 = TV - KP250000000 * T12;
    const T1g = T1a - T1f;
    const T1r = T1l - T1q;
    const T1s = KP951056516 * T1g + KP587785252 * T1r;
    const T1u = KP951056516 * T1r - KP587785252 * T1g;
    ro[10] = TV + T12;
    const T1t = T14 - T13;
    ro[7] = T1t - T1u;
    ro[13] = T1t + T1u;
    const T15 = T13 + T14;
    ro[4] = T15 - T1s;
    ro[1] = T15 + T1s;
  }

  return [ro, io];
}

module.exports = { n1_15 };
