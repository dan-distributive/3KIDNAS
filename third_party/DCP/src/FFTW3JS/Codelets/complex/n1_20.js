'use strict';

// =============================================================================
// n1_20.js -- faithful JS port of dft/scalar/codelets/n1_20.c (non-FMA
// branch), FFTW3's direct (base-case) radix-20 (4x5) complex DFT codelet.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;

function n1_20(ri, ii) {
  const T1 = ri[0], T2 = ri[10];
  const T3 = T1 + T2, T1Q = T1 - T2;
  const TL = ii[0], TM = ii[10];
  const TN = TL - TM, T2O = TL + TM;

  const T4 = ri[5], T5 = ri[15];
  const T6 = T4 + T5, TO = T4 - T5;
  const T1R = ii[5], T1S = ii[15];
  const T1T = T1R - T1S, T2P = T1R + T1S;

  const T7 = T3 - T6, T2Q = T2O - T2P, T3h = T2O + T2P, TD = T3 + T6;
  const TP = TN - TO, T1U = T1Q - T1T, T2l = T1Q + T1T, T1d = TO + TN;

  const Tn = ri[8], To = ri[18];
  const Tp = Tn + To, T1o = Tn - To;
  const T11 = ii[8], T12 = ii[18];
  const T13 = T11 - T12, T2u = T11 + T12;

  const Tq = ri[13], Tr = ri[3];
  const Ts = Tq + Tr, T14 = Tq - Tr;
  const T1p = ii[13], T1q = ii[3];
  const T1r = T1p - T1q, T2v = T1p + T1q;

  const Tu = ri[12], Tv = ri[2];
  const Tw = Tu + Tv, T1t = Tu - Tv;
  const T16 = ii[12], T17 = ii[2];
  const T18 = T16 - T17, T2x = T16 + T17;

  const Tx = ri[17], Ty = ri[7];
  const Tz = Tx + Ty, T19 = Tx - Ty;
  const T1u = ii[17], T1v = ii[7];
  const T1w = T1u - T1v, T2y = T1u + T1v;

  const Tt = Tp - Ts, TA = Tw - Tz, TB = Tt + TA;
  const T2w = T2u - T2v, T2z = T2x - T2y, T2S = T2w + T2z;
  const T35 = T2u + T2v, T36 = T2x + T2y, T3f = T35 + T36;
  const TH = Tp + Ts, TI = Tw + Tz, TJ = TH + TI;
  const T15 = T13 - T14, T1a = T18 - T19, T1b = T15 + T1a;
  const T1s = T1o - T1r, T1x = T1t - T1w, T1W = T1s + T1x;
  const T29 = T1o + T1r, T2a = T1t + T1w, T2j = T29 + T2a;
  const T1h = T14 + T13, T1i = T19 + T18, T1j = T1h + T1i;

  const T8 = ri[4], T9 = ri[14];
  const Ta = T8 + T9, T1z = T8 - T9;
  const TQ = ii[4], TR = ii[14];
  const TS = TQ - TR, T2B = TQ + TR;

  const Tb = ri[9], Tc = ri[19];
  const Td = Tb + Tc, TT = Tb - Tc;
  const T1A = ii[9], T1B = ii[19];
  const T1C = T1A - T1B, T2C = T1A + T1B;

  const Tf = ri[16], Tg = ri[6];
  const Th = Tf + Tg, T1E = Tf - Tg;
  const TV = ii[16], TW = ii[6];
  const TX = TV - TW, T2E = TV + TW;

  const Ti = ri[1], Tj = ri[11];
  const Tk = Ti + Tj, TY = Ti - Tj;
  const T1F = ii[1], T1G = ii[11];
  const T1H = T1F - T1G, T2F = T1F + T1G;

  const Te = Ta - Td, Tl = Th - Tk, Tm = Te + Tl;
  const T2D = T2B - T2C, T2G = T2E - T2F, T2R = T2D + T2G;
  const T32 = T2B + T2C, T33 = T2E + T2F, T3e = T32 + T33;
  const TE = Ta + Td, TF = Th + Tk, TG = TE + TF;
  const TU = TS - TT, TZ = TX - TY, T10 = TU + TZ;
  const T1D = T1z - T1C, T1I = T1E - T1H, T1V = T1D + T1I;
  const T26 = T1z + T1C, T27 = T1E + T1H, T2i = T26 + T27;
  const T1e = TT + TS, T1f = TY + TX, T1g = T1e + T1f;

  const ro = new Float64Array(20), io = new Float64Array(20);

  {
    const T2s = KP559016994 * (Tm - TB);
    const TC = Tm + TB;
    const T2r = T7 - KP250000000 * TC;
    const T2A = T2w - T2z;
    const T2H = T2D - T2G;
    const T2I = KP951056516 * T2A - KP587785252 * T2H;
    const T2K = KP951056516 * T2H + KP587785252 * T2A;
    ro[10] = T7 + TC;
    const T2J = T2s + T2r;
    ro[14] = T2J - T2K;
    ro[6] = T2J + T2K;
    const T2t = T2r - T2s;
    ro[2] = T2t - T2I;
    ro[18] = T2t + T2I;
  }

  {
    const T2V = KP559016994 * (T2R - T2S);
    const T2T = T2R + T2S;
    const T2U = T2Q - KP250000000 * T2T;
    const T2L = Tt - TA;
    const T2M = Te - Tl;
    const T2N = KP951056516 * T2L - KP587785252 * T2M;
    const T2Y = KP951056516 * T2M + KP587785252 * T2L;
    io[10] = T2Q + T2T;
    const T2X = T2V + T2U;
    io[6] = T2X - T2Y;
    io[14] = T2Y + T2X;
    const T2W = T2U - T2V;
    io[2] = T2N + T2W;
    io[18] = T2W - T2N;
  }

  {
    const T2Z = KP559016994 * (TG - TJ);
    const TK = TG + TJ;
    const T30 = TD - KP250000000 * TK;
    const T34 = T32 - T33;
    const T37 = T35 - T36;
    const T38 = KP951056516 * T34 + KP587785252 * T37;
    const T3a = KP951056516 * T37 - KP587785252 * T34;
    ro[0] = TD + TK;
    const T39 = T30 - T2Z;
    ro[12] = T39 - T3a;
    ro[8] = T39 + T3a;
    const T31 = T2Z + T30;
    ro[4] = T31 - T38;
    ro[16] = T31 + T38;
  }

  {
    const T3g = KP559016994 * (T3e - T3f);
    const T3i = T3e + T3f;
    const T3j = T3h - KP250000000 * T3i;
    const T3b = TE - TF;
    const T3c = TH - TI;
    const T3d = KP951056516 * T3b + KP587785252 * T3c;
    const T3m = KP951056516 * T3c - KP587785252 * T3b;
    io[0] = T3h + T3i;
    const T3l = T3j - T3g;
    io[8] = T3l - T3m;
    io[12] = T3m + T3l;
    const T3k = T3g + T3j;
    io[4] = T3d + T3k;
    io[16] = T3k - T3d;
  }

  {
    const T23 = KP559016994 * (T10 - T1b);
    const T1c = T10 + T1b;
    const T24 = TP - KP250000000 * T1c;
    const T28 = T26 - T27;
    const T2b = T29 - T2a;
    const T2c = KP951056516 * T28 + KP587785252 * T2b;
    const T2e = KP951056516 * T2b - KP587785252 * T28;
    io[5] = TP + T1c;
    const T2d = T24 - T23;
    io[13] = T2d - T2e;
    io[17] = T2d + T2e;
    const T25 = T23 + T24;
    io[1] = T25 - T2c;
    io[9] = T25 + T2c;
  }

  {
    const T2k = KP559016994 * (T2i - T2j);
    const T2m = T2i + T2j;
    const T2n = T2l - KP250000000 * T2m;
    const T2f = TU - TZ;
    const T2g = T15 - T1a;
    const T2h = KP951056516 * T2f + KP587785252 * T2g;
    const T2p = KP951056516 * T2g - KP587785252 * T2f;
    ro[5] = T2l + T2m;
    const T2q = T2n - T2k;
    ro[13] = T2p + T2q;
    ro[17] = T2q - T2p;
    const T2o = T2k + T2n;
    ro[1] = T2h + T2o;
    ro[9] = T2o - T2h;
  }

  {
    const T1m = KP559016994 * (T1g - T1j);
    const T1k = T1g + T1j;
    const T1l = T1d - KP250000000 * T1k;
    const T1y = T1s - T1x;
    const T1J = T1D - T1I;
    const T1K = KP951056516 * T1y - KP587785252 * T1J;
    const T1M = KP951056516 * T1J + KP587785252 * T1y;
    io[15] = T1d + T1k;
    const T1L = T1m + T1l;
    io[11] = T1L - T1M;
    io[19] = T1L + T1M;
    const T1n = T1l - T1m;
    io[3] = T1n - T1K;
    io[7] = T1n + T1K;
  }

  {
    const T1Z = KP559016994 * (T1V - T1W);
    const T1X = T1V + T1W;
    const T1Y = T1U - KP250000000 * T1X;
    const T1N = T1h - T1i;
    const T1O = T1e - T1f;
    const T1P = KP951056516 * T1N - KP587785252 * T1O;
    const T21 = KP951056516 * T1O + KP587785252 * T1N;
    ro[15] = T1U + T1X;
    const T22 = T1Z + T1Y;
    ro[11] = T21 + T22;
    ro[19] = T22 - T21;
    const T20 = T1Y - T1Z;
    ro[3] = T1P + T20;
    ro[7] = T20 - T1P;
  }

  return [ro, io];
}

module.exports = { n1_20 };
