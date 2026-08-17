'use strict';

// =============================================================================
// hf_15.js -- faithful JS port of rdft/scalar/r2cf/hf_15.c (non-FMA branch).
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hf_2/3/4/5/7/9.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP500000000 = 0.5;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function hf_15(cr, ci, Wc, Ws) {
  const T1 = cr[0], T2R = ci[0];
  const T3 = cr[5], T5 = ci[5];
  const T6 = Wc[5] * T3 + Ws[5] * T5;
  const T1o = Wc[5] * T5 - Ws[5] * T3;
  const T8 = cr[10], Ta = ci[10];
  const Tb = Wc[10] * T8 + Ws[10] * Ta;
  const T1p = Wc[10] * Ta - Ws[10] * T8;

  const T1q = KP866025403 * (T1o - T1p);
  const T2Q = KP866025403 * (Tb - T6);
  const Tc = T6 + Tb;
  const Td = T1 + Tc;
  const T1n = T1 - KP500000000 * Tc;
  const T2S = T1o + T1p;
  const T2T = T2R - KP500000000 * T2S;
  const T3l = T2S + T2R;

  const TO = cr[6], TQ = ci[6];
  const TR = Wc[6] * TO + Ws[6] * TQ;
  const T2c = Wc[6] * TQ - Ws[6] * TO;
  const T15 = cr[9], T17 = ci[9];
  const T18 = Wc[9] * T15 + Ws[9] * T17;
  const T2h = Wc[9] * T17 - Ws[9] * T15;
  const TT = cr[11], TV = ci[11];
  const TW = Wc[11] * TT + Ws[11] * TV;
  const T1E = Wc[11] * TV - Ws[11] * TT;
  const TY = cr[1], T10 = ci[1];
  const T11 = Wc[1] * TY + Ws[1] * T10;
  const T1F = Wc[1] * T10 - Ws[1] * TY;

  const T12 = TW + T11;
  const T2d = T1E + T1F;

  const T1a = cr[14], T1c = ci[14];
  const T1d = Wc[14] * T1a + Ws[14] * T1c;
  const T1J = Wc[14] * T1c - Ws[14] * T1a;
  const T1f = cr[4], T1h = ci[4];
  const T1i = Wc[4] * T1f + Ws[4] * T1h;
  const T1K = Wc[4] * T1h - Ws[4] * T1f;

  const T1j = T1d + T1i;
  const T2i = T1J + T1K;

  const T13 = TR + T12;
  const T1k = T18 + T1j;
  const T1l = T13 + T1k;
  const T2E = T2c + T2d;
  const T2F = T2h + T2i;
  const T3j = T2E + T2F;
  const T1D = TR - KP500000000 * T12;
  const T1G = KP866025403 * (T1E - T1F);
  const T1H = T1D - T1G, T1T = T1D + T1G;
  const T2g = KP866025403 * (T1d - T1i);
  const T2j = T2h - KP500000000 * T2i;
  const T2k = T2g - T2j, T2w = T2g + T2j;

  const T2b = KP866025403 * (T11 - TW);
  const T2e = T2c - KP500000000 * T2d;
  const T2f = T2b + T2e, T2v = T2e - T2b;
  const T1I = T18 - KP500000000 * T1j;
  const T1L = KP866025403 * (T1J - T1K);
  const T1M = T1I - T1L, T1U = T1I + T1L;

  const Tf = cr[3], Th = ci[3];
  const Ti = Wc[3] * Tf + Ws[3] * Th;
  const T21 = Wc[3] * Th - Ws[3] * Tf;
  const Tw = cr[12], Ty = ci[12];
  const Tz = Wc[12] * Tw + Ws[12] * Ty;
  const T26 = Wc[12] * Ty - Ws[12] * Tw;
  const Tk = cr[8], Tm = ci[8];
  const Tn = Wc[8] * Tk + Ws[8] * Tm;
  const T1t = Wc[8] * Tm - Ws[8] * Tk;
  const Tp = cr[13], Tr = ci[13];
  const Ts = Wc[13] * Tp + Ws[13] * Tr;
  const T1u = Wc[13] * Tr - Ws[13] * Tp;

  const Tt = Tn + Ts;
  const T22 = T1t + T1u;

  const TB = cr[2], TD = ci[2];
  const TE = Wc[2] * TB + Ws[2] * TD;
  const T1y = Wc[2] * TD - Ws[2] * TB;
  const TG = cr[7], TI = ci[7];
  const TJ = Wc[7] * TG + Ws[7] * TI;
  const T1z = Wc[7] * TI - Ws[7] * TG;

  const TK = TE + TJ;
  const T27 = T1y + T1z;

  const Tu = Ti + Tt;
  const TL = Tz + TK;
  const TM = Tu + TL;
  const T2H = T21 + T22;
  const T2I = T26 + T27;
  const T3i = T2H + T2I;
  const T1s = Ti - KP500000000 * Tt;
  const T1v = KP866025403 * (T1t - T1u);
  const T1w = T1s - T1v, T1Q = T1s + T1v;
  const T25 = KP866025403 * (TJ - TE);
  const T28 = T26 - KP500000000 * T27;
  const T29 = T25 + T28, T2t = T28 - T25;

  const T20 = KP866025403 * (Ts - Tn);
  const T23 = T21 - KP500000000 * T22;
  const T24 = T20 + T23, T2s = T23 - T20;
  const T1x = Tz - KP500000000 * TK;
  const T1A = KP866025403 * (T1y - T1z);
  const T1B = T1x - T1A, T1R = T1x + T1A;

  const outCr = new Float64Array(15), outCi = new Float64Array(15);

  const T2C = KP559016994 * (TM - T1l);
  const T1m = TM + T1l;
  const T2B = Td - KP250000000 * T1m;
  const T2G = T2E - T2F;
  const T2J = T2H - T2I;
  const T2K = KP951056516 * T2G - KP587785252 * T2J;
  const T2M = KP951056516 * T2J + KP587785252 * T2G;
  outCr[0] = Td + T1m;
  const T2L = T2C + T2B;
  outCi[5] = T2L - T2M;
  outCr[6] = T2L + T2M;
  const T2D = T2B - T2C;
  outCi[2] = T2D - T2K;
  outCr[3] = T2D + T2K;

  const T3k = KP559016994 * (T3i - T3j);
  const T3m = T3i + T3j;
  const T3n = T3l - KP250000000 * T3m;
  const T3f = T1k - T13;
  const T3g = Tu - TL;
  const T3h = KP587785252 * T3f - KP951056516 * T3g;
  const T3p = KP587785252 * T3g + KP951056516 * T3f;
  outCi[14] = T3m + T3l;
  const T3q = T3n - T3k;
  outCr[12] = T3p - T3q;
  outCi[11] = T3p + T3q;
  const T3o = T3k + T3n;
  outCr[9] = T3h - T3o;
  outCi[8] = T3h + T3o;

  const T2u = T2s - T2t;
  const T2x = T2v - T2w;
  const T2y = KP951056516 * T2u + KP587785252 * T2x;
  const T2A = KP951056516 * T2x - KP587785252 * T2u;
  const T1r = T1n - T1q;
  const T1C = T1w + T1B;
  const T1N = T1H + T1M;
  const T1O = T1C + T1N;
  const T2p = KP559016994 * (T1C - T1N);
  const T2q = T1r - KP250000000 * T1O;

  outCr[5] = T1r + T1O;
  const T2z = T2q - T2p;
  outCr[2] = T2z - T2A;
  outCi[6] = T2z + T2A;
  const T2r = T2p + T2q;
  outCi[0] = T2r - T2y;
  outCi[3] = T2r + T2y;

  const T33 = T1w - T1B;
  const T34 = T1H - T1M;
  const T35 = KP951056516 * T33 + KP587785252 * T34;
  const T3d = KP951056516 * T34 - KP587785252 * T33;
  const T39 = T2T - T2Q;
  const T36 = T2v + T2w;
  const T37 = T2s + T2t;
  const T3a = T37 + T36;
  const T38 = KP559016994 * (T36 - T37);
  const T3b = T39 - KP250000000 * T3a;

  outCi[9] = T3a + T39;
  const T3e = T38 + T3b;
  outCr[8] = T3d - T3e;
  outCi[12] = T3d + T3e;
  const T3c = T38 - T3b;
  outCr[11] = T35 + T3c;
  outCr[14] = T3c - T35;

  const T2V = T1T - T1U;
  const T2W = T1Q - T1R;
  const T2X = KP951056516 * T2V - KP587785252 * T2W;
  const T31 = KP951056516 * T2W + KP587785252 * T2V;
  const T2U = T2Q + T2T;
  const T2N = T2k - T2f;
  const T2O = T24 + T29;
  const T2P = T2N - T2O;
  const T2Y = KP250000000 * T2P + T2U;
  const T2Z = KP559016994 * (T2O + T2N);

  outCr[10] = T2P - T2U;
  const T32 = T2Z + T2Y;
  outCi[10] = T31 + T32;
  outCi[13] = T32 - T31;
  const T30 = T2Y - T2Z;
  outCr[13] = T2X - T30;
  outCi[7] = T2X + T30;

  const T2a = T24 - T29;
  const T2l = T2f + T2k;
  const T2m = KP951056516 * T2a + KP587785252 * T2l;
  const T2o = KP951056516 * T2l - KP587785252 * T2a;
  const T1P = T1n + T1q;
  const T1S = T1Q + T1R;
  const T1V = T1T + T1U;
  const T1W = T1S + T1V;
  const T1X = KP559016994 * (T1S - T1V);
  const T1Y = T1P - KP250000000 * T1W;

  outCi[4] = T1P + T1W;
  const T1Z = T1X + T1Y;
  outCr[4] = T1Z - T2m;
  outCr[1] = T1Z + T2m;
  const T2n = T1Y - T1X;
  outCr[7] = T2n - T2o;
  outCi[1] = T2n + T2o;

  return [outCr, outCi];
}

module.exports = { hf_15 };
