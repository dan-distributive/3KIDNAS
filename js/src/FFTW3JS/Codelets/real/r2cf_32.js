'use strict';

// =============================================================================
// r2cf_32.js -- faithful JS port of rdft/scalar/r2cf/r2cf_32.c (non-FMA
// branch). x[0..31] (real) -> O[0..31] packed halfcomplex: O[0..16]=Re0..
// Re16, O[17..31]=Im15..Im1 (O[n-k]=Im_k). R0/R1 (stride-2) input
// convention -- see r2cf_6.js's header.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP555570233 = 0.555570233019602224742830813948532874374937191;
const KP831469612 = 0.831469612302545237078788377617905756738560812;
const KP195090322 = 0.195090322016128267848284868477022240927691618;
const KP980785280 = 0.980785280403230449126182236134239036973933731;
const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function r2cf_32(x) {
  const T1 = x[0];
  const T2 = x[16];
  const T3 = T1 + T2;
  const T4 = x[8];
  const T5 = x[24];
  const T6 = T4 + T5;
  const T7 = T3 + T6;
  const T2b = T3 - T6;
  const Tv = T1 - T2;
  const T1l = T4 - T5;

  const T8 = x[4];
  const T9 = x[20];
  const Ta = T8 + T9;
  const Tw = T8 - T9;
  const Tb = x[28];
  const Tc = x[12];
  const Td = Tb + Tc;
  const Tx = Tb - Tc;
  const Te = Ta + Td;
  const T2o = Td - Ta;
  const Ty = KP707106781 * (Tw + Tx);
  const T1k = KP707106781 * (Tx - Tw);

  const Tn = x[30];
  const To = x[14];
  const Tp = Tn + To;
  const TD = Tn - To;
  const Tq = x[6];
  const Tr = x[22];
  const Ts = Tq + Tr;
  const TE = Tq - Tr;
  const Tt = Tp + Ts;
  const T2d = Tp - Ts;
  const TF = KP923879532 * TD + KP382683432 * TE; // FMA
  const T1h = KP382683432 * TD - KP923879532 * TE; // FNMS

  const Tg = x[2];
  const Th = x[18];
  const Ti = Tg + Th;
  const TA = Tg - Th;
  const Tj = x[10];
  const Tk = x[26];
  const Tl = Tj + Tk;
  const TB = Tj - Tk;
  const Tm = Ti + Tl;
  const T2c = Ti - Tl;
  const TC = KP923879532 * TA - KP382683432 * TB; // FNMS
  const T1i = KP382683432 * TA + KP923879532 * TB; // FMA

  const TZ = x[31];
  const T10 = x[15];
  const T11 = TZ - T10;
  const T1X = TZ + T10;
  const T1b = x[7];
  const T1c = x[23];
  const T1d = T1b - T1c;
  const T1Y = T1b + T1c;

  const T12 = x[3];
  const T13 = x[19];
  const T14 = T12 - T13;
  const T20 = T12 + T13;
  const T15 = x[27];
  const T16 = x[11];
  const T17 = T15 - T16;
  const T21 = T15 + T16;

  const T1Z = T1X + T1Y;
  const T22 = T20 + T21;
  const T2k = T21 - T20;
  const T2j = T1X - T1Y;
  const T1a = KP707106781 * (T17 - T14);
  const T1e = T1a - T1d;
  const T1C = T1d + T1a;
  const T18 = KP707106781 * (T14 + T17);
  const T19 = T11 + T18;
  const T1B = T11 - T18;

  const TI = x[1];
  const TJ = x[17];
  const TK = TI - TJ;
  const T1Q = TI + TJ;
  const TU = x[9];
  const TV = x[25];
  const TW = TU - TV;
  const T1R = TU + TV;

  const TL = x[5];
  const TM = x[21];
  const TN = TL - TM;
  const T1T = TL + TM;
  const TO = x[29];
  const TP = x[13];
  const TQ = TO - TP;
  const T1U = TO + TP;

  const T1S = T1Q + T1R;
  const T1V = T1T + T1U;
  const T2h = T1U - T1T;
  const T2g = T1Q - T1R;
  const TT = KP707106781 * (TQ - TN);
  const TX = TT - TW;
  const T1z = TW + TT;
  const TR = KP707106781 * (TN + TQ);
  const TS = TK + TR;
  const T1y = TK - TR;

  const O = new Float64Array(32);

  {
    const Tf = T7 + Te;
    const Tu = Tm + Tt;
    const T27 = Tf + Tu;
    const T28 = T1S + T1V;
    const T29 = T1Z + T22;
    const T2a = T28 + T29;
    O[8] = Tf - Tu;
    O[24] = T29 - T28;
    O[16] = T27 - T2a;
    O[0] = T27 + T2a;
  }
  {
    const T1P = T7 - Te;
    const T25 = Tt - Tm;
    const T1W = T1S - T1V;
    const T23 = T1Z - T22;
    const T24 = KP707106781 * (T1W + T23);
    const T26 = KP707106781 * (T23 - T1W);
    O[12] = T1P - T24;
    O[20] = T26 - T25;
    O[4] = T1P + T24;
    O[28] = T25 + T26;
  }
  {
    const T2e = KP707106781 * (T2c + T2d);
    const T2f = T2b + T2e;
    const T2v = T2b - T2e;
    const T2n = KP707106781 * (T2d - T2c);
    const T2p = T2n - T2o;
    const T2r = T2o + T2n;
    const T2i = KP923879532 * T2g + KP382683432 * T2h; // FMA
    const T2l = KP923879532 * T2j - KP382683432 * T2k; // FNMS
    const T2m = T2i + T2l;
    const T2q = T2l - T2i;
    const T2s = KP923879532 * T2h - KP382683432 * T2g; // FNMS
    const T2t = KP923879532 * T2k + KP382683432 * T2j; // FMA
    const T2u = T2s + T2t;
    const T2w = T2t - T2s;
    O[14] = T2f - T2m;
    O[18] = T2u - T2r;
    O[2] = T2f + T2m;
    O[30] = T2r + T2u;
    O[26] = T2p + T2q;
    O[6] = T2v + T2w;
    O[22] = T2q - T2p;
    O[10] = T2v - T2w;
  }
  {
    const Tz = Tv + Ty;
    const TG = TC + TF;
    const TH = Tz + TG;
    const T1t = Tz - TG;
    const T1q = KP980785280 * TX - KP195090322 * TS; // FNMS
    const T1r = KP195090322 * T19 + KP980785280 * T1e; // FMA
    const T1s = T1q + T1r;
    const T1u = T1r - T1q;

    const TY = KP980785280 * TS + KP195090322 * TX; // FMA
    const T1f = KP980785280 * T19 - KP195090322 * T1e; // FNMS
    const T1g = TY + T1f;
    const T1o = T1f - TY;
    const T1j = T1h - T1i;
    const T1m = T1k - T1l;
    const T1n = T1j - T1m;
    const T1p = T1m + T1j;

    O[15] = TH - T1g;
    O[17] = T1s - T1p;
    O[1] = TH + T1g;
    O[31] = T1p + T1s;
    O[25] = T1n + T1o;
    O[7] = T1t + T1u;
    O[23] = T1o - T1n;
    O[9] = T1t - T1u;
  }
  {
    const T1v = Tv - Ty;
    const T1w = T1i + T1h;
    const T1x = T1v + T1w;
    const T1N = T1v - T1w;
    const T1K = KP831469612 * T1z - KP555570233 * T1y; // FNMS
    const T1L = KP555570233 * T1B + KP831469612 * T1C; // FMA
    const T1M = T1K + T1L;
    const T1O = T1L - T1K;

    const T1A = KP831469612 * T1y + KP555570233 * T1z; // FMA
    const T1D = KP831469612 * T1B - KP555570233 * T1C; // FNMS
    const T1E = T1A + T1D;
    const T1I = T1D - T1A;
    const T1F = TF - TC;
    const T1G = T1l + T1k;
    const T1H = T1F - T1G;
    const T1J = T1G + T1F;

    O[13] = T1x - T1E;
    O[19] = T1M - T1J;
    O[3] = T1x + T1E;
    O[29] = T1J + T1M;
    O[27] = T1H + T1I;
    O[5] = T1N + T1O;
    O[21] = T1I - T1H;
    O[11] = T1N - T1O;
  }

  return O;
}

module.exports = { r2cf_32 };
