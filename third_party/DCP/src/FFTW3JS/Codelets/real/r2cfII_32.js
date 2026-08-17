'use strict';

// =============================================================================
// r2cfII_32.js -- faithful JS port of rdft/scalar/r2cf/r2cfII_32.c (non-FMA
// branch). R2HCII ("shifted") radix-32 direct codelet -- the "cldm"
// middle-column combine for an EVEN outer radix (see r2cfII_8.js's header
// for the general even-radix derivation: clean 16-and-16 Cr/Ci split,
// out[31-k]=Ci[k]).
// INPUT: ph[0..31], R0[k]=ph[2k] (k=0..15), R1[k]=ph[2k+1] (k=0..15).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP471396736 = 0.471396736825997648556387625905254377657460319;
const KP881921264 = 0.881921264348355029712756863660388349508442621;
const KP634393284 = 0.634393284163645498215171613225493370675687095;
const KP773010453 = 0.773010453362736960810906609758469800971041293;
const KP290284677 = 0.290284677254462367636192375817395274691476278;
const KP956940335 = 0.956940335732208864935797886980269969482849206;
const KP995184726 = 0.995184726672196886244836953109479921575474869;
const KP098017140 = 0.098017140329560601994195563888641845861136673;
const KP555570233 = 0.555570233019602224742830813948532874374937191;
const KP831469612 = 0.831469612302545237078788377617905756738560812;
const KP195090322 = 0.195090322016128267848284868477022240927691618;
const KP980785280 = 0.980785280403230449126182236134239036973933731;
const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function r2cfII_32(ph) {
  const T1 = ph[0];
  const T2p = ph[16];
  const T2 = ph[8];
  const T3 = ph[24];
  const T4 = KP707106781 * (T2 - T3);
  const T2o = KP707106781 * (T2 + T3);
  const T5 = T1 + T4;
  const T2D = T2p - T2o;
  const T1z = T1 - T4;
  const T2q = T2o + T2p;

  const T6 = ph[4];
  const T7 = ph[20];
  const T8 = KP923879532 * T6 - KP382683432 * T7; // FNMS
  const T1A = KP382683432 * T6 + KP923879532 * T7; // FMA
  const T9 = ph[12];
  const Ta = ph[28];
  const Tb = KP382683432 * T9 - KP923879532 * Ta; // FNMS
  const T1B = KP923879532 * T9 + KP382683432 * Ta; // FMA

  const Tc = T8 + Tb;
  const T2C = Tb - T8;
  const T1C = T1A - T1B;
  const T2n = T1A + T1B;

  const Te = ph[2];
  const Tk = ph[18];
  const Tf = ph[10];
  const Tg = ph[26];
  const Th = KP707106781 * (Tf - Tg);
  const Tj = KP707106781 * (Tf + Tg);
  const Ti = Te + Th;
  const Tl = Tj + Tk;
  const Tm = KP980785280 * Ti - KP195090322 * Tl; // FNMS
  const T1k = KP195090322 * Ti + KP980785280 * Tl; // FMA
  const T1H = Tk - Tj;
  const T1I = Te - Th;
  const T1J = KP831469612 * T1H - KP555570233 * T1I; // FNMS
  const T26 = KP831469612 * T1I + KP555570233 * T1H; // FMA

  const Tq = ph[30];
  const Tt = ph[14];
  const Tn = ph[6];
  const To = ph[22];
  const Tp = KP707106781 * (Tn - To);
  const Ts = KP707106781 * (Tn + To);
  const Tr = Tp - Tq;
  const Tu = Ts + Tt;
  const Tv = KP980785280 * Tr + KP195090322 * Tu; // FMA
  const T1l = KP195090322 * Tr - KP980785280 * Tu; // FNMS
  const T1E = Tt - Ts;
  const T1F = Tp + Tq;
  const T1G = KP831469612 * T1E - KP555570233 * T1F; // FNMS
  const T27 = KP831469612 * T1F + KP555570233 * T1E; // FMA

  const TW = ph[31];
  const T1a = ph[15];
  const TT = ph[7];
  const TU = ph[23];
  const TV = KP707106781 * (TT - TU);
  const T19 = KP707106781 * (TT + TU);

  const TY = ph[3];
  const TZ = ph[19];
  const T10 = KP923879532 * TY - KP382683432 * TZ; // FNMS
  const T16 = KP382683432 * TY + KP923879532 * TZ; // FMA
  const T11 = ph[11];
  const T12 = ph[27];
  const T13 = KP382683432 * T11 - KP923879532 * T12; // FNMS
  const T17 = KP923879532 * T11 + KP382683432 * T12; // FMA

  const TX = TV - TW;
  const T14 = T10 + T13;
  const T15 = TX + T14;
  const T1r = TX - T14;
  const T1W = T13 - T10;
  const T1X = T1a - T19;
  const T1Y = T1W - T1X;
  const T2e = T1W + T1X;

  const T18 = T16 + T17;
  const T1b = T19 + T1a;
  const T1c = T18 + T1b;
  const T1s = T1b - T18;
  const T1T = TV + TW;
  const T1U = T16 - T17;
  const T1V = T1T + T1U;
  const T2d = T1U - T1T;

  const Ty = ph[1];
  const TP = ph[17];
  const Tz = ph[9];
  const TA = ph[25];
  const TB = KP707106781 * (Tz - TA);
  const TO = KP707106781 * (Tz + TA);

  const TD = ph[5];
  const TE = ph[21];
  const TF = KP923879532 * TD - KP382683432 * TE; // FNMS
  const TL = KP382683432 * TD + KP923879532 * TE; // FMA
  const TG = ph[13];
  const TH = ph[29];
  const TI = KP382683432 * TG - KP923879532 * TH; // FNMS
  const TM = KP923879532 * TG + KP382683432 * TH; // FMA

  const TC = Ty + TB;
  const TJ = TF + TI;
  const TK = TC + TJ;
  const T1o = TC - TJ;
  const T1P = TI - TF;
  const T1Q = TP - TO;
  const T1R = T1P - T1Q;
  const T2b = T1P + T1Q;

  const TN = TL + TM;
  const TQ = TO + TP;
  const TR = TN + TQ;
  const T1p = TQ - TN;
  const T1M = Ty - TB;
  const T1N = TL - TM;
  const T1O = T1M - T1N;
  const T2a = T1M + T1N;

  const out = new Float64Array(32);

  {
    const Td = T5 + Tc;
    const Tw = Tm + Tv;
    const Tx = Td - Tw;
    const T1f = Td + Tw;
    const T2m = T1l - T1k;
    const T2r = T2n + T2q;
    const T2s = T2m - T2r;
    const T2u = T2m + T2r;

    const TS = KP098017140 * TK + KP995184726 * TR; // FMA
    const T1d = KP098017140 * T15 - KP995184726 * T1c; // FNMS
    const T1e = TS + T1d;
    const T2l = T1d - TS;
    const T1g = KP995184726 * TK - KP098017140 * TR; // FNMS
    const T1h = KP995184726 * T15 + KP098017140 * T1c; // FMA
    const T1i = T1g + T1h;
    const T2t = T1h - T1g;

    out[8] = Tx - T1e;
    out[23] = T2t - T2u;
    out[7] = Tx + T1e;
    out[24] = T2t + T2u;
    out[15] = T1f - T1i;
    out[16] = T2l - T2s;
    out[0] = T1f + T1i;
    out[31] = T2l + T2s;
  }
  {
    const T25 = T1z + T1C;
    const T28 = T26 - T27;
    const T29 = T25 + T28;
    const T2h = T25 - T28;
    const T2K = T1J + T1G;
    const T2L = T2C + T2D;
    const T2M = T2K - T2L;
    const T2O = T2K + T2L;

    const T2c = KP956940335 * T2a + KP290284677 * T2b; // FMA
    const T2f = KP956940335 * T2d - KP290284677 * T2e; // FNMS
    const T2g = T2c + T2f;
    const T2J = T2f - T2c;
    const T2i = KP290284677 * T2d + KP956940335 * T2e; // FMA
    const T2j = KP956940335 * T2b - KP290284677 * T2a; // FNMS
    const T2k = T2i - T2j;
    const T2N = T2j + T2i;

    out[14] = T29 - T2g;
    out[17] = T2N - T2O;
    out[1] = T29 + T2g;
    out[30] = T2N + T2O;
    out[9] = T2h - T2k;
    out[22] = T2J - T2M;
    out[6] = T2h + T2k;
    out[25] = T2J + T2M;
  }
  {
    const T1j = T5 - Tc;
    const T1m = T1k + T1l;
    const T1n = T1j + T1m;
    const T1v = T1j - T1m;
    const T2w = Tv - Tm;
    const T2x = T2q - T2n;
    const T2y = T2w - T2x;
    const T2A = T2w + T2x;

    const T1q = KP773010453 * T1o + KP634393284 * T1p; // FMA
    const T1t = KP773010453 * T1r - KP634393284 * T1s; // FNMS
    const T1u = T1q + T1t;
    const T2v = T1t - T1q;
    const T1w = KP634393284 * T1r + KP773010453 * T1s; // FMA
    const T1x = KP773010453 * T1p - KP634393284 * T1o; // FNMS
    const T1y = T1w - T1x;
    const T2z = T1x + T1w;

    out[12] = T1n - T1u;
    out[19] = T2z - T2A;
    out[3] = T1n + T1u;
    out[28] = T2z + T2A;
    out[11] = T1v - T1y;
    out[20] = T2v - T2y;
    out[4] = T1v + T1y;
    out[27] = T2v + T2y;
  }
  {
    const T1D = T1z - T1C;
    const T1K = T1G - T1J;
    const T1L = T1D + T1K;
    const T21 = T1D - T1K;
    const T2E = T2C - T2D;
    const T2F = T26 + T27;
    const T2G = T2E - T2F;
    const T2I = T2F + T2E;

    const T1S = KP881921264 * T1O + KP471396736 * T1R; // FMA
    const T1Z = KP881921264 * T1V + KP471396736 * T1Y; // FMA
    const T20 = T1S - T1Z;
    const T2H = T1S + T1Z;
    const T22 = KP881921264 * T1Y - KP471396736 * T1V; // FNMS
    const T23 = KP881921264 * T1R - KP471396736 * T1O; // FNMS
    const T24 = T22 - T23;
    const T2B = T23 + T22;

    out[13] = T1L - T20;
    out[18] = T2B - T2G;
    out[2] = T1L + T20;
    out[29] = T2B + T2G;
    out[10] = T21 - T24;
    out[21] = T2I - T2H;
    out[5] = T21 + T24;
    out[26] = -(T2H + T2I);
  }

  return out;
}

module.exports = { r2cfII_32 };
