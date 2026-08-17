'use strict';

// =============================================================================
// r2cbIII_32.js -- faithful JS port of rdft/scalar/r2cb/r2cbIII_32.c
// (non-FMA branch). HC2RIII ("shifted") radix-32 direct codelet -- the
// BACKWARD counterpart of r2cfII_32.js (see that file's header). INPUT
// in[0..31] uses r2cfII_32.js's OUTPUT convention (in[k]=Cr[k] for k=0..15,
// in[31-k]=Ci[k] for k=0..15); OUTPUT out[0..31] (phase p) uses
// r2cfII_32.js's INPUT convention (out[2k]=R0[k], out[2k+1]=R1[k]).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP1_913880671 = 1.913880671464417729871595773960539938965698411;
const KP580569354 = 0.580569354508924735272384751634790549382952557;
const KP942793473 = 0.942793473651995297112775251810508755314920638;
const KP1_763842528 = 1.763842528696710059425513727320776699016885241;
const KP1_546020906 = 1.546020906725473921621813219516939601942082586;
const KP1_268786568 = 1.268786568327290996430343226450986741351374190;
const KP196034280 = 0.196034280659121203988391127777283691722273346;
const KP1_990369453 = 1.990369453344393772489673906218959843150949737;
const KP765366864 = 0.765366864730179543456919968060797733522689125;
const KP1_847759065 = 1.847759065022573512256366378793576573644833252;
const KP1_961570560 = 1.961570560806460898252364472268478073947867462;
const KP390180644 = 0.390180644032256535696569736954044481855383236;
const KP1_111140466 = 1.111140466039204449485661627897065748749874382;
const KP1_662939224 = 1.662939224605090474157576755235811513477121624;
const KP1_414213562 = 1.414213562373095048801688724209698078569671875;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;
const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function r2cbIII_32(inArr) {
  const T1 = inArr[0];
  const T2 = inArr[15];
  const T3 = T1 + T2;
  const Tv = T1 - T2;
  const T1h = inArr[31];
  const T1i = inArr[16];
  const T1j = T1h + T1i;
  const T2h = T1i - T1h;

  const T4 = inArr[8];
  const T5 = inArr[7];
  const T6 = T4 + T5;
  const T1g = T4 - T5;
  const Tw = inArr[23];
  const Tx = inArr[24];
  const Ty = Tw + Tx;
  const T2g = Tw - Tx;

  const T7 = T3 + T6;
  const T2i = T2g + T2h;
  const T2F = T2h - T2g;
  const Tz = Tv - Ty;
  const T1k = T1g + T1j;
  const T1I = T1g - T1j;
  const T1Z = T3 - T6;
  const T1x = Tv + Ty;

  const T8 = inArr[4];
  const T9 = inArr[11];
  const Ta = T8 + T9;
  const TA = T8 - T9;
  const TB = inArr[27];
  const TC = inArr[20];
  const TD = TB + TC;
  const T21 = TB - TC;

  const Tb = inArr[3];
  const Tc = inArr[12];
  const Td = Tb + Tc;
  const TF = Tb - Tc;
  const TG = inArr[28];
  const TH = inArr[19];
  const TI = TG + TH;
  const T20 = TH - TG;

  const Te = Ta + Td;
  const T22 = T20 - T21;
  const T2E = T21 + T20;
  const T2j = Ta - Td;

  const T1d = TA + TD;
  const T1e = TF + TI;
  const T1f = KP707106781 * (T1d - T1e);
  const T1y = KP707106781 * (T1d + T1e);
  const TE = TA - TD;
  const TJ = TF - TI;
  const TK = KP707106781 * (TE + TJ);
  const T1J = KP707106781 * (TE - TJ);

  const Tg = inArr[2];
  const Th = inArr[13];
  const Ti = Tg + Th;
  const TM = Tg - Th;
  const TS = inArr[29];
  const TT = inArr[18];
  const TU = TS + TT;
  const T25 = TS - TT;

  const Tj = inArr[10];
  const Tk = inArr[5];
  const Tl = Tj + Tk;
  const TR = Tj - Tk;
  const TN = inArr[21];
  const TO = inArr[26];
  const TP = TN + TO;
  const T26 = TN - TO;

  const Tm = Ti + Tl;
  const T2B = T26 + T25;
  const TQ = TM - TP;
  const TV = TR + TU;
  const TW = KP923879532 * TQ - KP382683432 * TV; // FNMS
  const T1a = KP382683432 * TQ + KP923879532 * TV; // FMA

  const T1A = TM + TP;
  const T1B = TU - TR;
  const T1C = KP382683432 * T1A - KP923879532 * T1B; // FNMS
  const T1L = KP923879532 * T1A + KP382683432 * T1B; // FMA
  const T24 = Ti - Tl;
  const T27 = T25 - T26;
  const T28 = T24 - T27;
  const T2l = T24 + T27;

  const Tn = inArr[1];
  const To = inArr[14];
  const Tp = Tn + To;
  const TX = Tn - To;
  const T13 = inArr[30];
  const T14 = inArr[17];
  const T15 = T13 + T14;
  const T2a = T14 - T13;

  const Tq = inArr[6];
  const Tr = inArr[9];
  const Ts = Tq + Tr;
  const T12 = Tq - Tr;
  const TY = inArr[25];
  const TZ = inArr[22];
  const T10 = TY + TZ;
  const T2b = TY - TZ;

  const Tt = Tp + Ts;
  const T2A = T2b + T2a;
  const T11 = TX - T10;
  const T16 = T12 - T15;
  const T17 = KP923879532 * T11 + KP382683432 * T16; // FMA
  const T1b = KP923879532 * T16 - KP382683432 * T11; // FNMS

  const T1D = TX + T10;
  const T1E = T12 + T15;
  const T1F = KP382683432 * T1D - KP923879532 * T1E; // FNMS
  const T1M = KP923879532 * T1D + KP382683432 * T1E; // FMA
  const T29 = Tp - Ts;
  const T2c = T2a - T2b;
  const T2d = T29 + T2c;
  const T2m = T2c - T29;

  const out = new Float64Array(32);

  {
    const Tf = T7 + Te;
    const Tu = Tm + Tt;
    const T2L = Tf - Tu;
    const T2M = T2B + T2A;
    const T2N = T2F - T2E;
    const T2O = T2M + T2N;
    out[0] = KP2_000000000 * (Tf + Tu);
    out[16] = KP2_000000000 * (T2N - T2M);
    out[8] = KP1_414213562 * (T2L + T2O);
    out[24] = KP1_414213562 * (T2O - T2L);
  }
  {
    const T2r = T1Z - T22;
    const T2s = KP707106781 * (T2m - T2l);
    const T2t = T2r + T2s;
    const T2x = T2r - T2s;
    const T2u = T2j + T2i;
    const T2v = KP707106781 * (T28 - T2d);
    const T2w = T2u - T2v;
    const T2y = T2v + T2u;

    out[6] = KP1_662939224 * T2t + KP1_111140466 * T2w; // FMA
    out[30] = KP390180644 * T2y - KP1_961570560 * T2x; // FNMS
    out[22] = KP1_662939224 * T2w - KP1_111140466 * T2t; // FNMS
    out[14] = KP390180644 * T2x + KP1_961570560 * T2y; // FMA
  }
  {
    const T2z = T7 - Te;
    const T2C = T2A - T2B;
    const T2D = T2z + T2C;
    const T2J = T2z - T2C;
    const T2G = T2E + T2F;
    const T2H = Tm - Tt;
    const T2I = T2G - T2H;
    const T2K = T2H + T2G;

    out[4] = KP1_847759065 * T2D + KP765366864 * T2I; // FMA
    out[28] = KP765366864 * T2K - KP1_847759065 * T2J; // FNMS
    out[20] = KP1_847759065 * T2I - KP765366864 * T2D; // FNMS
    out[12] = KP765366864 * T2J + KP1_847759065 * T2K; // FMA
  }
  {
    const TL = Tz + TK;
    const T18 = TW + T17;
    const T19 = TL + T18;
    const T1n = TL - T18;
    const T1c = T1a + T1b;
    const T1l = T1f + T1k;
    const T1m = T1c + T1l;
    const T1o = T1c - T1l;

    out[1] = KP1_990369453 * T19 - KP196034280 * T1m; // FNMS
    out[25] = KP1_268786568 * T1o - KP1_546020906 * T1n; // FNMS
    out[17] = -(KP196034280 * T19 + KP1_990369453 * T1m); // -FMA
    out[9] = KP1_268786568 * T1n + KP1_546020906 * T1o; // FMA
  }
  {
    const T1p = Tz - TK;
    const T1q = T1b - T1a;
    const T1r = T1p + T1q;
    const T1v = T1p - T1q;
    const T1s = T1f - T1k;
    const T1t = TW - T17;
    const T1u = T1s - T1t;
    const T1w = T1t + T1s;

    out[5] = KP1_763842528 * T1r + KP942793473 * T1u; // FMA
    out[29] = KP580569354 * T1w - KP1_913880671 * T1v; // FNMS
    out[21] = KP1_763842528 * T1u - KP942793473 * T1r; // FNMS
    out[13] = KP580569354 * T1v + KP1_913880671 * T1w; // FMA
  }
  {
    const T1R = T1x + T1y;
    const T1S = T1L + T1M;
    const T1T = T1R - T1S;
    const T1X = T1R + T1S;
    const T1U = T1J + T1I;
    const T1V = T1C - T1F;
    const T1W = T1U - T1V;
    const T1Y = T1V + T1U;

    out[7] = KP1_546020906 * T1T + KP1_268786568 * T1W; // FMA
    out[31] = KP196034280 * T1Y - KP1_990369453 * T1X; // FNMS
    out[23] = KP1_546020906 * T1W - KP1_268786568 * T1T; // FNMS
    out[15] = KP196034280 * T1X + KP1_990369453 * T1Y; // FMA
  }
  {
    const T23 = T1Z + T22;
    const T2e = KP707106781 * (T28 + T2d);
    const T2f = T23 + T2e;
    const T2p = T23 - T2e;
    const T2k = T2i - T2j;
    const T2n = KP707106781 * (T2l + T2m);
    const T2o = T2k - T2n;
    const T2q = T2n + T2k;

    out[2] = KP1_961570560 * T2f + KP390180644 * T2o; // FMA
    out[26] = KP1_111140466 * T2q - KP1_662939224 * T2p; // FNMS
    out[18] = KP1_961570560 * T2o - KP390180644 * T2f; // FNMS
    out[10] = KP1_111140466 * T2p + KP1_662939224 * T2q; // FMA
  }
  {
    const T1z = T1x - T1y;
    const T1G = T1C + T1F;
    const T1H = T1z + T1G;
    const T1P = T1z - T1G;
    const T1K = T1I - T1J;
    const T1N = T1L - T1M;
    const T1O = T1K - T1N;
    const T1Q = T1N + T1K;

    out[3] = KP1_913880671 * T1H + KP580569354 * T1O; // FMA
    out[27] = KP942793473 * T1Q - KP1_763842528 * T1P; // FNMS
    out[19] = KP1_913880671 * T1O - KP580569354 * T1H; // FNMS
    out[11] = KP942793473 * T1P + KP1_763842528 * T1Q; // FMA
  }

  return out;
}

module.exports = { r2cbIII_32 };
