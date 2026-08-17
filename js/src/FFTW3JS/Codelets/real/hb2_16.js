'use strict';

// =============================================================================
// hb2_16.js -- faithful JS port of rdft/scalar/r2cb/hb2_16.c (non-FMA
// branch). Alternate-codegen EVEN-radix (r=16) backward twiddle codelet --
// same role as hb2_8.js, larger radix (see hf2_16.js's header for the
// forward-direction analogue). twinstr trig-generates W^1, W^3, W^9, W^15.
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hb2_8.js -- note the
// cross-indexed cr[k]/ci[k'] pairing pattern, transcribed literally.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP382683432 = 0.382683432365089771728459984030398866761344562;
const KP923879532 = 0.923879532511286756128183189396788286822416626;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function hb2_16(cr, ci, Wc, Ws) {
  const Tv = Wc[1], Ty = Ws[1], T1l = Wc[3], T1n = Ws[3];
  const T1m = Tv * T1l, T1s = Ty * T1l, T1o = Ty * T1n, T1r = Tv * T1n;
  const T1p = T1m + T1o, T1t = T1r - T1s, T27 = T1r + T1s, T25 = T1m - T1o;
  const Tz = Ws[9], Tw = Wc[9];
  const TA = Ty * Tz, T1J = T1l * Tz, T15 = Tv * Tz, T1G = T1n * Tz;
  const Tx = Tv * Tw, T1K = T1n * Tw, T16 = Ty * Tw, T1F = T1l * Tw;

  const TB = Tx - TA;
  const T21 = T1J + T1K;
  const T1P = T15 - T16;
  const T1H = T1F + T1G;
  const T1X = T1F - T1G;
  const T17 = T15 + T16;
  const T1L = T1J - T1K;
  const T1N = Tx + TA;
  const T1v = Wc[15], T1w = Ws[15];
  const T1x = Tv * T1v + Ty * T1w; // FMA
  const T1B = Tv * T1w - Ty * T1v; // FNMS

  const T2D = T25 * Tz, T2E = T27 * Tw;
  const T2F = T2D + T2E, T2T = T2D - T2E;
  const T29 = T25 * Tw, T2a = T27 * Tz;
  const T2b = T29 - T2a, T2R = T29 + T2a;

  const T3h = T1p * Tz, T3i = T1t * Tw;
  const T3j = T3h + T3i, T3x = T3h - T3i;
  const T33 = T1p * Tw, T34 = T1t * Tz;
  const T35 = T33 - T34, T3t = T33 + T34;

  const T1 = cr[0];
  const T2 = ci[7];
  const T3 = T1 + T2;
  const T2c = T1 - T2;
  const T1c = ci[11];
  const T1d = cr[12];
  const T1e = T1c - T1d;
  const T2d = T1c + T1d;

  const T4 = cr[4];
  const T5 = ci[3];
  const T6 = T4 + T5;
  const T2G = T4 - T5;
  const T19 = ci[15];
  const T1a = cr[8];
  const T1b = T19 - T1a;
  const T2H = T19 + T1a;

  const T7 = T3 + T6;
  const T36 = T2c + T2d;
  const T3k = T2H - T2G;
  const TC = T3 - T6;
  const T1f = T1b - T1e;
  const T2e = T2c - T2d;
  const T2I = T2G + T2H;
  const T1Q = T1b + T1e;

  const T8 = cr[2];
  const T9 = ci[5];
  const Ta = T8 + T9;
  const T2f = T8 - T9;
  const TG = ci[13];
  const TH = cr[10];
  const TI = TG - TH;
  const T2g = TG + TH;

  const Tb = ci[1];
  const Tc = cr[6];
  const Td = Tb + Tc;
  const T2i = Tb - Tc;
  const TD = ci[9];
  const TE = cr[14];
  const TF = TD - TE;
  const T2j = TD + TE;

  const Te = Ta + Td;
  const TJ = TF - TI;
  const T1R = TI + TF;
  const T18 = Ta - Td;

  const T2J = T2f + T2g;
  const T2K = T2i + T2j;
  const T2L = KP707106781 * (T2J - T2K);
  const T37 = KP707106781 * (T2J + T2K);
  const T2h = T2f - T2g;
  const T2k = T2i - T2j;
  const T2l = KP707106781 * (T2h + T2k);
  const T3l = KP707106781 * (T2h - T2k);

  const Tg = cr[1];
  const Th = ci[6];
  const Ti = Tg + Th;
  const T2x = Tg - Th;
  const TP = ci[10];
  const TQ = cr[13];
  const TR = TP - TQ;
  const T2y = TP + TQ;

  const Tj = cr[5];
  const Tk = ci[2];
  const Tl = Tj + Tk;
  const T2u = Tj - Tk;
  const TM = ci[14];
  const TN = cr[9];
  const TO = TM - TN;
  const T2v = TM + TN;

  const Tm = Ti + Tl;
  const T1T = TO + TR;
  const TL = Ti - Tl;
  const TS = TO - TR;
  const TT = TL - TS;
  const T1h = TL + TS;

  const T2w = T2u + T2v;
  const T2z = T2x - T2y;
  const T2A = KP923879532 * T2w + KP382683432 * T2z; // FMA
  const T2N = KP923879532 * T2z - KP382683432 * T2w; // FNMS
  const T39 = T2x + T2y;
  const T3a = T2v - T2u;
  const T3b = KP382683432 * T39 - KP923879532 * T3a; // FNMS
  const T3n = KP382683432 * T3a + KP923879532 * T39; // FMA

  const Tn = ci[0];
  const To = cr[7];
  const Tp = Tn + To;
  const T2q = Tn - To;
  const TY = ci[12];
  const TZ = cr[11];
  const T10 = TY - TZ;
  const T2r = TY + TZ;

  const Tq = cr[3];
  const Tr = ci[4];
  const Ts = Tq + Tr;
  const T2n = Tq - Tr;
  const TV = ci[8];
  const TW = cr[15];
  const TX = TV - TW;
  const T2o = TV + TW;

  const Tt = Tp + Ts;
  const T1U = TX + T10;
  const TU = Tp - Ts;
  const T11 = TX - T10;
  const T12 = TU + T11;
  const T1i = T11 - TU;

  const T2p = T2n - T2o;
  const T2s = T2q - T2r;
  const T2t = KP923879532 * T2p - KP382683432 * T2s; // FNMS
  const T2O = KP382683432 * T2p + KP923879532 * T2s; // FMA
  const T3c = T2q + T2r;
  const T3d = T2n + T2o;
  const T3e = KP382683432 * T3c - KP923879532 * T3d; // FNMS
  const T3o = KP382683432 * T3d + KP923879532 * T3c; // FMA

  const outCr = new Float64Array(16), outCi = new Float64Array(16);

  {
    const Tf = T7 + Te;
    const Tu = Tm + Tt;
    const T1O = Tf - Tu;
    const T1S = T1Q + T1R;
    const T1V = T1T + T1U;
    const T1W = T1S - T1V;
    outCr[0] = Tf + Tu;
    outCi[0] = T1S + T1V;
    outCr[8] = T1N * T1O - T1P * T1W; // FNMS
    outCi[8] = T1P * T1O + T1N * T1W; // FMA
  }
  {
    const T38 = T36 - T37;
    const T3f = T3b + T3e;
    const T3g = T38 - T3f;
    const T3r = T38 + T3f;
    const T3m = T3k + T3l;
    const T3p = T3n - T3o;
    const T3q = T3m - T3p;
    const T3s = T3m + T3p;

    outCr[11] = T35 * T3g - T3j * T3q; // FNMS
    outCi[11] = T3j * T3g + T35 * T3q; // FMA
    outCr[3] = T1l * T3r - T1n * T3s; // FNMS
    outCi[3] = T1n * T3r + T1l * T3s; // FMA
  }
  {
    const T3u = T36 + T37;
    const T3v = T3n + T3o;
    const T3w = T3u - T3v;
    const T3B = T3u + T3v;
    const T3y = T3k - T3l;
    const T3z = T3b - T3e;
    const T3A = T3y + T3z;
    const T3C = T3y - T3z;

    outCr[7] = T3t * T3w - T3x * T3A; // FNMS
    outCi[7] = T3t * T3A + T3x * T3w; // FMA
    outCr[15] = T1v * T3B - T1w * T3C; // FNMS
    outCi[15] = T1v * T3C + T1w * T3B; // FMA
  }
  {
    const TK = TC + TJ;
    const T13 = KP707106781 * (TT + T12);
    const T14 = TK - T13;
    const T1q = TK + T13;
    const T1g = T18 + T1f;
    const T1j = KP707106781 * (T1h + T1i);
    const T1k = T1g - T1j;
    const T1u = T1g + T1j;

    outCr[10] = TB * T14 - T17 * T1k; // FNMS
    outCi[10] = T17 * T14 + TB * T1k; // FMA
    outCr[2] = T1p * T1q - T1t * T1u; // FNMS
    outCi[2] = T1t * T1q + T1p * T1u; // FMA
  }
  {
    const T1y = TC - TJ;
    const T1z = KP707106781 * (T1i - T1h);
    const T1A = T1y - T1z;
    const T1I = T1y + T1z;
    const T1C = T1f - T18;
    const T1D = KP707106781 * (TT - T12);
    const T1E = T1C - T1D;
    const T1M = T1C + T1D;

    outCr[14] = T1x * T1A - T1B * T1E; // FNMS
    outCi[14] = T1x * T1E + T1B * T1A; // FMA
    outCr[6] = T1H * T1I - T1L * T1M; // FNMS
    outCi[6] = T1H * T1M + T1L * T1I; // FMA
  }
  {
    const T2m = T2e - T2l;
    const T2B = T2t - T2A;
    const T2C = T2m - T2B;
    const T2S = T2m + T2B;
    const T2M = T2I - T2L;
    const T2P = T2N - T2O;
    const T2Q = T2M - T2P;
    const T2U = T2M + T2P;

    outCr[13] = T2b * T2C - T2F * T2Q; // FNMS
    outCi[13] = T2F * T2C + T2b * T2Q; // FMA
    outCr[5] = T2R * T2S - T2T * T2U; // FNMS
    outCi[5] = T2T * T2S + T2R * T2U; // FMA
  }
  {
    const T2V = T2e + T2l;
    const T2W = T2N + T2O;
    const T2X = T2V - T2W;
    const T31 = T2V + T2W;
    const T2Y = T2I + T2L;
    const T2Z = T2A + T2t;
    const T30 = T2Y - T2Z;
    const T32 = T2Y + T2Z;

    outCr[9] = Tw * T2X - Tz * T30; // FNMS
    outCi[9] = Tw * T30 + Tz * T2X; // FMA
    outCr[1] = Tv * T31 - Ty * T32; // FNMS
    outCi[1] = Tv * T32 + Ty * T31; // FMA
  }
  {
    const T1Y = T7 - Te;
    const T1Z = T1U - T1T;
    const T20 = T1Y - T1Z;
    const T26 = T1Y + T1Z;
    const T22 = T1Q - T1R;
    const T23 = Tm - Tt;
    const T24 = T22 - T23;
    const T28 = T23 + T22;

    outCr[12] = T1X * T20 - T21 * T24; // FNMS
    outCi[12] = T21 * T20 + T1X * T24; // FMA
    outCr[4] = T25 * T26 - T27 * T28; // FNMS
    outCi[4] = T27 * T26 + T25 * T28; // FMA
  }

  return [outCr, outCi];
}

module.exports = { hb2_16 };
