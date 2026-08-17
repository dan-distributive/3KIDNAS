'use strict';

// =============================================================================
// hb2_20.js -- faithful JS port of rdft/scalar/r2cb/hb2_20.c (non-FMA
// branch). Alternate-codegen EVEN-radix (r=20) backward twiddle codelet --
// same role as hb2_16.js, radix 20 (see hf2_20.js's header for the
// forward-direction analogue). twinstr trig-generates W^1, W^3, W^9, W^19.
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hb2_16.js -- note the
// cross-indexed cr[k]/ci[k'] pairing pattern, transcribed literally.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.250000000000000000000000000000000000000000000;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;

function hb2_20(cr, ci, Wc, Ws) {
  const TD = Wc[1], TG = Ws[1], TE = Wc[3], TH = Ws[3];
  const TF = TD * TE, T1s = TG * TE, TI = TG * TH, T1r = TD * TH;
  const TJ = TF + TI, T1t = T1r - T1s, T27 = T1r + T1s, T25 = TF - TI;
  const T1T = Ws[9], T1R = Wc[9];
  const T1U = TH * T1T, T2l = TD * T1T, T1Z = TE * T1T, T2i = TG * T1T;
  const T1S = TE * T1R, T2m = TG * T1R, T20 = TH * T1R, T2h = TD * T1R;

  const T1V = T1S + T1U;
  const T2j = T2h - T2i;
  const T2Z = T1Z + T20;
  const T21 = T1Z - T20;
  const T2X = T1S - T1U;
  const T2T = T2l - T2m;
  const T2n = T2l + T2m;
  const T2P = T2h + T2i;

  const T3T = TJ * T1T, T3U = T1t * T1R;
  const T3V = T3T - T3U, T41 = T3T + T3U;
  const T3P = TJ * T1R, T3Q = T1t * T1T;
  const T3R = T3P + T3Q, T3X = T3P - T3Q;

  const T26 = T25 * T1R, T28 = T27 * T1T;
  const T29 = T26 + T28;
  const T2a = T25 * T1T, T2b = T27 * T1R;
  const T2c = T2a - T2b;
  const T4H = T26 - T28;
  const T4L = T2a + T2b;
  const T1L = Wc[19], T1M = Ws[19];
  const T1N = TD * T1L + TG * T1M; // FMA
  const T2d = T29 * T1L + T2c * T1M; // FMA
  const T4R = TJ * T1M - T1t * T1L; // FNMS
  const T1P = TD * T1M - TG * T1L; // FNMS
  const T4P = TJ * T1L + T1t * T1M; // FMA
  const T49 = T25 * T1M - T27 * T1L; // FNMS
  const T2N = TE * T1M - TH * T1L; // FNMS
  const T2f = T29 * T1M - T2c * T1L; // FNMS
  const T47 = T25 * T1L + T27 * T1M; // FMA
  const T2L = TE * T1L + TH * T1M; // FMA

  const T1 = cr[0];
  const T2 = ci[9];
  const T3 = T1 + T2;
  const T3g = T1 - T2;
  const T1A = ci[14];
  const T1B = cr[15];
  const T1C = T1A - T1B;
  const T3h = T1A + T1B;

  const T4 = cr[5];
  const T5 = ci[4];
  const T6 = T4 + T5;
  const T3D = T4 - T5;
  const T1x = ci[19];
  const T1y = cr[10];
  const T1z = T1x - T1y;
  const T3C = T1x + T1y;

  const T7 = T3 + T6;
  const T4i = T3g - T3h;
  const T4x = T3D + T3C;
  const TK = T3 - T6;
  const T1D = T1z - T1C;
  const T3i = T3g + T3h;
  const T3E = T3C - T3D;
  const T2D = T1z + T1C;

  const T8 = cr[4];
  const T9 = ci[5];
  const Ta = T8 + T9;
  const T3j = T8 - T9;
  const TY = ci[10];
  const TZ = cr[19];
  const T10 = TY - TZ;
  const T3k = TY + TZ;

  const Tb = cr[9];
  const Tc = ci[0];
  const Td = Tb + Tc;
  const T32 = Tb - Tc;
  const TV = ci[15];
  const TW = cr[14];
  const TX = TV - TW;
  const T31 = TV + TW;

  const Te = Ta + Td;
  const T4b = T3j - T3k;
  const T4m = T32 + T31;
  const TL = Ta - Td;
  const T11 = TX - T10;
  const T33 = T31 - T32;
  const T3l = T3j + T3k;
  const T2s = TX + T10;

  const Tu = ci[7];
  const Tv = cr[2];
  const Tw = Tu + Tv;
  const T3t = Tu - Tv;
  const Tx = ci[2];
  const Ty = cr[7];
  const Tz = Tx + Ty;
  const T3b = Tx - Ty;

  const T1h = ci[17];
  const T1i = cr[12];
  const T1j = T1h - T1i;
  const T3c = T1h + T1i;
  const T1k = ci[12];
  const T1l = cr[17];
  const T1m = T1k - T1l;
  const T3u = T1k + T1l;

  const TA = Tw + Tz;
  const T4f = T3t + T3u;
  const T4q = T3b - T3c;
  const TP = Tw - Tz;
  const T1n = T1j - T1m;
  const T3d = T3b + T3c;
  const T3v = T3t - T3u;
  const T2w = T1j + T1m;

  const Tf = ci[3];
  const Tg = cr[6];
  const Th = Tf + Tg;
  const T3m = Tf - Tg;
  const T15 = ci[18];
  const T16 = cr[11];
  const T17 = T15 - T16;
  const T3n = T15 + T16;

  const Ti = cr[1];
  const Tj = ci[8];
  const Tk = Ti + Tj;
  const T34 = Ti - Tj;
  const T12 = ci[13];
  const T13 = cr[16];
  const T14 = T12 - T13;
  const T35 = T12 + T13;

  const Tl = Th + Tk;
  const T4c = T3m - T3n;
  const T4n = T34 - T35;
  const TM = Th - Tk;
  const T18 = T14 - T17;
  const T36 = T34 + T35;
  const T3o = T3m + T3n;
  const T2t = T14 + T17;

  const Tn = cr[8];
  const To = ci[1];
  const Tp = Tn + To;
  const T3q = Tn - To;
  const T1d = ci[16];
  const T1e = cr[13];
  const T1f = T1d - T1e;
  const T3r = T1d + T1e;

  const Tq = ci[6];
  const Tr = cr[3];
  const Ts = Tq + Tr;
  const T39 = Tq - Tr;
  const T1a = ci[11];
  const T1b = cr[18];
  const T1c = T1a - T1b;
  const T38 = T1a + T1b;

  const Tt = Tp + Ts;
  const T4e = T3q + T3r;
  const T4p = T39 + T38;
  const TO = Tp - Ts;
  const T1g = T1c - T1f;
  const T3a = T38 - T39;
  const T3s = T3q - T3r;
  const T2v = T1c + T1f;

  const T19 = T11 - T18;
  const T3L = T3l - T3o;
  const T3M = T3s - T3v;
  const T1o = T1g - T1n;
  const T2x = T2v - T2w;
  const T4C = T4e - T4f;
  const T4B = T4b - T4c;
  const T2u = T2s - T2t;
  const T1v = TO - TP;
  const T4r = T4p - T4q;
  const T4o = T4m - T4n;
  const T1u = TL - TM;
  const T2H = Te - Tl;
  const T37 = T33 + T36;
  const T2I = Tt - TA;
  const T3e = T3a + T3d;
  const T3p = T3l + T3o;
  const T3w = T3s + T3v;
  const T3x = T3p + T3w;
  const Tm = Te + Tl;
  const TB = Tt + TA;
  const TC = Tm + TB;
  const T4u = T4m + T4n;
  const T4v = T4p + T4q;
  const T4y = T4u + T4v;
  const T2A = T2s + T2t;
  const T2B = T2v + T2w;
  const T2E = T2A + T2B;
  const T1E = T11 + T18;
  const T1F = T1g + T1n;
  const T1G = T1E + T1F;
  const T4d = T4b + T4c;
  const T4g = T4e + T4f;
  const T4j = T4d + T4g;
  const T3F = T33 - T36;
  const T3G = T3a - T3d;
  const T3H = T3F + T3G;
  const TN = TL + TM;
  const TQ = TO + TP;
  const TR = TN + TQ;

  const outCr = new Float64Array(20), outCi = new Float64Array(20);
  outCr[0] = T7 + TC;
  outCi[0] = T2D + T2E;
  {
    const T2k = TK + TR;
    const T2o = T1D + T1G;
    outCr[10] = T2j * T2k - T2n * T2o; // FNMS
    outCi[10] = T2n * T2k + T2j * T2o; // FMA
    const T4T = T4i + T4j;
    const T4U = T4x + T4y;
    outCr[5] = T29 * T4T - T2c * T4U; // FNMS
    outCi[5] = T2c * T4T + T29 * T4U; // FMA
  }
  {
    const T48 = T3i + T3x;
    const T4a = T3E + T3H;
    outCr[15] = T47 * T48 - T49 * T4a; // FNMS
    outCi[15] = T49 * T48 + T47 * T4a; // FMA
  }
  {
    const T2y = KP951056516 * T2u + KP587785252 * T2x; // FMA
    const T2J = KP951056516 * T2H + KP587785252 * T2I; // FMA
    const T2V = KP587785252 * T2H - KP951056516 * T2I; // FNMS
    const T2R = KP587785252 * T2u - KP951056516 * T2x; // FNMS

    const T2C = KP559016994 * (T2A - T2B);
    const T2F = T2D - KP250000000 * T2E; // FNMS
    const T2G = T2C + T2F;
    const T2U = T2F - T2C;
    const T2p = KP559016994 * (Tm - TB);
    const T2q = T7 - KP250000000 * TC; // FNMS
    const T2r = T2p + T2q;
    const T2Q = T2q - T2p;

    const T2z = T2r + T2y;
    const T2K = T2G - T2J;
    outCr[4] = T25 * T2z - T27 * T2K; // FNMS
    outCi[4] = T27 * T2z + T25 * T2K; // FMA
    const T2Y = T2Q - T2R;
    const T30 = T2V + T2U;
    outCr[12] = T2X * T2Y - T2Z * T30; // FNMS
    outCi[12] = T2Z * T2Y + T2X * T30; // FMA

    const T2M = T2r - T2y;
    const T2O = T2J + T2G;
    outCr[16] = T2L * T2M - T2N * T2O; // FNMS
    outCi[16] = T2N * T2M + T2L * T2O; // FMA
    const T2S = T2Q + T2R;
    const T2W = T2U - T2V;
    outCr[8] = T2P * T2S - T2T * T2W; // FNMS
    outCi[8] = T2T * T2S + T2P * T2W; // FMA
  }
  {
    const T4s = KP951056516 * T4o + KP587785252 * T4r; // FMA
    const T4D = KP951056516 * T4B + KP587785252 * T4C; // FMA
    const T4N = KP587785252 * T4B - KP951056516 * T4C; // FNMS
    const T4I = KP587785252 * T4o - KP951056516 * T4r; // FNMS

    const T4w = KP559016994 * (T4u - T4v);
    const T4z = T4x - KP250000000 * T4y; // FNMS
    const T4A = T4w + T4z;
    const T4M = T4z - T4w;
    const T4h = KP559016994 * (T4d - T4g);
    const T4k = T4i - KP250000000 * T4j; // FNMS
    const T4l = T4h + T4k;
    const T4J = T4k - T4h;

    const T4t = T4l - T4s;
    const T4E = T4A + T4D;
    outCr[1] = TD * T4t - TG * T4E; // FNMS
    outCi[1] = TG * T4t + TD * T4E; // FMA
    const T4Q = T4J - T4I;
    const T4S = T4M + T4N;
    outCr[17] = T4P * T4Q - T4R * T4S; // FNMS
    outCi[17] = T4R * T4Q + T4P * T4S; // FMA

    const T4F = T4s + T4l;
    const T4G = T4A - T4D;
    outCr[9] = T1R * T4F - T1T * T4G; // FNMS
    outCi[9] = T1T * T4F + T1R * T4G; // FMA
    const T4K = T4I + T4J;
    const T4O = T4M - T4N;
    outCr[13] = T4H * T4K - T4L * T4O; // FNMS
    outCi[13] = T4L * T4K + T4H * T4O; // FMA
  }
  {
    const T1p = KP587785252 * T19 - KP951056516 * T1o; // FNMS
    const T1w = KP587785252 * T1u - KP951056516 * T1v; // FNMS
    const T22 = KP951056516 * T1u + KP587785252 * T1v; // FMA
    const T1X = KP951056516 * T19 + KP587785252 * T1o; // FMA

    const T1H = T1D - KP250000000 * T1G; // FNMS
    const T1I = KP559016994 * (T1E - T1F);
    const T1J = T1H - T1I;
    const T23 = T1I + T1H;
    const TS = TK - KP250000000 * TR; // FNMS
    const TT = KP559016994 * (TN - TQ);
    const TU = TS - TT;
    const T1W = TT + TS;

    const T1q = TU - T1p;
    const T1K = T1w + T1J;
    outCr[2] = TJ * T1q - T1t * T1K; // FNMS
    outCi[2] = T1t * T1q + TJ * T1K; // FMA
    const T2e = T1W + T1X;
    const T2g = T23 - T22;
    outCr[14] = T2d * T2e - T2f * T2g; // FNMS
    outCi[14] = T2f * T2e + T2d * T2g; // FMA

    const T1O = TU + T1p;
    const T1Q = T1J - T1w;
    outCr[18] = T1N * T1O - T1P * T1Q; // FNMS
    outCi[18] = T1P * T1O + T1N * T1Q; // FMA
    const T1Y = T1W - T1X;
    const T24 = T22 + T23;
    outCr[6] = T1V * T1Y - T21 * T24; // FNMS
    outCi[6] = T21 * T1Y + T1V * T24; // FMA
  }
  {
    const T3f = KP587785252 * T37 - KP951056516 * T3e; // FNMS
    const T3N = KP587785252 * T3L - KP951056516 * T3M; // FNMS
    const T43 = KP951056516 * T3L + KP587785252 * T3M; // FMA
    const T3Z = KP951056516 * T37 + KP587785252 * T3e; // FMA

    const T3I = T3E - KP250000000 * T3H; // FNMS
    const T3J = KP559016994 * (T3F - T3G);
    const T3K = T3I - T3J;
    const T42 = T3J + T3I;
    const T3y = T3i - KP250000000 * T3x; // FNMS
    const T3z = KP559016994 * (T3p - T3w);
    const T3A = T3y - T3z;
    const T3Y = T3z + T3y;

    const T3B = T3f + T3A;
    const T3O = T3K - T3N;
    outCr[3] = TE * T3B - TH * T3O; // FNMS
    outCi[3] = TH * T3B + TE * T3O; // FMA
    const T45 = T3Z + T3Y;
    const T46 = T42 - T43;
    outCr[19] = T1L * T45 - T1M * T46; // FNMS
    outCi[19] = T1M * T45 + T1L * T46; // FMA

    const T3S = T3A - T3f;
    const T3W = T3K + T3N;
    outCr[7] = T3R * T3S - T3V * T3W; // FNMS
    outCi[7] = T3V * T3S + T3R * T3W; // FMA
    const T40 = T3Y - T3Z;
    const T44 = T42 + T43;
    outCr[11] = T3X * T40 - T41 * T44; // FNMS
    outCi[11] = T41 * T40 + T3X * T44; // FMA
  }

  return [outCr, outCi];
}

module.exports = { hb2_20 };
