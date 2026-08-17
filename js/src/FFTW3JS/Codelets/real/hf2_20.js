'use strict';

// =============================================================================
// hf2_20.js -- faithful JS port of rdft/scalar/r2cf/hf2_20.c (non-FMA
// branch). Alternate-codegen EVEN-radix (r=20) twiddle codelet -- same role
// as hf2_16.js, radix 20. twinstr trig-generates W^1, W^3, W^9, W^19; the
// rest are derived via complex products, exactly matching the C source's
// operation order. Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP250000000 = 0.250000000000000000000000000000000000000000000;
const KP559016994 = 0.559016994374947424102293417182819058860154590;

function hf2_20(cr, ci, Wc, Ws) {
  const T2 = Wc[1], T5 = Ws[1], Tg = Wc[3], Ti = Ws[3];
  const Th = T2 * Tg, Tn = T5 * Tg, Tj = T5 * Ti, Tm = T2 * Ti;
  const Tk = Th - Tj, To = Tm + Tn, T1h = Tm - Tn, T1f = Th + Tj;
  const T6 = Ws[9], T3 = Wc[9];
  const T7 = T5 * T6, T16 = Tg * T6, Ta = T2 * T6, T13 = Ti * T6;
  const T4 = T2 * T3, T17 = Ti * T3, Tb = T5 * T3, T12 = Tg * T3;

  const T8 = T4 - T7;
  const T14 = T12 + T13;
  const T1Q = T16 + T17;
  const Tc = Ta + Tb;
  const T1O = T12 - T13;
  const T1v = Ta - Tb;
  const T18 = T16 - T17;
  const T1t = T4 + T7;

  const T1l = T1f * T6, T1m = T1h * T3;
  const T1n = T1l + T1m, T24 = T1l - T1m;
  const T1g = T1f * T3, T1i = T1h * T6;
  const T1j = T1g - T1i, T22 = T1g + T1i;

  const Tl = Tk * T3, Tp = To * T6;
  const Tq = Tl + Tp;
  const Ts = Tk * T6, Tt = To * T3;
  const Tu = Ts - Tt;
  const T1E = Tl - Tp;
  const T1G = Ts + Tt;
  const Tx = Wc[19], Ty = Ws[19];
  const Tz = Tk * Tx + To * Ty; // FMA
  const TJ = Tq * Tx + Tu * Ty; // FMA
  const T1Z = T1f * Ty - T1h * Tx; // FNMS
  const TB = Tk * Ty - To * Tx; // FNMS
  const T1X = T1f * Tx + T1h * Ty; // FMA
  const T1A = T2 * Ty - T5 * Tx; // FNMS
  const TZ = Tg * Ty - Ti * Tx; // FNMS
  const TL = Tq * Ty - Tu * Tx; // FNMS
  const T1y = T2 * Tx + T5 * Ty; // FMA
  const TX = Tg * Tx + Ti * Ty; // FMA

  const T1 = cr[0];
  const T47 = ci[0];
  const T9 = cr[10], Td = ci[10];
  const Te = T8 * T9 + Tc * Td; // FMA
  const T46 = T8 * Td - Tc * T9; // FNMS

  const Tr = cr[5], Tv = ci[5];
  const Tw = Tq * Tr + Tu * Tv; // FMA
  const T2H = Tq * Tv - Tu * Tr; // FNMS
  const TA = cr[15], TC = ci[15];
  const TD = Tz * TA + TB * TC; // FMA
  const T2I = Tz * TC - TB * TA; // FNMS

  const Tf = T1 + Te;
  const TE = Tw + TD;
  const TF = Tf - TE;
  const T2b = Tf + TE;
  const T4B = T47 - T46;
  const T4C = Tw - TD;
  const T4D = T4B - T4C;
  const T4M = T4C + T4B;

  const T2G = T1 - Te;
  const T2J = T2H - T2I;
  const T2K = T2G - T2J;
  const T3r = T2G + T2J;
  const T48 = T46 + T47;
  const T49 = T2H + T2I;
  const T4a = T48 - T49;
  const T4m = T49 + T48;

  const T1u = cr[8], T1w = ci[8];
  const T1x = T1t * T1u + T1v * T1w; // FMA
  const T2Z = T1t * T1w - T1v * T1u; // FNMS
  const T1z = cr[18], T1B = ci[18];
  const T1C = T1y * T1z + T1A * T1B; // FMA
  const T30 = T1y * T1B - T1A * T1z; // FNMS

  const T1D = T1x + T1C;
  const T3A = T2Z + T30;
  const T2u = T1x - T1C;
  const T31 = T2Z - T30;

  const T1Y = cr[17], T20 = ci[17];
  const T21 = T1X * T1Y + T1Z * T20; // FMA
  const T2A = T1X * T20 - T1Z * T1Y; // FNMS
  const T23 = cr[7], T25 = ci[7];
  const T26 = T22 * T23 + T24 * T25; // FMA
  const T2B = T22 * T25 - T24 * T23; // FNMS

  const T27 = T21 + T26;
  const T3D = T2A + T2B;
  const T2C = T2A - T2B;
  const T37 = T21 - T26;

  const T1F = cr[13], T1H = ci[13];
  const T1I = T1E * T1F + T1G * T1H; // FMA
  const T2v = T1E * T1H - T1G * T1F; // FNMS
  const T1J = cr[3], T1K = ci[3];
  const T1L = Tg * T1J + Ti * T1K; // FMA
  const T2w = Tg * T1K - Ti * T1J; // FNMS

  const T1M = T1I + T1L;
  const T3B = T2v + T2w;
  const T2x = T2v - T2w;
  const T32 = T1I - T1L;

  const T1P = cr[12], T1R = ci[12];
  const T1S = T1O * T1P + T1Q * T1R; // FMA
  const T34 = T1O * T1R - T1Q * T1P; // FNMS
  const T1T = cr[2], T1U = ci[2];
  const T1V = T1f * T1T + T1h * T1U; // FMA
  const T35 = T1f * T1U - T1h * T1T; // FNMS

  const T1W = T1S + T1V;
  const T3E = T34 + T35;
  const T2z = T1S - T1V;
  const T36 = T34 - T35;

  const T1N = T1D - T1M;
  const T28 = T1W - T27;
  const T29 = T1N + T28;
  const T3C = T3A - T3B;
  const T3F = T3D - T3E;
  const T43 = T3F - T3C;
  const T3X = T3A + T3B;
  const T3Y = T3E + T3D;
  const T4o = T3X + T3Y;
  const T2f = T1D + T1M;
  const T2g = T1W + T27;
  const T2h = T2f + T2g;
  const T2y = T2u - T2x;
  const T2D = T2z - T2C;
  const T2E = T2y + T2D;
  const T3g = T31 - T32;
  const T3h = T36 - T37;
  const T4z = T3g + T3h;
  const T3n = T2u + T2x;
  const T3o = T2z + T2C;
  const T3p = T3n + T3o;
  const T33 = T31 + T32;
  const T38 = T36 + T37;
  const T4K = T33 + T38;

  const TG = cr[4], TH = ci[4];
  const TI = Tk * TG + To * TH; // FMA
  const T2O = Tk * TH - To * TG; // FNMS
  const TK = cr[14], TM = ci[14];
  const TN = TJ * TK + TL * TM; // FMA
  const T2P = TJ * TM - TL * TK; // FNMS

  const TO = TI + TN;
  const T3H = T2O + T2P;
  const T2j = TI - TN;
  const T2Q = T2O - T2P;

  const T1c = cr[1], T1d = ci[1];
  const T1e = T2 * T1c + T5 * T1d; // FMA
  const T2p = T2 * T1d - T5 * T1c; // FNMS
  const T1k = cr[11], T1o = ci[11];
  const T1p = T1j * T1k + T1n * T1o; // FMA
  const T2q = T1j * T1o - T1n * T1k; // FNMS

  const T1q = T1e + T1p;
  const T3L = T2p + T2q;
  const T2r = T2p - T2q;
  const T2T = T1p - T1e;

  const TP = cr[9], TQ = ci[9];
  const TR = T3 * TP + T6 * TQ; // FMA
  const T2k = T3 * TQ - T6 * TP; // FNMS
  const TS = cr[19], TT = ci[19];
  const TU = Tx * TS + Ty * TT; // FMA
  const T2l = Tx * TT - Ty * TS; // FNMS

  const TV = TR + TU;
  const T3I = T2k + T2l;
  const T2m = T2k - T2l;
  const T2R = TR - TU;

  const TY = cr[16], T10 = ci[16];
  const T11 = TX * TY + TZ * T10; // FMA
  const T2U = TX * T10 - TZ * TY; // FNMS
  const T15 = cr[6], T19 = ci[6];
  const T1a = T14 * T15 + T18 * T19; // FMA
  const T2V = T14 * T19 - T18 * T15; // FNMS

  const T1b = T11 + T1a;
  const T3K = T2U + T2V;
  const T2o = T11 - T1a;
  const T2W = T2U - T2V;

  const TW = TO - TV;
  const T1r = T1b - T1q;
  const T1s = TW + T1r;
  const T3J = T3H - T3I;
  const T3M = T3K - T3L;
  const T44 = T3J + T3M;
  const T3U = T3H + T3I;
  const T3V = T3K + T3L;
  const T4n = T3U + T3V;
  const T2c = TO + TV;
  const T2d = T1b + T1q;
  const T2e = T2c + T2d;
  const T2n = T2j - T2m;
  const T2s = T2o - T2r;
  const T2t = T2n + T2s;
  const T3d = T2Q - T2R;
  const T3e = T2W + T2T;
  const T4y = T3d + T3e;
  const T3k = T2j + T2m;
  const T3l = T2o + T2r;
  const T3m = T3k + T3l;
  const T2S = T2Q + T2R;
  const T2X = T2T - T2W;
  const T4J = T2X - T2S;

  const outCr = new Float64Array(20), outCi = new Float64Array(20);

  {
    const T3y = KP559016994 * (T1s - T29);
    const T2a = T1s + T29;
    const T3x = TF - KP250000000 * T2a; // FNMS
    const T3G = T3C + T3F;
    const T3N = T3J - T3M;
    const T3O = KP951056516 * T3G - KP587785252 * T3N; // FNMS
    const T3Q = KP951056516 * T3N + KP587785252 * T3G; // FMA

    outCi[9] = TF + T2a;
    const T3P = T3y + T3x;
    outCi[5] = T3P - T3Q;
    outCr[6] = T3P + T3Q;
    const T3z = T3x - T3y;
    outCr[2] = T3z - T3O;
    outCi[1] = T3z + T3O;
  }
  {
    const T3q = KP559016994 * (T3m - T3p);
    const T3s = T3m + T3p;
    const T3t = T3r - KP250000000 * T3s; // FNMS
    const T3f = T3d - T3e;
    const T3i = T3g - T3h;
    const T3j = KP951056516 * T3f + KP587785252 * T3i; // FMA
    const T3w = KP951056516 * T3i - KP587785252 * T3f; // FNMS

    outCr[5] = T3r + T3s;
    const T3v = T3t - T3q;
    outCi[2] = T3v - T3w;
    outCi[6] = T3w + T3v;
    const T3u = T3q + T3t;
    outCr[1] = T3j + T3u;
    outCr[9] = T3u - T3j;
  }
  {
    const T3R = KP559016994 * (T2e - T2h);
    const T2i = T2e + T2h;
    const T3S = T2b - KP250000000 * T2i; // FNMS
    const T3W = T3U - T3V;
    const T3Z = T3X - T3Y;
    const T40 = KP951056516 * T3W + KP587785252 * T3Z; // FMA
    const T42 = KP951056516 * T3Z - KP587785252 * T3W; // FNMS

    outCr[0] = T2b + T2i;
    const T41 = T3S - T3R;
    outCi[7] = T41 - T42;
    outCr[8] = T41 + T42;
    const T3T = T3R + T3S;
    outCr[4] = T3T - T40;
    outCi[3] = T3T + T40;
  }
  {
    const T2F = KP559016994 * (T2t - T2E);
    const T2L = T2t + T2E;
    const T2M = T2K - KP250000000 * T2L; // FNMS
    const T2Y = T2S + T2X;
    const T39 = T33 - T38;
    const T3a = KP951056516 * T2Y + KP587785252 * T39; // FMA
    const T3b = KP951056516 * T39 - KP587785252 * T2Y; // FNMS

    outCi[4] = T2K + T2L;
    const T3c = T2M - T2F;
    outCr[3] = T3b + T3c;
    outCr[7] = T3c - T3b;
    const T2N = T2F + T2M;
    outCi[0] = T2N - T3a;
    outCi[8] = T3a + T2N;
  }
  {
    const T4e = KP559016994 * (T44 + T43);
    const T45 = T43 - T44;
    const T4f = KP250000000 * T45 + T4a; // FMA
    const T4b = T1r - TW;
    const T4c = T1N - T28;
    const T4d = KP951056516 * T4b - KP587785252 * T4c; // FNMS
    const T4h = KP587785252 * T4b + KP951056516 * T4c; // FMA

    outCr[10] = T45 - T4a;
    const T4i = T4f - T4e;
    outCr[18] = T4h - T4i;
    outCi[17] = T4h + T4i;
    const T4g = T4e + T4f;
    outCr[14] = T4d - T4g;
    outCi[13] = T4d + T4g;
  }
  {
    const T4A = KP559016994 * (T4y - T4z);
    const T4E = T4y + T4z;
    const T4F = T4D - KP250000000 * T4E; // FNMS
    const T4v = T3n - T3o;
    const T4w = T3k - T3l;
    const T4x = KP951056516 * T4v - KP587785252 * T4w; // FNMS
    const T4H = KP951056516 * T4w + KP587785252 * T4v; // FMA

    outCi[14] = T4E + T4D;
    const T4I = T4A + T4F;
    outCi[10] = T4H + T4I;
    outCi[18] = T4I - T4H;
    const T4G = T4A - T4F;
    outCr[13] = T4x + T4G;
    outCr[17] = T4G - T4x;
  }
  {
    const T4r = KP559016994 * (T4n - T4o);
    const T4p = T4n + T4o;
    const T4q = T4m - KP250000000 * T4p; // FNMS
    const T4j = T2c - T2d;
    const T4k = T2f - T2g;
    const T4l = KP587785252 * T4j - KP951056516 * T4k; // FNMS
    const T4t = KP951056516 * T4j + KP587785252 * T4k; // FMA

    outCi[19] = T4p + T4m;
    const T4u = T4r + T4q;
    outCr[16] = T4t - T4u;
    outCi[15] = T4t + T4u;
    const T4s = T4q - T4r;
    outCr[12] = T4l - T4s;
    outCi[11] = T4l + T4s;
  }
  {
    const T4Q = KP559016994 * (T4J + T4K);
    const T4L = T4J - T4K;
    const T4R = KP250000000 * T4L + T4M; // FMA
    const T4N = T2n - T2s;
    const T4O = T2y - T2D;
    const T4P = KP951056516 * T4N + KP587785252 * T4O; // FMA
    const T4T = KP951056516 * T4O - KP587785252 * T4N; // FNMS

    outCr[15] = T4L - T4M;
    const T4U = T4Q + T4R;
    outCi[12] = T4T + T4U;
    outCi[16] = T4U - T4T;
    const T4S = T4Q - T4R;
    outCr[11] = T4P + T4S;
    outCr[19] = T4S - T4P;
  }

  return [outCr, outCi];
}

module.exports = { hf2_20 };
