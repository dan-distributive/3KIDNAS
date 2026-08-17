'use strict';

// =============================================================================
// t2_20.js -- faithful JS port of dft/scalar/codelets/t2_20.c (non-FMA
// branch), FFTW3's "twiddle-log3 / precompute-twiddles" radix-20 twiddle
// codelet. twinstr only trig-generates W^1, W^3, W^9, W^19 (four raw
// pairs); the other 15 needed multiples are DERIVED via one complex
// multiply each (same trick as t2_5/t2_8/t2_10/t2_16.js). Same calling
// convention as every other twiddle codelet here (br/bi/Wc/Ws with
// Composite1D.js's full r-1 = 19 pair table already built) -- only indices
// 1, 3, 9, 19 are actually read.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;

function t2_20(br, bi, Wc, Ws) {
  const T2 = Wc[1], T5 = Ws[1], Tg = Wc[3], Ti = Ws[3];
  const Th = T2 * Tg, Tn = T5 * Tg, Tj = T5 * Ti, Tm = T2 * Ti;
  const Tk = Th - Tj, To = Tm + Tn, T1h = Tm - Tn, T1f = Th + Tj;
  const T6 = Ws[9];
  const T7 = T5 * T6, T16 = Tg * T6, Ta = T2 * T6, T13 = Ti * T6;
  const T3 = Wc[9];
  const T4 = T2 * T3, T17 = Ti * T3, Tb = T5 * T3, T12 = Tg * T3;

  const T8 = T4 - T7, T14 = T12 + T13, T1Q = T16 + T17, Tc = Ta + Tb;
  const T1O = T12 - T13, T1v = Ta - Tb, T18 = T16 - T17, T1t = T4 + T7;

  const T1l = T1f * T6, T1m = T1h * T3, T1n = T1l + T1m, T24 = T1l - T1m;
  const T1g = T1f * T3, T1i = T1h * T6, T1j = T1g - T1i, T22 = T1g + T1i;

  const Tl = Tk * T3, Tp = To * T6, Tq = Tl + Tp;
  const Ts = Tk * T6, Tt = To * T3, Tu = Ts - Tt;
  const T1E = Tl - Tp, T1G = Ts + Tt;

  const Tx = Wc[19], Ty = Ws[19];
  const Tz = Tk * Tx + To * Ty;
  const TJ = Tq * Tx + Tu * Ty;
  const T1Z = T1f * Ty - T1h * Tx;
  const TB = Tk * Ty - To * Tx;
  const T1X = T1f * Tx + T1h * Ty;
  const T1A = T2 * Ty - T5 * Tx;
  const TZ = Tg * Ty - Ti * Tx;
  const TL = Tq * Ty - Tu * Tx;
  const T1y = T2 * Tx + T5 * Ty;
  const TX = Tg * Tx + Ti * Ty;

  const T1 = br[0], T48 = bi[0];
  const T9 = br[10], Td = bi[10];
  const Te = T8 * T9 + Tc * Td;
  const T47 = T8 * Td - Tc * T9;

  const Tr = br[5], Tv = bi[5];
  const Tw = Tq * Tr + Tu * Tv;
  const T2H = Tq * Tv - Tu * Tr;
  const TA = br[15], TC = bi[15];
  const TD = Tz * TA + TB * TC;
  const T2I = Tz * TC - TB * TA;

  const Tf = T1 + Te;
  const TE = Tw + TD;
  const TF = Tf - TE;
  const T2b = Tf + TE;
  const T4y = T48 - T47;
  const T4z = Tw - TD;
  const T4A = T4y - T4z;
  const T4J = T4z + T4y;

  const T2G = T1 - Te;
  const T2J = T2H - T2I;
  const T2K = T2G - T2J;
  const T3r = T2G + T2J;
  const T46 = T2H + T2I;
  const T49 = T47 + T48;
  const T4a = T46 + T49;
  const T4m = T49 - T46;

  const T1u = br[8], T1w = bi[8];
  const T1x = T1t * T1u + T1v * T1w;
  const T2j = T1t * T1w - T1v * T1u;
  const T1z = br[18], T1B = bi[18];
  const T1C = T1y * T1z + T1A * T1B;
  const T2k = T1y * T1B - T1A * T1z;

  const T1D = T1x + T1C;
  const T3A = T2j + T2k;
  const T2l = T2j - T2k;
  const T2W = T1x - T1C;

  const T1Y = br[17], T20 = bi[17];
  const T21 = T1X * T1Y + T1Z * T20;
  const T32 = T1X * T20 - T1Z * T1Y;
  const T23 = br[7], T25 = bi[7];
  const T26 = T22 * T23 + T24 * T25;
  const T33 = T22 * T25 - T24 * T23;

  const T27 = T21 + T26;
  const T3E = T32 + T33;
  const T2r = T21 - T26;
  const T34 = T32 - T33;

  const T1F = br[13], T1H = bi[13];
  const T1I = T1E * T1F + T1G * T1H;
  const T2X = T1E * T1H - T1G * T1F;
  const T1J = br[3], T1K = bi[3];
  const T1L = Tg * T1J + Ti * T1K;
  const T2Y = Tg * T1K - Ti * T1J;

  const T1M = T1I + T1L;
  const T3B = T2X + T2Y;
  const T2m = T1I - T1L;
  const T2Z = T2X - T2Y;

  const T1P = br[12], T1R = bi[12];
  const T1S = T1O * T1P + T1Q * T1R;
  const T2o = T1O * T1R - T1Q * T1P;
  const T1T = br[2], T1U = bi[2];
  const T1V = T1f * T1T + T1h * T1U;
  const T2p = T1f * T1U - T1h * T1T;

  const T1W = T1S + T1V;
  const T3D = T2o + T2p;
  const T2q = T2o - T2p;
  const T31 = T1S - T1V;

  const T1N = T1D - T1M;
  const T28 = T1W - T27;
  const T29 = T1N + T28;
  const T3C = T3A - T3B;
  const T3F = T3D - T3E;
  const T4o = T3C + T3F;
  const T3X = T3A + T3B;
  const T3Y = T3D + T3E;
  const T44 = T3X + T3Y;
  const T2f = T1D + T1M;
  const T2g = T1W + T27;
  const T2h = T2f + T2g;
  const T2n = T2l + T2m;
  const T2s = T2q + T2r;
  const T4L = T2n + T2s;
  const T3g = T2l - T2m;
  const T3h = T2q - T2r;
  const T4w = T3g + T3h;
  const T3n = T2W + T2Z;
  const T3o = T31 + T34;
  const T3p = T3n + T3o;
  const T30 = T2W - T2Z;
  const T35 = T31 - T34;
  const T36 = T30 + T35;

  const TG = br[4], TH = bi[4];
  const TI = Tk * TG + To * TH;
  const T2u = Tk * TH - To * TG;
  const TK = br[14], TM = bi[14];
  const TN = TJ * TK + TL * TM;
  const T2v = TJ * TM - TL * TK;

  const TO = TI + TN;
  const T3H = T2u + T2v;
  const T2w = T2u - T2v;
  const T2L = TI - TN;

  const T1c = br[1], T1d = bi[1];
  const T1e = T2 * T1c + T5 * T1d;
  const T2R = T2 * T1d - T5 * T1c;
  const T1k = br[11], T1o = bi[11];
  const T1p = T1j * T1k + T1n * T1o;
  const T2S = T1j * T1o - T1n * T1k;

  const T1q = T1e + T1p;
  const T3L = T2R + T2S;
  const T2C = T1e - T1p;
  const T2T = T2R - T2S;

  const TP = br[9], TQ = bi[9];
  const TR = T3 * TP + T6 * TQ;
  const T2M = T3 * TQ - T6 * TP;
  const TS = br[19], TT = bi[19];
  const TU = Tx * TS + Ty * TT;
  const T2N = Tx * TT - Ty * TS;

  const TV = TR + TU;
  const T3I = T2M + T2N;
  const T2x = TR - TU;
  const T2O = T2M - T2N;

  const TY = br[16], T10 = bi[16];
  const T11 = TX * TY + TZ * T10;
  const T2z = TX * T10 - TZ * TY;
  const T15 = br[6], T19 = bi[6];
  const T1a = T14 * T15 + T18 * T19;
  const T2A = T14 * T19 - T18 * T15;

  const T1b = T11 + T1a;
  const T3K = T2z + T2A;
  const T2B = T2z - T2A;
  const T2Q = T11 - T1a;

  const TW = TO - TV;
  const T1r = T1b - T1q;
  const T1s = TW + T1r;
  const T3J = T3H - T3I;
  const T3M = T3K - T3L;
  const T4n = T3J + T3M;
  const T3U = T3H + T3I;
  const T3V = T3K + T3L;
  const T43 = T3U + T3V;
  const T2c = TO + TV;
  const T2d = T1b + T1q;
  const T2e = T2c + T2d;
  const T2y = T2w + T2x;
  const T2D = T2B + T2C;
  const T4K = T2y + T2D;
  const T3d = T2w - T2x;
  const T3e = T2B - T2C;
  const T4v = T3d + T3e;
  const T3k = T2L + T2O;
  const T3l = T2Q + T2T;
  const T3m = T3k + T3l;
  const T2P = T2L - T2O;
  const T2U = T2Q - T2T;
  const T2V = T2P + T2U;

  const outR = new Float64Array(20), outI = new Float64Array(20);

  {
    const T3y = KP559016994 * (T1s - T29);
    const T2a = T1s + T29;
    const T3x = TF - KP250000000 * T2a;
    const T3G = T3C - T3F;
    const T3N = T3J - T3M;
    const T3O = KP951056516 * T3G - KP587785252 * T3N;
    const T3Q = KP951056516 * T3N + KP587785252 * T3G;
    outR[10] = TF + T2a;
    const T3P = T3y + T3x;
    outR[14] = T3P - T3Q;
    outR[6] = T3P + T3Q;
    const T3z = T3x - T3y;
    outR[2] = T3z - T3O;
    outR[18] = T3z + T3O;
  }

  {
    const T4r = KP559016994 * (T4n - T4o);
    const T4p = T4n + T4o;
    const T4q = T4m - KP250000000 * T4p;
    const T4j = T1N - T28;
    const T4k = TW - T1r;
    const T4l = KP951056516 * T4j - KP587785252 * T4k;
    const T4u = KP951056516 * T4k + KP587785252 * T4j;
    outI[10] = T4p + T4m;
    const T4t = T4r + T4q;
    outI[6] = T4t - T4u;
    outI[14] = T4u + T4t;
    const T4s = T4q - T4r;
    outI[2] = T4l + T4s;
    outI[18] = T4s - T4l;
  }

  {
    const T3R = KP559016994 * (T2e - T2h);
    const T2i = T2e + T2h;
    const T3S = T2b - KP250000000 * T2i;
    const T3W = T3U - T3V;
    const T3Z = T3X - T3Y;
    const T40 = KP951056516 * T3W + KP587785252 * T3Z;
    const T42 = KP951056516 * T3Z - KP587785252 * T3W;
    outR[0] = T2b + T2i;
    const T41 = T3S - T3R;
    outR[12] = T41 - T42;
    outR[8] = T41 + T42;
    const T3T = T3R + T3S;
    outR[4] = T3T - T40;
    outR[16] = T3T + T40;
  }

  {
    const T4e = KP559016994 * (T43 - T44);
    const T45 = T43 + T44;
    const T4f = T4a - KP250000000 * T45;
    const T4b = T2c - T2d;
    const T4c = T2f - T2g;
    const T4d = KP951056516 * T4b + KP587785252 * T4c;
    const T4i = KP951056516 * T4c - KP587785252 * T4b;
    outI[0] = T45 + T4a;
    const T4h = T4f - T4e;
    outI[8] = T4h - T4i;
    outI[12] = T4i + T4h;
    const T4g = T4e + T4f;
    outI[4] = T4d + T4g;
    outI[16] = T4g - T4d;
  }

  {
    const T39 = KP559016994 * (T2V - T36);
    const T37 = T2V + T36;
    const T38 = T2K - KP250000000 * T37;
    const T2t = T2n - T2s;
    const T2E = T2y - T2D;
    const T2F = KP951056516 * T2t - KP587785252 * T2E;
    const T3b = KP951056516 * T2E + KP587785252 * T2t;
    outR[15] = T2K + T37;
    const T3c = T39 + T38;
    outR[11] = T3b + T3c;
    outR[19] = T3c - T3b;
    const T3a = T38 - T39;
    outR[3] = T2F + T3a;
    outR[7] = T3a - T2F;
  }

  {
    const T4O = KP559016994 * (T4K - T4L);
    const T4M = T4K + T4L;
    const T4N = T4J - KP250000000 * T4M;
    const T4Q = T30 - T35;
    const T4R = T2P - T2U;
    const T4S = KP951056516 * T4Q - KP587785252 * T4R;
    const T4U = KP951056516 * T4R + KP587785252 * T4Q;
    outI[15] = T4M + T4J;
    const T4T = T4O + T4N;
    outI[11] = T4T - T4U;
    outI[19] = T4U + T4T;
    const T4P = T4N - T4O;
    outI[3] = T4P - T4S;
    outI[7] = T4S + T4P;
  }

  {
    const T3q = KP559016994 * (T3m - T3p);
    const T3s = T3m + T3p;
    const T3t = T3r - KP250000000 * T3s;
    const T3f = T3d - T3e;
    const T3i = T3g - T3h;
    const T3j = KP951056516 * T3f + KP587785252 * T3i;
    const T3v = KP951056516 * T3i - KP587785252 * T3f;
    outR[5] = T3r + T3s;
    const T3w = T3t - T3q;
    outR[13] = T3v + T3w;
    outR[17] = T3w - T3v;
    const T3u = T3q + T3t;
    outR[1] = T3j + T3u;
    outR[9] = T3u - T3j;
  }

  {
    const T4x = KP559016994 * (T4v - T4w);
    const T4B = T4v + T4w;
    const T4C = T4A - KP250000000 * T4B;
    const T4E = T3k - T3l;
    const T4F = T3n - T3o;
    const T4G = KP951056516 * T4E + KP587785252 * T4F;
    const T4I = KP951056516 * T4F - KP587785252 * T4E;
    outI[5] = T4B + T4A;
    const T4H = T4C - T4x;
    outI[13] = T4H - T4I;
    outI[17] = T4I + T4H;
    const T4D = T4x + T4C;
    outI[1] = T4D - T4G;
    outI[9] = T4G + T4D;
  }

  return [outR, outI];
}

module.exports = { t2_20 };
