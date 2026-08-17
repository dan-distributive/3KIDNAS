'use strict';

// =============================================================================
// t2_10.js -- faithful JS port of dft/scalar/codelets/t2_10.c (non-FMA
// branch), FFTW3's "twiddle-log3 / precompute-twiddles" radix-10 twiddle
// codelet. twinstr only trig-generates W^1, W^3, W^9 (three raw pairs);
// W^2, W^4, W^5, W^6, W^7, W^8 are DERIVED via one complex multiply each
// (same trick as t2_5.js/t2_8.js -- see t2_5.js's header for the general
// explanation):
//   (Tk,Tm) = conj(W^1)*W^3 = W^2
//   (T8,Tc) = W^1*W^3       = W^4
//   (Te,Tg) = conj(W^4)*W^9 = W^5
//   (TM,TO) = conj(W^3)*W^9 = W^6
//   (Tp,Tr) = conj(W^2)*W^9 = W^7
//   (Tv,Tx) = conj(W^1)*W^9 = W^8
// Same calling convention as every other twiddle codelet here (br/bi/Wc/Ws
// with Composite1D.js's full r-1 = 9 pair table already built) -- only
// indices 1, 3, 9 are actually read.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;

function t2_10(br, bi, Wc, Ws) {
  const T2 = Wc[1], T5 = Ws[1], T3 = Wc[3], T6 = Ws[3];
  const T4 = T2 * T3, Tb = T5 * T3, T7 = T5 * T6, Ta = T2 * T6;
  const T8 = T4 - T7, Tm = Ta - Tb, Tc = Ta + Tb, Tk = T4 + T7;
  const T9 = Wc[9], Td = Ws[9];
  const Te = T8 * T9 + Tc * Td;
  const TM = T3 * T9 + T6 * Td;
  const TO = T3 * Td - T6 * T9;
  const Tg = T8 * Td - Tc * T9;
  const Tp = Tk * T9 + Tm * Td;
  const Tv = T2 * T9 + T5 * Td;
  const Tx = T2 * Td - T5 * T9;
  const Tr = Tk * Td - Tm * T9;

  const T1 = br[0], T1F = bi[0];
  const Tf = br[5], Th = bi[5];
  const Ti = Te * Tf + Tg * Th;
  const T1E = Te * Th - Tg * Tf;
  const Tj = T1 - Ti;
  const T1S = T1F - T1E;
  const TX = T1 + Ti;
  const T1G = T1E + T1F;

  const TF = br[4], TG = bi[4];
  const TH = T8 * TF + Tc * TG;
  const T1f = T8 * TG - Tc * TF;
  const TR = br[1], TS = bi[1];
  const TT = T2 * TR + T5 * TS;
  const T1j = T2 * TS - T5 * TR;

  const TI = br[9], TJ = bi[9];
  const TK = T9 * TI + Td * TJ;
  const T1g = T9 * TJ - Td * TI;
  const TN = br[6], TP = bi[6];
  const TQ = TM * TN + TO * TP;
  const T1i = TM * TP - TO * TN;

  const TL = TH - TK;
  const TU = TQ - TT;
  const TV = TL + TU;
  const T1s = T1f + T1g;
  const T1t = T1i + T1j;
  const T1C = T1s + T1t;
  const T11 = TH + TK;
  const T12 = TQ + TT;
  const T13 = T11 + T12;
  const T1h = T1f - T1g;
  const T1k = T1i - T1j;
  const T1Q = T1h + T1k;

  const Tl = br[2], Tn = bi[2];
  const To = Tk * Tl + Tm * Tn;
  const T18 = Tk * Tn - Tm * Tl;
  const TA = br[3], TB = bi[3];
  const TC = T3 * TA + T6 * TB;
  const T1c = T3 * TB - T6 * TA;

  const Tq = br[7], Ts = bi[7];
  const Tt = Tp * Tq + Tr * Ts;
  const T19 = Tp * Ts - Tr * Tq;
  const Tw = br[8], Ty = bi[8];
  const Tz = Tv * Tw + Tx * Ty;
  const T1b = Tv * Ty - Tx * Tw;

  const Tu = To - Tt;
  const TD = Tz - TC;
  const TE = Tu + TD;
  const T1v = T18 + T19;
  const T1w = T1b + T1c;
  const T1B = T1v + T1w;
  const TY = To + Tt;
  const TZ = Tz + TC;
  const T10 = TY + TZ;
  const T1a = T18 - T19;
  const T1d = T1b - T1c;
  const T1P = T1a + T1d;

  const outR = new Float64Array(10), outI = new Float64Array(10);

  const T15 = KP559016994 * (TE - TV);
  const TW = TE + TV;
  const T16 = Tj - KP250000000 * TW;
  const T1e = T1a - T1d;
  const T1l = T1h - T1k;
  const T1m = KP951056516 * T1e + KP587785252 * T1l;
  const T1o = KP951056516 * T1l - KP587785252 * T1e;
  outR[5] = Tj + TW;
  const T1n = T16 - T15;
  outR[7] = T1n - T1o;
  outR[3] = T1n + T1o;
  const T17 = T15 + T16;
  outR[9] = T17 - T1m;
  outR[1] = T17 + T1m;

  const T1R = KP559016994 * (T1P - T1Q);
  const T1T = T1P + T1Q;
  const T1U = T1S - KP250000000 * T1T;
  const T1W = Tu - TD;
  const T1X = TL - TU;
  const T1Y = KP951056516 * T1W + KP587785252 * T1X;
  const T20 = KP951056516 * T1X - KP587785252 * T1W;
  outI[5] = T1T + T1S;
  const T1Z = T1U - T1R;
  outI[3] = T1Z - T20;
  outI[7] = T20 + T1Z;
  const T1V = T1R + T1U;
  outI[1] = T1V - T1Y;
  outI[9] = T1Y + T1V;

  const T1q = KP559016994 * (T10 - T13);
  const T14 = T10 + T13;
  const T1p = TX - KP250000000 * T14;
  const T1u = T1s - T1t;
  const T1x = T1v - T1w;
  const T1y = KP951056516 * T1u - KP587785252 * T1x;
  const T1A = KP951056516 * T1x + KP587785252 * T1u;
  outR[0] = TX + T14;
  const T1z = T1q + T1p;
  outR[4] = T1z - T1A;
  outR[6] = T1z + T1A;
  const T1r = T1p - T1q;
  outR[2] = T1r - T1y;
  outR[8] = T1r + T1y;

  const T1L = KP559016994 * (T1B - T1C);
  const T1D = T1B + T1C;
  const T1K = T1G - KP250000000 * T1D;
  const T1H = T11 - T12;
  const T1I = TY - TZ;
  const T1J = KP951056516 * T1H - KP587785252 * T1I;
  const T1N = KP951056516 * T1I + KP587785252 * T1H;
  outI[0] = T1D + T1G;
  const T1O = T1L + T1K;
  outI[4] = T1N + T1O;
  outI[6] = T1O - T1N;
  const T1M = T1K - T1L;
  outI[2] = T1J + T1M;
  outI[8] = T1M - T1J;

  return [outR, outI];
}

module.exports = { t2_10 };
