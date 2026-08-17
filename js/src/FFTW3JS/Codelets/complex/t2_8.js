'use strict';

// =============================================================================
// t2_8.js -- faithful JS port of dft/scalar/codelets/t2_8.c (non-FMA branch),
// FFTW3's "twiddle-log3 / precompute-twiddles" alternate-codegen sibling of
// a plain t1_8 (which doesn't exist as a separate ported file -- t2_8 is the
// only radix-8 twiddle codelet in this port). twinstr only trig-generates
// W^1, W^3, W^7 (three raw pairs); W^2, W^4, W^5, W^6 are DERIVED via one
// complex multiply each (same "twiddle-log3" trick as t2_5.js -- see that
// file's header for the general explanation):
//   (Tg,Ti) = conj(W^1)*W^3 = W^2
//   (T8,Tc) = W^1*W^3       = W^4
//   (Tx,Tz) = conj(W^2)*W^7 = W^5
//   (Tn,Tp) = conj(W^1)*W^7 = W^6
// Same calling convention as every other twiddle codelet here (br/bi/Wc/Ws
// with Composite1D.js's full r-1 = 7 pair table already built) -- only
// indices 1, 3, 7 are actually read.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP707106781 = 0.707106781186547524400844362104849039284835938;

function t2_8(br, bi, Wc, Ws) {
  const T2 = Wc[1], T5 = Ws[1], T3 = Wc[3], T6 = Ws[3];
  const T4 = T2 * T3, Tb = T5 * T3, T7 = T5 * T6, Ta = T2 * T6;
  const T8 = T4 - T7, Tc = Ta + Tb, Tg = T4 + T7, Ti = Ta - Tb;
  const Tl = Wc[7], Tm = Ws[7];
  const Tn = T2 * Tl + T5 * Tm;
  const Tz = Tg * Tm - Ti * Tl;
  const Tp = T2 * Tm - T5 * Tl;
  const Tx = Tg * Tl + Ti * Tm;

  const T1 = br[0], T1c = bi[0];
  const T9 = br[4], Td = bi[4];
  const Te = T8 * T9 + Tc * Td;
  const T1b = T8 * Td - Tc * T9;
  const Tf = T1 + Te;
  const T1i = T1c - T1b;
  const TL = T1 - Te;
  const T1d = T1b + T1c;

  const TD = br[7], TE = bi[7];
  const TF = Tl * TD + Tm * TE;
  const TW = Tl * TE - Tm * TD;
  const TG = br[3], TH = bi[3];
  const TI = T3 * TG + T6 * TH;
  const TX = T3 * TH - T6 * TG;
  const TJ = TF + TI;
  const T17 = TW + TX;
  const TV = TF - TI;
  const TY = TW - TX;

  const Th = br[2], Tj = bi[2];
  const Tk = Tg * Th + Ti * Tj;
  const TM = Tg * Tj - Ti * Th;
  const To = br[6], Tq = bi[6];
  const Tr = Tn * To + Tp * Tq;
  const TN = Tn * Tq - Tp * To;
  const Ts = Tk + Tr;
  const T1j = Tk - Tr;
  const TO = TM - TN;
  const T1a = TM + TN;

  const Tu = br[1], Tv = bi[1];
  const Tw = T2 * Tu + T5 * Tv;
  const TR = T2 * Tv - T5 * Tu;
  const Ty = br[5], TA = bi[5];
  const TB = Tx * Ty + Tz * TA;
  const TS = Tx * TA - Tz * Ty;
  const TC = Tw + TB;
  const T16 = TR + TS;
  const TQ = Tw - TB;
  const TT = TR - TS;

  const outR = new Float64Array(8), outI = new Float64Array(8);

  const Tt = Tf + Ts;
  const TK = TC + TJ;
  outR[4] = Tt - TK;
  outR[0] = Tt + TK;

  const T19 = T16 + T17;
  const T1e = T1a + T1d;
  outI[0] = T19 + T1e;
  outI[4] = T1e - T19;
  const T15 = Tf - Ts;
  const T18 = T16 - T17;
  outR[6] = T15 - T18;
  outR[2] = T15 + T18;

  const T1f = TJ - TC;
  const T1g = T1d - T1a;
  outI[2] = T1f + T1g;
  outI[6] = T1g - T1f;

  const T11 = TL - TO;
  const T1k = T1i - T1j;
  const T12 = TT - TQ;
  const T13 = TV + TY;
  const T14 = KP707106781 * (T12 - T13);
  const T1h = KP707106781 * (T12 + T13);
  outR[7] = T11 - T14;
  outI[5] = T1k - T1h;
  outR[3] = T11 + T14;
  outI[1] = T1h + T1k;

  const TP = TL + TO;
  const T1m = T1j + T1i;
  const TU = TQ + TT;
  const TZ = TV - TY;
  const T10 = KP707106781 * (TU + TZ);
  const T1l = KP707106781 * (TZ - TU);
  outR[5] = TP - T10;
  outI[7] = T1m - T1l;
  outR[1] = TP + T10;
  outI[3] = T1l + T1m;

  return [outR, outI];
}

module.exports = { t2_8 };
