'use strict';

// =============================================================================
// t2_5.js -- faithful JS port of dft/scalar/codelets/t2_5.c (non-FMA branch),
// FFTW3's "twiddle-log3 / precompute-twiddles" alternate-codegen sibling of
// t1_5 (same math, different instruction schedule, different bit-exact
// rounding). Real FFTW's twinstr for this codelet is {TW_CEXP,0,1},
// {TW_CEXP,0,3} -- it only trig-generates the col=1 and col=3 twiddle pairs
// (W^j, W^3j) and DERIVES W^2j/W^4j from them via one complex multiply each
// (Tb/Tf = W^j*W^3j = W^4j; Tj/Tl = conj(W^j*conj(W^3j)) = W^2j), rather than
// the plain "all r-1 pairs from the trig table directly" convention t1_5
// uses. To keep this file's calling convention IDENTICAL to every other
// twiddle codelet here (so Composite1D.js's dispatch needs no special case),
// this still takes the full Wc[1..4]/Ws[1..4] table Composite1D already
// builds and simply reads only indices 1 and 3 out of it -- mathematically
// and bit-for-bit equivalent to the C source reading W[0..3], since
// Wc[1]/Ws[1]=W[0]/W[1] and Wc[3]/Ws[3]=W[2]/W[3] there.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;

function t2_5(br, bi, Wc, Ws) {
  const T2 = Wc[1], T4 = Ws[1], T7 = Wc[3], T9 = Ws[3];
  const T8 = T2 * T7, Te = T4 * T7, Ta = T4 * T9, Td = T2 * T9;
  const Tb = T8 - Ta, Tl = Td - Te, Tf = Td + Te, Tj = T8 + Ta;

  const T1 = br[0], TI = bi[0];

  const T3 = br[1], T5 = bi[1];
  const T6 = T2 * T3 + T4 * T5;
  const Tw = T2 * T5 - T4 * T3;

  const To = br[3], Tp = bi[3];
  const Tq = T7 * To + T9 * Tp;
  const TA = T7 * Tp - T9 * To;

  const Tc = br[4], Tg = bi[4];
  const Th = Tb * Tc + Tf * Tg;
  const Tx = Tb * Tg - Tf * Tc;

  const Tk = br[2], Tm = bi[2];
  const Tn = Tj * Tk + Tl * Tm;
  const Tz = Tj * Tm - Tl * Tk;

  const Ty = Tw - Tx;
  const TB = Tz - TA;
  const TN = Tn - Tq;
  const TM = T6 - Th;
  const TF = Tw + Tx;
  const TG = Tz + TA;
  const TH = TF + TG;
  const Ti = T6 + Th;
  const Tr = Tn + Tq;
  const Ts = Ti + Tr;

  const outR = new Float64Array(5), outI = new Float64Array(5);
  outR[0] = T1 + Ts;
  outI[0] = TH + TI;

  const TC = KP951056516 * Ty + KP587785252 * TB;
  const TE = KP951056516 * TB - KP587785252 * Ty;
  const Tt = KP559016994 * (Ti - Tr);
  const Tu = T1 - KP250000000 * Ts;
  const Tv = Tt + Tu;
  const TD = Tu - Tt;
  outR[4] = Tv - TC;
  outR[3] = TD + TE;
  outR[1] = Tv + TC;
  outR[2] = TD - TE;

  const TO = KP951056516 * TM + KP587785252 * TN;
  const TP = KP951056516 * TN - KP587785252 * TM;
  const TJ = KP559016994 * (TF - TG);
  const TK = TI - KP250000000 * TH;
  const TL = TJ + TK;
  const TQ = TK - TJ;
  outI[1] = TL - TO;
  outI[3] = TQ - TP;
  outI[4] = TO + TL;
  outI[2] = TP + TQ;

  return [outR, outI];
}

module.exports = { t2_5 };
