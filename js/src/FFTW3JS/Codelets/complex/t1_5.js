'use strict';

// =============================================================================
// t1_5.js -- faithful JS port of dft/scalar/codelets/t1_5.c (non-FMA branch),
// FFTW3's radix-5 Cooley-Tukey twiddle/combine codelet (decimation-in-time).
// br/bi[0] is the untwiddled DC block; phases 1-4 are twiddled by
// Wc[p]/Ws[p] = cos/sin(2*pi*col*p/n) before combining. Same radix-5
// butterfly structure as n1_5.js, applied to the twiddled phase values.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;

function t1_5(br, bi, Wc, Ws) {
  const T1 = br[0], TE = bi[0];

  const T6 = Wc[1] * br[1] + Ws[1] * bi[1];
  const Ts = Wc[1] * bi[1] - Ws[1] * br[1];

  const Tm = Wc[3] * br[3] + Ws[3] * bi[3];
  const Tw = Wc[3] * bi[3] - Ws[3] * br[3];

  const Tb = Wc[4] * br[4] + Ws[4] * bi[4];
  const Tt = Wc[4] * bi[4] - Ws[4] * br[4];

  const Th = Wc[2] * br[2] + Ws[2] * bi[2];
  const Tv = Wc[2] * bi[2] - Ws[2] * br[2];

  const Tu = Ts - Tt;
  const Tx = Tv - Tw;
  const TJ = Th - Tm;
  const TI = T6 - Tb;
  const TB = Ts + Tt;
  const TC = Tv + Tw;
  const TD = TB + TC;
  const Tc = T6 + Tb;
  const Tn = Th + Tm;
  const To = Tc + Tn;

  const outR = new Float64Array(5), outI = new Float64Array(5);
  outR[0] = T1 + To;
  outI[0] = TD + TE;

  const Ty = KP951056516 * Tu + KP587785252 * Tx;
  const TA = KP951056516 * Tx - KP587785252 * Tu;
  const Tp = KP559016994 * (Tc - Tn);
  const Tq = T1 - KP250000000 * To;
  const Tr = Tp + Tq;
  const Tz = Tq - Tp;
  outR[4] = Tr - Ty;
  outR[3] = Tz + TA;
  outR[1] = Tr + Ty;
  outR[2] = Tz - TA;

  const TK = KP951056516 * TI + KP587785252 * TJ;
  const TL = KP951056516 * TJ - KP587785252 * TI;
  const TF = KP559016994 * (TB - TC);
  const TG = TE - KP250000000 * TD;
  const TH = TF + TG;
  const TM = TG - TF;
  outI[1] = TH - TK;
  outI[3] = TM - TL;
  outI[4] = TK + TH;
  outI[2] = TL + TM;

  return [outR, outI];
}

module.exports = { t1_5 };
