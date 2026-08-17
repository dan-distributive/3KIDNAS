'use strict';

// =============================================================================
// n1_8.js -- faithful JS port of dft/scalar/codelets/n1_8.c (non-FMA branch),
// FFTW3's direct (base-case) radix-8 complex DFT codelet.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP707106781 = 0.707106781186547524400844362104849039284835938;

function n1_8(ri, ii) {
  const T1 = ri[0], T2 = ri[4];
  const T3 = T1 + T2;
  const Tn = T1 - T2;
  const Tg = ii[0], Th = ii[4];
  const Ti = Tg + Th;
  const TC = Tg - Th;
  const T4 = ri[2], T5 = ri[6];
  const T6 = T4 + T5;
  const TB = T4 - T5;
  const Tj = ii[2], Tk = ii[6];
  const Tl = Tj + Tk;
  const To = Tj - Tk;

  const Tb = ri[7], Tc = ri[3];
  const Tv = Tb - Tc;
  const Tw = ii[7], Tx = ii[3];
  const Ty = Tw - Tx;
  const Td = Tb + Tc;
  const TN = Tw + Tx;
  const Tz = Tv - Ty;
  const TH = Tv + Ty;

  const T8 = ri[1], T9 = ri[5];
  const Tq = T8 - T9;
  const Tr = ii[1], Ts = ii[5];
  const Tt = Tr - Ts;
  const Ta = T8 + T9;
  const TM = Tr + Ts;
  const Tu = Tq + Tt;
  const TG = Tt - Tq;

  const ro = new Float64Array(8), io = new Float64Array(8);

  const T7 = T3 + T6, Te = Ta + Td;
  ro[4] = T7 - Te;
  ro[0] = T7 + Te;
  const TP = Ti + Tl, TQ = TM + TN;
  io[4] = TP - TQ;
  io[0] = TP + TQ;

  const Tf = Td - Ta, Tm = Ti - Tl;
  io[2] = Tf + Tm;
  io[6] = Tm - Tf;
  const TL = T3 - T6, TO = TM - TN;
  ro[6] = TL - TO;
  ro[2] = TL + TO;

  const Tp = Tn + To, TA = KP707106781 * (Tu + Tz);
  ro[5] = Tp - TA;
  ro[1] = Tp + TA;
  const TJ = TC - TB, TK = KP707106781 * (TG + TH);
  io[5] = TJ - TK;
  io[1] = TJ + TK;

  const TD = TB + TC, TE = KP707106781 * (Tz - Tu);
  io[7] = TD - TE;
  io[3] = TD + TE;
  const TF = Tn - To, TI = KP707106781 * (TG - TH);
  ro[7] = TF - TI;
  ro[3] = TF + TI;

  return [ro, io];
}

module.exports = { n1_8 };
