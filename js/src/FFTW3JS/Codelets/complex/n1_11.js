'use strict';

// =============================================================================
// n1_11.js -- faithful JS port of dft/scalar/codelets/n1_11.c (non-FMA
// branch), FFTW3's direct (base-case) radix-11 complex DFT codelet.
//
// IMPORTANT: several outputs are sums of 2-3 separate FMA/FNMS/FNMA macro
// results (e.g. "Th = FMA(...) + FNMS(...) + FNMA(...)"). Each macro call's
// OWN internal addition must be computed as its own intermediate BEFORE
// being summed with the others -- flattening straight into one left-to-right
// JS expression changes the floating-point rounding (addition isn't
// associative) and was caught by strict bit-exact verification against real
// FFTW3 (a ~1-2 ULP mismatch on exactly the outputs built from 3-way sums).
// So every such macro result gets its own named Mn variable here, combined
// in the same left-to-right order the C source shows.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP654860733 = 0.654860733945285064056925072466293553183791199;
const KP142314838 = 0.142314838273285140443792668616369668791051361;
const KP959492973 = 0.959492973614497389890368057066327699062454848;
const KP415415013 = 0.415415013001886425529274149229623203524004910;
const KP841253532 = 0.841253532831181168861811648919367717513292498;
const KP989821441 = 0.989821441880932732376092037776718787376519372;
const KP909631995 = 0.909631995354518371411715383079028460060241051;
const KP281732556 = 0.281732556841429697711417915346616899035777899;
const KP540640817 = 0.540640817455597582107635954318691695431770608;
const KP755749574 = 0.755749574354258283774035843972344420179717445;

function n1_11(ri, ii) {
  const T1 = ri[0], TM = ii[0];
  const T2 = ri[1], T3 = ri[10];
  const T4 = T2 + T3, TG = T3 - T2;
  const Ti = ii[1], Tj = ii[10];
  const Tk = Ti - Tj, TR = Ti + Tj;
  const Tu = ii[2], Tv = ii[9];
  const Tw = Tu - Tv, TN = Tu + Tv;
  const T5 = ri[2], T6 = ri[9];
  const T7 = T5 + T6, TK = T6 - T5;

  const T8 = ri[3], T9 = ri[8];
  const Ta = T8 + T9, TH = T9 - T8;
  const Tl = ii[3], Tm = ii[8];
  const Tn = Tl - Tm, TQ = Tl + Tm;
  const Tb = ri[4], Tc = ri[7];
  const Td = Tb + Tc, TJ = Tc - Tb;
  const To = ii[4], Tp = ii[7];
  const Tq = To - Tp, TO = To + Tp;
  const Tr = ii[5], Ts = ii[6];
  const Tt = Tr - Ts, TP = Tr + Ts;
  const Te = ri[5], Tf = ri[6];
  const Tg = Te + Tf, TI = Tf - Te;

  const ro = new Float64Array(11), io = new Float64Array(11);

  ro[0] = T1 + T4 + T7 + Ta + Td + Tg;
  io[0] = TM + TR + TN + TQ + TO + TP;

  // Tx, Th -> ro[7], ro[4]
  {
    const M1 = KP755749574 * Tk + KP540640817 * Tn;
    const M2 = KP281732556 * Tq - KP909631995 * Tt;
    const Tx = M1 + M2 - KP989821441 * Tw;
    const N1 = KP841253532 * Ta + T1;
    const N2 = KP415415013 * Tg - KP959492973 * Td;
    const N3 = -(KP142314838 * T7 + KP654860733 * T4);
    const Th = N1 + N2 + N3;
    ro[7] = Th - Tx;
    ro[4] = Th + Tx;
  }

  // TZ, T10 -> io[4], io[7]
  {
    const M1 = KP755749574 * TG + KP540640817 * TH;
    const M2 = KP281732556 * TJ - KP909631995 * TI;
    const TZ = M1 + M2 - KP989821441 * TK;
    const N1 = KP841253532 * TQ + TM;
    const N2 = KP415415013 * TP - KP959492973 * TO;
    const N3 = -(KP142314838 * TN + KP654860733 * TR);
    const T10 = N1 + N2 + N3;
    io[4] = TZ + T10;
    io[7] = T10 - TZ;
  }

  // TX, TY -> io[2], io[9]
  {
    const M1 = KP909631995 * TG + KP755749574 * TK;
    const M2 = -(KP540640817 * TI + KP989821441 * TJ);
    const TX = M1 + M2 - KP281732556 * TH;
    const N1 = KP415415013 * TR + TM;
    const N2 = KP841253532 * TP - KP142314838 * TO;
    const N3 = -(KP959492973 * TQ + KP654860733 * TN);
    const TY = N1 + N2 + N3;
    io[2] = TX + TY;
    io[9] = TY - TX;
  }

  // Tz, Ty -> ro[9], ro[2]
  {
    const M1 = KP909631995 * Tk + KP755749574 * Tw;
    const M2 = -(KP540640817 * Tt + KP989821441 * Tq);
    const Tz = M1 + M2 - KP281732556 * Tn;
    const N1 = KP415415013 * T4 + T1;
    const N2 = KP841253532 * Tg - KP142314838 * Td;
    const N3 = -(KP959492973 * Ta + KP654860733 * T7);
    const Ty = N1 + N2 + N3;
    ro[9] = Ty - Tz;
    ro[2] = Ty + Tz;
  }

  // TB, TA -> ro[10], ro[1]
  {
    const M1 = KP540640817 * Tk + KP909631995 * Tw;
    const M2 = KP989821441 * Tn + KP755749574 * Tq;
    const TB = M1 + M2 + KP281732556 * Tt;
    const N1 = KP841253532 * T4 + T1;
    const N2 = KP415415013 * T7 - KP959492973 * Tg;
    const N3 = -(KP654860733 * Td + KP142314838 * Ta);
    const TA = N1 + N2 + N3;
    ro[10] = TA - TB;
    ro[1] = TA + TB;
  }

  // TV, TW -> io[1], io[10]
  {
    const M1 = KP540640817 * TG + KP909631995 * TK;
    const M2 = KP989821441 * TH + KP755749574 * TJ;
    const TV = M1 + M2 + KP281732556 * TI;
    const N1 = KP841253532 * TR + TM;
    const N2 = KP415415013 * TN - KP959492973 * TP;
    const N3 = -(KP654860733 * TO + KP142314838 * TQ);
    const TW = N1 + N2 + N3;
    io[1] = TV + TW;
    io[10] = TW - TV;
  }

  // TD, TC -> ro[8], ro[3]
  {
    const M1 = KP989821441 * Tk + KP540640817 * Tq;
    const M2 = KP755749574 * Tt - KP909631995 * Tn;
    const TD = M1 + M2 - KP281732556 * Tw;
    const N1 = KP415415013 * Ta + T1;
    const N2 = KP841253532 * Td - KP654860733 * Tg;
    const N3 = -(KP959492973 * T7 + KP142314838 * T4);
    const TC = N1 + N2 + N3;
    ro[8] = TC - TD;
    ro[3] = TC + TD;
  }

  // TT, TU -> io[3], io[8]
  {
    const M1 = KP989821441 * TG + KP540640817 * TJ;
    const M2 = KP755749574 * TI - KP909631995 * TH;
    const TT = M1 + M2 - KP281732556 * TK;
    const N1 = KP415415013 * TQ + TM;
    const N2 = KP841253532 * TO - KP654860733 * TP;
    const N3 = -(KP959492973 * TN + KP142314838 * TR);
    const TU = N1 + N2 + N3;
    io[3] = TT + TU;
    io[8] = TU - TT;
  }

  // TL, TS -> io[5], io[6]
  {
    const M1 = KP281732556 * TG + KP755749574 * TH;
    const M2 = KP989821441 * TI - KP909631995 * TJ;
    const TL = M1 + M2 - KP540640817 * TK;
    const N1 = KP841253532 * TN + TM;
    const N2 = KP415415013 * TO - KP142314838 * TP;
    const N3 = -(KP654860733 * TQ + KP959492973 * TR);
    const TS = N1 + N2 + N3;
    io[5] = TL + TS;
    io[6] = TS - TL;
  }

  // TF, TE -> ro[6], ro[5]
  {
    const M1 = KP281732556 * Tk + KP755749574 * Tn;
    const M2 = KP989821441 * Tt - KP909631995 * Tq;
    const TF = M1 + M2 - KP540640817 * Tw;
    const N1 = KP841253532 * T7 + T1;
    const N2 = KP415415013 * Td - KP142314838 * Tg;
    const N3 = -(KP654860733 * Ta + KP959492973 * T4);
    const TE = N1 + N2 + N3;
    ro[6] = TE - TF;
    ro[5] = TE + TF;
  }

  return [ro, io];
}

module.exports = { n1_11 };
