'use strict';

// =============================================================================
// hf_6.js -- faithful JS port of rdft/scalar/r2cf/hf_6.c (non-FMA branch).
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hf_2.js/hf_3.js/hf_4.js
// -- the PLAIN (non-alt-codegen) family, TW_FULL twinstr (all r-1 twiddle
// phases read from Wc[]/Ws[], unlike hf2_r's restricted subset). Radix 6
// has no hf2_6 in real FFTW's own registered codelet set, so this plain
// family is the only way to cover it.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function hf_6(cr, ci, Wc, Ws) {
  const T1 = cr[0];
  const TM = ci[0];
  const T3 = cr[3], T5 = ci[3];
  const T6 = Wc[3] * T3 + Ws[3] * T5; // FMA
  const TN = Wc[3] * T5 - Ws[3] * T3; // FNMS
  const T7 = T1 - T6;
  const TS = TN + TM;
  const Tv = T1 + T6;
  const TO = TM - TN;

  const Tk = cr[4], Tm = ci[4];
  const Tn = Wc[4] * Tk + Ws[4] * Tm; // FMA
  const TE = Wc[4] * Tm - Ws[4] * Tk; // FNMS
  const Tp = cr[1], Tr = ci[1];
  const Ts = Wc[1] * Tp + Ws[1] * Tr; // FMA
  const TD = Wc[1] * Tr - Ws[1] * Tp; // FNMS
  const Tt = Tn - Ts;
  const TJ = TE + TD;
  const Tx = Tn + Ts;
  const TF = TD - TE;

  const T9 = cr[2], Tb = ci[2];
  const Tc = Wc[2] * T9 + Ws[2] * Tb; // FMA
  const TA = Wc[2] * Tb - Ws[2] * T9; // FNMS
  const Te = cr[5], Tg = ci[5];
  const Th = Wc[5] * Te + Ws[5] * Tg; // FMA
  const TB = Wc[5] * Tg - Ws[5] * Te; // FNMS
  const Ti = Tc - Th;
  const TI = TA + TB;
  const Tw = Tc + Th;
  const TC = TA - TB;

  const outCr = new Float64Array(6), outCi = new Float64Array(6);

  {
    const TG = KP866025403 * (TC + TF);
    const Tu = Ti + Tt;
    const Tz = T7 - KP500000000 * Tu; // FNMS
    outCi[2] = T7 + Tu;
    outCr[1] = Tz + TG;
    outCi[0] = Tz - TG;
  }
  {
    const TK = KP866025403 * (TI - TJ);
    const Ty = Tw + Tx;
    const TH = Tv - KP500000000 * Ty; // FNMS
    outCr[0] = Tv + Ty;
    outCi[1] = TH + TK;
    outCr[2] = TH - TK;
  }
  {
    const TP = KP866025403 * (Tt - Ti);
    const TL = TF - TC;
    const TQ = KP500000000 * TL + TO; // FMA
    outCr[3] = TL - TO;
    outCi[4] = TP + TQ;
    outCr[5] = TP - TQ;
  }
  {
    const TR = KP866025403 * (Tw - Tx);
    const TT = TI + TJ;
    const TU = TS - KP500000000 * TT; // FNMS
    outCr[4] = TR - TU;
    outCi[5] = TT + TS;
    outCi[3] = TR + TU;
  }

  return [outCr, outCi];
}

module.exports = { hf_6 };
