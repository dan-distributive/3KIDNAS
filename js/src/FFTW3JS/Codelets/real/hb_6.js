'use strict';

// =============================================================================
// hb_6.js -- faithful JS port of rdft/scalar/r2cb/hb_6.c (non-FMA branch).
// Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hb_4.js -- the PLAIN
// (non-alt-codegen) backward family, BACKWARD counterpart of hf_6.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP500000000 = 0.500000000000000000000000000000000000000000000;
const KP866025403 = 0.866025403784438646763723170752936183471402627;

function hb_6(cr, ci, Wc, Ws) {
  const T1 = cr[0];
  const T2 = ci[2];
  const T3 = T1 + T2;
  const Ty = T1 - T2;

  const T4 = cr[2];
  const T5 = ci[0];
  const T6 = T4 + T5;
  const Tz = T4 - T5;
  const T7 = ci[1];
  const T8 = cr[1];
  const T9 = T7 + T8;
  const TA = T7 - T8;

  const Ta = T6 + T9;
  const TO = KP866025403 * (Tz - TA);
  const Tr = KP866025403 * (T6 - T9);
  const TB = Tz + TA;

  const Tb = ci[5];
  const Tc = cr[3];
  const Td = Tb - Tc;
  const TE = Tb + Tc;

  const Te = ci[3];
  const Tf = cr[5];
  const Tg = Te - Tf;
  const TG = Te + Tf;
  const Th = ci[4];
  const Ti = cr[4];
  const Tj = Th - Ti;
  const TF = Th + Ti;

  const Tk = Tg + Tj;
  const TL = KP866025403 * (TG + TF);
  const Tn = KP866025403 * (Tj - Tg);
  const TH = TF - TG;

  const outCr = new Float64Array(6), outCi = new Float64Array(6);
  outCr[0] = T3 + Ta;
  outCi[0] = Td + Tk;

  {
    const TC = Ty + TB;
    const TI = TE - TH;
    const Tx = Wc[3], TD = Ws[3];
    outCr[3] = Tx * TC - TD * TI; // FNMS
    outCi[3] = TD * TC + Tx * TI; // FMA
  }
  {
    const Tm = T3 - KP500000000 * Ta; // FNMS
    const To = Tm - Tn;
    const Tu = Tm + Tn;
    const Tq = Td - KP500000000 * Tk; // FNMS
    const Ts = Tq - Tr;
    const Tw = Tr + Tq;

    const Tl = Wc[2], Tp = Ws[2];
    outCr[2] = Tl * To - Tp * Ts; // FNMS
    outCi[2] = Tl * Ts + Tp * To; // FMA
    const Tt = Wc[4], Tv = Ws[4];
    outCr[4] = Tt * Tu - Tv * Tw; // FNMS
    outCi[4] = Tt * Tw + Tv * Tu; // FMA
  }
  {
    const TK = Ty - KP500000000 * TB; // FNMS
    const TM = TK - TL;
    const TS = TK + TL;
    const TP = KP500000000 * TH + TE; // FMA
    const TQ = TO + TP;
    const TU = TP - TO;

    const TJ = Wc[1], TN = Ws[1];
    outCr[1] = TJ * TM - TN * TQ; // FNMS
    outCi[1] = TN * TM + TJ * TQ; // FMA
    const TR = Wc[5], TT = Ws[5];
    outCr[5] = TR * TS - TT * TU; // FNMS
    outCi[5] = TR * TU + TT * TS; // FMA
  }

  return [outCr, outCi];
}

module.exports = { hb_6 };
