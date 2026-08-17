'use strict';

// =============================================================================
// hb_10.js -- faithful JS port of rdft/scalar/r2cb/hb_10.c (non-FMA
// branch). Same (cr,ci,Wc,Ws) -> (outCr,outCi) convention as hb_6.js/
// hb_12.js -- the PLAIN (non-alt-codegen) backward family, BACKWARD
// counterpart of hf_10.js.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.250000000000000000000000000000000000000000000;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP559016994 = 0.559016994374947424102293417182819058860154590;

function hb_10(cr, ci, Wc, Ws) {
  const T1 = cr[0];
  const T2 = ci[4];
  const T3 = T1 + T2;
  const T18 = T1 - T2;

  const T4 = cr[2], T5 = ci[2];
  const T6 = T4 + T5;
  const T19 = T4 - T5;
  const Te = ci[3], Tf = cr[1];
  const Tg = Te + Tf;
  const T1d = Te - Tf;

  const T7 = ci[1], T8 = cr[3];
  const T9 = T7 + T8;
  const T1a = T7 - T8;
  const Tb = cr[4], Tc = ci[0];
  const Td = Tb + Tc;
  const T1c = Tb - Tc;

  const TE = T6 - T9;
  const TF = Td - Tg;
  const T1B = T1c - T1d;
  const T1A = T19 - T1a;

  const T1b = T19 + T1a;
  const T1e = T1c + T1d;
  const T1f = T1b + T1e;
  const T1t = KP559016994 * (T1b - T1e);
  const Ta = T6 + T9;
  const Th = Td + Tg;
  const Ti = Ta + Th;
  const Tl = KP559016994 * (Ta - Th);

  const TH = ci[9], TI = cr[5];
  const TJ = TH - TI;
  const T1i = TH + TI;

  const Tn = ci[7], To = cr[7];
  const Tp = Tn - To;
  const T1j = Tn + To;
  const Tx = ci[8], Ty = cr[6];
  const Tz = Tx - Ty;
  const T1n = Tx + Ty;

  const Tq = ci[6], Tr = cr[8];
  const Ts = Tq - Tr;
  const T1k = Tq + Tr;
  const Tu = ci[5], Tv = cr[9];
  const Tw = Tu - Tv;
  const T1m = Tu + Tv;

  const Tt = Tp - Ts;
  const TA = Tw - Tz;
  const T1w = T1m + T1n;
  const T1v = T1j + T1k;

  const T1l = T1j - T1k;
  const T1o = T1m - T1n;
  const T1p = T1l + T1o;
  const T1E = KP559016994 * (T1l - T1o);
  const TK = Tp + Ts;
  const TL = Tw + Tz;
  const TM = TK + TL;
  const TO = KP559016994 * (TK - TL);

  const outCr = new Float64Array(10), outCi = new Float64Array(10);
  outCr[0] = T3 + Ti;
  outCi[0] = TJ + TM;

  {
    const T1g = T18 + T1f;
    const T1q = T1i + T1p;
    const T17 = Wc[5], T1h = Ws[5];
    outCr[5] = T17 * T1g - T1h * T1q; // FNMS
    outCi[5] = T1h * T1g + T17 * T1q; // FMA
  }
  {
    const TB = KP587785252 * Tt - KP951056516 * TA; // FNMS
    const TG = KP587785252 * TE - KP951056516 * TF; // FNMS
    const T11 = KP951056516 * TE + KP587785252 * TF; // FMA
    const TX = KP951056516 * Tt + KP587785252 * TA; // FMA
    const TN = TJ - KP250000000 * TM; // FNMS
    const TP = TN - TO;
    const T10 = TO + TN;
    const Tk = T3 - KP250000000 * Ti; // FNMS
    const Tm = Tk - Tl;
    const TW = Tl + Tk;

    const TC = Tm - TB;
    const TQ = TG + TP;
    const Tj = Wc[2], TD = Ws[2];
    outCr[2] = Tj * TC - TD * TQ; // FNMS
    outCi[2] = TD * TC + Tj * TQ; // FMA

    const T14 = TW - TX;
    const T16 = T11 + T10;
    const T13 = Wc[6], T15 = Ws[6];
    outCr[6] = T13 * T14 - T15 * T16; // FNMS
    outCi[6] = T15 * T14 + T13 * T16; // FMA

    const TS = Tm + TB;
    const TU = TP - TG;
    const TR = Wc[8], TT = Ws[8];
    outCr[8] = TR * TS - TT * TU; // FNMS
    outCi[8] = TT * TS + TR * TU; // FMA

    const TY = TW + TX;
    const T12 = T10 - T11;
    const TV = Wc[4], TZ = Ws[4];
    outCr[4] = TV * TY - TZ * T12; // FNMS
    outCi[4] = TZ * TY + TV * T12; // FMA
  }
  {
    const T1x = KP587785252 * T1v - KP951056516 * T1w; // FNMS
    const T1C = KP587785252 * T1A - KP951056516 * T1B; // FNMS
    const T1Q = KP951056516 * T1A + KP587785252 * T1B; // FMA
    const T1N = KP951056516 * T1v + KP587785252 * T1w; // FMA
    const T1D = T1i - KP250000000 * T1p; // FNMS
    const T1F = T1D - T1E;
    const T1R = T1E + T1D;
    const T1s = T18 - KP250000000 * T1f; // FNMS
    const T1u = T1s - T1t;
    const T1M = T1t + T1s;

    const T1y = T1u - T1x;
    const T1G = T1C + T1F;
    const T1r = Wc[7], T1z = Ws[7];
    outCr[7] = T1r * T1y - T1z * T1G; // FNMS
    outCi[7] = T1r * T1G + T1z * T1y; // FMA

    const T1U = T1M + T1N;
    const T1W = T1R - T1Q;
    const T1T = Wc[9], T1V = Ws[9];
    outCr[9] = T1T * T1U - T1V * T1W; // FNMS
    outCi[9] = T1T * T1W + T1V * T1U; // FMA

    const T1I = T1u + T1x;
    const T1K = T1F - T1C;
    const T1H = Wc[3], T1J = Ws[3];
    outCr[3] = T1H * T1I - T1J * T1K; // FNMS
    outCi[3] = T1H * T1K + T1J * T1I; // FMA

    const T1O = T1M - T1N;
    const T1S = T1Q + T1R;
    const T1L = Wc[1], T1P = Ws[1];
    outCr[1] = T1L * T1O - T1P * T1S; // FNMS
    outCi[1] = T1L * T1S + T1P * T1O; // FMA
  }

  return [outCr, outCi];
}

module.exports = { hb_10 };
