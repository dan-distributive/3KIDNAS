'use strict';

// =============================================================================
// r2cfII_20.js -- faithful JS port of rdft/scalar/r2cf/r2cfII_20.c (non-FMA
// branch). R2HCII ("shifted") radix-20 direct codelet -- the "cldm"
// middle-column combine for an EVEN outer radix (see r2cfII_8.js's header
// for the general even-radix derivation: clean 10-and-10 Cr/Ci split,
// out[19-k]=Ci[k]).
// INPUT: ph[0..19], R0[k]=ph[2k] (k=0..9), R1[k]=ph[2k+1] (k=0..9).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP572061402 = 0.572061402817684297600072783580302076536153377;
const KP218508012 = 0.218508012224410535399650602527877556893735408;
const KP309016994 = 0.309016994374947424102293417182819058860154590;
const KP809016994 = 0.809016994374947424102293417182819058860154590;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP250000000 = 0.250000000000000000000000000000000000000000000;
const KP176776695 = 0.176776695296636881100211090526212259821208984;
const KP395284707 = 0.395284707521047416499861693054089816714944392;
const KP672498511 = 0.672498511963957326960058968885748755876783111;
const KP415626937 = 0.415626937777453428589967464113135184222253485;
const KP707106781 = 0.707106781186547524400844362104849039284835938;

function r2cfII_20(ph) {
  const T8 = ph[5];
  const TD = KP707106781 * T8;
  const Tm = ph[15];
  const TN = KP707106781 * Tm;

  const T9 = ph[13];
  const Ta = ph[17];
  const TA = T9 + Ta;
  const Tb = ph[1];
  const Tc = ph[9];
  const Td = Tb + Tc;
  const TB = Tb - Tc;
  const TC = KP415626937 * TA + KP672498511 * TB; // FMA
  const TY = KP672498511 * TA - KP415626937 * TB; // FNMS
  const TE = KP395284707 * (Ta - Td);
  const Te = Ta + Td;
  const TF = KP176776695 * Te;

  const Tg = ph[3];
  const Tl = ph[7];
  const TJ = Tg + Tl;
  const Th = ph[11];
  const Ti = ph[19];
  const Tj = Th + Ti;
  const TI = Th - Ti;
  const TK = KP672498511 * TI - KP415626937 * TJ; // FNMS
  const T12 = KP415626937 * TI + KP672498511 * TJ; // FMA
  const TL = KP395284707 * (Tg - Tj);
  const Tk = Tg + Tj;
  const TM = KP176776695 * Tk;

  const T2 = ph[12];
  const T5 = ph[16];
  const T3 = ph[4];
  const T4 = ph[8];
  const T1a = T4 + T2;
  const T1b = T5 + T3;
  const T6 = T2 + T3 - (T4 + T5);
  const T1 = ph[0];
  const Tq = KP250000000 * T6 + T1; // FMA
  const T1l = KP587785252 * T1a - KP951056516 * T1b; // FNMS
  const T1c = KP951056516 * T1a + KP587785252 * T1b; // FMA
  const Tp = KP559016994 * (T5 + T2 - (T4 + T3));

  const T1f = ph[10];

  const Tv = ph[18];
  const Tw = ph[2];
  const Tx = Tv - Tw;
  const T1e = Tv + Tw;
  const Ts = ph[6];
  const Tt = ph[14];
  const Tu = Ts - Tt;
  const T1d = Ts + Tt;

  const Ty = KP951056516 * Tu + KP587785252 * Tx; // FMA
  const TW = KP587785252 * Tu - KP951056516 * Tx; // FNMS
  const T1g = KP809016994 * T1d + KP309016994 * T1e + T1f; // FMA + add
  const T1m = (T1f - KP809016994 * T1e) - KP309016994 * T1d; // FNMS - term

  const out = new Float64Array(20);

  {
    const T7 = T1 - T6;
    const T1r = T1e + T1f - T1d;
    const Tf = T8 + (T9 - Te);
    const Tn = (Tk - Tl) - Tm;
    const To = KP707106781 * (Tf + Tn);
    const T1q = KP707106781 * (Tf - Tn);
    out[2] = T7 - To;
    out[17] = T1q - T1r;
    out[7] = T7 + To;
    out[12] = T1q + T1r;
  }
  {
    const T1h = T1c - T1g;
    const T1j = T1c + T1g;
    const TV = Tq - Tp;
    const TX = TV - TW;
    const T15 = TV + TW;
    const TZ = KP218508012 * T9 + TD + TF - TE; // FMA + adds
    const T10 = TY + TZ;
    const T16 = TZ - TY;
    const T11 = (TL - KP218508012 * Tl) - (TM + TN); // FNMS - term
    const T13 = T11 - T12;
    const T17 = T11 + T12;

    const T14 = T10 + T13;
    out[5] = TX - T14;
    out[4] = TX + T14;
    const T19 = T17 - T16;
    out[14] = T19 - T1h;
    out[15] = T19 + T1h;
    const T18 = T16 + T17;
    out[9] = T15 - T18;
    out[0] = T15 + T18;
    const T1i = T13 - T10;
    out[19] = T1i - T1j;
    out[10] = T1i + T1j;
  }
  {
    const T1n = T1l + T1m;
    const T1p = T1m - T1l;
    const Tr = Tp + Tq;
    const Tz = Tr + Ty;
    const TR = Tr - Ty;
    const TG = TD + TE + (TF - KP572061402 * T9); // FNMS + adds
    const TH = TC + TG;
    const TS = TC - TG;
    const TO = TL + TM + (TN - KP572061402 * Tl); // FNMS + adds
    const TP = TK - TO;
    const TT = TK + TO;

    const TQ = TH + TP;
    out[6] = Tz - TQ;
    out[3] = Tz + TQ;
    const T1o = TT - TS;
    out[13] = T1o - T1p;
    out[16] = T1o + T1p;
    const TU = TS + TT;
    out[8] = TR - TU;
    out[1] = TR + TU;
    const T1k = TP - TH;
    out[11] = T1k - T1n;
    out[18] = T1k + T1n;
  }

  return out;
}

module.exports = { r2cfII_20 };
