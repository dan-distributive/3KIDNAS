'use strict';

// =============================================================================
// r2cf_20.js -- faithful JS port of rdft/scalar/r2cf/r2cf_20.c (non-FMA
// branch). x[0..19] (real) -> O[0..19] packed halfcomplex: O[0..10]=Re0..
// Re10, O[11..19]=Im9..Im1 (O[n-k]=Im_k). R0/R1 (stride-2) input convention
// -- see r2cf_6.js's header.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.250000000000000000000000000000000000000000000;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;

function r2cf_20(x) {
  const T1 = x[0];
  const T2 = x[10];
  const T15 = T1 + T2;
  const TD = x[15];
  const TE = x[5];
  const T16 = TE + TD;
  const T3 = T1 - T2;
  const T1m = T15 + T16;
  const TF = TD - TE;
  const T17 = T15 - T16;

  const T4 = x[4];
  const T5 = x[14];
  const T6 = T4 - T5;
  const TU = T4 + T5;
  const Tt = x[17];
  const Tu = x[7];
  const Tv = Tt - Tu;
  const T12 = Tt + Tu;

  const Tw = x[13];
  const Tx = x[3];
  const Ty = Tw - Tx;
  const TZ = Tw + Tx;
  const T7 = x[16];
  const T8 = x[6];
  const T9 = T7 - T8;
  const TR = T7 + T8;

  const Tb = x[8];
  const Tc = x[18];
  const Td = Tb - Tc;
  const TY = Tb + Tc;
  const Tm = x[1];
  const Tn = x[11];
  const To = Tm - Tn;
  const TS = Tm + Tn;

  const Tp = x[9];
  const Tq = x[19];
  const Tr = Tp - Tq;
  const TV = Tp + Tq;
  const Te = x[12];
  const Tf = x[2];
  const Tg = Te - Tf;
  const T11 = Te + Tf;

  const Ts = To - Tr;
  const TM = T6 - T9;
  const TN = Td - Tg;
  const Tz = Tv - Ty;
  const Ta = T6 + T9;
  const Th = Td + Tg;
  const Ti = Ta + Th;
  const T1g = TY + TZ;
  const T1h = T11 + T12;
  const T1k = T1g + T1h;
  const T10 = TY - TZ;
  const T13 = T11 - T12;
  const T19 = T10 + T13;
  const TG = Tr + To;
  const TH = Ty + Tv;
  const TI = TG + TH;
  const T1d = TU + TV;
  const T1e = TR + TS;
  const T1j = T1d + T1e;
  const TT = TR - TS;
  const TW = TU - TV;
  const T18 = TW + TT;

  const O = new Float64Array(20);
  O[5] = T3 + Ti;
  O[15] = TF - TI;

  const TX = TT - TW;
  const T14 = T10 - T13;
  O[14] = KP951056516 * TX - KP587785252 * T14; // FNMS
  O[18] = KP587785252 * TX + KP951056516 * T14; // FMA
  const T1f = T1d - T1e;
  const T1i = T1g - T1h;
  O[12] = KP587785252 * T1f - KP951056516 * T1i; // FNMS
  O[16] = KP951056516 * T1f + KP587785252 * T1i; // FMA

  const T1l = KP559016994 * (T1j - T1k);
  const T1n = T1j + T1k;
  const T1o = T1m - KP250000000 * T1n; // FNMS
  O[4] = T1l + T1o;
  O[0] = T1m + T1n;
  O[8] = T1o - T1l;
  const T1c = KP559016994 * (T18 - T19);
  const T1a = T18 + T19;
  const T1b = T17 - KP250000000 * T1a; // FNMS
  O[2] = T1b - T1c;
  O[10] = T17 + T1a;
  O[6] = T1c + T1b;

  const TA = KP951056516 * Ts + KP587785252 * Tz; // FMA
  const TC = KP951056516 * Tz - KP587785252 * Ts; // FNMS
  const Tj = KP559016994 * (Ta - Th);
  const Tk = T3 - KP250000000 * Ti; // FNMS
  const Tl = Tj + Tk;
  const TB = Tk - Tj;
  O[9] = Tl - TA;
  O[7] = TB + TC;
  O[1] = Tl + TA;
  O[3] = TB - TC;

  const TO = KP951056516 * TM + KP587785252 * TN; // FMA
  const TQ = KP951056516 * TN - KP587785252 * TM; // FNMS
  const TJ = KP250000000 * TI + TF; // FMA
  const TK = KP559016994 * (TH - TG);
  const TL = TJ + TK;
  const TP = TK - TJ;
  O[19] = TL - TO;
  O[13] = TQ + TP;
  O[11] = TO + TL;
  O[17] = TP - TQ;

  return O;
}

module.exports = { r2cf_20 };
