'use strict';

// =============================================================================
// n1_10.js -- faithful JS port of dft/scalar/codelets/n1_10.c (non-FMA
// branch), FFTW3's direct (base-case) radix-10 complex DFT codelet.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP951056516 = 0.951056516295153572116439333379382143405698634;

function n1_10(ri, ii) {
  const T1 = ri[0], T2 = ri[5];
  const T3 = T1 - T2;
  const Tj = T1 + T2;
  const TO = ii[0], TP = ii[5];
  const TQ = TO - TP;
  const T1e = TO + TP;

  const T4 = ri[2], T5 = ri[7];
  const T6 = T4 - T5;
  const Tk = T4 + T5;
  const Te = ri[6], Tf = ri[1];
  const Tg = Te - Tf;
  const To = Te + Tf;
  const T7 = ri[8], T8 = ri[3];
  const T9 = T7 - T8;
  const Tl = T7 + T8;
  const Tb = ri[4], Tc = ri[9];
  const Td = Tb - Tc;
  const Tn = Tb + Tc;

  const TU = T6 - T9;
  const TV = Td - Tg;
  const T1c = Tk - Tl;
  const T1b = Tn - To;
  const Tm = Tk + Tl;
  const Tp = Tn + To;
  const Tq = Tm + Tp;
  const Ta = T6 + T9;
  const Th = Td + Tg;
  const Ti = Ta + Th;

  const Tu = ii[2], Tv = ii[7];
  const Tw = Tu - Tv;
  const T15 = Tu + Tv;
  const TE = ii[6], TF = ii[1];
  const TG = TE - TF;
  const T13 = TE + TF;
  const Tx = ii[8], Ty = ii[3];
  const Tz = Tx - Ty;
  const T16 = Tx + Ty;
  const TB = ii[4], TC = ii[9];
  const TD = TB - TC;
  const T12 = TB + TC;

  const TA = Tw - Tz;
  const TH = TD - TG;
  const T17 = T15 - T16;
  const T14 = T12 - T13;
  const T1f = T15 + T16;
  const T1g = T12 + T13;
  const T1h = T1f + T1g;
  const TL = Tw + Tz;
  const TM = TD + TG;
  const TR = TL + TM;

  const ro = new Float64Array(10), io = new Float64Array(10);
  ro[5] = T3 + Ti;
  io[5] = TQ + TR;
  ro[0] = Tj + Tq;
  io[0] = T1e + T1h;

  const TI = KP951056516 * TA + KP587785252 * TH;
  const TK = KP951056516 * TH - KP587785252 * TA;
  const Tr = KP559016994 * (Ta - Th);
  const Ts = T3 - KP250000000 * Ti;
  const Tt = Tr + Ts;
  const TJ = Ts - Tr;
  ro[9] = Tt - TI;
  ro[3] = TJ + TK;
  ro[1] = Tt + TI;
  ro[7] = TJ - TK;

  const TW = KP951056516 * TU + KP587785252 * TV;
  const TY = KP951056516 * TV - KP587785252 * TU;
  const TN = KP559016994 * (TL - TM);
  const TS = TQ - KP250000000 * TR;
  const TT = TN + TS;
  const TX = TS - TN;
  io[1] = TT - TW;
  io[7] = TY + TX;
  io[9] = TW + TT;
  io[3] = TX - TY;

  const T18 = KP951056516 * T14 - KP587785252 * T17;
  const T1a = KP951056516 * T17 + KP587785252 * T14;
  const TZ = Tj - KP250000000 * Tq;
  const T10 = KP559016994 * (Tm - Tp);
  const T11 = TZ - T10;
  const T19 = T10 + TZ;
  ro[2] = T11 - T18;
  ro[6] = T19 + T1a;
  ro[8] = T11 + T18;
  ro[4] = T19 - T1a;

  const T1d = KP951056516 * T1b - KP587785252 * T1c;
  const T1l = KP951056516 * T1c + KP587785252 * T1b;
  const T1i = T1e - KP250000000 * T1h;
  const T1j = KP559016994 * (T1f - T1g);
  const T1k = T1i - T1j;
  const T1m = T1j + T1i;
  io[2] = T1d + T1k;
  io[6] = T1m - T1l;
  io[8] = T1k - T1d;
  io[4] = T1l + T1m;

  return [ro, io];
}

module.exports = { n1_10 };
