'use strict';

// =============================================================================
// r2cf_15.js -- faithful JS port of rdft/scalar/r2cf/r2cf_15.c (non-FMA
// branch). x[0..14] (real) -> O[0..14] packed halfcomplex: O[0..7]=Re0..
// Re7, O[8..14]=Im7..Im1 (O[n-k]=Im_k, same convention as prior r2cf_*.js).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP484122918 = 0.484122918275927110647408174972799951354115213;
const KP216506350 = 0.216506350946109661690930792688234045867850657;
const KP951056516 = 0.951056516295153572116439333379382143405698634;
const KP587785252 = 0.587785252292473129168705954639072768597652438;
const KP250000000 = 0.25;
const KP559016994 = 0.559016994374947424102293417182819058860154590;
const KP509036960 = 0.509036960455127183450980863393907648510733164;
const KP823639103 = 0.823639103546331925877420039278190003029660514;
const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP500000000 = 0.5;

function r2cf_15(x) {
  const TJ = x[0];
  const Tg = x[10], Th = x[5];
  const TK = Th + Tg;
  const Ti = Tg - Th;
  const TR = TJ + TK;
  const TL = TJ - KP500000000 * TK;

  const Tm = x[3], Tt = x[6], Tw = x[9], Tp = x[12];

  const T1 = x[14], T2 = x[4];
  const T3 = T1 - T2, Tx = T1 + T2;
  const T8 = x[13], T9 = x[8];
  const Ta = T8 - T9, Tn = T9 + T8;

  const Tb = x[7], Tc = x[2];
  const Td = Tb - Tc, Tq = Tc + Tb;
  const T4 = x[1], T5 = x[11];
  const T6 = T4 - T5, Tu = T5 + T4;

  const TD = Ta - Td;
  const TE = T6 + T3;
  const T7 = T3 - T6;
  const Te = Ta + Td;
  const Tf = T7 - Te;
  const TV = Tt + Tu;
  const TW = Tw + Tx;
  const TX = TV + TW;
  const Tv = Tt - KP500000000 * Tu;
  const Ty = Tw - KP500000000 * Tx;
  const TH = Tv + Ty;
  const To = Tm - KP500000000 * Tn;
  const Tr = Tp - KP500000000 * Tq;
  const TG = To + Tr;
  const TS = Tm + Tn;
  const TT = Tp + Tq;
  const TU = TS + TT;

  const O = new Float64Array(15);
  O[10] = KP866025403 * (Tf - Ti);

  const TF = KP823639103 * TD + KP509036960 * TE;
  const TP = KP823639103 * TE - KP509036960 * TD;
  const TI = KP559016994 * (TG - TH);
  const TM = TG + TH;
  const TN = TL - KP250000000 * TM;
  O[5] = TL + TM;
  const TQ = TN - TI;
  O[2] = TP + TQ;
  O[7] = TQ - TP;
  const TO = TI + TN;
  O[1] = TF + TO;
  O[4] = TO - TF;

  const T11 = TS - TT;
  const T12 = TW - TV;
  O[12] = KP587785252 * T11 + KP951056516 * T12;
  O[9] = KP587785252 * T12 - KP951056516 * T11;
  const T10 = KP559016994 * (TU - TX);
  const TY = TU + TX;
  const TZ = TR - KP250000000 * TY;
  O[3] = TZ - T10;
  O[0] = TR + TY;
  O[6] = T10 + TZ;

  const Tj = KP866025403 * Ti + KP216506350 * Tf;
  const Tk = KP484122918 * (Te + T7);
  const Tl = Tj + Tk;
  const TB = Tk - Tj;
  const Ts = To - Tr;
  const Tz = Tv - Ty;
  const TA = KP951056516 * Ts + KP587785252 * Tz;
  const TC = KP951056516 * Tz - KP587785252 * Ts;

  O[14] = Tl - TA;
  O[8] = TC - TB;
  O[11] = Tl + TA;
  O[13] = TB + TC;

  return O;
}

module.exports = { r2cf_15 };
