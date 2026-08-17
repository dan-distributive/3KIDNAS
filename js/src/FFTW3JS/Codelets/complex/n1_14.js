'use strict';

// =============================================================================
// n1_14.js -- faithful JS port of dft/scalar/codelets/n1_14.c (non-FMA
// branch), FFTW3's direct (base-case) radix-14 complex DFT codelet.
// Every combining expression here sums at most 2 top-level chunks (one
// macro result + one plain term, or two macro results), so -- unlike
// n1_11.js -- there's no risk of the "3-way flattened sum changes rounding"
// bug; still, each macro call gets its own named intermediate for clarity
// and to keep the pattern consistent across files.
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP222520933 = 0.222520933956314404288902564496794759466355569;
const KP900968867 = 0.900968867902419126236102319507445051165919162;
const KP623489801 = 0.623489801858733530525004884004239810632274731;
const KP433883739 = 0.433883739117558120475768332848358754609990728;
const KP781831482 = 0.781831482468029808708444526674057750232334519;
const KP974927912 = 0.974927912181823607018131682993931217232785801;

function n1_14(ri, ii) {
  const T1 = ri[0], T2 = ri[7];
  const T3 = T1 - T2, Tp = T1 + T2;
  const T14 = ii[0], T15 = ii[7];
  const T16 = T14 - T15, T1f = T14 + T15;

  const T4 = ri[2], T5 = ri[9];
  const T6 = T4 - T5, Tq = T4 + T5;
  const T7 = ri[12], T8 = ri[5];
  const T9 = T7 - T8, Tr = T7 + T8;
  const Ta = T6 + T9, T1q = Tr - Tq, Ts = Tq + Tr, T10 = T9 - T6;

  const TA = ii[2], TB = ii[9];
  const TC = TA - TB, T1g = TA + TB;
  const TD = ii[12], TE = ii[5];
  const TF = TD - TE, T1h = TD + TE;
  const TG = TC - TF, T1z = T1g - T1h, T19 = TC + TF, T1i = T1g + T1h;

  const Tb = ri[4], Tc = ri[11];
  const Td = Tb - Tc, Tt = Tb + Tc;
  const Te = ri[10], Tf = ri[3];
  const Tg = Te - Tf, Tu = Te + Tf;
  const Th = Td + Tg, T1s = Tt - Tu, Tv = Tt + Tu, T12 = Tg - Td;

  const TO = ii[4], TP = ii[11];
  const TQ = TO - TP, T1m = TO + TP;
  const TR = ii[10], TS = ii[3];
  const TT = TR - TS, T1n = TR + TS;
  const TU = TQ - TT, T1B = T1n - T1m, T17 = TQ + TT, T1o = T1m + T1n;

  const Ti = ri[6], Tj = ri[13];
  const Tk = Ti - Tj, Tw = Ti + Tj;
  const Tl = ri[8], Tm = ri[1];
  const Tn = Tl - Tm, Tx = Tl + Tm;
  const To = Tk + Tn, T1r = Tw - Tx, Ty = Tw + Tx, T11 = Tn - Tk;

  const TH = ii[6], TI = ii[13];
  const TJ = TH - TI, T1j = TH + TI;
  const TK = ii[8], TL = ii[1];
  const TM = TK - TL, T1k = TK + TL;
  const TN = TJ - TM, T1A = T1k - T1j, T18 = TJ + TM, T1l = T1j + T1k;

  const ro = new Float64Array(14), io = new Float64Array(14);

  ro[7] = T3 + Ta + Th + To;
  io[7] = T16 + T19 + T17 + T18;
  ro[0] = Tp + Ts + Tv + Ty;
  io[0] = T1f + T1i + T1o + T1l;

  {
    const M1 = KP974927912 * TG - KP781831482 * TN;
    const TV = M1 - KP433883739 * TU;
    const N1 = T3 + KP623489801 * To;
    const N2 = -(KP900968867 * Th + KP222520933 * Ta);
    const Tz = N1 + N2;
    ro[5] = Tz - TV;
    ro[9] = Tz + TV;
    const M2 = KP974927912 * T10 - KP781831482 * T11;
    const T1e = M2 - KP433883739 * T12;
    const N3 = T16 + KP623489801 * T18;
    const N4 = -(KP900968867 * T17 + KP222520933 * T19);
    const T1d = N3 + N4;
    io[5] = T1d - T1e;
    io[9] = T1e + T1d;
  }

  {
    const M1 = KP781831482 * TG + KP974927912 * TU;
    const TX = M1 + KP433883739 * TN;
    const N1 = T3 + KP623489801 * Ta;
    const N2 = -(KP900968867 * To + KP222520933 * Th);
    const TW = N1 + N2;
    ro[13] = TW - TX;
    ro[1] = TW + TX;
    const M2 = KP781831482 * T10 + KP974927912 * T12;
    const T1b = M2 + KP433883739 * T11;
    const N3 = T16 + KP623489801 * T19;
    const N4 = -(KP900968867 * T18 + KP222520933 * T17);
    const T1c = N3 + N4;
    io[1] = T1b + T1c;
    io[13] = T1c - T1b;
  }

  {
    const M1 = KP433883739 * TG + KP974927912 * TN;
    const TZ = M1 - KP781831482 * TU;
    const N1 = T3 + KP623489801 * Th;
    const N2 = -(KP222520933 * To + KP900968867 * Ta);
    const TY = N1 + N2;
    ro[11] = TY - TZ;
    ro[3] = TY + TZ;
    const M2 = KP433883739 * T10 + KP974927912 * T11;
    const T13 = M2 - KP781831482 * T12;
    const N3 = T16 + KP623489801 * T17;
    const N4 = -(KP222520933 * T18 + KP900968867 * T19);
    const T1a = N3 + N4;
    io[3] = T13 + T1a;
    io[11] = T1a - T13;
  }

  {
    const M1 = KP781831482 * T1q - KP433883739 * T1r;
    const T1t = M1 - KP974927912 * T1s;
    const N1 = T1f + KP623489801 * T1i;
    const N2 = -(KP900968867 * T1l + KP222520933 * T1o);
    const T1p = N1 + N2;
    io[6] = T1p - T1t;
    io[8] = T1t + T1p;
    const M2 = KP781831482 * T1z - KP433883739 * T1A;
    const T1C = M2 - KP974927912 * T1B;
    const N3 = Tp + KP623489801 * Ts;
    const N4 = -(KP900968867 * Ty + KP222520933 * Tv);
    const T1y = N3 + N4;
    ro[6] = T1y - T1C;
    ro[8] = T1y + T1C;
  }

  {
    const M1 = KP433883739 * T1q + KP781831482 * T1s;
    const T1v = M1 - KP974927912 * T1r;
    const N1 = T1f + KP623489801 * T1o;
    const N2 = -(KP222520933 * T1l + KP900968867 * T1i);
    const T1u = N1 + N2;
    io[4] = T1u - T1v;
    io[10] = T1v + T1u;
    const M2 = KP433883739 * T1z + KP781831482 * T1B;
    const T1E = M2 - KP974927912 * T1A;
    const N3 = Tp + KP623489801 * Tv;
    const N4 = -(KP222520933 * Ty + KP900968867 * Ts);
    const T1D = N3 + N4;
    ro[4] = T1D - T1E;
    ro[10] = T1D + T1E;
  }

  {
    const M1 = KP974927912 * T1q + KP433883739 * T1s;
    const T1w = M1 + KP781831482 * T1r;
    const N1 = T1f + KP623489801 * T1l;
    const N2 = -(KP900968867 * T1o + KP222520933 * T1i);
    const T1x = N1 + N2;
    io[2] = T1w + T1x;
    io[12] = T1x - T1w;
    const M2 = KP974927912 * T1z + KP433883739 * T1B;
    const T1G = M2 + KP781831482 * T1A;
    const N3 = Tp + KP623489801 * Ty;
    const N4 = -(KP900968867 * Tv + KP222520933 * Ts);
    const T1F = N3 + N4;
    ro[12] = T1F - T1G;
    ro[2] = T1F + T1G;
  }

  return [ro, io];
}

module.exports = { n1_14 };
