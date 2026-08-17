'use strict';

// =============================================================================
// r2cf_14.js -- faithful JS port of rdft/scalar/r2cf/r2cf_14.c (non-FMA
// branch). x[0..13] (real) -> O[0..13] packed halfcomplex: O[0..7]=Re0..Re7,
// O[8..13]=Im6..Im1 (O[n-k]=Im_k). R0/R1 (stride-2) input convention -- see
// r2cf_6.js's header. FMA(a,b,c)=a*b+c, FNMA(a,b,c)=-(a*b+c),
// FNMS(a,b,c)=c-a*b -- named intermediates used wherever two such macro
// results are combined via a bare +/- (the "n1_11 lesson": flattening would
// silently regroup the floating-point sum).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP900968867 = 0.900968867902419126236102319507445051165919162;
const KP222520933 = 0.222520933956314404288902564496794759466355569;
const KP623489801 = 0.623489801858733530525004884004239810632274731;
const KP433883739 = 0.433883739117558120475768332848358754609990728;
const KP974927912 = 0.974927912181823607018131682993931217232785801;
const KP781831482 = 0.781831482468029808708444526674057750232334519;

function r2cf_14(x) {
  const T1 = x[0];
  const T2 = x[7];
  const T3 = T1 - T2;
  const TB = T1 + T2;
  const T4 = x[4];
  const T5 = x[11];
  const T6 = T4 - T5;
  const Tv = T4 + T5;
  const Tl = x[12];
  const Tm = x[5];
  const Tn = Tl - Tm;
  const Ts = Tl + Tm;
  const Ti = x[2];
  const Tj = x[9];
  const Tk = Ti - Tj;
  const Tt = Ti + Tj;
  const Tb = x[6];
  const Tc = x[13];
  const Td = Tb - Tc;
  const Ty = Tb + Tc;
  const T7 = x[10];
  const T8 = x[3];
  const T9 = T7 - T8;
  const Tw = T7 + T8;
  const Te = x[8];
  const Tf = x[1];
  const Tg = Te - Tf;
  const Tz = Te + Tf;

  const O = new Float64Array(14);

  const Tp = Tn - Tk;
  const Tr = Tg - Td;
  const Tq = T9 - T6;
  O[13] = KP781831482 * Tp + KP974927912 * Tq + KP433883739 * Tr;
  O[9] = KP433883739 * Tq + KP781831482 * Tr - KP974927912 * Tp;
  O[11] = KP433883739 * Tp + KP974927912 * Tr - KP781831482 * Tq;
  const Ta = T6 + T9;
  const To = Tk + Tn;
  const Th = Td + Tg;
  const M1 = KP623489801 * Ta + T3;
  const M2 = -(KP222520933 * Th + KP900968867 * To);
  O[3] = M1 + M2;
  O[7] = T3 + To + Ta + Th;
  const M3 = KP623489801 * To + T3;
  const M4 = -(KP900968867 * Th + KP222520933 * Ta);
  O[1] = M3 + M4;
  const M5 = KP623489801 * Th + T3;
  const M6 = -(KP900968867 * Ta + KP222520933 * To);
  O[5] = M5 + M6;

  const Tu = Ts - Tt;
  const TA = Ty - Tz;
  const Tx = Tv - Tw;
  O[12] = KP974927912 * Tu + KP433883739 * Tx + KP781831482 * TA;
  O[8] = KP974927912 * Tx + KP433883739 * TA - KP781831482 * Tu;
  const Tfnms4 = KP974927912 * TA - KP781831482 * Tx; // FNMS(KP781831482,Tx,KP974927912*TA)
  O[10] = Tfnms4 - KP433883739 * Tu;
  const TC = Tt + Ts;
  const TE = Tv + Tw;
  const TD = Ty + Tz;
  const M7 = KP623489801 * TC + TB;
  const M8 = -(KP900968867 * TD + KP222520933 * TE);
  O[6] = M7 + M8;
  const M9 = KP623489801 * TD + TB;
  const M10 = -(KP900968867 * TE + KP222520933 * TC);
  O[2] = M9 + M10;
  const M11 = KP623489801 * TE + TB;
  const M12 = -(KP222520933 * TD + KP900968867 * TC);
  O[4] = M11 + M12;
  O[0] = TB + TC + TE + TD;

  return O;
}

module.exports = { r2cf_14 };
