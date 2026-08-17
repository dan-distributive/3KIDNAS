'use strict';

// =============================================================================
// n1_7.js -- faithful JS port of dft/scalar/codelets/n1_7.c (non-FMA branch),
// FFTW3's direct (base-case) radix-7 complex DFT codelet.
// FMA(a,b,c)=a*b+c, FNMS(a,b,c)=c-a*b, FNMA(a,b,c)=-(a*b+c) macro-expand to
// plain arithmetic here (kernel/ifftw.h's exact definitions, confirmed by
// direct source read -- FNMA is new at this radix, not seen in radix 2-5).
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP222520933 = 0.222520933956314404288902564496794759466355569;
const KP900968867 = 0.900968867902419126236102319507445051165919162;
const KP623489801 = 0.623489801858733530525004884004239810632274731;
const KP433883739 = 0.433883739117558120475768332848358754609990728;
const KP781831482 = 0.781831482468029808708444526674057750232334519;
const KP974927912 = 0.974927912181823607018131682993931217232785801;

function n1_7(ri, ii) {
  const T1 = ri[0], Tu = ii[0];

  const T2 = ri[1], T3 = ri[6];
  const T4 = T2 + T3;
  const Tq = T3 - T2;
  const Tc = ii[1], Td = ii[6];
  const Te = Tc - Td;
  const Tx = Tc + Td;

  const T5 = ri[2], T6 = ri[5];
  const T7 = T5 + T6;
  const Ts = T6 - T5;
  const Ti = ii[2], Tj = ii[5];
  const Tk = Ti - Tj;
  const Tv = Ti + Tj;

  const T8 = ri[3], T9 = ri[4];
  const Ta = T8 + T9;
  const Tr = T9 - T8;
  const Tf = ii[3], Tg = ii[4];
  const Th = Tf - Tg;
  const Tw = Tf + Tg;

  const ro = new Float64Array(7), io = new Float64Array(7);
  ro[0] = T1 + T4 + T7 + Ta;
  io[0] = Tu + Tx + Tv + Tw;

  const Tl = (KP974927912 * Te - KP781831482 * Th) - KP433883739 * Tk;
  const Tb = (KP623489801 * Ta + T1) - (KP900968867 * T7 + KP222520933 * T4);
  ro[5] = Tb - Tl;
  ro[2] = Tb + Tl;
  const TB = (KP974927912 * Tq - KP781831482 * Tr) - KP433883739 * Ts;
  const TC = (KP623489801 * Tw + Tu) - (KP900968867 * Tv + KP222520933 * Tx);
  io[2] = TB + TC;
  io[5] = TC - TB;

  const Tn = (KP781831482 * Te + KP974927912 * Tk) + KP433883739 * Th;
  const Tm = (KP623489801 * T4 + T1) - (KP900968867 * Ta + KP222520933 * T7);
  ro[6] = Tm - Tn;
  ro[1] = Tm + Tn;
  const Tz = (KP781831482 * Tq + KP974927912 * Ts) + KP433883739 * Tr;
  const TA = (KP623489801 * Tx + Tu) - (KP900968867 * Tw + KP222520933 * Tv);
  io[1] = Tz + TA;
  io[6] = TA - Tz;

  const Tp = (KP433883739 * Te + KP974927912 * Th) - KP781831482 * Tk;
  const To = (KP623489801 * T7 + T1) - (KP222520933 * Ta + KP900968867 * T4);
  ro[4] = To - Tp;
  ro[3] = To + Tp;
  const Tt = (KP433883739 * Tq + KP974927912 * Tr) - KP781831482 * Ts;
  const Ty = (KP623489801 * Tv + Tu) - (KP222520933 * Tw + KP900968867 * Tx);
  io[3] = Tt + Ty;
  io[4] = Ty - Tt;

  return [ro, io];
}

module.exports = { n1_7 };
