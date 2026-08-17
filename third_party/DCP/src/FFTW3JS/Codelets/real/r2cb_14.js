'use strict';

// =============================================================================
// r2cb_14.js -- faithful JS port of rdft/scalar/r2cb/r2cb_14.c (non-FMA
// branch). O[0..13] packed halfcomplex -> x[0..13] real, UNNORMALIZED.
// Output uses the R0/R1 (stride-2) convention -- see r2cb_6.js's header.
// FMA(a,b,c)=a*b+c, FNMA(a,b,c)=-(a*b+c), FNMS(a,b,c)=c-a*b -- named
// intermediates used wherever two such macro results combine via a bare
// +/- (the "n1_11 lesson").
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const KP1_801937735 = 1.801937735804838252472204639014890102331838324;
const KP445041867 = 0.445041867912628808577805128993589518932711138;
const KP1_246979603 = 1.246979603717467061050009768008479621264549462;
const KP867767478 = 0.867767478235116240951536665696717509219981456;
const KP1_949855824 = 1.949855824363647214036263365987862434465571601;
const KP1_563662964 = 1.563662964936059617416889053348115500464669037;
const KP2_000000000 = 2.000000000000000000000000000000000000000000000;

function r2cb_14(O) {
  const T1 = O[0];
  const T2 = O[7];
  const T3 = T1 - T2;
  const Td = T1 + T2;
  const T4 = O[2];
  const T5 = O[5];
  const T6 = T4 - T5;
  const Te = T4 + T5;
  const To = O[12];
  const Tp = O[9];
  const Tq = To - Tp;
  const Tz = To + Tp;
  const Tl = O[8];
  const Tm = O[13];
  const Tn = Tl - Tm;
  const Ty = Tl + Tm;
  const Ta = O[6];
  const Tb = O[1];
  const Tc = Ta - Tb;
  const Tg = Ta + Tb;
  const Ti = O[10];
  const Tj = O[11];
  const Tk = Ti - Tj;
  const Tx = Ti + Tj;
  const T7 = O[4];
  const T8 = O[3];
  const T9 = T7 - T8;
  const Tf = T7 + T8;

  const x = new Float64Array(14);
  x[7] = KP2_000000000 * (T6 + T9 + Tc) + T3; // R1[WS(rs,3)]
  x[0] = KP2_000000000 * (Te + Tf + Tg) + Td; // R0[0]

  {
    const Y1 = KP1_563662964 * Tk - KP1_949855824 * Tn; // FNMS
    const Tr = Y1 - KP867767478 * Tq;
    const X1 = KP1_246979603 * Tf + Td; // FMA
    const X2 = -(KP445041867 * Tg + KP1_801937735 * Te); // FNMA
    const Th = X1 + X2;
    x[4] = Th - Tr; // R0[WS(rs,2)]
    x[10] = Th + Tr; // R0[WS(rs,5)]
    const Y2 = KP867767478 * Tx + KP1_563662964 * Ty; // FMA
    const TE = Y2 - KP1_949855824 * Tz;
    const X3 = KP1_246979603 * Tc + T3; // FMA
    const X4 = -(KP1_801937735 * T9 + KP445041867 * T6); // FNMA
    const TD = X3 + X4;
    x[5] = TD - TE; // R1[WS(rs,2)]
    x[9] = TD + TE; // R1[WS(rs,4)]
  }
  {
    const Y3 = KP867767478 * Tk + KP1_563662964 * Tn; // FMA
    const Tt = Y3 - KP1_949855824 * Tq;
    const X5 = KP1_246979603 * Tg + Td; // FMA
    const X6 = -(KP1_801937735 * Tf + KP445041867 * Te); // FNMA
    const Ts = X5 + X6;
    x[12] = Ts - Tt; // R0[WS(rs,6)]
    x[2] = Ts + Tt; // R0[WS(rs,1)]
    const Y4 = KP1_563662964 * Tx - KP1_949855824 * Ty; // FNMS
    const TA = Y4 - KP867767478 * Tz;
    const X7 = KP1_246979603 * T9 + T3; // FMA
    const X8 = -(KP445041867 * Tc + KP1_801937735 * T6); // FNMA
    const Tw = X7 + X8;
    x[11] = Tw - TA; // R1[WS(rs,5)]
    x[3] = Tw + TA; // R1[WS(rs,1)]
  }
  {
    const Y5 = KP1_563662964 * Tz + KP1_949855824 * Tx; // FMA
    const TC = Y5 + KP867767478 * Ty;
    const X9 = KP1_246979603 * T6 + T3; // FMA
    const X10 = -(KP1_801937735 * Tc + KP445041867 * T9); // FNMA
    const TB = X9 + X10;
    x[1] = TB - TC; // R1[0]
    x[13] = TB + TC; // R1[WS(rs,6)]
    const Y6 = KP1_563662964 * Tq + KP1_949855824 * Tk; // FMA
    const Tv = Y6 + KP867767478 * Tn;
    const X11 = KP1_246979603 * Te + Td; // FMA
    const X12 = -(KP1_801937735 * Tg + KP445041867 * Tf); // FNMA
    const Tu = X11 + X12;
    x[8] = Tu - Tv; // R0[WS(rs,4)]
    x[6] = Tu + Tv; // R0[WS(rs,3)]
  }

  return x;
}

module.exports = { r2cb_14 };
