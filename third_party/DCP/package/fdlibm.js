module.declare(["./fma.js"], function (require, exports, module) {
'use strict';

// =============================================================================
// fdlibm.js
// Faithful JS port of fdlibm (Sun Microsystems 1993, via FreeBSD's lib/msun),
// double-precision sin/cos/log/atanh. Ported line-for-line, preserving exact
// operation order and coefficients, so that Fortran -- linking the SAME
// fdlibm C source instead of the platform's system libm -- and this JS port
// produce bit-identical results for the same double-precision input.
//
// WHY: IEEE 754 only mandates correctly-rounded results for +,-,*,/,sqrt.
// Transcendental functions (sin, cos, log, atanh, ...) are explicitly allowed
// up to ~1ULP of implementation-defined slop, and different platforms' libm
// (Apple's on macOS vs V8's own internal Math implementation) are NOT
// required to agree. That gap was traced this session to being the likely
// source of JS's amoeba optimizer following a different search trajectory
// than Fortran's on the same noisy objective function, landing in a
// different (but equally valid, per the objective function itself) local
// optimum for high-inclination bootstrap realizations.
//
// SCOPE: only the fast + medium argument-reduction paths from fdlibm's
// __ieee754_rem_pio2 are ported (covering |x| up to ~2^20*(pi/2), i.e. any
// physically meaningful angle by a huge margin). The large-argument
// Payne-Hanek path (__kernel_rem_pio2, for |x| approaching double's dynamic
// range) is intentionally NOT ported -- nothing in this pipeline ever
// produces an angle anywhere near that scale -- and calling into that regime
// throws loudly rather than silently returning a wrong answer.
//
// Source: https://github.com/freebsd/freebsd-src lib/msun/src/
//   s_sin.c, s_cos.c, k_sin.c, k_cos.c, e_rem_pio2.c, e_log.c, e_atanh.c
//   Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//   "Permission to use, copy, modify, and distribute this software is
//   freely granted, provided that this notice is preserved."
// =============================================================================

// ---------------------------------------------------------------------------
// IEEE754 double bit-word access. fdlibm's GET_HIGH_WORD/GET_LOW_WORD/
// INSERT_WORDS macros are endian-normalized in the C source (msw/lsw always
// mean "most/least significant 32 bits", regardless of host byte order) --
// a big-endian DataView read/write reproduces that directly.
// ---------------------------------------------------------------------------
const _buf = new ArrayBuffer(8);
const _dv  = new DataView(_buf);

function getHighWord(x) {
  _dv.setFloat64(0, x, false);
  return _dv.getInt32(0, false);
}
function getLowWordU(x) {
  _dv.setFloat64(0, x, false);
  return _dv.getUint32(4, false);
}
function insertWords(hi, lo) {
  _dv.setInt32(0, hi | 0, false);
  _dv.setUint32(4, lo >>> 0, false);
  return _dv.getFloat64(0, false);
}

// ---------------------------------------------------------------------------
// rnint(x): round-to-nearest-integer via the classic double-precision
// magic-constant trick (e_rem_pio2.c / math_private.h). Exact in JS doubles.
// ---------------------------------------------------------------------------
const RNINT_MAGIC = 6755399441055744.0; // 0x1.8p52 = 1.5 * 2^52
function rnint(x) {
  return (x + RNINT_MAGIC) - RNINT_MAGIC;
}

// =============================================================================
// __kernel_sin / __kernel_cos  (k_sin.c / k_cos.c)
// Valid for |x| <~ pi/4. Pure polynomial evaluation -- direct port.
// =============================================================================
const K_half = 5.00000000000000000000e-01;
const K_S1 = -1.66666666666666324348e-01;
const K_S2 =  8.33333333332248946124e-03;
const K_S3 = -1.98412698298579493134e-04;
const K_S4 =  2.75573137070700676789e-06;
const K_S5 = -2.50507602534068634195e-08;
const K_S6 =  1.58969099521155010221e-10;

// FMA NOTE: this project's actual C build (src/makeflags CFLAGS=-O, no
// -ffp-contract=off) compiles fdlibm_k_sin.c/fdlibm_k_cos.c WITH hardware
// FMA for these polynomial evaluations (verified via -S disassembly on the
// real object files) -- confirmed to change results by 1 ULP for real
// inputs (theta=0.18265073567382517084 differs between the FMA and
// non-FMA compiled versions). Since the goal is matching THIS project's
// actual, current Fortran build -- not a hypothetical portable one --
// these functions replicate the exact fusion graph the compiler chose
// (identified by reading the disassembly instruction-by-instruction),
// using fma.js's software FMA emulation (JS has no native FMA). See
// fma.js's header for the full rationale.
const { fma } = require('./fma.js');

function __kernel_sin(x, y, iy) {
  const z = x * x;
  const w = z * z;
  // r = S2+z*(S3+z*S4) + z*w*(S5+z*S6), fused as:
  const t1 = fma(z, K_S4, K_S3);        // S3+z*S4
  const t2 = fma(z, t1, K_S2);           // S2+z*t1
  const t3 = fma(z, K_S6, K_S5);         // S5+z*S6
  const zw = z * w;                      // plain multiply, not fused
  const r = fma(zw, t3, t2);             // t2+zw*t3
  const v = z * x;
  if (iy === 0) return fma(v, fma(z, r, K_S1), x); // x+v*(S1+z*r)
  // iy!=0: x-((z*(half*y-v*r)-y)-v*S1), fused per real -S -g disassembly as:
  const p1 = -(v * r);                // fnmul, plain (single rounding)
  const tb = fma(y, K_half, p1);      // fmadd: half*y + p1 == half*y-v*r
  const tc = fma(z, tb, -y);          // fnmsub: z*tb - y
  const td = fma(v, -K_S1, tc);       // fmadd: tc + v*(-S1) == tc-v*S1
  return x - td;                       // fsub, plain
}

const K_one =  1.00000000000000000000e+00;
const K_C1  =  4.16666666666666019037e-02;
const K_C2  = -1.38888888888741095749e-03;
const K_C3  =  2.48015872894767294178e-05;
const K_C4  = -2.75573143513906633035e-07;
const K_C5  =  2.08757232129817482790e-09;
const K_C6  = -1.13596475577881948265e-11;

function __kernel_cos(x, y) {
  const z = x * x;
  const w = z * z;
  // r = z*(C1+z*(C2+z*C3)) + w*w*(C4+z*(C5+z*C6)), fused as (see fdlibm.js's
  // __kernel_sin comment on why -- same real, verified FMA usage):
  const c1 = fma(z, K_C3, K_C2);         // C2+z*C3
  const c2 = fma(z, c1, K_C1);           // C1+z*c1
  const c3 = fma(z, K_C6, K_C5);         // C5+z*C6
  const c4 = fma(z, c3, K_C4);           // C4+z*c3
  const w2sq = w * w;                    // plain multiply, not fused
  const r = fma(z, c2, w2sq * c4);       // z*c2 + w2sq*c4
  const hz = 0.5 * z;
  const w2 = K_one - hz;
  // return w2+(((one-w2)-hz)+(z*r-x*y)), fused per real -S -g disassembly as:
  const oneMinusW2 = K_one - w2;         // fsub, plain
  const t = oneMinusW2 - hz;             // fsub, plain
  const negxy = -(x * y);                // fnmul, plain (single rounding)
  const zrMinusXy = fma(z, r, negxy);    // fmadd: z*r + negxy == z*r-x*y
  return w2 + (t + zrMinusXy);
}

// =============================================================================
// __ieee754_rem_pio2  (e_rem_pio2.c) -- fast + medium paths only, see SCOPE.
// Returns { n, y0, y1 } matching the C function's return value + y[0]/y[1].
// =============================================================================
const R_invpio2 =  6.36619772367581382433e-01;
const R_pio2_1  =  1.57079632673412561417e+00;
const R_pio2_1t =  6.07710050650619224932e-11;
const R_pio2_2  =  6.07710050630396597660e-11;
const R_pio2_2t =  2.02226624879595063154e-21;
const R_pio2_3  =  2.02226624871116645580e-21;
const R_pio2_3t =  8.47842766036889956997e-32;

function __ieee754_rem_pio2(x) {
  const hx = getHighWord(x);
  const ix = hx & 0x7fffffff;

  if (ix <= 0x400f6a7a) {           // |x| ~<= 5pi/4
    if ((ix & 0xfffff) === 0x921fb) return _remPio2Medium(x, ix);
    if (ix <= 0x4002d97c) {         // |x| ~<= 3pi/4
      if (hx > 0) {
        const z = x - R_pio2_1;
        const y0 = z - R_pio2_1t;
        const y1 = (z - y0) - R_pio2_1t;
        return { n: 1, y0, y1 };
      } else {
        const z = x + R_pio2_1;
        const y0 = z + R_pio2_1t;
        const y1 = (z - y0) + R_pio2_1t;
        return { n: -1, y0, y1 };
      }
    } else {
      if (hx > 0) {
        const z = x - 2 * R_pio2_1;
        const y0 = z - 2 * R_pio2_1t;
        const y1 = (z - y0) - 2 * R_pio2_1t;
        return { n: 2, y0, y1 };
      } else {
        const z = x + 2 * R_pio2_1;
        const y0 = z + 2 * R_pio2_1t;
        const y1 = (z - y0) + 2 * R_pio2_1t;
        return { n: -2, y0, y1 };
      }
    }
  }
  if (ix <= 0x401c463b) {           // |x| ~<= 9pi/4
    if (ix <= 0x4015fdbc) {         // |x| ~<= 7pi/4
      if (ix === 0x4012d97c) return _remPio2Medium(x, ix);
      if (hx > 0) {
        const z = x - 3 * R_pio2_1;
        const y0 = z - 3 * R_pio2_1t;
        const y1 = (z - y0) - 3 * R_pio2_1t;
        return { n: 3, y0, y1 };
      } else {
        const z = x + 3 * R_pio2_1;
        const y0 = z + 3 * R_pio2_1t;
        const y1 = (z - y0) + 3 * R_pio2_1t;
        return { n: -3, y0, y1 };
      }
    } else {
      if (ix === 0x401921fb) return _remPio2Medium(x, ix);
      if (hx > 0) {
        const z = x - 4 * R_pio2_1;
        const y0 = z - 4 * R_pio2_1t;
        const y1 = (z - y0) - 4 * R_pio2_1t;
        return { n: 4, y0, y1 };
      } else {
        const z = x + 4 * R_pio2_1;
        const y0 = z + 4 * R_pio2_1t;
        const y1 = (z - y0) + 4 * R_pio2_1t;
        return { n: -4, y0, y1 };
      }
    }
  }
  if (ix < 0x413921fb) {            // |x| ~< 2^20*(pi/2), medium size
    return _remPio2Medium(x, ix);
  }
  throw new Error(
    'fdlibm.__ieee754_rem_pio2: |x| out of the ported fast/medium range ' +
    '(large-argument Payne-Hanek path intentionally not ported -- see ' +
    'fdlibm.js SCOPE comment). x=' + x);
}

function _remPio2Medium(x, ix) {
  const fn = rnint(x * R_invpio2);
  const n  = Math.trunc(fn);
  let r  = x - fn * R_pio2_1;
  let w  = fn * R_pio2_1t;
  const j = ix >> 20;
  let y0 = r - w;
  let high = getHighWord(y0);
  let i = j - ((high >> 20) & 0x7ff);
  if (i > 16) {
    const t = r;
    w = fn * R_pio2_2;
    r = t - w;
    w = fn * R_pio2_2t - ((t - r) - w);
    y0 = r - w;
    high = getHighWord(y0);
    i = j - ((high >> 20) & 0x7ff);
    if (i > 49) {
      const t2 = r;
      w = fn * R_pio2_3;
      r = t2 - w;
      w = fn * R_pio2_3t - ((t2 - r) - w);
      y0 = r - w;
    }
  }
  const y1 = (r - y0) - w;
  return { n, y0, y1 };
}

// =============================================================================
// sin / cos  (s_sin.c / s_cos.c)
// =============================================================================
function fdSin(x) {
  const ix = getHighWord(x) & 0x7fffffff;
  if (ix <= 0x3fe921fb) {              // |x| ~< pi/4
    if (ix < 0x3e500000) return x;     // |x| < 2**-26 -> sin(x) ~ x
    return __kernel_sin(x, 0.0, 0);
  }
  if (ix >= 0x7ff00000) return NaN;    // sin(Inf or NaN) is NaN
  const { n, y0, y1 } = __ieee754_rem_pio2(x);
  switch (n & 3) {
    case 0: return  __kernel_sin(y0, y1, 1);
    case 1: return  __kernel_cos(y0, y1);
    case 2: return -__kernel_sin(y0, y1, 1);
    default: return -__kernel_cos(y0, y1);
  }
}

function fdCos(x) {
  const ix = getHighWord(x) & 0x7fffffff;
  if (ix <= 0x3fe921fb) {              // |x| ~< pi/4
    if (ix < 0x3e400000) return 1.0;   // |x| < 2**-27 -> cos(x) ~ 1
    return __kernel_cos(x, 0.0);
  }
  if (ix >= 0x7ff00000) return NaN;    // cos(Inf or NaN) is NaN
  const { n, y0, y1 } = __ieee754_rem_pio2(x);
  switch (n & 3) {
    case 0: return  __kernel_cos(y0, y1);
    case 1: return -__kernel_sin(y0, y1, 1);
    case 2: return -__kernel_cos(y0, y1);
    default: return  __kernel_sin(y0, y1, 1);
  }
}

// ---------------------------------------------------------------------------
// Extra word-access helpers needed by log/log1p/atanh (e_log.c, s_log1p.c,
// e_atanh.c). EXTRACT_WORDS(hx,lx,x) / SET_HIGH_WORD(x,v) equivalents.
// ---------------------------------------------------------------------------
function getLowWord(x) {         // signed int32 read, matches C's u_int32_t
  _dv.setFloat64(0, x, false);   // bit pattern -- JS's bitwise ops correctly
  return _dv.getInt32(4, false); // ToInt32-wrap it either way (see fdlibm.js
}                                 // log/atanh port notes on signedness).
function setHighWord(x, hi) {
  const lo = getLowWordU(x);
  return insertWords(hi, lo);
}

// =============================================================================
// log(x)  (e_log.c) -- direct port, including the subnormal/zero/inf/nan
// special cases even though this pipeline's actual inputs (gasdev's rsq,
// always in (0,1]) never reach them.
// =============================================================================
const L_ln2_hi = 6.93147180369123816490e-01;
const L_ln2_lo = 1.90821492927058770002e-10;
const L_two54  = 1.80143985094819840000e+16;
const L_Lg1 = 6.666666666666735130e-01;
const L_Lg2 = 3.999999999940941908e-01;
const L_Lg3 = 2.857142874366239149e-01;
const L_Lg4 = 2.222219843214978396e-01;
const L_Lg5 = 1.818357216161805012e-01;
const L_Lg6 = 1.531383769920937332e-01;
const L_Lg7 = 1.479819860511658591e-01;

function fdLog(x) {
  let hx = getHighWord(x);
  const lx = getLowWordU(x);
  let k = 0;

  if (hx < 0x00100000) {                       // x < 2**-1022
    if (((hx & 0x7fffffff) | lx) === 0) return -Infinity;  // log(+-0)
    if (hx < 0) return NaN;                    // log(-#)
    k -= 54; x *= L_two54;
    hx = getHighWord(x);
  }
  if (hx >= 0x7ff00000) return x + x;
  k += (hx >> 20) - 1023;
  hx &= 0x000fffff;
  const i = (hx + 0x95f64) & 0x100000;
  x = setHighWord(x, hx | (i ^ 0x3ff00000));
  k += (i >> 20);
  const f = x - 1.0;

  if ((0x000fffff & (2 + hx)) < 3) {           // -2**-20 <= f < 2**-20
    if (f === 0) {
      if (k === 0) return 0;
      const dk = k;
      return dk * L_ln2_hi + dk * L_ln2_lo;
    }
    const R = f * f * (0.5 - 0.33333333333333333 * f);
    if (k === 0) return f - R;
    const dk = k;
    return dk * L_ln2_hi - ((R - dk * L_ln2_lo) - f);
  }
  const s = f / (2.0 + f);
  const dk = k;
  const z = s * s;
  let i2 = hx - 0x6147a;
  const w = z * z;
  const j2 = 0x6b851 - hx;
  // t1 = w*(Lg2+w*(Lg4+w*Lg6)), t2 = z*(Lg1+w*(Lg3+w*(Lg5+w*Lg7))): this
  // project's actual C build (src/makeflags CFLAGS=-O, no -ffp-contract=off)
  // compiles fdlibm_log.c's Horner-chain inner additions with hardware FMA
  // (verified via -S disassembly and cross-checked against a debug-
  // instrumented build's own printed t1/t2/R intermediates for representative
  // gasdev-range inputs) -- same rationale as __kernel_sin/__kernel_cos
  // above. The outer w*(...)/z*(...) multiplies and the final R=t2+t1 are
  // NOT fused (confirmed: matches bit-exact as plain ops against the same
  // debug trace).
  const t1 = w * fma(w, fma(w, L_Lg6, L_Lg4), L_Lg2);
  const t2 = z * fma(w, fma(w, fma(w, L_Lg7, L_Lg5), L_Lg3), L_Lg1);
  i2 |= j2;
  const R = t2 + t1;
  if (i2 > 0) {
    const hfsq = 0.5 * f * f;
    // hfsq-s*(hfsq+R) fuses to one fmsub (a*b-c form => fma(-s,hfsq+R,hfsq)).
    if (k === 0) return f - fma(-s, hfsq + R, hfsq);
    return dk * L_ln2_hi - ((hfsq - fma(s, hfsq + R, dk * L_ln2_lo)) - f);
  } else {
    // f-s*(f-R) fuses to a single fmsub covering the WHOLE expression.
    if (k === 0) return fma(-s, f - R, f);
    return dk * L_ln2_hi - (fma(s, f - R, -(dk * L_ln2_lo)) - f);
  }
}

// =============================================================================
// log1p(x)  (s_log1p.c) -- needed by atanh.
// =============================================================================
const P_ln2_hi = 6.93147180369123816490e-01;
const P_ln2_lo = 1.90821492927058770002e-10;
const P_two54  = 1.80143985094819840000e+16;
const P_Lp1 = 6.666666666666735130e-01;
const P_Lp2 = 3.999999999940941908e-01;
const P_Lp3 = 2.857142874366239149e-01;
const P_Lp4 = 2.222219843214978396e-01;
const P_Lp5 = 1.818357216161805012e-01;
const P_Lp6 = 1.531383769920937332e-01;
const P_Lp7 = 1.479819860511658591e-01;

function fdLog1p(x) {
  const hxIn = getHighWord(x);
  const ax = hxIn & 0x7fffffff;

  let k = 1, f, hu = 0, c = 0;
  let hx = hxIn;

  if (hx < 0x3FDA827A) {                       // 1+x < sqrt(2)
    if (ax >= 0x3ff00000) {                    // x <= -1.0
      if (x === -1.0) return -Infinity;
      return NaN;
    }
    if (ax < 0x3e200000) {                     // |x| < 2**-29
      if ((P_two54 + x) > 0 && ax < 0x3c900000) return x;
      return x - x * x * 0.5;
    }
    if (hx > 0 || hx <= (0xbfd2bec4 | 0)) { k = 0; f = x; hu = 1; }
  }
  if (hx >= 0x7ff00000) return x + x;

  let u;
  if (k !== 0) {
    if (hx < 0x43400000) {
      u = 1.0 + x;
      hu = getHighWord(u);
      k = (hu >> 20) - 1023;
      c = (k > 0) ? 1.0 - (u - x) : x - (u - 1.0);
      c /= u;
    } else {
      u = x;
      hu = getHighWord(u);
      k = (hu >> 20) - 1023;
      c = 0;
    }
    hu &= 0x000fffff;
    if (hu < 0x6a09e) {
      u = setHighWord(u, hu | 0x3ff00000);
    } else {
      k += 1;
      u = setHighWord(u, hu | 0x3fe00000);
      hu = (0x00100000 - hu) >> 2;
    }
    f = u - 1.0;
  }
  const hfsq = 0.5 * f * f;
  if (hu === 0) {                              // |f| < 2**-20
    if (f === 0) {
      if (k === 0) return 0;
      c += k * P_ln2_lo;
      return k * P_ln2_hi + c;
    }
    const R = hfsq * (1.0 - 0.66666666666666666 * f);
    if (k === 0) return f - R;
    return k * P_ln2_hi - ((R - (k * P_ln2_lo + c)) - f);
  }
  const s = f / (2.0 + f);
  const z = s * s;
  // Same real, verified FMA usage as fdLog's t1/t2 Horner chains and final
  // hfsq combination (this project's actual -O build fuses these) -- see
  // fdLog's comment above for the full rationale. Outer z* multiply is
  // plain (matches fdLog's t1/t2 pattern: inner Horner steps fuse, the
  // final outer multiply by z does not).
  const h1 = fma(z, P_Lp7, P_Lp6);
  const h2 = fma(z, h1, P_Lp5);
  const h3 = fma(z, h2, P_Lp4);
  const h4 = fma(z, h3, P_Lp3);
  const h5 = fma(z, h4, P_Lp2);
  const h6 = fma(z, h5, P_Lp1);
  const R = z * h6;
  if (k === 0) return f - fma(-s, hfsq + R, hfsq);
  return k * P_ln2_hi - ((hfsq - fma(s, hfsq + R, k * P_ln2_lo + c)) - f);
}

// =============================================================================
// atanh(x)  (e_atanh.c)
// =============================================================================
function fdAtanh(x) {
  const hx = getHighWord(x);
  const lx = getLowWordU(x);
  const ix = hx & 0x7fffffff;

  // (lx|(-lx))>>>31 in C is a logical (unsigned) shift -- use >>> in JS to match.
  const lxNonZero = ((lx | (-lx)) >>> 31);
  if ((ix | lxNonZero) > 0x3ff00000) return NaN;    // |x| > 1
  if (ix === 0x3ff00000) return x > 0 ? Infinity : -Infinity; // |x| == 1
  if (ix < 0x3e300000 && (1e300 + x) > 0) return x; // |x| < 2**-28

  const xn = setHighWord(x, ix);                    // x = |x|, keep magnitude
  let t;
  if (ix < 0x3fe00000) {                             // |x| < 0.5
    const t2 = xn + xn;
    t = 0.5 * fdLog1p(t2 + t2 * xn / (1.0 - xn));
  } else {
    t = 0.5 * fdLog1p((xn + xn) / (1.0 - xn));
  }
  return hx >= 0 ? t : -t;
}

// =============================================================================
// exp(x)  (e_exp.c) -- direct port. Newly added (not part of the original
// sin/cos/log/atanh set) after tracing a beam-kernel discrepancy
// (CalculateBeamKernel.f/.js) to Fortran's native EXP intrinsic and JS's
// native Math.exp disagreeing for the same input -- neither side had ever
// been routed through fdlibm for this function. Verified via
// fdlibm_log(fdlibm_exp(x)) round-trip (machine-epsilon-level error) and
// direct comparison against the same fdlibm_exp.c compiled and linked into
// the real Fortran binary. ln2HI/ln2LO reuse the SAME split-ln2 constants as
// fdLog above (fdlibm shares this double-double ln(2) representation across
// log/log1p/exp).
// =============================================================================
const E_one = 1.0;
const E_halF = [0.5, -0.5];
const E_huge = 1.0e300;
const E_twom1000 = 9.33263618503218878990e-302;
const E_o_threshold = 7.09782712893383973096e+02;
const E_u_threshold = -7.45133219101941108420e+02;
const E_ln2HI = [6.93147180369123816490e-01, -6.93147180369123816490e-01];
const E_ln2LO = [1.90821492927058770002e-10, -1.90821492927058770002e-10];
const E_invln2 = 1.44269504088896338700e+00;
const E_P1 = 1.66666666666666019037e-01;
const E_P2 = -2.77777777770155933842e-03;
const E_P3 = 6.61375632143793436117e-05;
const E_P4 = -1.65339022054652515390e-06;
const E_P5 = 4.13813679705723846039e-08;

function fdExp(x) {
  let hx = getHighWord(x);
  const xsb = (hx >>> 31) & 1;
  hx &= 0x7fffffff;

  // filter out non-finite argument
  if (hx >= 0x40862E42) {
    if (hx >= 0x7ff00000) {
      const lx = getLowWordU(x);
      if (((hx & 0xfffff) | lx) !== 0) return NaN;
      return xsb === 0 ? x : 0.0;
    }
    if (x > E_o_threshold) return E_huge * E_huge;
    if (x < E_u_threshold) return E_twom1000 * E_twom1000;
  }

  // argument reduction. FMA NOTE: this project's actual C build (-O, no
  // -ffp-contract=off) compiles fdlibm_exp.c's `invln2*x+halF[xsb]` and
  // `x-t*ln2HI[0]` with hardware FMA (verified via -S disassembly) -- same
  // real, verified FMA usage as __kernel_sin/__kernel_cos above. The
  // "x-t*ln2HI[0]" case fuses as fma(t,-ln2HI[0],x) (compiler negates the
  // constant at compile time, same trick seen in __kernel_sin's iy!=0
  // branch), not a literal x-(t*ln2HI[0]) two-step.
  let k = 0, hi = 0, lo = 0;
  if (hx > 0x3fd62e42) {
    if (hx < 0x3FF0A2B2) {
      hi = x - E_ln2HI[xsb]; lo = E_ln2LO[xsb]; k = 1 - xsb - xsb;
    } else {
      k = Math.trunc(fma(x, E_invln2, E_halF[xsb])) | 0;
      const t = k;
      hi = fma(t, -E_ln2HI[0], x);
      lo = t * E_ln2LO[0];
    }
    x = hi - lo;
  } else if (hx < 0x3e300000) {
    if (E_huge + x > E_one) return E_one + x;
  }

  // x is now in primary range. FMA NOTE: the P1..P5 Horner chain and the
  // final x-t*(...) combine all fuse (5 total fmadd/fmsub instructions,
  // verified via disassembly); the k==0/k!=0 y-computation branches below
  // do NOT fuse (division breaks the contraction opportunity there).
  const t = x * x;
  const h1 = fma(t, E_P5, E_P4);
  const h2 = fma(t, h1, E_P3);
  const h3 = fma(t, h2, E_P2);
  const h4 = fma(t, h3, E_P1);
  const c = fma(-t, h4, x);
  let y;
  if (k === 0) return E_one - ((x * c) / (c - 2.0) - x);
  y = E_one - ((lo - (x * c) / (2.0 - c)) - hi);
  let hy = getHighWord(y);
  if (k >= -1021) {
    hy = hy + (k << 20);
    return setHighWord(y, hy);
  }
  hy = hy + ((k + 1000) << 20);
  return setHighWord(y, hy) * E_twom1000;
}

// =============================================================================
// atan(x)  (s_atan.c) -- direct port. Newly added (not part of the original
// sin/cos/log/atanh/exp set) after tracing a real Fortran/JS shape/
// inclination estimate divergence to FullCircTrig.f's/.js's native atan()/
// Math.atan() -- unlike sin/cos, V8's Math.atan is NOT already a bit-exact
// fdlibm derivative (confirmed by direct comparison: mismatches occur for
// ordinary inputs, e.g. atan(-0.3), atan(0.618)), and neither had ever been
// routed through fdlibm before. Source: FreeBSD lib/msun/src/s_atan.c.
//
// FMA: this project's fdlibm_atan.c, compiled with the actual project
// CFLAGS (-O, no -ffp-contract=off), gets hardware-FMA-fused for every
// coefficient-chain multiply-add in the two Horner polynomials AND for the
// two reduced-argument numerator/denominator terms that contain a literal
// multiply (id==0's 2*x-1 and id==2's 1+1.5*x) -- verified via -S
// disassembly, same methodology as __kernel_sin/__kernel_cos above. Ported
// with fma() at each such site; id==1/id==3's reductions have no multiply
// in their numerator/denominator and are plain ops, matching the
// disassembly showing no fmadd/fmsub there.
// =============================================================================
const AT_atanhi = [
  4.63647609000806093515e-01, // atan(0.5)hi
  7.85398163397448278999e-01, // atan(1.0)hi
  9.82793723247329054082e-01, // atan(1.5)hi
  1.57079632679489655800e+00, // atan(inf)hi
];
const AT_atanlo = [
  2.26987774529616870924e-17, // atan(0.5)lo
  3.06161699786838301793e-17, // atan(1.0)lo
  1.39033110312309984516e-17, // atan(1.5)lo
  6.12323399573676603587e-17, // atan(inf)lo
];
const AT0  =  3.33333333333329318027e-01;
const AT1  = -1.99999999998764832476e-01;
const AT2  =  1.42857142725034663711e-01;
const AT3  = -1.11111104054623557880e-01;
const AT4  =  9.09090142995028394523e-02;
const AT5  = -7.69187620504482999495e-02;
const AT6  =  6.66107313738753120669e-02;
const AT7  = -5.83357013379057348645e-02;
const AT8  =  4.97687799461593236017e-02;
const AT9  = -3.65315727442169155270e-02;
const AT10 =  1.62858201153657823623e-02;
const AT_huge = 1.0e300;

function fdAtan(x) {
  const hx = getHighWord(x);
  const ix = hx & 0x7fffffff;
  const lx = getLowWordU(x);

  if (ix >= 0x44100000) {                    // |x| >= 2^66
    if (ix > 0x7ff00000 || (ix === 0x7ff00000 && lx !== 0)) return x + x; // NaN
    return hx > 0 ? (AT_atanhi[3] + AT_atanlo[3]) : -(AT_atanhi[3] + AT_atanlo[3]);
  }

  let id;
  if (ix < 0x3fdc0000) {                     // |x| < 0.4375
    if (ix < 0x3e400000) {                   // |x| < 2^-27
      if (AT_huge + x > 1.0) return x;       // raise inexact
    }
    id = -1;
  } else {
    x = Math.abs(x);
    if (ix < 0x3ff30000) {                   // |x| < 1.1875
      if (ix < 0x3fe60000) {                 // 7/16 <= |x| < 11/16
        id = 0; x = fma(x, 2.0, -1.0) / (2.0 + x);      // (2x-1)/(2+x)
      } else {                               // 11/16 <= |x| < 19/16
        id = 1; x = (x - 1.0) / (x + 1.0);
      }
    } else {
      if (ix < 0x40038000) {                 // |x| < 2.4375
        id = 2; x = (x - 1.5) / fma(x, 1.5, 1.0);       // (x-1.5)/(1+1.5x)
      } else {                               // 2.4375 <= |x| < 2^66
        id = 3; x = -1.0 / x;
      }
    }
  }
  // end of argument reduction
  const z = x * x;
  const w = z * z;
  // s1 = z*(aT0+w*(aT2+w*(aT4+w*(aT6+w*(aT8+w*aT10))))), fused per-level:
  let t = fma(w, AT10, AT8);
  t = fma(w, t, AT6);
  t = fma(w, t, AT4);
  t = fma(w, t, AT2);
  t = fma(w, t, AT0);
  const s1 = z * t;
  // s2 = w*(aT1+w*(aT3+w*(aT5+w*(aT7+w*aT9)))), fused per-level:
  let u = fma(w, AT9, AT7);
  u = fma(w, u, AT5);
  u = fma(w, u, AT3);
  u = fma(w, u, AT1);
  const s2 = w * u;

  if (id < 0) {
    return fma(x, -(s1 + s2), x);            // x - x*(s1+s2)
  }
  const sum = s1 + s2;
  const t2 = fma(x, sum, -AT_atanlo[id]);    // x*(s1+s2) - atanlo[id]
  const t3 = t2 - x;
  const z2 = AT_atanhi[id] - t3;
  return hx < 0 ? -z2 : z2;
}

module.exports = {
  fdSin, fdCos, fdLog, fdLog1p, fdAtanh, fdExp, fdAtan,
  // exposed for testing
  __kernel_sin, __kernel_cos, __ieee754_rem_pio2,
  getHighWord, getLowWordU, insertWords, rnint,
};

// ---------------------------------------------------------------------------
// Self-test (node fdlibm.js): sanity check against native Math.sin/cos --
// NOT a substitute for the real cross-language verification against a
// compiled fdlibm C harness, just a quick smoke test.
// ---------------------------------------------------------------------------
if (require.main === module) {
  const tests = [0, 1e-10, 0.5, 1.0, Math.PI / 4, Math.PI / 2, 1.0, 2.0, 3.0,
                  Math.PI, 4.0, 5.0, 2 * Math.PI, -1.5, 6.28];
  console.log('=== sin/cos vs native Math ===');
  for (const x of tests) {
    const s = fdSin(x), c = fdCos(x);
    console.log(x, s, Math.sin(x), s - Math.sin(x), c, Math.cos(x), c - Math.cos(x));
  }
  console.log('=== log vs native Math.log (positive test values) ===');
  const logTests = [1e-10, 0.001, 0.01, 0.1, 0.3333333333333333, 0.5, 0.9,
                     0.9999999, 1.0, 1.5, 2.0, 10.0, 1e10];
  for (const x of logTests) {
    const v = fdLog(x);
    console.log(x, v, Math.log(x), v - Math.log(x));
  }
  console.log('=== atanh vs native Math.atanh ===');
  const atanhTests = [1e-10, 0.001, 0.1, 0.3, 0.49, 0.5, 0.51, 0.7, 0.9,
                       0.9999999, -0.3, -0.7];
  for (const x of atanhTests) {
    const v = fdAtanh(x);
    console.log(x, v, Math.atanh(x), v - Math.atanh(x));
  }
}

});
