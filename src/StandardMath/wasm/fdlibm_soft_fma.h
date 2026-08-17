/*
 * fdlibm_soft_fma.h -- software-emulated IEEE754 fused multiply-add
 * (a*b+c, single rounding), for the wasm build ONLY.
 *
 * WHY THIS EXISTS: WebAssembly's scalar f64 arithmetic has no fused
 * multiply-add instruction (the relaxed-simd proposal's fma variants only
 * cover f32x4/f64x2 SIMD lanes, not plain scalars). The canonical native
 * build of fdlibm_k_sin.c/fdlibm_k_cos.c relies on implicit compiler
 * contraction of `a*b+c` into a real hardware FMA (mandatory baseline on
 * arm64, this project's dev machine -- confirmed via disassembly, see
 * ../fma.js's header). That contraction cannot survive the trip to wasm:
 * LLVM's wasm backend legalizes any contracted fmuladd back to separate
 * fmul+fadd (double rounding), which is a REAL, measured divergence from
 * the native/Fortran-matching behavior (confirmed empirically: building
 * fdlibm_k_sin.c/fdlibm_k_cos.c as-is for wasm and comparing bit-exact
 * against the verified JS port produced 1271/200000 sin mismatches,
 * 1107/200000 cos mismatches -- same order of magnitude as the native
 * platform-libm divergence problem this whole fdlibm port exists to fix).
 *
 * THE FIX: this header provides soft_fma(a,b,c), a straight C port of
 * ../fdlibm.js's fma.js (Dekker/Veltkamp TwoProduct/TwoSum error-free
 * transformations -- see that file's header for the algorithm reference
 * and correctness discussion). fdlibm_k_sin_wasm.c/fdlibm_k_cos_wasm.c
 * call this explicitly wherever the JS port's __kernel_sin/__kernel_cos
 * calls its own fma(), instead of relying on `a*b+c` contraction -- so
 * the wasm build now performs the exact same explicit computation the
 * verified JS port does, immune to what the target's FMA hardware
 * support (or lack thereof) would otherwise do.
 *
 * NOT used by the canonical fdlibm_k_sin.c/fdlibm_k_cos.c (those stay
 * exactly as verified against the real native/Fortran-linked build).
 */
#ifndef FDLIBM_SOFT_FMA_H
#define FDLIBM_SOFT_FMA_H

static const double FDLIBM_FMA_SPLITTER = 134217729.0; /* 2^27 + 1, Veltkamp's algorithm */

static inline void fdlibm_fma_split(double a, double *hi, double *lo) {
  double c = FDLIBM_FMA_SPLITTER * a;
  *hi = c - (c - a);
  *lo = a - *hi;
}

/* TwoProduct(a,b): p+e == a*b exactly, p = round(a*b). */
static inline void fdlibm_fma_two_product(double a, double b, double *p, double *e) {
  double aHi, aLo, bHi, bLo;
  *p = a * b;
  fdlibm_fma_split(a, &aHi, &aLo);
  fdlibm_fma_split(b, &bHi, &bLo);
  *e = ((aHi * bHi - *p) + aHi * bLo + aLo * bHi) + aLo * bLo;
}

/* TwoSum(a,b): s+t == a+b exactly, s = round(a+b). */
static inline void fdlibm_fma_two_sum(double a, double b, double *s, double *t) {
  double v;
  *s = a + b;
  v = *s - a;
  *t = (a - (*s - v)) + (b - v);
}

/* soft_fma(a,b,c) -- correctly-rounded a*b+c (single rounding). */
static inline double soft_fma(double a, double b, double c) {
  double p, e, s, t;
  fdlibm_fma_two_product(a, b, &p, &e);
  fdlibm_fma_two_sum(p, c, &s, &t);
  return s + (t + e);
}

#endif /* FDLIBM_SOFT_FMA_H */
