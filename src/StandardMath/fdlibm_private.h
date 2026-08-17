/* Minimal replacement for FreeBSD's math_private.h, providing only what
 * s_sin.c / s_cos.c / k_sin.c / k_cos.c / e_rem_pio2.c actually use.
 * The original is riddled with long-double portability macros (LDBL_MANT_DIG
 * concatenation tricks) that don't compile cleanly against macOS's SDK and
 * are irrelevant here -- we only need double precision.
 */
#ifndef MATH_PRIVATE_H
#define MATH_PRIVATE_H

#include <stdint.h>
typedef uint32_t u_int32_t;
typedef double double_t;

typedef union {
  double value;
  struct { uint32_t lsw; uint32_t msw; } parts; /* little-endian host */
} ieee_double_shape_type;

#define GET_HIGH_WORD(i,d) do { \
  ieee_double_shape_type gh_u; gh_u.value = (d); (i) = gh_u.parts.msw; \
} while (0)

#define GET_LOW_WORD(i,d) do { \
  ieee_double_shape_type gl_u; gl_u.value = (d); (i) = gl_u.parts.lsw; \
} while (0)

#define INSERT_WORDS(d,ix0,ix1) do { \
  ieee_double_shape_type iw_u; iw_u.parts.msw = (ix0); iw_u.parts.lsw = (ix1); \
  (d) = iw_u.value; \
} while (0)

#define EXTRACT_WORDS(ix0,ix1,d) do { \
  ieee_double_shape_type ew_u; ew_u.value = (d); \
  (ix0) = ew_u.parts.msw; (ix1) = ew_u.parts.lsw; \
} while (0)

#define SET_HIGH_WORD(d,v) do { \
  ieee_double_shape_type sh_u; sh_u.value = (d); sh_u.parts.msw = (v); \
  (d) = sh_u.value; \
} while (0)

#define STRICT_ASSIGN(type,lval,rval) ((lval) = (rval))

static inline double rnint(double_t x) {
  return ((double)(x + 0x1.8p52) - 0x1.8p52);
}
#define irint(x) ((int)(x))

#define __always_inline inline __attribute__((always_inline))
#define __weak_reference(sym,alias) /* not needed: double precision only */

double __kernel_sin(double x, double y, int iy);
double __kernel_cos(double x, double y);
/* Large-argument Payne-Hanek path -- intentionally not implemented, see
 * fdlibm.js SCOPE comment. Nothing in this pipeline's angle range reaches
 * it; abort loudly if that assumption is ever wrong. */
int __kernel_rem_pio2(double *x, double *y, int e0, int nx, int prec);
double fdlibm_log1p(double x);

#endif
