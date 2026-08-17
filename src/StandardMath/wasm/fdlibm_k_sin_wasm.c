/*
 * fdlibm_k_sin_wasm.c -- wasm-only variant of ../fdlibm_k_sin.c.
 *
 * Identical polynomial/algorithm to the canonical file (same Sun
 * Microsystems fdlibm __kernel_sin, same coefficients) EXCEPT every
 * multiply-add the native build gets via hardware FMA contraction is
 * here an explicit soft_fma() call instead -- see fdlibm_soft_fma.h's
 * header for why. Fusion structure mirrors ../fdlibm.js's __kernel_sin
 * exactly (that JS port's fma() call sites are the verified reference
 * for which operations are fused vs. left plain).
 *
 * Do not edit ../fdlibm_k_sin.c to match this file, or vice versa --
 * they intentionally diverge in HOW the fusion is achieved (compiler
 * contraction vs. explicit software emulation) because they target
 * different backends with different native FMA availability.
 */
#include "math.h"
#include "fdlibm_private.h"
#include "fdlibm_soft_fma.h"

static const double
half =  5.00000000000000000000e-01, /* 0x3FE00000, 0x00000000 */
S1  = -1.66666666666666324348e-01, /* 0xBFC55555, 0x55555549 */
S2  =  8.33333333332248946124e-03, /* 0x3F811111, 0x1110F8A6 */
S3  = -1.98412698298579493134e-04, /* 0xBF2A01A0, 0x19C161D5 */
S4  =  2.75573137070700676789e-06, /* 0x3EC71DE3, 0x57B1FE7D */
S5  = -2.50507602534068634195e-08, /* 0xBE5AE5E6, 0x8A2B9CEB */
S6  =  1.58969099521155010221e-10; /* 0x3DE5D93A, 0x5ACFD57C */

double
__kernel_sin(double x, double y, int iy)
{
	double z,w,r,v;
	double t1,t2,t3,zw;

	z = x*x;
	w = z*z;
	t1 = soft_fma(z, S4, S3);       /* S3+z*S4 */
	t2 = soft_fma(z, t1, S2);        /* S2+z*t1 */
	t3 = soft_fma(z, S6, S5);        /* S5+z*S6 */
	zw = z*w;                        /* plain multiply, not fused */
	r  = soft_fma(zw, t3, t2);       /* t2+zw*t3 */
	v  = z*x;
	if (iy==0) return soft_fma(v, soft_fma(z, r, S1), x); /* x+v*(S1+z*r) */
	else {
		double p1 = -(v*r);
		double tb = soft_fma(y, half, p1);   /* half*y-v*r */
		double tc = soft_fma(z, tb, -y);      /* z*(half*y-v*r)-y */
		double td = soft_fma(v, -S1, tc);     /* tc-v*S1 */
		return x - td;
	}
}
