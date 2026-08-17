/*
 * fdlibm_k_cos_wasm.c -- wasm-only variant of ../fdlibm_k_cos.c.
 * See fdlibm_k_sin_wasm.c's header for why this exists and why it must
 * not be merged back into the canonical native file. Fusion structure
 * mirrors ../fdlibm.js's __kernel_cos exactly.
 */
#include "math.h"
#include "fdlibm_private.h"
#include "fdlibm_soft_fma.h"

static const double
one =  1.00000000000000000000e+00, /* 0x3FF00000, 0x00000000 */
C1  =  4.16666666666666019037e-02, /* 0x3FA55555, 0x5555554C */
C2  = -1.38888888888741095749e-03, /* 0xBF56C16C, 0x16C15177 */
C3  =  2.48015872894767294178e-05, /* 0x3EFA01A0, 0x19CB1590 */
C4  = -2.75573143513906633035e-07, /* 0xBE927E4F, 0x809C52AD */
C5  =  2.08757232129817482790e-09, /* 0x3E21EE9E, 0xBDB4B1C4 */
C6  = -1.13596475577881948265e-11; /* 0xBDA8FAE9, 0xBE8838D4 */

double
__kernel_cos(double x, double y)
{
	double hz,z,r,w;
	double c1,c2,c3,c4,w2sq,hz2,w2,oneMinusW2,t,negxy,zrMinusXy;

	z  = x*x;
	w  = z*z;
	c1 = soft_fma(z, C3, C2);        /* C2+z*C3 */
	c2 = soft_fma(z, c1, C1);         /* C1+z*c1 */
	c3 = soft_fma(z, C6, C5);         /* C5+z*C6 */
	c4 = soft_fma(z, c3, C4);         /* C4+z*c3 */
	w2sq = w*w;                       /* plain multiply, not fused */
	r  = soft_fma(z, c2, w2sq*c4);    /* z*c2 + w2sq*c4 */
	hz2 = 0.5*z;
	w2 = one-hz2;
	oneMinusW2 = one-w2;
	t = oneMinusW2-hz2;
	negxy = -(x*y);
	zrMinusXy = soft_fma(z, r, negxy); /* z*r-x*y */
	return w2 + (t + zrMinusXy);
}
