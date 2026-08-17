/* Fortran-callable wrappers for the fdlibm port (fdlibm_sin.c, fdlibm_cos.c,
 * fdlibm_log.c, fdlibm_atanh.c). gfortran's default (non-bind(c)) calling
 * convention: arguments passed by reference, external symbol names
 * lowercased with a trailing underscore.
 *
 * WHY: Fortran's SIN/COS/LOG/ATANH intrinsics call the platform's system
 * libm (Apple's on macOS), which is a different implementation from
 * fdlibm/FreeBSD's msun -- and, confirmed by direct comparison against a
 * compiled fdlibm reference, occasionally diverges from it by 1 ULP
 * (e.g. sin(4.0), atanh(0.5)). The JS port's V8 runtime appears to already
 * implement fdlibm (or a bit-identical derivative) for these functions.
 * Routing Fortran's calls through this same fdlibm source closes that gap
 * so both languages compute bit-identical results for the same input.
 *
 * The pipeline's Fortran `real` is single precision (4 bytes) throughout --
 * matching the JS port's f32(Math.sin(x)) convention (compute in double,
 * round to float32 at the end), these wrappers widen the real(4) input to
 * double, call the double-precision fdlibm function, and let the C->Fortran
 * return-value narrowing do the same float32 rounding JS's f32() does.
 */
double fdlibm_sin(double x);
double fdlibm_cos(double x);
double fdlibm_log(double x);
double fdlibm_atanh(double x);
double fdlibm_exp(double x);
double fdlibm_atan(double x);

float fd_sin_(float *x)   { return (float)fdlibm_sin((double)(*x)); }
float fd_cos_(float *x)   { return (float)fdlibm_cos((double)(*x)); }
float fd_log_(float *x)   { return (float)fdlibm_log((double)(*x)); }
float fd_atanh_(float *x) { return (float)fdlibm_atanh((double)(*x)); }
float fd_exp_(float *x)   { return (float)fdlibm_exp((double)(*x)); }
float fd_atan_(float *x)  { return (float)fdlibm_atan((double)(*x)); }
