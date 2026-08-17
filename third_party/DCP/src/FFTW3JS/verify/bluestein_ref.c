/* bluestein_ref.c -- independent C reference reimplementation of
 * dft/bluestein.c's apply()/mktwiddle()/bluestein_sequence(), using the
 * PUBLIC fftw3.h API (fftw_plan_dft_1d, FFTW_FORWARD, FFTW_ESTIMATE) for the
 * two internal size-nb DFT calls, and a hand-copied real_cexp (same
 * algorithm as kernel/trig.c, same fdlibm sin/cos) for the chirp sequence.
 * Purpose: isolate whether a JS Bluestein transcription bug is in chirp
 * generation or in the convolution combine steps, by diffing this
 * independent C path against both real FFTW's direct N-point DFT and the
 * JS BluesteinSolver.js output for the same N.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

extern double fdlibm_sin(double x);
extern double fdlibm_cos(double x);

static const double K2PI = 6.2831853071795864769252867665590057683943388;

static void real_cexp(long m, long n, double *out) {
    double theta, c, s, t;
    unsigned octant = 0;
    long quarter_n = n;
    n += n; n += n;
    m += m; m += m;
    if (m < 0) m += n;
    if (m > n - m) { m = n - m; octant |= 4; }
    if (m - quarter_n > 0) { m = m - quarter_n; octant |= 2; }
    if (m > quarter_n - m) { m = quarter_n - m; octant |= 1; }
    theta = (K2PI * m) / n;
    c = fdlibm_cos(theta); s = fdlibm_sin(theta);
    if (octant & 1) { t = c; c = s; s = t; }
    if (octant & 2) { t = c; c = -s; s = t; }
    if (octant & 4) { s = -s; }
    out[0] = c; out[1] = s;
}

static int factors_into_small_primes(long n) {
    long p;
    int ps[3] = {2,3,5};
    for (int i = 0; i < 3; i++) { p = ps[i]; while (n % p == 0) n /= p; }
    return n == 1;
}
static long choose_transform_size(long minsz) {
    while (!factors_into_small_primes(minsz)) minsz++;
    return minsz;
}

/* run fftw_plan_dft_1d(nb, FFTW_FORWARD, FFTW_ESTIMATE) in place on buf
 * (interleaved re,im), matching cldf->apply(cldf, buf,buf+1,buf,buf+1). */
static void dft_forward_inplace(int nb, double *buf) {
    fftw_complex *in = fftw_malloc(sizeof(fftw_complex) * nb);
    fftw_complex *out = fftw_malloc(sizeof(fftw_complex) * nb);
    for (int i = 0; i < nb; i++) { in[i][0] = buf[2*i]; in[i][1] = buf[2*i+1]; }
    fftw_plan p = fftw_plan_dft_1d(nb, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    for (int i = 0; i < nb; i++) { buf[2*i] = out[i][0]; buf[2*i+1] = out[i][1]; }
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s N\n", argv[0]); return 1; }
    int n = atoi(argv[1]);
    int n2 = 2 * n;
    int nb = (int)choose_transform_size(2*n - 1);

    double *w = malloc(sizeof(double) * 2 * n);
    long ksq = 0;
    for (int k = 0; k < n; k++) {
        real_cexp(ksq, n2, w + 2*k);
        ksq += 2*k + 1;
        while (ksq > n2) ksq -= n2;
    }

    double *W = calloc(2 * nb, sizeof(double));
    W[0] = w[0] / nb; W[1] = w[1] / nb;
    for (int i = 1; i < n; i++) {
        W[2*i] = W[2*(nb-i)] = w[2*i] / nb;
        W[2*i+1] = W[2*(nb-i)+1] = w[2*i+1] / nb;
    }
    dft_forward_inplace(nb, W);

    double *ri = malloc(sizeof(double)*n), *ii = malloc(sizeof(double)*n);
    for (int i = 0; i < n; i++) {
        ri[i] = sin(i * 1.7 + 0.3);
        ii[i] = cos(i * 0.9 + 1.1);
    }

    double *b = calloc(2 * nb, sizeof(double));
    for (int i = 0; i < n; i++) {
        double xr = ri[i], xi = ii[i];
        double wr = w[2*i], wi = w[2*i+1];
        b[2*i] = xr*wr + xi*wi;
        b[2*i+1] = xi*wr - xr*wi;
    }
    dft_forward_inplace(nb, b);
    for (int i = 0; i < nb; i++) {
        double xr = b[2*i], xi = b[2*i+1];
        double wr = W[2*i], wi = W[2*i+1];
        b[2*i] = xi*wr + xr*wi;
        b[2*i+1] = xr*wr - xi*wi;
    }
    dft_forward_inplace(nb, b);

    double *ro = malloc(sizeof(double)*n), *io = malloc(sizeof(double)*n);
    for (int i = 0; i < n; i++) {
        double xi = b[2*i], xr = b[2*i+1];
        double wr = w[2*i], wi = w[2*i+1];
        ro[i] = xr*wr + xi*wi;
        io[i] = xi*wr - xr*wi;
    }

    for (int i = 0; i < n; i++) printf("%.17g %.17g\n", ri[i], ii[i]);
    for (int i = 0; i < n; i++) printf("%.17g %.17g\n", ro[i], io[i]);

    return 0;
}
