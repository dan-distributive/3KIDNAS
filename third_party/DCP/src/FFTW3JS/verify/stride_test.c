/* stride_test.c -- does fftw_plan_many_dft with istride=ostride=2 (matching
 * dft/bluestein.c's internal mktensor_1d(nb,2,2) sub-problem shape) produce
 * bit-different output than a plain stride-1 fftw_plan_dft_1d call, for the
 * same logical input? Tests the hypothesis that bluestein.c's cldf sub-plan
 * picks a different solver (hence different rounding) than a top-level
 * stride-1 call of the same size. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

int main(int argc, char **argv) {
    int n = atoi(argv[1]);

    /* stride-1 */
    fftw_complex *in1 = fftw_malloc(sizeof(fftw_complex) * n);
    fftw_complex *out1 = fftw_malloc(sizeof(fftw_complex) * n);
    for (int i = 0; i < n; i++) { in1[i][0] = sin(i*1.7+0.3); in1[i][1] = cos(i*0.9+1.1); }
    fftw_plan p1 = fftw_plan_dft_1d(n, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p1);

    /* stride-2, in-place, interleaved -- same shape as bluestein.c's cldf:
     * mktensor_1d(nb,2,2) for BOTH sz (with istride=ostride=2). */
    double *buf = fftw_malloc(sizeof(double) * 2 * n);
    for (int i = 0; i < n; i++) { buf[2*i] = sin(i*1.7+0.3); buf[2*i+1] = cos(i*0.9+1.1); }
    fftw_iodim dims[1];
    dims[0].n = n; dims[0].is = 2; dims[0].os = 2;
    fftw_plan p2 = fftw_plan_guru_dft(1, dims, 0, NULL,
                                       (fftw_complex*)buf, (fftw_complex*)buf,
                                       FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p2);

    double maxDiff = 0;
    for (int i = 0; i < n; i++) {
        double dr = fabs(out1[i][0] - buf[2*i]);
        double di = fabs(out1[i][1] - buf[2*i+1]);
        if (dr > maxDiff) maxDiff = dr;
        if (di > maxDiff) maxDiff = di;
    }
    printf("n=%d maxDiff(stride1 vs stride2)=%.17g\n", n, maxDiff);

    fftw_fprint_plan(p1, stdout);
    printf("\n");
    fftw_fprint_plan(p2, stdout);
    printf("\n");
    return 0;
}
