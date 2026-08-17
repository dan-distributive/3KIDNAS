/* inplace_test.c -- disentangle stride vs in-place-ness as the trigger for
 * q1 (twiddle-squared) selection. Tests all 4 combinations: {stride1,
 * stride2} x {in-place, out-of-place}. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

static int has_q1(fftw_plan p) {
    char buf[8192];
    FILE *f = fmemopen(buf, sizeof(buf), "w");
    fftw_fprint_plan(p, f);
    fclose(f);
    return strstr(buf, "\"q1_") != NULL;
}

int main(int argc, char **argv) {
    int n = atoi(argv[1]);

    /* stride1, out-of-place */
    fftw_complex *in1 = fftw_malloc(sizeof(fftw_complex) * n);
    fftw_complex *out1 = fftw_malloc(sizeof(fftw_complex) * n);
    fftw_plan p1 = fftw_plan_dft_1d(n, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);

    /* stride1, in-place */
    fftw_complex *io2 = fftw_malloc(sizeof(fftw_complex) * n);
    fftw_plan p2 = fftw_plan_dft_1d(n, io2, io2, FFTW_FORWARD, FFTW_ESTIMATE);

    /* stride2, out-of-place */
    double *in3 = fftw_malloc(sizeof(double) * 2 * n);
    double *out3 = fftw_malloc(sizeof(double) * 2 * n);
    fftw_iodim dims3[1]; dims3[0].n = n; dims3[0].is = 2; dims3[0].os = 2;
    fftw_plan p3 = fftw_plan_guru_dft(1, dims3, 0, NULL,
                                       (fftw_complex*)in3, (fftw_complex*)out3,
                                       FFTW_FORWARD, FFTW_ESTIMATE);

    /* stride2, in-place */
    double *io4 = fftw_malloc(sizeof(double) * 2 * n);
    fftw_iodim dims4[1]; dims4[0].n = n; dims4[0].is = 2; dims4[0].os = 2;
    fftw_plan p4 = fftw_plan_guru_dft(1, dims4, 0, NULL,
                                       (fftw_complex*)io4, (fftw_complex*)io4,
                                       FFTW_FORWARD, FFTW_ESTIMATE);

    printf("n=%d stride1_outofplace_q1=%d stride1_inplace_q1=%d stride2_outofplace_q1=%d stride2_inplace_q1=%d\n",
           n, has_q1(p1), has_q1(p2), has_q1(p3), has_q1(p4));
    return 0;
}
