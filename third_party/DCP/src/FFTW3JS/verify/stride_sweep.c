/* stride_sweep.c -- for a given N, print the fftw_print_plan output for
 * BOTH a stride-1 (plain fftw_plan_dft_1d) and stride-2 (matching dft/
 * bluestein.c's / dft/rader.c's internal cldf/cld_omega tensor shape:
 * mktensor_1d(N,2,2)) complex DFT of the same size, to characterize when
 * FFTW's planner picks q1 (twiddle-squared) under non-unit stride but not
 * under unit stride for the SAME logical decomposition. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

int main(int argc, char **argv) {
    int n = atoi(argv[1]);

    fftw_complex *in1 = fftw_malloc(sizeof(fftw_complex) * n);
    fftw_complex *out1 = fftw_malloc(sizeof(fftw_complex) * n);
    fftw_plan p1 = fftw_plan_dft_1d(n, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);

    double *buf = fftw_malloc(sizeof(double) * 2 * n);
    fftw_iodim dims[1];
    dims[0].n = n; dims[0].is = 2; dims[0].os = 2;
    fftw_plan p2 = fftw_plan_guru_dft(1, dims, 0, NULL,
                                       (fftw_complex*)buf, (fftw_complex*)buf,
                                       FFTW_FORWARD, FFTW_ESTIMATE);

    char buf1[4096], buf2[4096];
    FILE *f1 = fmemopen(buf1, sizeof(buf1), "w");
    FILE *f2 = fmemopen(buf2, sizeof(buf2), "w");
    fftw_fprint_plan(p1, f1);
    fftw_fprint_plan(p2, f2);
    fclose(f1); fclose(f2);

    int hasQ1_stride1 = strstr(buf1, "\"q1_") != NULL;
    int hasQ1_stride2 = strstr(buf2, "\"q1_") != NULL;

    printf("n=%d stride1_q1=%d stride2_q1=%d %s\n", n, hasQ1_stride1, hasQ1_stride2,
           (hasQ1_stride1 != hasQ1_stride2) ? "*** DIVERGES ***" : "");
    if (hasQ1_stride1 != hasQ1_stride2) {
        printf("  stride1: %s\n", buf1);
        printf("  stride2: %s\n", buf2);
    }
    return 0;
}
