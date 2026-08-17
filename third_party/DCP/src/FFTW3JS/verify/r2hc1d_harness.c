/* r2hc1d_harness.c -- minimal 1D real<->halfcomplex ground-truth generator,
 * same vendored/patched libfftw3.a as dft1d_harness.c (see build.sh).
 * Usage: ./r2hc1d_harness r2hc N   -- prints N lines "x[i]" of real input,
 *                                     then N lines "O[i]" of
 *                                     fftw_plan_r2r_1d(FFTW_R2HC,
 *                                     FFTW_ESTIMATE) packed-halfcomplex
 *                                     output.
 *        ./r2hc1d_harness hc2r N   -- prints N lines "O[i]" of packed-
 *                                     halfcomplex input (synthesized from a
 *                                     smooth signal's R2HC transform, so it's
 *                                     a valid halfcomplex vector), then N
 *                                     lines "x[i]" of
 *                                     fftw_plan_r2r_1d(FFTW_HC2R,
 *                                     FFTW_ESTIMATE) real output.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s r2hc|hc2r N\n", argv[0]); return 1; }
    int n = atoi(argv[2]);

    if (strcmp(argv[1], "r2hc") == 0) {
        double *in = fftw_malloc(sizeof(double) * n);
        double *out = fftw_malloc(sizeof(double) * n);
        for (int i = 0; i < n; i++) in[i] = sin(i * 1.7 + 0.3) + 0.5 * cos(i * 0.9 + 1.1);
        fftw_plan p = fftw_plan_r2r_1d(n, in, out, FFTW_R2HC, FFTW_ESTIMATE);
        fftw_execute(p);
        for (int i = 0; i < n; i++) printf("%.17g\n", in[i]);
        for (int i = 0; i < n; i++) printf("%.17g\n", out[i]);
        fftw_destroy_plan(p);
        fftw_free(in); fftw_free(out);
    } else if (strcmp(argv[1], "hc2r") == 0) {
        /* Build a valid halfcomplex input by running R2HC on a smooth
         * signal first (guarantees valid decay/structure, not required by
         * HC2R but keeps values well-scaled). */
        double *sig = fftw_malloc(sizeof(double) * n);
        double *hc = fftw_malloc(sizeof(double) * n);
        for (int i = 0; i < n; i++) sig[i] = sin(i * 2.3 + 0.7) - 0.3 * cos(i * 1.3 + 0.2);
        fftw_plan pf = fftw_plan_r2r_1d(n, sig, hc, FFTW_R2HC, FFTW_ESTIMATE);
        fftw_execute(pf);
        fftw_destroy_plan(pf);

        double *hc_in = fftw_malloc(sizeof(double) * n);
        double *out = fftw_malloc(sizeof(double) * n);
        memcpy(hc_in, hc, sizeof(double) * n);
        fftw_plan pb = fftw_plan_r2r_1d(n, hc_in, out, FFTW_HC2R, FFTW_ESTIMATE);
        fftw_execute(pb);
        for (int i = 0; i < n; i++) printf("%.17g\n", hc[i]);
        for (int i = 0; i < n; i++) printf("%.17g\n", out[i]);
        fftw_destroy_plan(pb);
        fftw_free(sig); fftw_free(hc); fftw_free(hc_in); fftw_free(out);
    } else {
        fprintf(stderr, "unknown mode %s\n", argv[1]);
        return 1;
    }
    return 0;
}
