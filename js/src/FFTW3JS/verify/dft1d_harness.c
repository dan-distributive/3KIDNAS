/* dft1d_harness.c -- minimal 1D complex-DFT ground-truth generator, same
 * vendored/patched libfftw3.a as ground_truth_harness.c (see build.sh).
 * Usage: ./dft1d_harness N seed  -- prints N lines "re im" of input, then
 * N lines "re im" of fftw_plan_dft_1d(FFTW_FORWARD, FFTW_ESTIMATE) output.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s N\n", argv[0]); return 1; }
    int n = atoi(argv[1]);

    fftw_complex *in = fftw_malloc(sizeof(fftw_complex) * n);
    fftw_complex *out = fftw_malloc(sizeof(fftw_complex) * n);
    for (int i = 0; i < n; i++) {
        in[i][0] = sin(i * 1.7 + 0.3);
        in[i][1] = cos(i * 0.9 + 1.1);
    }
    fftw_plan p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    for (int i = 0; i < n; i++) printf("%.17g %.17g\n", in[i][0], in[i][1]);
    for (int i = 0; i < n; i++) printf("%.17g %.17g\n", out[i][0], out[i][1]);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
    return 0;
}
