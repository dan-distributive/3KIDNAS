/* q1_probe.c -- probe a real, in-place FFTW complex DFT of size N with each
 * standard basis vector (delta at k, value 1+0i), dumping the full N x N
 * complex transformation matrix. Used to empirically recover the q1-based
 * inner sub-transform as a pure black-box linear map (see chooseDecomposition
 * .js's chooseComplexEmbedded header for context) without needing to decode
 * dft/ct.c's tensor construction directly. */
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>

int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    fftw_complex *io = fftw_malloc(sizeof(fftw_complex) * n);
    fftw_plan p = fftw_plan_dft_1d(n, io, io, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) { io[i][0] = 0; io[i][1] = 0; }
        io[k][0] = 1.0; io[k][1] = 0.0;
        fftw_execute(p);
        for (int i = 0; i < n; i++) printf("%.17g %.17g\n", io[i][0], io[i][1]);
    }
    return 0;
}
