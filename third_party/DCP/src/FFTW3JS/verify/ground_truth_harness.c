/* ground_truth_harness.c
 * ==========================================================================
 * Permanent ground-truth generator for the FFTW3-to-JS port
 * (third_party/DCP/src/FFTW3JS/). Links the project's own vendored,
 * -ffp-contract=off, fdlibm-patched libfftw3.a (see
 * third_party/fftw-3.3.8/kernel/trig.c and src/makeflags for why those two
 * patches exist) and drives it exactly the way Fortran's TwoDConvolution.f
 * does: FFTW_ESTIMATE, real-to-complex / complex-to-real, no runtime
 * planning/benchmarking, so results are reproducible run to run.
 *
 * Dumps RAW double-precision bytes, not decimal text -- decimal printf can
 * silently mask a genuine 1-ULP mismatch (this bit us once already this
 * session: default Fortran list-directed output looked identical across a
 * real numeric change until printed at full width). The companion
 * compare.js does a byte-exact diff against the JS port's output.
 *
 * Usage:
 *   ground_truth_harness fwd  <N0> <N1> <input_type> [input_file] <out_file>
 *   ground_truth_harness inv  <N0> <N1> <input_type> [input_file] <out_file>
 *   ground_truth_harness plan <N0> <N1>
 *
 * input_type: delta   -- impulse at index 0 (real input for fwd; for inv,
 *                         a unit impulse at the DC complex bin)
 *             random   -- deterministic pseudo-random input (srand(42))
 *             file     -- read raw double bytes from input_file
 *                         (N0*N1 reals for fwd, N0*(N1/2+1) interleaved
 *                         re/im complex doubles for inv)
 *
 * fwd:  real N0xN1 array -> complex N0 x (N1/2+1) array (interleaved re/im),
 *       written raw to out_file. UNNORMALIZED (matches dfftw_execute_dft_r2c).
 * inv:  complex N0 x (N1/2+1) array (interleaved re/im) -> real N0xN1 array,
 *       written raw to out_file. UNNORMALIZED (matches dfftw_execute_dft_c2r
 *       -- Fortran applies its own explicit /SizePad(1)/SizePad(2) after;
 *       replicate that on whichever side consumes this, not here).
 * plan: prints fftw_print_plan output to stdout for documentation --
 *       this is how every (N0,N1) size's actual FFTW3 decomposition gets
 *       recorded before porting it (see PlanTable.js).
 *
 * Source: this harness's structure (link the vendored lib, dump ground
 * truth, diff byte-for-byte in a companion script) mirrors the existing
 * fdlibm verification pattern already used in this project
 * (src/StandardMath/fdlibm_fortran_wrappers.c and its test harnesses).
 *
 * FFTW_PRESERVE_INPUT: intentionally NOT used, despite TwoDConvolution.f's
 * source text appearing to request it. Traced and confirmed: Fortran's
 * legacy dfftw_plan_dft_c2r_2d wrapper (api/f77funcs.h) has a C signature
 * taking exactly one `flags` int; the Fortran call site passes it as TWO
 * trailing arguments (FFTW_ESTIMATE, FFTW_PRESERVE_INPUT), and the classic
 * F77 calling convention silently drops the second, unread one. So
 * Fortran's actual plans have always been FFTW_ESTIMATE-only (default
 * DESTROY_INPUT semantics for c2r) -- confirmed harmless in the real
 * pipeline (TwoDConvolution.f's c2r input array, ComplexConvolve, is never
 * read again after the c2r call), but it does mean requesting
 * FFTW_PRESERVE_INPUT here would both mismatch Fortran's actual behavior
 * AND hit a separate, real FFTW3 limitation: multi-dimensional c2r
 * transforms have no input-preserving algorithm at all, and the planner
 * returns NULL if PRESERVE_INPUT is requested for one (confirmed via
 * FFTW3's own doc/tutorial.texi and empirically -- direct C-API calls with
 * that flag crash on a NULL plan; matching Fortran's actual ESTIMATE-only
 * calls do not).
 * ==========================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

static double *read_doubles(const char *path, size_t count) {
    FILE *f = fopen(path, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    double *buf = malloc(count * sizeof(double));
    size_t got = fread(buf, sizeof(double), count, f);
    if (got != count) {
        fprintf(stderr, "expected %zu doubles, read %zu from %s\n", count, got, path);
        exit(1);
    }
    fclose(f);
    return buf;
}

static void write_doubles(const char *path, const double *buf, size_t count) {
    FILE *f = fopen(path, "wb");
    if (!f) { fprintf(stderr, "cannot open %s for write\n", path); exit(1); }
    fwrite(buf, sizeof(double), count, f);
    fclose(f);
}

static void fill_delta(double *buf, size_t count) {
    memset(buf, 0, count * sizeof(double));
    buf[0] = 1.0;
}

static void fill_random(double *buf, size_t count) {
    srand(42);
    for (size_t i = 0; i < count; i++) {
        buf[i] = ((double)rand() / (double)RAND_MAX) * 2.0 - 1.0;
    }
}

static void do_fwd(int N0, int N1, const char *input_type, const char *input_file, const char *out_file) {
    size_t nreal = (size_t)N0 * (size_t)N1;
    size_t ncomplex = (size_t)N0 * (size_t)(N1 / 2 + 1);

    double *in = fftw_malloc(sizeof(double) * nreal);
    fftw_complex *out = fftw_malloc(sizeof(fftw_complex) * ncomplex);

    if (strcmp(input_type, "delta") == 0) {
        fill_delta(in, nreal);
    } else if (strcmp(input_type, "random") == 0) {
        fill_random(in, nreal);
    } else if (strcmp(input_type, "file") == 0) {
        double *fromfile = read_doubles(input_file, nreal);
        memcpy(in, fromfile, nreal * sizeof(double));
        free(fromfile);
    } else {
        fprintf(stderr, "unknown input_type %s\n", input_type);
        exit(1);
    }

    fftw_plan p = fftw_plan_dft_r2c_2d(N0, N1, in, out, FFTW_ESTIMATE);
    fftw_execute(p);

    /* out is fftw_complex = double[2] (re, im) per element -- already the
     * exact interleaved layout we want to compare against. */
    write_doubles(out_file, (double *)out, ncomplex * 2);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}

static void do_inv(int N0, int N1, const char *input_type, const char *input_file, const char *out_file) {
    size_t nreal = (size_t)N0 * (size_t)N1;
    size_t ncomplex = (size_t)N0 * (size_t)(N1 / 2 + 1);

    fftw_complex *in = fftw_malloc(sizeof(fftw_complex) * ncomplex);
    double *out = fftw_malloc(sizeof(double) * nreal);

    if (strcmp(input_type, "delta") == 0) {
        memset(in, 0, ncomplex * sizeof(fftw_complex));
        in[0][0] = 1.0; /* unit DC bin */
    } else if (strcmp(input_type, "random") == 0) {
        fill_random((double *)in, ncomplex * 2);
    } else if (strcmp(input_type, "file") == 0) {
        double *fromfile = read_doubles(input_file, ncomplex * 2);
        memcpy(in, fromfile, ncomplex * 2 * sizeof(double));
        free(fromfile);
    } else {
        fprintf(stderr, "unknown input_type %s\n", input_type);
        exit(1);
    }

    fftw_plan p = fftw_plan_dft_c2r_2d(N0, N1, in, out, FFTW_ESTIMATE);
    fftw_execute(p);

    write_doubles(out_file, out, nreal);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}

static void do_plan(int N0, int N1) {
    double *in = fftw_malloc(sizeof(double) * N0 * N1);
    fftw_complex *out = fftw_malloc(sizeof(fftw_complex) * N0 * (N1 / 2 + 1));

    fftw_plan pf = fftw_plan_dft_r2c_2d(N0, N1, in, out, FFTW_ESTIMATE);
    printf("=== r2c plan for %dx%d ===\n", N0, N1);
    fftw_print_plan(pf);
    printf("\n\n");

    fftw_plan pi = fftw_plan_dft_c2r_2d(N0, N1, out, in, FFTW_ESTIMATE);
    printf("=== c2r plan for %dx%d ===\n", N0, N1);
    fftw_print_plan(pi);
    printf("\n\n");

    fftw_destroy_plan(pf);
    fftw_destroy_plan(pi);
    fftw_free(in);
    fftw_free(out);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s fwd|inv <N0> <N1> <input_type> [input_file] <out_file>\n", argv[0]);
        fprintf(stderr, "       %s plan <N0> <N1>\n", argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "plan") == 0) {
        if (argc != 4) { fprintf(stderr, "plan needs N0 N1\n"); return 1; }
        do_plan(atoi(argv[2]), atoi(argv[3]));
        return 0;
    }

    if (strcmp(argv[1], "fwd") == 0 || strcmp(argv[1], "inv") == 0) {
        int N0 = atoi(argv[2]);
        int N1 = atoi(argv[3]);
        const char *input_type = argv[4];
        const char *input_file = NULL;
        const char *out_file;
        if (strcmp(input_type, "file") == 0) {
            if (argc != 7) { fprintf(stderr, "file mode needs input_file and out_file\n"); return 1; }
            input_file = argv[5];
            out_file = argv[6];
        } else {
            if (argc != 6) { fprintf(stderr, "%s mode needs out_file\n", input_type); return 1; }
            out_file = argv[5];
        }

        if (strcmp(argv[1], "fwd") == 0) {
            do_fwd(N0, N1, input_type, input_file, out_file);
        } else {
            do_inv(N0, N1, input_type, input_file, out_file);
        }
        return 0;
    }

    fprintf(stderr, "unknown mode %s\n", argv[1]);
    return 1;
}
