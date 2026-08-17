/*
 * fftw-driver.c -- thin C entry points around libfftw3, exported for
 * Emscripten's ccall/cwrap. Unlike SoFiA-2 (a CLI tool with a natural
 * callMain() entry point), FFTW is a library: there's no main() to call,
 * so the wasm/JS boundary here is a couple of exported C functions
 * operating directly on caller-allocated heap buffers, not file I/O.
 *
 * Both functions accept plain double* buffers (not fftw_malloc'd) --
 * safe because this build has no SIMD codelets enabled (see build.sh:
 * plain `./configure`, no --enable-sse2/--enable-avx/etc.), so there's
 * no alignment requirement beyond what ordinary malloc already gives.
 */
#include <fftw3.h>
#include <emscripten.h>
#include <string.h>

/* Complex-to-complex 1D DFT. sign is FFTW_FORWARD (-1) or FFTW_BACKWARD (+1).
 * re_in/im_in/re_out/im_out are each n-element double arrays in wasm
 * memory (caller allocates via Module._malloc, frees the same way). */
EMSCRIPTEN_KEEPALIVE
int fftw_dft_1d_wasm(double *re_in, double *im_in, double *re_out, double *im_out, int n, int sign) {
  if (n <= 0) return 1;

  fftw_complex *in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (size_t) n);
  fftw_complex *out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (size_t) n);
  if (!in || !out) { fftw_free(in); fftw_free(out); return 2; }

  for (int i = 0; i < n; i++) {
    in[i][0] = re_in[i];
    in[i][1] = im_in[i];
  }

  fftw_plan p = fftw_plan_dft_1d(n, in, out, sign, FFTW_ESTIMATE);
  if (!p) { fftw_free(in); fftw_free(out); return 3; }
  fftw_execute(p);

  for (int i = 0; i < n; i++) {
    re_out[i] = out[i][0];
    im_out[i] = out[i][1];
  }

  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
  return 0;
}

/* Real-to-complex 1D DFT (forward only -- r2c is inherently one direction).
 * re_in is n reals; re_out/im_out are each (n/2+1)-element arrays (the
 * non-redundant half of the spectrum, FFTW's standard r2c output size). */
EMSCRIPTEN_KEEPALIVE
int fftw_r2c_1d_wasm(double *re_in, double *re_out, double *im_out, int n) {
  if (n <= 0) return 1;
  int nc = n / 2 + 1;

  double *in = (double *) fftw_malloc(sizeof(double) * (size_t) n);
  fftw_complex *out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (size_t) nc);
  if (!in || !out) { fftw_free(in); fftw_free(out); return 2; }

  memcpy(in, re_in, sizeof(double) * (size_t) n);

  fftw_plan p = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
  if (!p) { fftw_free(in); fftw_free(out); return 3; }
  fftw_execute(p);

  for (int i = 0; i < nc; i++) {
    re_out[i] = out[i][0];
    im_out[i] = out[i][1];
  }

  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
  return 0;
}

/* Real-to-complex 2D DFT (forward), using FFTW's OWN native 2D planner
 * directly -- NOT composed from row-then-column 1D calls the way
 * FFTW3WasmRank2.js's rdft2R2cSync does (see that file's header for why
 * it was ever built that way: no native 2D entry point existed before
 * this). That composition costs ~150-170 individual 1D FFTW calls per
 * channel (N0 rows + NC columns, forward and inverse combined) -- found
 * to be the dominant cost of the whole objective-function evaluation
 * (~85-90% of eval time) via direct per-eval instrumentation. This
 * function does the real work in ONE fftw_execute call instead.
 *
 * re_in is n0*n1 reals (row-major). re_out/im_out are each n0*nc
 * doubles, nc=n1/2+1 (FFTW's r2c_2d non-redundant-spectrum output size,
 * same shape the 1D-composed path already produces). */
EMSCRIPTEN_KEEPALIVE
int fftw_r2c_2d_wasm(double *re_in, double *re_out, double *im_out, int n0, int n1) {
  if (n0 <= 0 || n1 <= 0) return 1;
  int nc = n1 / 2 + 1;

  double *in = (double *) fftw_malloc(sizeof(double) * (size_t) n0 * n1);
  fftw_complex *out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (size_t) n0 * nc);
  if (!in || !out) { fftw_free(in); fftw_free(out); return 2; }

  memcpy(in, re_in, sizeof(double) * (size_t) n0 * n1);

  // FFTW_PRESERVE_INPUT matches TwoDConvolution.f's own
  // dfftw_plan_dft_r2c_2d call exactly (FFTW_ESTIMATE,FFTW_PRESERVE_INPUT)
  // -- omitting it here was a real, confirmed mismatch: this plan and
  // Fortran's were asking FFTW for different guarantees on the same
  // transform, letting the planner pick a different internal algorithm
  // even though both ultimately call fftw_plan_dft_r2c_2d. Found via a
  // bit-exact stage checksum (sum of the model cube's flux, printed as a
  // hex double) taken right before and right after CubeBeamConvolution:
  // the two implementations agreed to ~1e-9 relative before convolution
  // and only ~1e-7 after -- convolution was the only step where the
  // agreement measurably degraded.
  fftw_plan p2 = fftw_plan_dft_r2c_2d(n0, n1, in, out,
      FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  if (!p2) { fftw_free(in); fftw_free(out); return 3; }
  fftw_execute(p2);

  for (int i = 0; i < n0 * nc; i++) {
    re_out[i] = out[i][0];
    im_out[i] = out[i][1];
  }

  fftw_destroy_plan(p2);
  fftw_free(in);
  fftw_free(out);
  return 0;
}

/* Complex-to-real 2D DFT (inverse), FFTW's native 2D planner. UNNORMALIZED
 * (matching this project's existing convention -- forward+inverse without
 * dividing by n0*n1 recovers n0*n1 * original; CubeKernelConvolution.js's
 * caller already does that division itself).
 *
 * re_in/im_in are each n0*nc doubles (nc=n1/2+1, half-complex spectrum,
 * same layout fftw_r2c_2d_wasm produces). re_out is n0*n1 reals
 * (row-major).
 *
 * FFTW_PRESERVE_INPUT: TRIED here to match TwoDConvolution.f's own
 * dfftw_plan_dft_c2r_2d call (FFTW_ESTIMATE,FFTW_PRESERVE_INPUT) exactly --
 * confirmed EMPIRICALLY NOT SUPPORTED by this build: fftw_plan_dft_c2r_2d
 * returns NULL (this function's own return-3 path) the moment the flag is
 * added, for this exact transform size. Fortran's native, fully-configured
 * FFTW build has the codelets to honor it for a multi-dimensional c2r
 * transform; this WASM build's deliberately reduced "generic, no-SIMD"
 * codelet set (see build.sh) apparently doesn't. Reverted to plain
 * FFTW_ESTIMATE here -- the forward r2c plan above keeps
 * FFTW_PRESERVE_INPUT (that one DOES succeed), so the two sides are now
 * matched on the forward transform and knowingly mismatched on the
 * inverse one, for a build-configuration reason, not an oversight. See
 * fftw_r2c_2d_wasm's own comment for the checksum trace that found the
 * original discrepancy this was chasing. */
EMSCRIPTEN_KEEPALIVE
int fftw_c2r_2d_wasm(double *re_in, double *im_in, double *re_out, int n0, int n1) {
  if (n0 <= 0 || n1 <= 0) return 1;
  int nc = n1 / 2 + 1;

  fftw_complex *in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (size_t) n0 * nc);
  double *out = (double *) fftw_malloc(sizeof(double) * (size_t) n0 * n1);
  if (!in || !out) { fftw_free(in); fftw_free(out); return 2; }

  for (int i = 0; i < n0 * nc; i++) {
    in[i][0] = re_in[i];
    in[i][1] = im_in[i];
  }

  fftw_plan p3 = fftw_plan_dft_c2r_2d(n0, n1, in, out, FFTW_ESTIMATE);
  if (!p3) { fftw_free(in); fftw_free(out); return 3; }
  fftw_execute(p3);

  memcpy(re_out, out, sizeof(double) * (size_t) n0 * n1);

  fftw_destroy_plan(p3);
  fftw_free(in);
  fftw_free(out);
  return 0;
}
