/*
 * fdlibm-driver.c -- thin scalar entry points around this project's
 * vendored fdlibm sin/cos/atanh, exported for Emscripten's ccall/cwrap.
 * Unlike fftw-driver.c/cfitsio-driver.c, these are plain double-in/
 * double-out calls -- no caller-allocated heap buffers, no malloc/free,
 * no HEAPF64 bookkeeping. cwrap's 'number' arg/return type marshals a
 * JS double directly, so the JS wrapper (fdlibm-wasm.js) needs nothing
 * beyond a cwrap'd function reference per entry point.
 */
#include <emscripten.h>

double fdlibm_sin(double x);
double fdlibm_cos(double x);
double fdlibm_atanh(double x);

EMSCRIPTEN_KEEPALIVE
double fd_sin_wasm(double x) { return fdlibm_sin(x); }

EMSCRIPTEN_KEEPALIVE
double fd_cos_wasm(double x) { return fdlibm_cos(x); }

EMSCRIPTEN_KEEPALIVE
double fd_atanh_wasm(double x) { return fdlibm_atanh(x); }
