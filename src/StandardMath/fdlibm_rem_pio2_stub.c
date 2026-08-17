#include <stdio.h>
#include <stdlib.h>
/* Stub: the large-argument Payne-Hanek path is intentionally not ported
 * (see fdlibm.js SCOPE comment) -- nothing in this test range should ever
 * reach it. Abort loudly rather than silently doing the wrong thing. */
int __kernel_rem_pio2(double *x, double *y, int e0, int nx, int prec) {
  fprintf(stderr, "__kernel_rem_pio2 stub called -- large-argument path "
                  "not implemented, this should never happen for our test range\n");
  abort();
  return 0;
}
