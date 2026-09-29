/* One-off diagnostic (Dan probe, 2026-09-27): exposes fftw_export_wisdom_to_string()
 * to Fortran via a plain char-buffer signature (Fortran's F77 wisdom API
 * uses an awkward write-callback; this sidesteps it). Not currently wired
 * into the build system Makefiles -- compiled and linked manually for
 * this investigation. See TwoDConvolution.f's matching caller. */
#include <fftw3.h>
#include <string.h>
#include <stdlib.h>

void fftw_export_wisdom_cstr_(char *buf, int *buflen) {
  char *s = fftw_export_wisdom_to_string();
  if (s) {
    strncpy(buf, s, (size_t)(*buflen - 1));
    buf[*buflen - 1] = '\0';
    free(s);
  } else {
    buf[0] = '\0';
  }
}
