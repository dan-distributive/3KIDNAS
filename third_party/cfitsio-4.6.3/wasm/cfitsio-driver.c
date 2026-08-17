/*
 * cfitsio-driver.c -- thin C entry points around libcfitsio, exported for
 * Emscripten's ccall/cwrap. Like SoFiA-2 and unlike FFTW, cfitsio is
 * fundamentally file-I/O shaped (fits_open_file/fits_create_file operate
 * on paths), so these entry points work against Emscripten's in-memory
 * filesystem (MEMFS) -- preload a FITS file's bytes, call one of these,
 * read the result back out -- the same shape as SoFiA's MEMFS story, just
 * without a callMain()/CLI layer since cfitsio is a library, not a tool.
 */
#include <fitsio.h>
#include <emscripten.h>
#include <string.h>
#include <stdlib.h>

/* Fills naxis (number of axes) and up to 8 axis lengths (0-padded) for
 * the current HDU's image, plus bitpix (FITS pixel data type code).
 * Returns a cfitsio status code (0 = success). */
EMSCRIPTEN_KEEPALIVE
int cfits_read_image_info_wasm(const char *path, int *naxis, long *naxes /* [8] */, int *bitpix) {
  fitsfile *fptr;
  int status = 0;

  if (fits_open_file(&fptr, path, READONLY, &status)) return status;

  int nfound = 0;
  if (fits_get_img_dim(fptr, naxis, &status)) { fits_close_file(fptr, &status); return status; }
  if (fits_get_img_type(fptr, bitpix, &status)) { fits_close_file(fptr, &status); return status; }

  long dims[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int n = *naxis > 8 ? 8 : *naxis;
  if (fits_get_img_size(fptr, n, dims, &status)) { fits_close_file(fptr, &status); return status; }
  for (int i = 0; i < 8; i++) naxes[i] = dims[i];
  (void) nfound;

  fits_close_file(fptr, &status);
  return status;
}

/* Reads the whole image array (all axes flattened, row-major/FITS order)
 * as doubles into caller-allocated outbuf (must be >= product of naxes
 * elements, from cfits_read_image_info_wasm). Returns a cfitsio status
 * code. */
EMSCRIPTEN_KEEPALIVE
int cfits_read_image_double_wasm(const char *path, double *outbuf, long nelements) {
  fitsfile *fptr;
  int status = 0;
  int anynul = 0;

  if (fits_open_file(&fptr, path, READONLY, &status)) return status;
  fits_read_img(fptr, TDOUBLE, 1, nelements, NULL, outbuf, &anynul, &status);
  fits_close_file(fptr, &status);
  return status;
}

/* Creates a new single-HDU FITS file at path (MEMFS) containing a
 * naxis-dimensional double-precision image, written from data (row-major,
 * FITS axis order, product-of-naxes elements). Overwrites any existing
 * file at path. Returns a cfitsio status code. */
EMSCRIPTEN_KEEPALIVE
int cfits_write_image_double_wasm(const char *path, int naxis, long *naxes, double *data) {
  fitsfile *fptr;
  int status = 0;
  long nelements = 1;
  for (int i = 0; i < naxis; i++) nelements *= naxes[i];

  /* cfitsio's own "!" overwrite-if-exists prefix, since MEMFS files from
   * a prior run in the same worker instance would otherwise collide. */
  char clobberPath[1024];
  snprintf(clobberPath, sizeof(clobberPath), "!%s", path);

  if (fits_create_file(&fptr, clobberPath, &status)) return status;
  if (fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status)) { fits_close_file(fptr, &status); return status; }
  fits_write_img(fptr, TDOUBLE, 1, nelements, data, &status);
  fits_close_file(fptr, &status);
  return status;
}

/* Reads a string-valued header keyword (e.g. "OBJECT", "TELESCOP") into
 * caller-allocated outbuf (>= 71 bytes, FITS's own max keyword-value
 * length). Returns a cfitsio status code (202 = keyword not found). */
EMSCRIPTEN_KEEPALIVE
int cfits_read_key_str_wasm(const char *path, const char *keyname, char *outbuf) {
  fitsfile *fptr;
  int status = 0;
  char comment[128];

  if (fits_open_file(&fptr, path, READONLY, &status)) return status;
  fits_read_key(fptr, TSTRING, keyname, outbuf, comment, &status);
  fits_close_file(fptr, &status);
  return status;
}

/* Reads a double-valued header keyword (e.g. "CRVAL1", "BSCALE"). */
EMSCRIPTEN_KEEPALIVE
int cfits_read_key_dbl_wasm(const char *path, const char *keyname, double *outval) {
  fitsfile *fptr;
  int status = 0;
  char comment[128];

  if (fits_open_file(&fptr, path, READONLY, &status)) return status;
  fits_read_key(fptr, TDOUBLE, keyname, outval, comment, &status);
  fits_close_file(fptr, &status);
  return status;
}

/* Human-readable text for a cfitsio status code, into caller-allocated
 * outbuf (>= 31 bytes, FITS_ERRMSG_LEN). */
EMSCRIPTEN_KEEPALIVE
void cfits_status_message_wasm(int status, char *outbuf) {
  fits_get_errstatus(status, outbuf);
}

/*
 * ---------------------------------------------------------------------
 * Stateful create/write-header/write-data/close primitives.
 *
 * cfits_write_image_double_wasm above is one-shot (create + write data,
 * no header keywords beyond what fits_create_img sets automatically) --
 * fine for a bare data dump, not enough for a file another tool (e.g.
 * SoFiA) needs to parse a WCS/beam out of. These four functions expose
 * cfitsio's real create/write-key/write-data/close sequence as separate
 * calls so a JS caller can write an arbitrary number of header keywords
 * in between create and data-write.
 *
 * Holds one fitsfile* in a static global rather than threading a handle
 * back through JS: safe because each DCP work function call gets a
 * fresh Emscripten Module instance (see sofia-wasm.js's own doc comment
 * for the same reasoning), so there's never more than one file open at
 * a time within a single instance, and JS itself is single-threaded, so
 * these calls can't interleave.
 * --------------------------------------------------------------------- */
static fitsfile *g_wfptr = NULL;

EMSCRIPTEN_KEEPALIVE
int cfits_create_image_wasm(const char *path, int naxis, long *naxes) {
  int status = 0;
  if (g_wfptr) { fits_close_file(g_wfptr, &status); g_wfptr = NULL; }

  /* cfitsio's own "!" overwrite-if-exists prefix, since MEMFS files from
   * a prior run in the same worker instance would otherwise collide. */
  char clobberPath[1024];
  snprintf(clobberPath, sizeof(clobberPath), "!%s", path);

  if (fits_create_file(&g_wfptr, clobberPath, &status)) return status;
  if (fits_create_img(g_wfptr, DOUBLE_IMG, naxis, naxes, &status)) {
    fits_close_file(g_wfptr, &status);
    g_wfptr = NULL;
    return status;
  }
  return status;
}

/* value passed by address (matches cfitsio's own fits_write_key contract
 * for numeric datatypes). */
EMSCRIPTEN_KEEPALIVE
int cfits_write_key_dbl_wasm(const char *keyname, double value, const char *comment) {
  if (!g_wfptr) return -1;
  int status = 0;
  fits_write_key(g_wfptr, TDOUBLE, (char *) keyname, &value, (char *) comment, &status);
  return status;
}

EMSCRIPTEN_KEEPALIVE
int cfits_write_key_str_wasm(const char *keyname, const char *value, const char *comment) {
  if (!g_wfptr) return -1;
  int status = 0;
  fits_write_key(g_wfptr, TSTRING, (char *) keyname, (void *) value, (char *) comment, &status);
  return status;
}

EMSCRIPTEN_KEEPALIVE
int cfits_write_image_data_wasm(double *data, long nelements) {
  if (!g_wfptr) return -1;
  int status = 0;
  fits_write_img(g_wfptr, TDOUBLE, 1, nelements, data, &status);
  return status;
}

/* Closes the file opened by cfits_create_image_wasm and returns its
 * bytes have already been written to path by cfitsio itself (readable
 * afterward via Module.FS.readFile(path), same as every other output
 * file in this driver) -- this just finalizes and releases the handle. */
EMSCRIPTEN_KEEPALIVE
int cfits_close_image_wasm(void) {
  if (!g_wfptr) return -1;
  int status = 0;
  fits_close_file(g_wfptr, &status);
  g_wfptr = NULL;
  return status;
}
