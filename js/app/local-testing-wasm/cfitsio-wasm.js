/**
 * @file        cfitsio-wasm.js
 * @description Single entry point for running cfitsio 4.6.3 (real
 *              standards-compliant FITS I/O, see build.sh) inside a DCP
 *              worker. job.requires(['./cfitsio-wasm']) is enough -- this
 *              file require()s cfitsio-module.js, the Emscripten build of
 *              libcfitsio, with the wasm binary embedded directly in the
 *              JS (SINGLE_FILE=1).
 *
 *              Like SoFiA-2 and unlike FFTW, cfitsio is file-I/O shaped
 *              (fits_open_file/fits_create_file operate on paths), so
 *              this wraps Emscripten's in-memory filesystem (MEMFS):
 *              write bytes in, call a driver function, read bytes back
 *              out. Each exported function here gets a fresh Module
 *              instance (fresh MEMFS), so concurrent/repeated calls in
 *              the same worker never see each other's files.
 *
 * @usage       const cfitsio = require('./cfitsio-wasm');
 *              const { naxes, bitpix } = await cfitsio.readImageInfo(fitsBytes);
 *              const { data, naxes } = await cfitsio.readImageDouble(fitsBytes);
 *              const value = await cfitsio.readKeyDouble(fitsBytes, 'CRPIX1');
 *              const fitsBytes = await cfitsio.writeImageDouble(naxes, data);
 */
'use strict';

// require() exists in Node and inside a DCP sandbox (via job.requires()'s
// own module bundle) but not in a raw browser tab -- cfitsio-module.js
// itself is already Emscripten-UMD (a bare `var CfitsioModule = ...` at
// module scope, so loading it via a plain <script> tag already puts it on
// globalThis for free; only THIS hand-written wrapper needs the guard.
const createCfitsioModule = (typeof require === 'function')
  ? require('./cfitsio-module')
  : globalThis.CfitsioModule;

const WORK_PATH = '/work/cube.fits';
const OUT_PATH = '/work/out.fits';
const MAX_NAXIS = 8;

async function withModule(fn) {
  const Module = await createCfitsioModule();
  return fn(Module);
}

function checkStatus(Module, rc, label) {
  if (rc === 0) return;
  const p = Module._malloc(32); // FLEN_ERRMSG
  Module.ccall('cfits_status_message_wasm', null, ['number', 'number'], [rc, p]);
  const msg = Module.UTF8ToString(p);
  Module._free(p);
  throw new Error(`${label} failed: cfitsio status ${rc} (${msg})`);
}

/**
 * Reads basic image dimension info from a FITS file's primary/current HDU.
 * @param {Uint8Array} fitsBytes
 * @returns {Promise<{naxis: number, naxes: number[], bitpix: number}>}
 */
async function readImageInfo(fitsBytes) {
  return withModule((Module) => {
    Module.FS.mkdir('/work');
    Module.FS.writeFile(WORK_PATH, fitsBytes);

    const naxisP = Module._malloc(4);
    const naxesP = Module._malloc(MAX_NAXIS * 4);
    const bitpixP = Module._malloc(4);
    try {
      const rc = Module.ccall(
        'cfits_read_image_info_wasm', 'number',
        ['string', 'number', 'number', 'number'],
        [WORK_PATH, naxisP, naxesP, bitpixP]
      );
      checkStatus(Module, rc, 'readImageInfo');
      const naxis = Module.HEAP32[naxisP / 4];
      const naxes = Array.from({ length: naxis }, (_, i) => Module.HEAP32[naxesP / 4 + i]);
      const bitpix = Module.HEAP32[bitpixP / 4];
      return { naxis, naxes, bitpix };
    } finally {
      [naxisP, naxesP, bitpixP].forEach((p) => Module._free(p));
    }
  });
}

/**
 * Reads a FITS image's full pixel array as doubles.
 * @param {Uint8Array} fitsBytes
 * @returns {Promise<{naxes: number[], data: number[]}>}
 */
async function readImageDouble(fitsBytes) {
  return withModule(async (Module) => {
    Module.FS.mkdir('/work');
    Module.FS.writeFile(WORK_PATH, fitsBytes);

    const naxisP = Module._malloc(4);
    const naxesP = Module._malloc(MAX_NAXIS * 4);
    const bitpixP = Module._malloc(4);
    let naxes;
    try {
      const rc = Module.ccall(
        'cfits_read_image_info_wasm', 'number',
        ['string', 'number', 'number', 'number'],
        [WORK_PATH, naxisP, naxesP, bitpixP]
      );
      checkStatus(Module, rc, 'readImageDouble (info)');
      const naxis = Module.HEAP32[naxisP / 4];
      naxes = Array.from({ length: naxis }, (_, i) => Module.HEAP32[naxesP / 4 + i]);
    } finally {
      [naxisP, naxesP, bitpixP].forEach((p) => Module._free(p));
    }

    const nelements = naxes.reduce((a, b) => a * b, 1);
    const dataP = Module._malloc(nelements * 8);
    try {
      const rc = Module.ccall(
        'cfits_read_image_double_wasm', 'number',
        ['string', 'number', 'number'],
        [WORK_PATH, dataP, nelements]
      );
      checkStatus(Module, rc, 'readImageDouble');
      return { naxes, data: Array.from(Module.HEAPF64.subarray(dataP / 8, dataP / 8 + nelements)) };
    } finally {
      Module._free(dataP);
    }
  });
}

/**
 * Reads a double-valued header keyword (e.g. CRVAL1, BSCALE).
 * @param {Uint8Array} fitsBytes
 * @param {string} keyName
 * @returns {Promise<number>}
 */
async function readKeyDouble(fitsBytes, keyName) {
  return withModule((Module) => {
    Module.FS.mkdir('/work');
    Module.FS.writeFile(WORK_PATH, fitsBytes);

    const outP = Module._malloc(8);
    try {
      const rc = Module.ccall(
        'cfits_read_key_dbl_wasm', 'number',
        ['string', 'string', 'number'],
        [WORK_PATH, keyName, outP]
      );
      checkStatus(Module, rc, `readKeyDouble(${keyName})`);
      return Module.HEAPF64[outP / 8];
    } finally {
      Module._free(outP);
    }
  });
}

/**
 * Reads a string-valued header keyword (e.g. OBJECT, TELESCOP).
 * @param {Uint8Array} fitsBytes
 * @param {string} keyName
 * @returns {Promise<string>}
 */
async function readKeyString(fitsBytes, keyName) {
  return withModule((Module) => {
    Module.FS.mkdir('/work');
    Module.FS.writeFile(WORK_PATH, fitsBytes);

    const outP = Module._malloc(72); // FLEN_VALUE
    try {
      const rc = Module.ccall(
        'cfits_read_key_str_wasm', 'number',
        ['string', 'string', 'number'],
        [WORK_PATH, keyName, outP]
      );
      checkStatus(Module, rc, `readKeyString(${keyName})`);
      return Module.UTF8ToString(outP);
    } finally {
      Module._free(outP);
    }
  });
}

/**
 * Creates a new single-HDU FITS file containing a double-precision image.
 * @param {number[]} naxes - axis lengths, FITS order (naxes[0] is the fastest-varying axis)
 * @param {number[]} data - row-major pixel data, product-of-naxes elements
 * @returns {Promise<Uint8Array>} the resulting FITS file's bytes
 */
async function writeImageDouble(naxes, data) {
  return withModule((Module) => {
    Module.FS.mkdir('/work');

    const naxesP = Module._malloc(naxes.length * 4);
    Module.HEAP32.set(naxes, naxesP / 4);
    const dataP = Module._malloc(data.length * 8);
    Module.HEAPF64.set(data, dataP / 8);
    try {
      const rc = Module.ccall(
        'cfits_write_image_double_wasm', 'number',
        ['string', 'number', 'number', 'number'],
        [OUT_PATH, naxes.length, naxesP, dataP]
      );
      checkStatus(Module, rc, 'writeImageDouble');
      return Module.FS.readFile(OUT_PATH);
    } finally {
      [naxesP, dataP].forEach((p) => Module._free(p));
    }
  });
}

/**
 * Creates a new single-HDU FITS file containing a double-precision image,
 * with arbitrary header keywords written after the image is created but
 * before the pixel data -- e.g. a WCS (CRPIX/CRVAL/CDELT/CTYPE/CUNIT) and
 * beam (BMAJ/BMIN/BPA) block, which plain writeImageDouble() has no way
 * to attach.
 * @param {number[]} naxes - axis lengths, FITS order (naxes[0] is the fastest-varying axis)
 * @param {number[]} data - row-major pixel data, product-of-naxes elements
 * @param {Object<string, number|string>} headerEntries - keyword -> value;
 *   type (TDOUBLE vs TSTRING) is inferred from typeof value.
 * @returns {Promise<Uint8Array>} the resulting FITS file's bytes
 */
async function writeImageDoubleWithHeader(naxes, data, headerEntries = {}) {
  return withModule((Module) => {
    Module.FS.mkdir('/work');

    const createImage = Module.cwrap('cfits_create_image_wasm', 'number', ['string', 'number', 'number']);
    const writeKeyDbl = Module.cwrap('cfits_write_key_dbl_wasm', 'number', ['string', 'number', 'string']);
    const writeKeyStr = Module.cwrap('cfits_write_key_str_wasm', 'number', ['string', 'string', 'string']);
    const writeImageData = Module.cwrap('cfits_write_image_data_wasm', 'number', ['number', 'number']);
    const closeImage = Module.cwrap('cfits_close_image_wasm', 'number', []);

    const naxesP = Module._malloc(naxes.length * 4);
    Module.HEAP32.set(naxes, naxesP / 4);
    const dataP = Module._malloc(data.length * 8);
    Module.HEAPF64.set(data, dataP / 8);
    try {
      checkStatus(Module, createImage(OUT_PATH, naxes.length, naxesP), 'writeImageDoubleWithHeader (create)');
      for (const [key, value] of Object.entries(headerEntries)) {
        const rc = typeof value === 'string'
          ? writeKeyStr(key, value, '')
          : writeKeyDbl(key, value, '');
        checkStatus(Module, rc, `writeImageDoubleWithHeader (key ${key})`);
      }
      checkStatus(Module, writeImageData(dataP, data.length), 'writeImageDoubleWithHeader (data)');
      checkStatus(Module, closeImage(), 'writeImageDoubleWithHeader (close)');
      return Module.FS.readFile(OUT_PATH);
    } finally {
      [naxesP, dataP].forEach((p) => Module._free(p));
    }
  });
}

const cfitsioWasmApi = {
  readImageInfo, readImageDouble, readKeyDouble, readKeyString,
  writeImageDouble, writeImageDoubleWithHeader,
};
if (typeof module !== 'undefined' && module.exports) {
  module.exports = cfitsioWasmApi;
} else if (typeof globalThis !== 'undefined') {
  globalThis.CfitsioWasm = cfitsioWasmApi;
}
