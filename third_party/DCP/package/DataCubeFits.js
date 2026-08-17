module.declare(["./DataCube.js"], function (require, exports, module) {
'use strict';

// =============================================================================
// DataCubeFits.js
// Converts a DataCube (ObjectDefinitions/DataCube.js) to real FITS bytes via
// cfitsio-wasm, so it can be handed to SoFiA (which needs an actual FITS
// file, not a JS array) -- and back: fitsBytesToMaskDataCube reads SoFiA's
// returned mask FITS bytes directly into a DataCube inside the JS worker.
//
// CORRECTION: an earlier version of this file's header claimed no
// FITS->DataCube direction was needed, since the mask used to be handed to
// local Python (astropy.io.fits) for decoding. That local step is exactly
// what the merged per-realization DCP launcher removes -- so this direction
// is needed now. cfitsio-wasm already exposes readImageDouble (see its own
// header), unused until now.
//
// UNIT CONVENTIONS -- derived from DataCube.js's own coordinate formula,
// not guessed, and cross-checked against a real WALLABY cube's actual FITS
// header (WALLABY_J103554-475245_cube.fits) plus a real diskfit_fixture.json
// on disk in this repo:
//
//   allocateDataCube's per-pixel formula (DataCube.js:114-123) is
//     world(i) = refVal[j] + (i - refLocation[j]) * pixelSize[j]     (i: 0-based)
//   The standard FITS WCS formula is
//     world(p) = CRVAL + (p - CRPIX) * CDELT                        (p: 1-based)
//   Substituting p = i+1 and matching terms:
//     CRPIX = refLocation + 1,  CRVAL = refVal,  CDELT = pixelSize
//
//   pixelSize/refVal for the two spatial axes are in ARCSEC, not degrees --
//   confirmed empirically: a real fixture's refValX (581523.125 arcsec)
//   matches that same cube's real CRVAL1 (161.534191667 deg) to within
//   f32 rounding (161.534191667*3600 = 581522.99). FITS wants CRVAL/CDELT
//   in degrees for RA---SIN/DEC--SIN, hence the /3600 below.
//
//   channelSize/refVal[2] are already km/s (confirmed: refValV/startV in a
//   real fixture are ~5700-5900, the right order of magnitude for a WALLABY
//   source's LSR velocity in km/s, not Hz) -- written as a VELO-LSR axis in
//   km/s directly, no conversion, rather than forcing a FREQ-axis
//   conversion this pipeline has no rest-frequency input for.
//
//   Beam major/minor axis (Beam2D.beamMajorAxis/beamMinorAxis) are in
//   PIXELS (confirmed: used directly as a pixel count in
//   GenerateBootstrap.js's block-resampling size calculation) -- converted
//   to degrees for BMAJ/BMIN via pixelSize (arcsec/pixel) / 3600.
// =============================================================================

const ARCSEC_PER_DEG = 3600;

/**
 * Serializes a DataCube (dh + flux) to real FITS bytes via cfitsio-wasm,
 * with a WCS + beam header sufficient for SoFiA to run against.
 *
 * @param {Object} cfitsio - an already-loaded cfitsio-wasm module (or the
 *   published cfitsio4wasm package's require('cfitsio-wasm.js')) -- passed
 *   in rather than required here, so this file doesn't hardcode whether
 *   the caller is using the local dev build or the published package.
 * @param {import('../ObjectDefinitions/DataCube.js').DataCube} dataCube
 * @param {import('../ObjectDefinitions/Beam.js').Beam2D} beam
 * @returns {Promise<Uint8Array>}
 */
async function dataCubeToFitsBytes(cfitsio, dataCube, beam) {
  const dh = dataCube.dh;
  const naxes = [dh.nPixels[0], dh.nPixels[1], dh.nChannels];

  const header = {
    CRPIX1: dh.refLocation[0] + 1, CRVAL1: dh.refVal[0] / ARCSEC_PER_DEG,
    CDELT1: dh.pixelSize[0] / ARCSEC_PER_DEG, CTYPE1: 'RA---SIN', CUNIT1: 'deg',

    CRPIX2: dh.refLocation[1] + 1, CRVAL2: dh.refVal[1] / ARCSEC_PER_DEG,
    CDELT2: dh.pixelSize[1] / ARCSEC_PER_DEG, CTYPE2: 'DEC--SIN', CUNIT2: 'deg',

    CRPIX3: dh.refLocation[2] + 1, CRVAL3: dh.refVal[2],
    CDELT3: dh.channelSize, CTYPE3: 'VELO-LSR', CUNIT3: 'km/s',

    BMAJ: (beam.beamMajorAxis * Math.abs(dh.pixelSize[0])) / ARCSEC_PER_DEG,
    BMIN: (beam.beamMinorAxis * Math.abs(dh.pixelSize[1])) / ARCSEC_PER_DEG,
    BPA: beam.beamPositionAngle || 0,
    BUNIT: 'Jy/beam',
  };

  const data = flattenToFitsOrder(dataCube.flux, dh.nPixels[0], dh.nPixels[1], dh.nChannels);

  return cfitsio.writeImageDoubleWithHeader(naxes, data, header);
}

// FITS storage order varies axis[0] (NAXIS1=x) fastest: m = i + j*nPixX +
// k*nPixX*nPixY. DataCube.flux varies channel fastest instead (flatIndxCalc,
// DataCube.js:162-164): l = k + j*nChan + i*nChan*nPixY. These are genuinely
// different orderings, not the same order under a different name -- checked
// numerically (a 3x2x4 test cube has 22/24 cells at different flat offsets
// between the two), not just asserted. This is the inverse of
// make_bootstrap_delta.py's `cube.transpose(2,1,0)` (FITS/astropy order ->
// DataCube order): same transpose, opposite direction, done directly on the
// flat array here instead of via a real ndarray transpose.
function flattenToFitsOrder(flux, nPixX, nPixY, nChan) {
  const out = new Array(nPixX * nPixY * nChan);
  for (let i = 0; i < nPixX; i++) {
    for (let j = 0; j < nPixY; j++) {
      for (let k = 0; k < nChan; k++) {
        out[i + j * nPixX + k * nPixX * nPixY] = flux[k + j * nChan + i * nChan * nPixY];
      }
    }
  }
  return out;
}

// Inverse of flattenToFitsOrder: FITS order (x fastest) -> DataCube's
// flatIndxCalc order (channel fastest).
function unflattenFromFitsOrder(data, nPixX, nPixY, nChan) {
  const out = new Float32Array(nPixX * nPixY * nChan);
  for (let i = 0; i < nPixX; i++) {
    for (let j = 0; j < nPixY; j++) {
      for (let k = 0; k < nChan; k++) {
        out[k + j * nChan + i * nChan * nPixY] = Math.fround(data[i + j * nPixX + k * nPixX * nPixY]);
      }
    }
  }
  return out;
}

/**
 * Reads FITS bytes (e.g. SoFiA's returned mask) into a DataCube, reusing
 * refDataCube's header (pixel size, channel size, ref location/value) --
 * valid because the mask shares the exact same pixel grid as the cube it
 * was produced from, so there's no need to reconstruct WCS from the FITS
 * keywords themselves.
 *
 * @param {Object} cfitsio - an already-loaded cfitsio-wasm module.
 * @param {Uint8Array} fitsBytes
 * @param {import('../ObjectDefinitions/DataCube.js').DataCube} refDataCube
 * @returns {Promise<import('../ObjectDefinitions/DataCube.js').DataCube>}
 */
async function fitsBytesToDataCube(cfitsio, fitsBytes, refDataCube) {
  const { DataCube, allocateDataCube } = require('./DataCube.js');
  const { naxes, data } = await cfitsio.readImageDouble(fitsBytes);
  const [nPixX, nPixY, nChan] = naxes;

  const dc = new DataCube();
  const dh = dc.dh;
  const refDh = refDataCube.dh;
  dh.nPixels[0] = nPixX;
  dh.nPixels[1] = nPixY;
  dh.nChannels  = nChan;
  dh.pixelSize[0]   = refDh.pixelSize[0];
  dh.pixelSize[1]   = refDh.pixelSize[1];
  dh.channelSize    = refDh.channelSize;
  dh.refLocation[0] = refDh.refLocation[0];
  dh.refLocation[1] = refDh.refLocation[1];
  dh.refLocation[2] = refDh.refLocation[2];
  dh.refVal[0]      = refDh.refVal[0];
  dh.refVal[1]      = refDh.refVal[1];
  dh.refVal[2]      = refDh.refVal[2];
  dh.uncertainty    = refDh.uncertainty;

  allocateDataCube(dc);
  dh.start[0] = refDh.start[0];
  dh.start[1] = refDh.start[1];
  dh.start[2] = refDh.start[2];

  dc.flux = unflattenFromFitsOrder(data, nPixX, nPixY, nChan);
  return dc;
}

module.exports = { dataCubeToFitsBytes, fitsBytesToDataCube };


// ---------------------------------------------------------------------------
// Self-test (node DataCubeFits.js) -- real round trip via the local cfitsio
// wasm build (third_party/cfitsio-4.6.3/wasm/cfitsio-wasm.js), not a mock.
// ---------------------------------------------------------------------------
if (require.main === module) {
  (async () => {
    const cfitsio = require('../../../cfitsio-4.6.3/wasm/cfitsio-wasm.js');
    const { DataCube, allocateDataCube, flatIndxCalc } = require('./DataCube.js');
    const { Beam2D } = require('./Beam.js');
    const f32 = Math.fround;

    const dc = new DataCube();
    const dh = dc.dh;
    dh.nPixels[0] = 5; dh.nPixels[1] = 4; dh.nChannels = 3;
    dh.pixelSize[0] = f32(6.0); dh.pixelSize[1] = f32(6.0);
    dh.channelSize = f32(4.0);
    dh.refLocation[0] = f32(2.0); dh.refLocation[1] = f32(1.5); dh.refLocation[2] = f32(1.0);
    dh.refVal[0] = f32(150.0); dh.refVal[1] = f32(30.0); dh.refVal[2] = f32(1000.0);
    allocateDataCube(dc);
    // Distinct value per cell so any axis-order mistake shows up clearly.
    for (let i = 0; i < 5; i++)
      for (let j = 0; j < 4; j++)
        for (let k = 0; k < 3; k++)
          dc.flux[flatIndxCalc(i, j, k, dh)] = f32(i * 100 + j * 10 + k);

    const beam = new Beam2D();
    beam.beamMajorAxis = f32(3.0); beam.beamMinorAxis = f32(2.0);

    console.log('=== DataCubeFits round trip (real cfitsio wasm) ===');
    const bytes = await dataCubeToFitsBytes(cfitsio, dc, beam);
    console.log('  wrote', bytes.length, 'bytes');

    const back = await fitsBytesToDataCube(cfitsio, bytes, dc);
    console.log('  read back nPixels:', Array.from(back.dh.nPixels), 'nChannels:', back.dh.nChannels);

    let allMatch = true;
    for (let i = 0; i < 5; i++)
      for (let j = 0; j < 4; j++)
        for (let k = 0; k < 3; k++) {
          const expected = f32(i * 100 + j * 10 + k);
          const actual = back.flux[flatIndxCalc(i, j, k, back.dh)];
          if (actual !== expected) {
            allMatch = false;
            console.log(`  MISMATCH at (${i},${j},${k}): expected ${expected}, got ${actual}`);
          }
        }
    console.log('  all', 5 * 4 * 3, 'cells match:', allMatch ? 'OK' : 'FAIL');
  })();
}

});
