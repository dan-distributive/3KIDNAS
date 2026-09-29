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

    // BUG FIX (2026-09-29, Dan): hardcoded CTYPE3='VELO-LSR'/CUNIT3='km/s'
    // regardless of what the real cube actually uses -- every WALLABY test
    // cube in this pipeline is actually CTYPE3='VOPT'/CUNIT3='m/s' (confirmed
    // directly against WALLABY_J100336-262923_VelCube.fits's real header).
    // dh.refVal[2]/dh.channelSize are this DataCube's own internal km/s
    // working units, so converting to m/s here (*1000) reproduces the same
    // physical values Fortran's own resampled-cube FITS output has. Found
    // via a header diff between this function's own output and Fortran's
    // real BootStrapSampler output for the same resample -- SoFiA reads
    // this file for its own source-finding, and a wrong spectral axis
    // convention can change which pixels it segments into the mask.
    CRPIX3: dh.refLocation[2] + 1, CRVAL3: dh.refVal[2] * 1000,
    CDELT3: dh.channelSize * 1000, CTYPE3: 'VOPT', CUNIT3: 'm/s',

    BMAJ: (beam.beamMajorAxis * Math.abs(dh.pixelSize[0])) / ARCSEC_PER_DEG,
    BMIN: (beam.beamMinorAxis * Math.abs(dh.pixelSize[1])) / ARCSEC_PER_DEG,
    // BUG FIX (2026-09-29, Dan): beam.beamPositionAngle is stored in RADIANS
    // (Beam.js's own convention -- fitBeam.beamPositionAngle is built as
    // bpaDeg*Pi/180 everywhere it's constructed), but the FITS BPA keyword
    // is degrees. Writing the raw radian value straight through was silently
    // wrong by a factor of ~180/pi -- masked for a long time because
    // runBootstrapRealization's resampleBeam never set beamPositionAngle at
    // all (stayed at Beam2D's default 0, and 0 rad happens to equal 0 deg),
    // so this only became visible once that separate missing-field bug
    // (same file, runBootstrapRealization) was fixed and started passing a
    // real nonzero radian value through. SoFiA reads this BPA to orient its
    // own segmentation, so a wrong value here gives a genuinely different
    // (not just ULP-off) source mask.
    BPA: (beam.beamPositionAngle || 0) * 180 / Math.PI,
    BUNIT: 'Jy/beam',
  };

  const data = flattenToFitsOrder(dataCube.flux, dh.nPixels[0], dh.nPixels[1], dh.nChannels);

  // BUG FIX (2026-09-29, Dan): wrote BITPIX=-64 (double precision) for
  // every FITS file this pipeline ever handed to SoFiA, while Fortran's
  // own resampled-cube output is BITPIX=-32 (real4/single precision) --
  // confirmed via a direct header diff. dataCube.flux is already a
  // Float32Array (every value already float32-rounded), so this was pure
  // precision INFLATION on the way out, not real extra precision -- but
  // SoFiA's own noise/threshold statistics could still read differently
  // off a file declaring double- vs single-precision pixels. Switched to
  // writeImageFloatWithHeader (BITPIX=-32) to match Fortran's real output
  // exactly.
  return cfitsio.writeImageFloatWithHeader(naxes, data, header);
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

  // BUG FIX (Dan, 2026): FITS blank/undefined pixels round-trip through
  // cfitsio's readImageDouble as real IEEE NaN, but this port previously
  // left them as NaN in dc.flux -- unlike Fortran's ReadFullDataCube
  // (DataCubeInput.f:317-382), which passes a finite nullval (-1010) to
  // ftgpve, explicitly detects it, zeroes the pixel, and builds
  // FlattendValidIndices/nValid to EXCLUDE those cells from every
  // downstream sum (critically, CubeComparison.f's chi2 likelihood, which
  // sums only over Cube1%FlattendValidIndices -- see
  // FullModelComparison.js's own header comment, which already assumed
  // this was "handled upstream" and was never actually done). Left
  // unfixed, a NaN pixel poisons every sum it touches (observed directly:
  // MaskCube's flux*=maskFlux still yields NaN*0=NaN, which propagated all
  // the way into EstimateShape's flux-weighted center, producing NaN
  // geometry and crashing pre-analysis entirely on a real, faint test
  // galaxy -- WALLABY_J100336-262923 -- that happens to have blanked edge
  // pixels; never caught before because earlier bit-exactness test
  // galaxies apparently didn't hit this path). Fixed here, once, for every
  // caller (both the observed cube and the mask go through this same
  // function) rather than patching each consumer.
  const flux = dc.flux;
  const validIndices = new Int32Array(flux.length);
  let nValid = 0;
  for (let idx = 0; idx < flux.length; idx++) {
    if (Number.isNaN(flux[idx])) {
      flux[idx] = 0;
    } else {
      validIndices[nValid++] = idx;
    }
  }
  dh.nValid = nValid;
  dc.flattendValidIndices = validIndices.slice(0, nValid);

  return dc;
}

module.exports = { dataCubeToFitsBytes, fitsBytesToDataCube };


// ---------------------------------------------------------------------------
// Self-test (node DataCubeFits.js) -- real round trip via the local cfitsio
// wasm build (third_party/cfitsio-4.6.3/wasm/cfitsio-wasm.js), not a mock.
// ---------------------------------------------------------------------------
if (require.main === module) {
  (async () => {
    const cfitsio = require('../../../third_party/cfitsio-4.6.3/wasm/cfitsio-wasm.js');
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
