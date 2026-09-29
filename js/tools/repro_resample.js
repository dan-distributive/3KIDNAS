'use strict';
// One-off standalone repro (2026-09-29, Dan) -- kept as a template for the
// next time this exact question comes up, NOT wired into any test runner.
// Hardcoded to WALLABY_J100336-262923's realization 10 of an 11-bootstrap,
// cloudDensity=100 run (BootstrapSeed=42 -> idum=-(42+10+1)=-53) -- edit
// cubePath/modelPath/bsCent/idum for a different galaxy/realization.
//
// Runs JS's genFlipBootstrapSample directly (bypassing the whole fit/
// bootstrap-loop machinery) with the exact geometry Fortran's own
// MakeBootstrapSample.WriteBootstrapFile would have written for this
// realization (see /tmp/resample_repro/bs_repro.txt, built by hand from
// this run's own IniEstimate/AvgModel/RawGeom_v1.txt output), to test
// whether the resampled cube ITSELF already diverges from Fortran's
// Programs/BootStrapSampler output for this realization, independent of
// the amoeba optimizer entirely.
//
// RESULT (2026-09-29): proved bit-exact. A sorted-value comparison against
// Fortran's own ReproBS10.fits (via Programs/BootStrapSampler on the same
// bs_repro.txt) showed zero differences -- the two arrays contain the
// identical multiset of flux values (a direct index-aligned diff differs
// only because this script's flat array and astropy's FITS-order read use
// different axis orderings, not because any value actually differs). This
// ruled out bootstrap resampling as the source of a real, larger
// (non-ULP) divergence isolated to this one realization -- see
// repro_sofia.js for the next stage this pointed to (SoFiA wasm vs
// native), and JS_FORTRAN_PARITY_SESSION_2026-09-28.md for the full story.
const fs = require('fs');

const { DataCube, allocateDataCube } = require('../src/ObjectDefinitions/DataCube');
const { dataCubeToFitsBytes, fitsBytesToDataCube } = require('../src/BootstrapSampler/DataCubeFits');
const { genFlipBootstrapSample } = require('../src/BootstrapSampler/FlipBootstrap');
const { makeRng } = require('../src/StandardMath/random');

async function main() {
  const cfitsio = require('../../third_party/cfitsio-4.6.3/wasm/cfitsio-wasm.js');
  await cfitsio.ready;

  const cubePath = '/Users/dandesjardins/DCP/3KIDNAS/3KIDNASTests/TestData/WALLABY_Test_sources/WALLABY_J100336-262923/WALLABY_J100336-262923_VelCube.fits';
  const modelPath = '/Users/dandesjardins/DCP/3KIDNAS/3KIDNASTests/SingleGalaxyTest/TestFits_RunAllThree_FortranLocal/WALLABY_J100336-262923/WALLABY_J100336-262923_AverageModel_v1.fits';

  // Header fields, matching DataCubeInput.f's exact FITS-keyword mapping
  // (ChannelSize=CDELT3, RefLocation=CRPIX-1 (0-indexed), RefVal=CRVAL).
  const hdrObserved = new DataCube();
  hdrObserved.dh.nPixels[0] = 43; hdrObserved.dh.nPixels[1] = 44; hdrObserved.dh.nChannels = 54;
  hdrObserved.dh.pixelSize[0] = -0.00166666666667; hdrObserved.dh.pixelSize[1] = 0.00166666666667;
  hdrObserved.dh.channelSize = -3931.7167548543075;
  hdrObserved.dh.refLocation[0] = -1611.0 - 1; hdrObserved.dh.refLocation[1] = -488.0 - 1; hdrObserved.dh.refLocation[2] = 27 - 1;
  hdrObserved.dh.refVal[0] = 153.949308333; hdrObserved.dh.refVal[1] = -27.3743805556; hdrObserved.dh.refVal[2] = 889520.5453497537;
  hdrObserved.dh.uncertainty = 0;
  allocateDataCube(hdrObserved);

  const observedDC = await fitsBytesToDataCube(cfitsio, fs.readFileSync(cubePath), hdrObserved);
  const modelDC = await fitsBytesToDataCube(cfitsio, fs.readFileSync(modelPath), hdrObserved);

  const bsCent = {
    centX: 22.432947158813477,
    centY: 22.443210601806641,
    centV: 19.909175662526987,
    pa: 3.11680603,
    inc: 0.534571469,
  };

  const rng = makeRng(-53);
  const bootstrapCube = genFlipBootstrapSample(observedDC, modelDC, bsCent, 1.0, rng);

  const stats = { sum: 0, min: Infinity, max: -Infinity, n: bootstrapCube.flux.length };
  for (let i = 0; i < bootstrapCube.flux.length; i++) {
    const v = bootstrapCube.flux[i];
    stats.sum += v;
    if (v < stats.min) stats.min = v;
    if (v > stats.max) stats.max = v;
  }
  console.log('JS resample stats', JSON.stringify(stats));

  fs.writeFileSync('/tmp/resample_repro/js_flux.f32', Buffer.from(bootstrapCube.flux.buffer));
  console.log('wrote /tmp/resample_repro/js_flux.f32');
}

main().catch((e) => { console.error(e); process.exit(1); });
