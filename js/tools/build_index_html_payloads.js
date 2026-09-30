'use strict';
// One-off: build the SAME initial-fit payload index.html would build for
// WALLABY_J100336-262923, cloudDensity=100, matching run_both.js's earlier
// seed=42/cloudDensity=100 test exactly, using the SAME isomorphic
// buildFitPayloads.js index.html itself uses (not a hand-rolled payload).
const fs = require('fs');
const path = require('path');
const { defaultFittingOptions } = require('../src/PipelineConfig/defaultFittingOptions.js');
const { buildInitialFitPayload } = require('../src/PayloadBuilder/buildFitPayloads.js');

const TEST_ROOT = path.join(__dirname, '..', '..', '3KIDNASTests', 'TestData', 'WALLABY_Test_sources', 'WALLABY_J100336-262923');
const cubeBytes = new Uint8Array(fs.readFileSync(path.join(TEST_ROOT, 'WALLABY_J100336-262923_VelCube.fits')));
const maskBytes = new Uint8Array(fs.readFileSync(path.join(TEST_ROOT, 'SoFiA_J100336-262923_mask.fits')));

const options = defaultFittingOptions();
options.cloudBaseSurfDens = 100.0; // matches run_both.js --cloudDensity 100

const payload = buildInitialFitPayload({
  cubeBytes, maskBytes,
  paEstDeg: 81.2713489724864,
  incEstDeg: 31.490766615048877,
  options,
});

fs.writeFileSync(path.join(__dirname, 'index_html_initial_fit_payload.json'), JSON.stringify(payload));
console.log('wrote index_html_initial_fit_payload.json, cloudBaseSurfDens=', options.cloudBaseSurfDens);
