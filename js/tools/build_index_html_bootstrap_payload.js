'use strict';
const fs = require('fs');
const path = require('path');
const { buildBootstrapPayload } = require('../src/PayloadBuilder/buildFitPayloads.js');
const { defaultFittingOptions } = require('../src/PipelineConfig/defaultFittingOptions.js');

const TEST_ROOT = path.join(__dirname, '..', '..', '3KIDNASTests', 'TestData', 'WALLABY_Test_sources', 'WALLABY_J100336-262923');
const cubeBytes = new Uint8Array(fs.readFileSync(path.join(TEST_ROOT, 'WALLABY_J100336-262923_VelCube.fits')));
const sofiaParTemplateText = fs.readFileSync(path.join(__dirname, '..', 'sofia-template-par-file.par'), 'utf8');
const initialFitResult = JSON.parse(fs.readFileSync(path.join(__dirname, 'index_html_initial_fit_result5.json'), 'utf8'));

const options = defaultFittingOptions();
options.cloudBaseSurfDens = 100.0;

const payload = buildBootstrapPayload({
  initialFitResult, cubeBytes, sofiaParTemplateText,
  bootstrapSeed: 42, options,
});

fs.writeFileSync(path.join(__dirname, 'index_html_bootstrap_payload.json'), JSON.stringify(payload));
console.log('wrote index_html_bootstrap_payload.json');
