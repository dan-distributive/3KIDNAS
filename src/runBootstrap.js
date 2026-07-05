'use strict';

// =============================================================================
// runBootstrap.js
// End-to-end bootstrap pipeline test for the 3KIDNAS JS port.
//
// Loads diskfit_fixture.json (observed cube + fitting state) and
// model_cube_bestfit.json (Fortran best-fit model cube), generates N
// bootstrap realizations, runs GalaxyFit_Simple on each, and reports
// the resulting parameter distributions.
//
// Usage:
//   node runBootstrap.js <fixture.json> <bestfit_model.json> [nRealizations]
//
// Example:
//   node runBootstrap.js diskfit_fixture.json model_cube_bestfit.json 10
// =============================================================================

const f32 = Math.fround;
const fs  = require('fs');

const { makeRng }            = require('./StandardMath/random.js');
const { DataCube, allocateDataCube, flatIndxCalc }
                             = require('./ObjectDefinitions/DataCube.js');
const { Beam2D }             = require('./ObjectDefinitions/Beam.js');
const { TiltedRingModel, TiltedRingFittingOptions,
        tiltRing_Allocate, tiltRingFittingOptions_Allocate,
        logicalTiltedRingIndexing }
                             = require('./ObjectDefinitions/TiltedRing.js');
const { ParameterVector, allocateParamVector }
                             = require('./ObjectDefinitions/ParameterVector.js');
const { calculate2DBeamKernel }
                             = require('./ConvolveCube/CalculateBeamKernel.js');
const { generalizedParamVectorToTiltedRing }
                             = require('./ParameterToTiltedRingInterface/ParameterToTiltedRingVector.js');
const { galaxyFit_Simple }   = require('./GalaxyAnalysis/GalaxyFit.js');
const { genBootstrapSample } = require('./BootstrapSampler/GenerateBootstrap.js');


// ---------------------------------------------------------------------------
// loadFixture — reuse from runFixture.js logic
// ---------------------------------------------------------------------------
function loadFixture(fixturePath) {
  const d = JSON.parse(fs.readFileSync(fixturePath, 'utf8'));

  const pvIni = new ParameterVector();
  pvIni.nParams = d.nParams;
  allocateParamVector(pvIni);
  for (let i = 0; i < d.nParams; i++) {
    pvIni.param[i]          = f32(d.param[i]);
    pvIni.paramLowerLims[i] = f32(d.paramLowerLims[i]);
    pvIni.paramUpperLims[i] = f32(d.paramUpperLims[i]);
    pvIni.paramRange[i]     = f32(d.paramRange[i]);
    pvIni.cyclicSwitch[i]   = d.cyclicSwitch[i];
  }

  const trfo = new TiltedRingFittingOptions();
  trfo.nRings = d.nRings;
  trfo.nRingsPerBeam = 1;
  trfo.nTargRings = -1;
  for (let i = 0; i <= 12; i++) {
    trfo.constParams[i]    = d.constParams[i];
    trfo.fixedParams[i]    = d.fixedParams[i];
    trfo.paramLowerLims[i] = f32(d.paramLowerLims13[i]);
    trfo.paramUpperLims[i] = f32(d.paramUpperLims13[i]);
    trfo.paramRange[i]     = f32(d.paramRange13[i]);
    trfo.cyclicSwitch[i]   = d.cyclicSwitch13[i];
  }
  tiltRingFittingOptions_Allocate(trfo);
  logicalTiltedRingIndexing(trfo);
  for (let r = 0; r < d.nRings; r++) {
    const rp = trfo.radialProfiles[r];
    const dp = d.radialProfiles[r];
    rp.rmid = f32(dp.rmid); rp.rwidth = f32(dp.rwidth);
    rp.centPos[0] = f32(dp.centPos0); rp.centPos[1] = f32(dp.centPos1);
    rp.inclination = f32(dp.inclination); rp.positionAngle = f32(dp.positionAngle);
    rp.vSys = f32(dp.vSys); rp.vRot = f32(dp.vRot); rp.vRad = f32(dp.vRad);
    rp.vDisp = f32(dp.vDisp); rp.vvert = f32(dp.vvert); rp.dvdz = f32(dp.dvdz);
    rp.sigUse = f32(dp.sigUse); rp.z0 = f32(dp.z0);
    rp.zGradiantStart = f32(dp.zGradiantStart);
  }

  // Observed cube
  const obs = new DataCube();
  const dh  = obs.dh;
  dh.nPixels[0] = d.observedDC.nPixelsX;
  dh.nPixels[1] = d.observedDC.nPixelsY;
  dh.nChannels  = d.observedDC.nChannels;
  dh.pixelSize[0] = f32(d.observedDC.pixelSizeX);
  dh.pixelSize[1] = f32(d.observedDC.pixelSizeY);
  dh.channelSize  = f32(d.observedDC.channelSize);
  dh.refLocation[0] = f32(d.observedDC.refLocX);
  dh.refLocation[1] = f32(d.observedDC.refLocY);
  dh.refLocation[2] = f32(d.observedDC.refLocV);
  dh.refVal[0] = f32(d.observedDC.refValX);
  dh.refVal[1] = f32(d.observedDC.refValY);
  dh.refVal[2] = f32(d.observedDC.refValV);
  dh.uncertainty = f32(d.observedDC.uncertainty);
  dh.nValid = d.observedDC.nValid;
  const nCells = dh.nPixels[0] * dh.nPixels[1] * dh.nChannels;
  obs.pixels  = new Float32Array(2 * Math.max(dh.nPixels[0], dh.nPixels[1]));
  obs.channels = new Float32Array(dh.nChannels);
  obs.flux     = new Float32Array(nCells);
  obs.flattendValidIndices = new Int32Array(dh.nValid);
  for (let l = 0; l < dh.nValid; l++) {
    obs.flattendValidIndices[l] = d.observedDC.flattendValidIndices[l];
  }
  for (let l = 0; l < dh.nValid; l++) {
    obs.flux[obs.flattendValidIndices[l]] = f32(d.observedDC.validFlux[l]);
  }
  obs.dh.start[0] = f32(d.observedDC.startX);
  obs.dh.start[1] = f32(d.observedDC.startY);
  obs.dh.start[2] = f32(d.observedDC.startV);

  // Beam
  const beam = new Beam2D();
  beam.beamSigmaVector[0] = f32(d.observedBeam.beamSigma0);
  beam.beamSigmaVector[1] = f32(d.observedBeam.beamSigma1);
  beam.beamSigmaVector[2] = f32(d.observedBeam.beamSigma2);
  beam.nRadialCells  = d.observedBeam.nRadialCells;
  beam.sigmaLengths  = f32(d.observedBeam.sigmaLengths);
  beam.pixelSize[0]  = f32(d.observedBeam.pixelSizeX);
  beam.pixelSize[1]  = f32(d.observedBeam.pixelSizeY);
  beam.beamMajorAxis = f32(d.observedBeam.beamMajorAxis);
  beam.beamMinorAxis = f32(d.observedBeam.beamMinorAxis);
  const n = beam.nRadialCells;
  const kSz = 2 * n + 1;
  beam.kernel = new Float32Array(kSz * kSz);
  beam.paddedSize[0] = 2*n+1 + dh.nPixels[0];
  beam.paddedSize[1] = 2*n+1 + dh.nPixels[1];
  beam.complexSize[0] = Math.trunc(beam.paddedSize[0]/2)+1;
  beam.complexSize[1] = beam.paddedSize[1];
  beam.complexKernelCreated = false;
  calculate2DBeamKernel(beam, beam.pixelSize);

  return { pvIni, trfo, obs, beam,
           ftol: f32(d.ftol), idum: d.idum,
           linearLogSDSwitch: d.linearLogSDSwitch };
}


// ---------------------------------------------------------------------------
// loadModelCube — load Fortran best-fit model cube JSON
// ---------------------------------------------------------------------------
function loadModelCube(modelPath, obs) {
  const d = JSON.parse(fs.readFileSync(modelPath, 'utf8'));

  const mdl = new DataCube();
  mdl.dh.nPixels[0] = d.nPixelsX;
  mdl.dh.nPixels[1] = d.nPixelsY;
  mdl.dh.nChannels  = d.nChannels;
  mdl.dh.pixelSize[0] = obs.dh.pixelSize[0];
  mdl.dh.pixelSize[1] = obs.dh.pixelSize[1];
  mdl.dh.channelSize  = obs.dh.channelSize;
  mdl.dh.refLocation[0] = obs.dh.refLocation[0];
  mdl.dh.refLocation[1] = obs.dh.refLocation[1];
  mdl.dh.refLocation[2] = obs.dh.refLocation[2];
  mdl.dh.refVal[0] = obs.dh.refVal[0];
  mdl.dh.refVal[1] = obs.dh.refVal[1];
  mdl.dh.refVal[2] = obs.dh.refVal[2];
  mdl.dh.uncertainty = obs.dh.uncertainty;
  mdl.dh.nValid = obs.dh.nValid;
  allocateDataCube(mdl);
  mdl.dh.start[0] = obs.dh.start[0];
  mdl.dh.start[1] = obs.dh.start[1];
  mdl.dh.start[2] = obs.dh.start[2];
  mdl.flattendValidIndices.set(obs.flattendValidIndices);

  // Load flux in Fortran loop order (i,j,k)
  let idx = 0;
  for (let i = 0; i < d.nPixelsX; i++)
    for (let j = 0; j < d.nPixelsY; j++)
      for (let k = 0; k < d.nChannels; k++)
        mdl.flux[flatIndxCalc(i, j, k, mdl.dh)] = f32(d.flux[idx++]);

  return mdl;
}


// ---------------------------------------------------------------------------
// buildFitState — build the state object for galaxyFit_Simple
// ---------------------------------------------------------------------------
function buildFitState(pvIni, trfo, obs, beam, ftol, idum,
                       linearLogSDSwitch, bootstrapCube) {
  const pvModel    = new ParameterVector();
  const pvFirstFit = new ParameterVector();
  pvModel.nParams = pvIni.nParams;
  allocateParamVector(pvModel);
  pvModel.param.set(pvIni.param);
  pvModel.paramLowerLims.set(pvIni.paramLowerLims);
  pvModel.paramUpperLims.set(pvIni.paramUpperLims);
  pvModel.cyclicSwitch.set(pvIni.cyclicSwitch);
  pvModel.paramRange.set(pvIni.paramRange);

  const modelTR = new TiltedRingModel();
  modelTR.nRings = trfo.nRings;
  modelTR.cmode  = 0;
  modelTR.cloudBaseSurfDens = f32(100.0);
  tiltRing_Allocate(modelTR);
  for (let r = 0; r < trfo.nRings; r++) {
    modelTR.r[r].rmid   = f32(trfo.radialProfiles[r].rmid);
    modelTR.r[r].rwidth = f32(trfo.radialProfiles[r].rwidth);
  }

  // Use bootstrap cube as observed cube for fitting
  const fitObs = bootstrapCube !== null ? bootstrapCube : obs;

  const modelDC = new DataCube();
  modelDC.dh.nPixels[0]    = fitObs.dh.nPixels[0];
  modelDC.dh.nPixels[1]    = fitObs.dh.nPixels[1];
  modelDC.dh.nChannels     = fitObs.dh.nChannels;
  modelDC.dh.pixelSize[0]  = fitObs.dh.pixelSize[0];
  modelDC.dh.pixelSize[1]  = fitObs.dh.pixelSize[1];
  modelDC.dh.channelSize   = fitObs.dh.channelSize;
  modelDC.dh.refLocation[0]= fitObs.dh.refLocation[0];
  modelDC.dh.refLocation[1]= fitObs.dh.refLocation[1];
  modelDC.dh.refLocation[2]= fitObs.dh.refLocation[2];
  modelDC.dh.refVal[0]     = fitObs.dh.refVal[0];
  modelDC.dh.refVal[1]     = fitObs.dh.refVal[1];
  modelDC.dh.refVal[2]     = fitObs.dh.refVal[2];
  modelDC.dh.uncertainty   = fitObs.dh.uncertainty;
  modelDC.dh.nValid        = fitObs.dh.nValid;
  allocateDataCube(modelDC);
  modelDC.dh.start[0] = fitObs.dh.start[0];
  modelDC.dh.start[1] = fitObs.dh.start[1];
  modelDC.dh.start[2] = fitObs.dh.start[2];
  modelDC.flattendValidIndices.set(fitObs.flattendValidIndices);

  return {
    pvIni, pvModel, pvFirstFit,
    modelTiltedRing:  modelTR,
    modelDC,
    observedDC:       fitObs,
    observedBeam:     beam,
    trFittingOptions: trfo,
    rng:              makeRng(idum),
    linearLogSDSwitch,
    ftol,
    iniGuessWidth:    f32(1.0),
    paramToTiltedRing: generalizedParamVectorToTiltedRing
  };
}


// ---------------------------------------------------------------------------
// collate — compute per-parameter statistics from N best-fit PV arrays
// ---------------------------------------------------------------------------
function collate(results, nParams) {
  console.log('\n=== Bootstrap Results ===');
  console.log(`N realizations: ${results.length}`);
  console.log(`\n${'Param'.padEnd(6)} ${'Mean'.padEnd(12)} ${'Std'.padEnd(12)} ${'p16'.padEnd(12)} ${'p50'.padEnd(12)} ${'p84'.padEnd(12)}`);
  console.log('-'.repeat(66));

  const stats = [];
  for (let p = 0; p < nParams; p++) {
    const vals = results.map(r => r.params[p]).sort((a,b) => a-b);
    const mean = vals.reduce((a,v) => a+v, 0) / vals.length;
    const std  = Math.sqrt(vals.reduce((a,v) => a+(v-mean)**2, 0) / vals.length);
    const p16  = vals[Math.floor(0.16 * vals.length)];
    const p50  = vals[Math.floor(0.50 * vals.length)];
    const p84  = vals[Math.floor(0.84 * vals.length)];
    stats.push({ mean, std, p16, p50, p84 });
    console.log(`[${String(p).padEnd(3)}]  ${mean.toFixed(4).padEnd(12)} ${std.toFixed(4).padEnd(12)} ${p16.toFixed(4).padEnd(12)} ${p50.toFixed(4).padEnd(12)} ${p84.toFixed(4).padEnd(12)}`);
  }

  const chi2vals = results.map(r => r.chi2).sort((a,b) => a-b);
  const chi2mean = chi2vals.reduce((a,v) => a+v, 0) / chi2vals.length;
  console.log(`\nchi² mean: ${chi2mean.toExponential(4)}`);
  console.log(`chi² range: [${chi2vals[0].toExponential(4)}, ${chi2vals[chi2vals.length-1].toExponential(4)}]`);

  return stats;
}


// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
async function main() {
  const fixturePath  = process.argv[2];
  const modelPath    = process.argv[3];
  const nReal        = parseInt(process.argv[4] || '5');

  if (!fixturePath || !modelPath) {
    console.error('Usage: node runBootstrap.js <fixture.json> <bestfit_model.json> [nRealizations]');
    process.exit(1);
  }

  console.log(`Loading fixture: ${fixturePath}`);
  const { pvIni, trfo, obs, beam, ftol, idum, linearLogSDSwitch }
    = loadFixture(fixturePath);

  console.log(`Loading best-fit model cube: ${modelPath}`);
  const modelCube = loadModelCube(modelPath, obs);

  console.log(`\nGalaxy: ${obs.dh.nPixels[0]}x${obs.dh.nPixels[1]}x${obs.dh.nChannels}`);
  console.log(`nParams: ${pvIni.nParams}, nRings: ${trfo.nRings}`);
  console.log(`beamMajorAxis: ${beam.beamMajorAxis.toFixed(3)} pixels`);
  console.log(`Observed flux sum: ${obs.flux.reduce((a,v)=>a+v,0).toFixed(6)}`);
  console.log(`Model flux sum:    ${modelCube.flux.reduce((a,v)=>a+v,0).toFixed(6)}`);

  // Bootstrap geometry — use best-fit centre from radial profiles ring 0
  const bsCent = {
    centX: f32(trfo.radialProfiles[0].centPos[0]),
    centY: f32(trfo.radialProfiles[0].centPos[1]),
    centV: f32(trfo.radialProfiles[0].vSys),
    pa:    f32(trfo.radialProfiles[0].positionAngle),
    inc:   f32(trfo.radialProfiles[0].inclination),
  };
  const spatialBlockSize = f32(1.0);
  const velBlockSize     = f32(1.0);

  console.log(`\nBootstrap centre: (${bsCent.centX.toFixed(2)}, ${bsCent.centY.toFixed(2)})`);
  console.log(`PA: ${bsCent.pa.toFixed(4)}, Inc: ${bsCent.inc.toFixed(4)}`);
  console.log(`Block size: ${Math.round(spatialBlockSize * beam.beamMajorAxis)} x ${Math.round(spatialBlockSize * beam.beamMajorAxis)} x ${Math.round(velBlockSize)} channels`);
  console.log(`\nRunning ${nReal} bootstrap realizations...`);

  const results = [];
  const t0 = Date.now();

  for (let n = 0; n < nReal; n++) {
    const tN = Date.now();
    process.stdout.write(`  Realization ${n+1}/${nReal}... `);

    // Generate bootstrap cube
    const bootstrapCube = genBootstrapSample(
      obs, modelCube, beam, bsCent,
      spatialBlockSize, velBlockSize
    );

    // Run fit on bootstrap cube with unique seed
    const bsIdum = idum - (n + 1);
    const state = buildFitState(
      pvIni, trfo, obs, beam, ftol, bsIdum,
      linearLogSDSwitch, bootstrapCube
    );

    const pvBest = galaxyFit_Simple(state);
    const elapsed = ((Date.now() - tN) / 1000).toFixed(1);

    results.push({
      params: Array.from(pvBest.param),
      chi2:   pvBest.bestLike
    });

    console.log(`chi²=${pvBest.bestLike.toExponential(4)} (${elapsed}s)`);
  }

  const totalElapsed = ((Date.now() - t0) / 1000).toFixed(1);
  console.log(`\nTotal time: ${totalElapsed}s`);

  collate(results, pvIni.nParams);

  // Save results to JSON
  const outPath = 'bootstrap_results.json';
  fs.writeFileSync(outPath, JSON.stringify({
    nRealizations: nReal,
    nParams: pvIni.nParams,
    results
  }, null, 2));
  console.log(`\nResults saved to ${outPath}`);
}

main();
