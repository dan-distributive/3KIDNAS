'use strict';

// =============================================================================
// runFixture.js
// End-to-end integration test for the 3KIDNAS JS port.
// Loads the JSON fixture dumped by DumpFittingFixture() from Fortran,
// reconstructs all state objects, runs GalaxyFit_Simple, and reports
// the final chi² against the Fortran reference value.
//
// Usage:
//   node runFixture.js diskfit_fixture.json [fortran_chi2]
//
// Example:
//   node runFixture.js diskfit_fixture.json 1234.56
//
// If fortran_chi2 is supplied, the final chi² is compared against it.
// =============================================================================

const f32 = Math.fround;
const fs  = require('fs');
const path = require('path');

// Adjust these paths to match your folder structure
const { makeRng }            = require('./StandardMath/random.js');
const { DataCube, DataCubeHeader, allocateDataCube, flatIndxCalc, threeDIndxCalc }
                             = require('./ObjectDefinitions/DataCube.js');
const { Beam2D, allocate_Beam2D }
                             = require('./ObjectDefinitions/Beam.js');
const { TiltedRingModel, TiltedRingFittingOptions, Ring,
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
const { tiltedRingModelComparison }
                             = require('./CompareCubes/FullModelComparison.js');


// ---------------------------------------------------------------------------
// loadFixture
// Parses the JSON fixture and reconstructs all state objects.
// ---------------------------------------------------------------------------
function loadFixture(fixturePath) {
  console.log(`Loading fixture: ${fixturePath}`);
  const raw = fs.readFileSync(fixturePath, 'utf8');
  const d   = JSON.parse(raw);

  // ---- PVIni ----
  const pvIni    = new ParameterVector();
  pvIni.nParams  = d.nParams;
  allocateParamVector(pvIni);
  for (let i = 0; i < d.nParams; i++) {
    pvIni.param[i]          = f32(d.param[i]);
    pvIni.paramLowerLims[i] = f32(d.paramLowerLims[i]);
    pvIni.paramUpperLims[i] = f32(d.paramUpperLims[i]);
    pvIni.paramRange[i]     = f32(d.paramRange[i]);
    pvIni.cyclicSwitch[i]   = d.cyclicSwitch[i];
  }

  // ---- TiltedRingFittingOptions ----
  const trfo        = new TiltedRingFittingOptions();
  trfo.nRings       = d.nRings;
  trfo.nRingsPerBeam = 1; // not in fixture, safe default
  trfo.nTargRings   = -1;
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

  // Populate radial profiles
  for (let r = 0; r < d.nRings; r++) {
    const rp = trfo.radialProfiles[r];
    const dp = d.radialProfiles[r];
    rp.rmid           = f32(dp.rmid);
    rp.rwidth         = f32(dp.rwidth);
    rp.centPos[0]     = f32(dp.centPos0);
    rp.centPos[1]     = f32(dp.centPos1);
    rp.inclination    = f32(dp.inclination);
    rp.positionAngle  = f32(dp.positionAngle);
    rp.vSys           = f32(dp.vSys);
    rp.vRot           = f32(dp.vRot);
    rp.vRad           = f32(dp.vRad);
    rp.vDisp          = f32(dp.vDisp);
    rp.vvert          = f32(dp.vvert);
    rp.dvdz           = f32(dp.dvdz);
    rp.sigUse         = f32(dp.sigUse);
    rp.z0             = f32(dp.z0);
    rp.zGradiantStart = f32(dp.zGradiantStart);
  }

  // ---- ObservedDC ----
  const obs = new DataCube();
  const dh  = obs.dh;
  dh.nPixels[0]    = d.observedDC.nPixelsX;
  dh.nPixels[1]    = d.observedDC.nPixelsY;
  dh.nChannels     = d.observedDC.nChannels;
  dh.pixelSize[0]  = f32(d.observedDC.pixelSizeX);
  dh.pixelSize[1]  = f32(d.observedDC.pixelSizeY);
  dh.channelSize   = f32(d.observedDC.channelSize);
  dh.start[0]      = f32(d.observedDC.startX);
  dh.start[1]      = f32(d.observedDC.startY);
  dh.start[2]      = f32(d.observedDC.startV);
  dh.refLocation[0]= f32(d.observedDC.refLocX);
  dh.refLocation[1]= f32(d.observedDC.refLocY);
  dh.refLocation[2]= f32(d.observedDC.refLocV);
  dh.refVal[0]     = f32(d.observedDC.refValX);
  dh.refVal[1]     = f32(d.observedDC.refValY);
  dh.refVal[2]     = f32(d.observedDC.refValV);
  dh.uncertainty   = f32(d.observedDC.uncertainty);
  dh.nValid        = d.observedDC.nValid;

  // Allocate but override flux and valid indices from fixture
  const nCells     = dh.nPixels[0] * dh.nPixels[1] * dh.nChannels;
  obs.pixels       = new Float32Array(2 * Math.max(dh.nPixels[0], dh.nPixels[1]));
  obs.channels     = new Float32Array(dh.nChannels);
  obs.flux         = new Float32Array(nCells);
  obs.flattendValidIndices = new Int32Array(dh.nValid);

  // Populate valid indices
  for (let l = 0; l < dh.nValid; l++) {
    obs.flattendValidIndices[l] = d.observedDC.flattendValidIndices[l];
  }

  // Populate flux at valid voxels
  for (let l = 0; l < dh.nValid; l++) {
    const idx = obs.flattendValidIndices[l];
    obs.flux[idx] = f32(d.observedDC.validFlux[l]);
  }

  // Override start after allocation in case allocateDataCube was called
  obs.dh.start[0] = f32(d.observedDC.startX);
  obs.dh.start[1] = f32(d.observedDC.startY);
  obs.dh.start[2] = f32(d.observedDC.startV);

  // ---- ObservedBeam ----
  const beam = new Beam2D();
  beam.beamSigmaVector[0] = f32(d.observedBeam.beamSigma0);
  beam.beamSigmaVector[1] = f32(d.observedBeam.beamSigma1);
  beam.beamSigmaVector[2] = f32(d.observedBeam.beamSigma2);
  beam.nRadialCells        = d.observedBeam.nRadialCells;
  beam.sigmaLengths        = f32(d.observedBeam.sigmaLengths);
  beam.pixelSize[0]        = f32(d.observedBeam.pixelSizeX);
  beam.pixelSize[1]        = f32(d.observedBeam.pixelSizeY);

  // Allocate beam arrays without recomputing sigmas (already set from fixture)
  const n    = beam.nRadialCells;
  const kSz  = 2 * n + 1;
  beam.kernel = new Float32Array(kSz * kSz);
  for (let i = 0; i < kSz * kSz; i++) {
    beam.kernel[i] = f32(d.observedBeam.kernel[i]);
  }

  // Set paddedSize and complexSize matching Fortran Allocate_Beam2D
  beam.paddedSize[0]  = 2 * n + 1 + dh.nPixels[0];
  beam.paddedSize[1]  = 2 * n + 1 + dh.nPixels[1];
  beam.complexSize[0] = Math.trunc(beam.paddedSize[0] / 2) + 1;
  beam.complexSize[1] = beam.paddedSize[1];
  beam.complexKernelCreated = false;

  // Recompute kernel — fixture kernel is zero because calculate2DBeamKernel()
  // runs inside GalaxyFit_Simple() after the dump point
  calculate2DBeamKernel(beam, beam.pixelSize);

  return {
    pvIni,
    trfo,
    obs,
    beam,
    ftol:              f32(d.ftol),
    idum:              d.idum,
    linearLogSDSwitch: d.linearLogSDSwitch
  };
}


// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
function main() {
  const fixturePath    = process.argv[2];
  const fortranChi2    = process.argv[3] ? parseFloat(process.argv[3]) : null;

  if (!fixturePath) {
    console.error('Usage: node runFixture.js <fixture.json> [fortran_chi2]');
    process.exit(1);
  }

  const { pvIni, trfo, obs, beam, ftol, idum, linearLogSDSwitch }
    = loadFixture(fixturePath);

  console.log('\n=== Fixture loaded ===');
  console.log(`nParams:         ${pvIni.nParams}`);
  console.log(`nRings:          ${trfo.nRings}`);
  console.log(`nValid voxels:   ${obs.dh.nValid}`);
  console.log(`cube dims:       ${obs.dh.nPixels[0]} x ${obs.dh.nPixels[1]} x ${obs.dh.nChannels}`);
  console.log(`beam nRadial:    ${beam.nRadialCells}`);
  console.log(`ftol:            ${ftol}`);
  console.log(`idum:            ${idum}`);
  console.log(`linearLogSD:     ${linearLogSDSwitch}`);

  // Print initial PV
  console.log('\nInitial PV params:');
  for (let i = 0; i < pvIni.nParams; i++) {
    console.log(`  [${i}] ${pvIni.param[i].toFixed(6)}  [${pvIni.paramLowerLims[i].toFixed(4)}, ${pvIni.paramUpperLims[i].toFixed(4)}]  range=${pvIni.paramRange[i].toFixed(4)}  cyclic=${pvIni.cyclicSwitch[i]}`);
  }

  // Build working copies
  const pvModel    = new ParameterVector();
  const pvFirstFit = new ParameterVector();

  // Allocate pvModel so CP5 chi2 call works directly
  pvModel.nParams = pvIni.nParams;
  allocateParamVector(pvModel);
  pvModel.param.set(pvIni.param);
  pvModel.paramLowerLims.set(pvIni.paramLowerLims);
  pvModel.paramUpperLims.set(pvIni.paramUpperLims);
  pvModel.cyclicSwitch.set(pvIni.cyclicSwitch);
  pvModel.paramRange.set(pvIni.paramRange);

  const modelTR = new TiltedRingModel();
  modelTR.nRings            = trfo.nRings;
  modelTR.cmode             = 0;
  modelTR.cloudBaseSurfDens = f32(100.0);
  tiltRing_Allocate(modelTR);

  // Copy ring geometry from radial profiles
  for (let r = 0; r < trfo.nRings; r++) {
    modelTR.r[r].rmid   = f32(trfo.radialProfiles[r].rmid);
    modelTR.r[r].rwidth = f32(trfo.radialProfiles[r].rwidth);
  }

  const modelDC = new DataCube();

  // Allocate modelDC matching observedDC dimensions for CP5
  modelDC.dh.nPixels[0]    = obs.dh.nPixels[0];
  modelDC.dh.nPixels[1]    = obs.dh.nPixels[1];
  modelDC.dh.nChannels     = obs.dh.nChannels;
  modelDC.dh.pixelSize[0]  = obs.dh.pixelSize[0];
  modelDC.dh.pixelSize[1]  = obs.dh.pixelSize[1];
  modelDC.dh.channelSize   = obs.dh.channelSize;
  modelDC.dh.start[0]      = obs.dh.start[0];
  modelDC.dh.start[1]      = obs.dh.start[1];
  modelDC.dh.start[2]      = obs.dh.start[2];
  modelDC.dh.refLocation[0]= obs.dh.refLocation[0];
  modelDC.dh.refLocation[1]= obs.dh.refLocation[1];
  modelDC.dh.refLocation[2]= obs.dh.refLocation[2];
  modelDC.dh.refVal[0]     = obs.dh.refVal[0];
  modelDC.dh.refVal[1]     = obs.dh.refVal[1];
  modelDC.dh.refVal[2]     = obs.dh.refVal[2];
  modelDC.dh.uncertainty   = obs.dh.uncertainty;
  modelDC.dh.nValid        = obs.dh.nValid;
  allocateDataCube(modelDC);
  // allocateDataCube recomputes start from refVal/refLocation — override with fixture values
  modelDC.dh.start[0] = obs.dh.start[0];
  modelDC.dh.start[1] = obs.dh.start[1];
  modelDC.dh.start[2] = obs.dh.start[2];
  // Copy valid indices from observed
  modelDC.flattendValidIndices.set(obs.flattendValidIndices);

  const rng = makeRng(idum);

  const state = {
    pvIni,
    pvModel,
    pvFirstFit,
    modelTiltedRing:   modelTR,
    modelDC,
    observedDC:        obs,
    observedBeam:      beam,
    trFittingOptions:  trfo,
    rng,
    linearLogSDSwitch,
    ftol,
    iniGuessWidth:     f32(1.0),
    paramToTiltedRing: generalizedParamVectorToTiltedRing
  };

  // ---- CP5: evaluate chi² at initial params ----
  console.log('\n=== CP5: chi² at initial PV ===');
  const chi2Ini = tiltedRingModelComparison(Array.from(pvIni.param), state);
  console.log(`chi² (initial):  ${chi2Ini.toExponential(6)}`);
  if (fortranChi2 !== null) {
    const pct = Math.abs(chi2Ini - fortranChi2) / fortranChi2 * 100;
    console.log(`Fortran chi²:    ${fortranChi2.toExponential(6)}`);
    console.log(`Difference:      ${pct.toFixed(3)}%`);
    console.log(`CP5 status:      ${pct < 5.0 ? 'PASS (<5%)' : 'FAIL (>5%)'}`);
  }

  // Debug: test a single perturbation
  console.log('=== Debug: test perturbed params ===');
  const testParams = Array.from(pvIni.param);
  testParams[2] = 0.5;
  testParams[5] = 150.0;
  const chi2Test = tiltedRingModelComparison(testParams, state);
  console.log('chi2 at perturbed params:', chi2Test.toExponential(6));
  const { badModelCheck } = require('./CompareCubes/FullModelComparison.js');
  console.log('badModel check:', badModelCheck(state.modelTiltedRing, state.observedDC));

  // ---- CP6: run full optimizer ----
  console.log('\n=== CP6: galaxyFit_Simple ===');
  const t0     = Date.now();
  const pvBest = galaxyFit_Simple(state);
  const elapsed = ((Date.now() - t0) / 1000).toFixed(1);

  console.log(`\nOptimizer completed in ${elapsed}s`);
  console.log(`Best chi²: ${pvBest.bestLike.toExponential(6)}`);
  console.log('\nBest-fit parameters:');
  for (let i = 0; i < pvBest.nParams; i++) {
    console.log(`  [${i}] ${pvBest.param[i].toFixed(6)}`);
  }

  if (fortranChi2 !== null) {
    const pct = Math.abs(pvBest.bestLike - fortranChi2) / fortranChi2 * 100;
    console.log(`\nFortran best chi²: ${fortranChi2.toExponential(6)}`);
    console.log(`JS best chi²:      ${pvBest.bestLike.toExponential(6)}`);
    console.log(`Difference:        ${pct.toFixed(3)}%`);
    console.log(`CP6 status:        ${pct < 10.0 ? 'PASS (<10%)' : 'CHECK (>10%) — compare ring params manually'}`);
  }
}

main();
