'use strict';

// =============================================================================
// FullModelComparison.js
// High-fidelity port of src/CompareCubes/FullModelComparison.f
// (FullModelComparisonMod)
//
// PORTING NOTES
// -------------
// The Fortran uses PipelineGlobals for shared state (PVModel, ModelTiltedRing,
// ModelDC, ObservedDC, ObservedBeam, TR_FittingOptions, idum, PFlags).
// In JS all state is passed explicitly via a `state` context object,
// making the function pure and suitable for DCP worker use.
//
// state object fields:
//   pvModel           — ParameterVector (current candidate)
//   modelTiltedRing   — TiltedRingModel (modified in place each call)
//   modelDC           — DataCube (model cube, reset each call)
//   observedDC        — DataCube (observed data, read-only)
//   observedBeam      — Beam2D (beam kernel)
//   trFittingOptions  — TiltedRingFittingOptions
//   rng               — makeRng() object from random.js
//   linearLogSDSwitch — 0=linear, 1=log10 surface density
//   likePoint         — chi² function (default chi2Calc)
//
// Pi from CommonConsts: f32(Math.PI)
//
// MaskCubeMod is used in the Fortran but only for masking support which
// is handled upstream (ObservedDC.flattendValidIndices already set).
// No explicit masking needed here.
// =============================================================================

const f32 = Math.fround;
const Pi  = f32(Math.PI);

const { ring_ParticleDeAllocation }    = require('../ObjectDefinitions/TiltedRing.js');
const { buildTiltedRingModel }          = require('../TiltedModelGeneration/TiltedRingModelGeneration.js');
const { fillDataCubeWithTiltedRing }    = require('../TiltedRingToDataCube/FillDataCubeByTiltedRing.js');
const { cubeBeamConvolution }           = require('../ConvolveCube/CubeKernelConvolution.js');
const { cubeCompare }                   = require('./CubeComparison.js');
const { chi2Calc }                      = require('./LikelihoodFunctions.js');


// ---------------------------------------------------------------------------
// badModelCheck
// Fortran: BadModelCheck(BadModelFlag)
//
// Returns true if any ring parameter is unphysical:
//   - inclination outside [0, Pi/2]
//   - VRot < 0
//   - Sigma < 0
//   - CentPos outside cube bounds [0, nPixels-1]
//
// This is a fast early exit — if bad, chi2=1e20 and we skip cube synthesis.
// ---------------------------------------------------------------------------
function badModelCheck(modelTiltedRing, observedDC) {
  const dh = observedDC.dh;
  for (let i = 0; i < modelTiltedRing.nRings; i++) {
    const r = modelTiltedRing.r[i];
    if (r.inclination < f32(0.0) || r.inclination > f32(Pi / 2.0)) return true;
    if (r.vRot < f32(0.0))  return true;
    if (r.sigma < f32(0.0)) return true;
    for (let j = 0; j <= 1; j++) {
      if (r.centPos[j] < f32(0.0) ||
          r.centPos[j] > f32(dh.nPixels[j] - 1)) return true;
    }
  }
  return false;
}


// ---------------------------------------------------------------------------
// tiltedRingModelComparison
// Fortran: TiltedRingModelComparison(TestParams, chi2)
//
// The function called on every chi² evaluation by the optimizer.
// Steps:
//   1. Load candidate params into pvModel
//   2. Deserialize PV → TiltedRing via paramToTiltedRing
//   3. Apply SD switch (linear or log10)
//   4. BadModelCheck — return 1e20 if unphysical
//   5. BuildTiltedRingModel — generate particles
//   6. FillDataCubeWithTiltedRing — splat to model cube
//   7. CubeBeamConvolution — convolve with beam
//   8. CubeCompare — compute chi²
//   9. Deallocate particle arrays
//
// paramToTiltedRing is a function pointer in Fortran — pass explicitly.
// likePoint defaults to chi2Calc.
//
// Returns chi2 as f32.
// ---------------------------------------------------------------------------
function tiltedRingModelComparison(testParams, state) {
  const {
    pvModel,
    modelTiltedRing,
    modelDC,
    observedDC,
    observedBeam,
    trFittingOptions,
    rng,
    linearLogSDSwitch = 0,
    paramToTiltedRing,
    likePoint = chi2Calc
  } = state;

  // Step 1: load candidate params
  for (let i = 0; i < pvModel.nParams; i++) {
    pvModel.param[i] = f32(testParams[i]);
  }

  // Step 2: deserialize PV → TiltedRing
  paramToTiltedRing(pvModel, modelTiltedRing, trFittingOptions);

  // Step 3: apply SD switch
  for (let i = 0; i < modelTiltedRing.nRings; i++) {
    const r = modelTiltedRing.r[i];
    if (linearLogSDSwitch === 0) {
      r.sigma = f32(r.sigUse);
    } else if (linearLogSDSwitch === 1) {
      r.sigma = f32(Math.pow(10.0, r.sigUse));
    }
  }

  // Step 4: physicality check
  if (badModelCheck(modelTiltedRing, observedDC)) {
    return f32(1.0e20);
  }

  // Step 5: build particle model
  buildTiltedRingModel(modelTiltedRing, rng);

  // Step 6: fill model cube
  fillDataCubeWithTiltedRing(modelDC, modelTiltedRing);

  // Step 7: beam convolution
  cubeBeamConvolution(modelDC, observedBeam);

  // Step 8: compare cubes
  const chi2 = cubeCompare(
    observedDC,
    modelDC,
    observedDC.dh.uncertainty,
    likePoint
  );

  // Step 9: deallocate particle arrays
  for (let i = 0; i < modelTiltedRing.nRings; i++) {
    ring_ParticleDeAllocation(modelTiltedRing.r[i]);
  }

  return chi2;
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = { tiltedRingModelComparison, badModelCheck };


// ---------------------------------------------------------------------------
// Self-test (node FullModelComparison.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { makeRng }              = require('../StandardMath/random.js');
  const { DataCube, allocateDataCube } = require('../ObjectDefinitions/DataCube.js');
  const { Beam2D, allocate_Beam2D }    = require('../ObjectDefinitions/Beam.js');
  const { TiltedRingModel, TiltedRingFittingOptions,
          tiltRing_Allocate, tiltRingFittingOptions_Allocate,
          logicalTiltedRingIndexing }  = require('../ObjectDefinitions/TiltedRing.js');
  const { ParameterVector, allocateParamVector } = require('../ObjectDefinitions/ParameterVector.js');
  const { calculate2DBeamKernel }      = require('../ConvolveCube/CalculateBeamKernel.js');
  const { generalizedParamVectorToTiltedRing } = require('../ParameterToTiltedRingInterface/ParameterToTiltedRingVector.js');

  // --- Build observed cube (point source) ---
  const obs = new DataCube();
  const dh  = obs.dh;
  dh.nPixels[0] = 64; dh.nPixels[1] = 64; dh.nChannels = 64;
  dh.pixelSize[0] = f32(-1.0); dh.pixelSize[1] = f32(1.0);
  dh.channelSize  = f32(10.0);
  dh.refLocation[0] = f32(32.0); dh.refLocation[1] = f32(32.0);
  dh.refLocation[2] = f32(32.0);
  dh.refVal[0] = f32(0.0); dh.refVal[1] = f32(0.0);
  dh.refVal[2] = f32(1000.0);
  dh.uncertainty = f32(0.001);
  allocateDataCube(obs);

  // --- Build model cube (same dimensions) ---
  const mdl = new DataCube();
  mdl.dh = Object.assign(Object.create(Object.getPrototypeOf(dh)), dh);
  mdl.dh.nPixels    = new Int32Array(dh.nPixels);
  mdl.dh.pixelSize  = new Float32Array(dh.pixelSize);
  mdl.dh.refLocation = new Float32Array(dh.refLocation);
  mdl.dh.refVal      = new Float32Array(dh.refVal);
  allocateDataCube(mdl);

  // --- Build beam ---
  const beam = new Beam2D();
  beam.beamFWHM    = f32(3.0);
  beam.sigmaLengths = f32(3.0);
  beam.pixelSize[0] = f32(-1.0);
  beam.pixelSize[1] = f32(1.0);
  allocate_Beam2D(beam, new Int32Array([64, 64]));
  calculate2DBeamKernel(beam, beam.pixelSize);

  // --- Build a simple 2-ring TiltedRingModel ---
  const tr = new TiltedRingModel();
  tr.nRings = 2; tr.cmode = 1; tr.cloudBaseSurfDens = f32(1.0);
  tiltRing_Allocate(tr);

  const trfo = new TiltedRingFittingOptions();
  trfo.nRings = 2; trfo.nRingsPerBeam = 2; trfo.nTargRings = -1;
  // All params constant across rings, all fitted
  for (let i = 0; i <= 12; i++) {
    trfo.constParams[i] = true;
    trfo.fixedParams[i] = false;
  }
  // Fix params we don't want to fit (just keep it simple)
  trfo.fixedParams[6]  = true;  // VRad
  trfo.fixedParams[8]  = true;  // Vvert
  trfo.fixedParams[9]  = true;  // dvdz
  trfo.fixedParams[11] = true;  // z0
  trfo.fixedParams[12] = true;  // zGradiantStart
  tiltRingFittingOptions_Allocate(trfo);
  logicalTiltedRingIndexing(trfo);

  // Set up fitting options radial profiles
  const r0 = trfo.radialProfiles[0];
  const r1 = trfo.radialProfiles[1];
  for (const r of [r0, r1]) {
    r.centPos[0] = f32(32.0); r.centPos[1] = f32(32.0);
    r.inclination = f32(45.0 * Math.PI / 180.0);
    r.positionAngle = f32(30.0 * Math.PI / 180.0);
    r.vSys = f32(1000.0); r.vRot = f32(150.0);
    r.vDisp = f32(8.0); r.sigma = f32(0.01); r.sigUse = f32(0.01);
    r.z0 = f32(0.0); r.zGradiantStart = f32(0.0);
  }
  r0.rmid = f32(5.0); r0.rwidth = f32(2.0);
  r1.rmid = f32(10.0); r1.rwidth = f32(2.0);

  // Set param limits
  trfo.paramLowerLims[0] = f32(0.0);   trfo.paramUpperLims[0] = f32(64.0);
  trfo.paramLowerLims[1] = f32(0.0);   trfo.paramUpperLims[1] = f32(64.0);
  trfo.paramLowerLims[2] = f32(0.0);   trfo.paramUpperLims[2] = f32(Math.PI/2);
  trfo.paramLowerLims[3] = f32(0.0);   trfo.paramUpperLims[3] = f32(2*Math.PI);
  trfo.paramLowerLims[4] = f32(900.0); trfo.paramUpperLims[4] = f32(1100.0);
  trfo.paramLowerLims[5] = f32(50.0);  trfo.paramUpperLims[5] = f32(300.0);
  trfo.paramLowerLims[7] = f32(0.0);   trfo.paramUpperLims[7] = f32(20.0);
  trfo.paramLowerLims[10]= f32(0.0);   trfo.paramUpperLims[10]= f32(1.0);

  // Build PVModel
  const pv = new ParameterVector();
  pv.nParams = trfo.nFittedParamsTotal;
  allocateParamVector(pv);

  // Fill PV with initial values from trfo
  const { tiltedRingOptionsToPV } = require('../ParameterToTiltedRingInterface/ParameterToTiltedRingVector.js');
  tiltedRingOptionsToPV(pv, trfo);

  const state = {
    pvModel:           pv,
    modelTiltedRing:   tr,
    modelDC:           mdl,
    observedDC:        obs,
    observedBeam:      beam,
    trFittingOptions:  trfo,
    rng:               makeRng(-1),
    linearLogSDSwitch: 0,
    paramToTiltedRing: generalizedParamVectorToTiltedRing
  };

  console.log('=== tiltedRingModelComparison ===');
  console.log('nParams:', pv.nParams);
  const chi2 = tiltedRingModelComparison(Array.from(pv.param), state);
  console.log('chi2:', chi2.toExponential(4));
  console.log('chi2 is finite:', isFinite(chi2) ? 'OK' : 'FAIL (bad model)');

  // Bad model test
  const badParams = Array.from(pv.param);
  badParams[2] = -0.1; // negative inclination → bad
  const chi2Bad = tiltedRingModelComparison(badParams, state);
  console.log('chi2 (bad model):', chi2Bad.toExponential(4), '(expect 1e20)');
}
