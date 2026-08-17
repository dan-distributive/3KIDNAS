'use strict';

// =============================================================================
// ModellingInitializations.js
// Port of src/PreAnalysis/ModellingInitializations.f (ModellingInitializationsMod)
//
// SCOPE: SetUpFittingParameters (the top-level Fortran driver) only ever
// takes the PFlags%CoreCodeSwitch==1 branch in this project (the inbuilt
// downhill-simplex path -- the only optimizer this whole port implements),
// so that branch-check itself isn't ported; callers here just run the three
// steps directly: setupTRFittingOptionsFromModelTR ->
// setupTRParamLimitsFromModelTR -> setInitialParamVector.
//
// PRECISION
// ---------
// Fortran `real` -> Math.fround() per step, matching the -O0 build.
// `ApproxSize=nRings/TR_FittingOptions%nRingsPerBeam`: BOTH operands are
// INTEGER, so this is integer (truncating) division, THEN converted to real
// -- NOT floating-point division. Ported as Math.trunc(nRings/nRingsPerBeam)
// before the f32 conversion.
// =============================================================================

const f32 = Math.fround;
const Pi = f32(Math.PI);

const { TiltedRingModel, TiltedRingFittingOptions, tiltRing_Allocate,
        tiltRingFittingOptions_Allocate, logicalTiltedRingIndexing }
  = require('../ObjectDefinitions/TiltedRing.js');
const { ParameterVector, allocateParamVector }
  = require('../ObjectDefinitions/ParameterVector.js');
const { tiltedRingOptionsToPV }
  = require('../ParameterToTiltedRingInterface/ParameterToTiltedRingVector.js');
const { msolPc2ToJyAS2 } = require('../UnitConversions/UnitConversions.js');

// ---------------------------------------------------------------------------
// setupTiltedRingStructures
// Fortran: SetupTiltedRingStructures(nRings)
//
// Allocates a fresh ModelTiltedRing + TiltedRingFittingOptions at nRings,
// sets every ring's width to beamMajorAxis/nRingsPerBeam. sizeFlag is a
// diagnostic only (Fortran: GalaxyDict%Flags%SizeFlag) -- computed for
// fidelity but not consumed by the fit itself.
// ---------------------------------------------------------------------------
function setupTiltedRingStructures(nRings, nRingsPerBeam, beamMajorAxis, sizeLims) {
  const approxSize = f32(Math.trunc(nRings / nRingsPerBeam));
  let sizeFlag;
  if (approxSize < sizeLims[0]) sizeFlag = 1;
  else if (approxSize > sizeLims[1]) sizeFlag = 2;
  else sizeFlag = 0;

  const modelTiltedRing = new TiltedRingModel();
  modelTiltedRing.nRings = nRings;
  tiltRing_Allocate(modelTiltedRing);

  const trFittingOptions = new TiltedRingFittingOptions();
  trFittingOptions.nRingsPerBeam = nRingsPerBeam;
  trFittingOptions.nRings = nRings;
  tiltRingFittingOptions_Allocate(trFittingOptions);

  const rWidth = f32(beamMajorAxis / f32(nRingsPerBeam));
  for (let i = 0; i < nRings; i++) {
    modelTiltedRing.r[i].rwidth = rWidth;
  }

  return { modelTiltedRing, trFittingOptions, sizeFlag };
}

// ---------------------------------------------------------------------------
// setupTRFittingOptionsFromModelTR
// Fortran: SetupTRFittingOptionsFromModelTR(SDSwitch)
//
// Copies ModelTiltedRing.R(i) (full struct) into TR_FittingOptions'
// radialProfiles(i), sets sigUse per sdSwitch (0=linear Sigma, 1=logSigma),
// then runs logicalTiltedRingIndexing to compute the free-parameter counts.
// ---------------------------------------------------------------------------
function setupTRFittingOptionsFromModelTR(sdSwitch, modelTiltedRing, trFittingOptions) {
  for (let i = 0; i < trFittingOptions.nRings; i++) {
    const src = modelTiltedRing.r[i];
    const dst = trFittingOptions.radialProfiles[i];
    dst.rmid           = src.rmid;
    dst.rwidth         = src.rwidth;
    dst.centPos[0]     = src.centPos[0];
    dst.centPos[1]     = src.centPos[1];
    dst.inclination    = src.inclination;
    dst.positionAngle  = src.positionAngle;
    dst.vSys           = src.vSys;
    dst.vRot           = src.vRot;
    dst.vRad           = src.vRad;
    dst.vDisp          = src.vDisp;
    dst.vvert          = src.vvert;
    dst.dvdz           = src.dvdz;
    dst.sigma          = src.sigma;
    dst.logSigma       = src.logSigma;
    dst.z0             = src.z0;
    dst.zGradiantStart = src.zGradiantStart;

    if (sdSwitch === 0) dst.sigUse = dst.sigma;
    else if (sdSwitch === 1) dst.sigUse = dst.logSigma;
  }

  logicalTiltedRingIndexing(trFittingOptions);
}

// ---------------------------------------------------------------------------
// setupTRParamLimitsFromModelTR
// Fortran: SetupTRParamLimitsFromModelTR(SDSwitch, CentSize)
//
// Sets the 13 parameter lower/upper limits, ranges, and cyclic switches.
// Most of these are fixed constants (identical for the initial fit and
// every bootstrap realization); center (0,1) depends on THIS cube's own
// center estimate (modelTiltedRing.r[0].centPos), VSys (4) on the channel
// axis range, and Sigma (10, sdSwitch=0 branch) on this cube's pixel/channel
// geometry via a Jy/beam<->Msol/pc^2-style conversion.
//
// ctx: { minChannel, maxChannel, vRotLims: [lo,hi], linSDLims: [lo,hi],
//        logSDLims: [lo,hi], pixelSize0, pixelSize1, channelSize }
// ---------------------------------------------------------------------------
function setupTRParamLimitsFromModelTR(sdSwitch, centSize, modelTiltedRing, trFittingOptions, ctx) {
  const trfo = trFittingOptions;
  for (let i = 0; i <= 12; i++) trfo.cyclicSwitch[i] = 0;

  const r0 = modelTiltedRing.r[0];
  for (let i = 0; i <= 1; i++) {
    trfo.paramLowerLims[i] = f32(r0.centPos[i] - centSize);
    trfo.paramUpperLims[i] = f32(r0.centPos[i] + centSize);
    trfo.paramRange[i] = f32(5.0);
  }

  trfo.paramLowerLims[2] = f32(0.0);
  trfo.paramUpperLims[2] = f32(Pi / f32(2.0));
  trfo.paramRange[2] = f32(f32(10.0) * Pi / f32(180.0));

  trfo.paramLowerLims[3] = f32(0.0);
  trfo.paramUpperLims[3] = f32(f32(2.0) * Pi);
  trfo.paramRange[3] = f32(f32(10.0) * Pi / f32(180.0));
  trfo.cyclicSwitch[3] = 1;

  trfo.paramLowerLims[4] = ctx.minChannel;
  trfo.paramUpperLims[4] = ctx.maxChannel;
  trfo.paramRange[4] = f32(20.0);

  trfo.paramLowerLims[5] = ctx.vRotLims[0];
  trfo.paramUpperLims[5] = ctx.vRotLims[1];
  trfo.paramRange[5] = f32(50.0);

  trfo.paramLowerLims[6] = f32(-50.0);
  trfo.paramUpperLims[6] = f32(50.0);
  trfo.paramRange[6] = f32(50.0);

  trfo.paramLowerLims[7] = f32(0.0);
  trfo.paramUpperLims[7] = f32(20.0);
  trfo.paramRange[7] = f32(5.0);

  trfo.paramLowerLims[8] = f32(-20.0);
  trfo.paramUpperLims[8] = f32(20.0);
  trfo.paramRange[8] = f32(5.0);

  trfo.paramLowerLims[9] = f32(0.0);
  trfo.paramUpperLims[9] = f32(5.0);
  trfo.paramRange[9] = f32(1.0);

  if (sdSwitch === 0) {
    trfo.paramLowerLims[10] = ctx.linSDLims[0];
    trfo.paramUpperLims[10] = ctx.linSDLims[1];

    const prod    = f32(ctx.pixelSize0 * ctx.pixelSize1);
    const absProd = f32(Math.abs(prod));
    const sdConv1 = f32(1.0 / absProd);
    let sdTemp = msolPc2ToJyAS2(f32(2.0));
    sdTemp = f32(sdTemp / sdConv1);
    sdTemp = f32(sdTemp / f32(Math.abs(ctx.channelSize)));
    trfo.paramRange[10] = sdTemp;
  } else if (sdSwitch === 1) {
    trfo.paramLowerLims[10] = ctx.logSDLims[0];
    trfo.paramUpperLims[10] = ctx.logSDLims[1];
    trfo.paramRange[10] = f32(1.0);
  }

  trfo.paramLowerLims[11] = f32(0.0);
  trfo.paramUpperLims[11] = f32(5.0);
  trfo.paramRange[11] = f32(1.0);

  trfo.paramLowerLims[12] = f32(0.5);
  trfo.paramUpperLims[12] = f32(10.0);
  trfo.paramRange[12] = f32(1.0);
}

// ---------------------------------------------------------------------------
// setInitialParamVector
// Fortran: SetInitialParamVector()
//
// Allocates PVIni at trFittingOptions.nFittedParamsTotal and serializes the
// fitting options' radial profiles/limits into it via the already-ported
// tiltedRingOptionsToPV.
// ---------------------------------------------------------------------------
function setInitialParamVector(trFittingOptions) {
  const pvIni = new ParameterVector();
  pvIni.nParams = trFittingOptions.nFittedParamsTotal;
  allocateParamVector(pvIni);
  tiltedRingOptionsToPV(pvIni, trFittingOptions);
  return pvIni;
}

// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  setupTiltedRingStructures,
  setupTRFittingOptionsFromModelTR,
  setupTRParamLimitsFromModelTR,
  setInitialParamVector,
};


// ---------------------------------------------------------------------------
// Self-test (node ModellingInitializations.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  console.log('=== setupTiltedRingStructures ===');
  const nRings = 5, nRingsPerBeam = 2, beamMajorAxis = f32(10.0);
  const { modelTiltedRing, trFittingOptions, sizeFlag } =
    setupTiltedRingStructures(nRings, nRingsPerBeam, beamMajorAxis, [1, 100]);
  console.log('  nRings:', modelTiltedRing.nRings, '(expect 5)');
  console.log('  ring[0].rwidth:', modelTiltedRing.r[0].rwidth, '(expect 5)');
  console.log('  sizeFlag:', sizeFlag, '(expect 0, approxSize=trunc(5/2)=2, within [1,100])');

  // Populate the model with plausible values before deriving fitting options.
  for (let i = 0; i < nRings; i++) {
    const r = modelTiltedRing.r[i];
    r.rmid = f32(i * 5 + 2.5);
    r.centPos[0] = f32(20.0); r.centPos[1] = f32(18.0);
    r.inclination = f32(0.5); r.positionAngle = f32(1.0);
    r.vSys = f32(1000.0); r.vRot = f32(100.0 + i * 5);
    r.vRad = f32(0.0); r.vDisp = f32(8.0);
    r.vvert = f32(0.0); r.dvdz = f32(0.0);
    r.sigma = f32(2.0); r.logSigma = f32(Math.log10(2.0));
    r.z0 = f32(0.0); r.zGradiantStart = f32(0.0);
  }

  console.log('\n=== setupTRFittingOptionsFromModelTR (sdSwitch=0) ===');
  setupTRFittingOptionsFromModelTR(0, modelTiltedRing, trFittingOptions);
  console.log('  radialProfiles[2].sigUse (expect 2.0):', trFittingOptions.radialProfiles[2].sigUse);
  console.log('  nFittedParamsTotal (all radial, free, 13*5=65):', trFittingOptions.nFittedParamsTotal);

  console.log('\n=== setupTRParamLimitsFromModelTR ===');
  const ctx = {
    minChannel: f32(500.0), maxChannel: f32(1500.0),
    vRotLims: [f32(0.0), f32(300.0)],
    linSDLims: [f32(0.001), f32(10.0)],
    logSDLims: [f32(-3.0), f32(1.0)],
    pixelSize0: f32(6.0), pixelSize1: f32(6.0),
    channelSize: f32(4.0),
  };
  setupTRParamLimitsFromModelTR(0, f32(2.0 * beamMajorAxis), modelTiltedRing, trFittingOptions, ctx);
  console.log('  centerLowerLims[0] (expect 20-20=0):', trFittingOptions.paramLowerLims[0]);
  console.log('  inclUpperLims[2] (expect Pi/2):', trFittingOptions.paramUpperLims[2], f32(Math.PI / 2));
  console.log('  paCyclicSwitch[3] (expect 1):', trFittingOptions.cyclicSwitch[3]);
  console.log('  vSysLowerLims[4] (expect 500):', trFittingOptions.paramLowerLims[4]);

  console.log('\n=== setInitialParamVector ===');
  const pvIni = setInitialParamVector(trFittingOptions);
  console.log('  nParams (expect 65):', pvIni.nParams);
  console.log('  param.length (expect 65):', pvIni.param.length);
  console.log('  param[0] (X center, expect 20):', pvIni.param[0]);
}
