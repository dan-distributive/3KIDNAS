'use strict';

// =============================================================================
// InitialAnalysis.js
// Port of src/PreAnalysis/InitialAnalysis.f (InitialAnalysisMod) +
// src/PreAnalysis/InitialAnalysis_PreProcessing.f
// (InitialAnalysisPreProcessingMod) -- consolidated into one file since
// they're tightly coupled (InitialAnalysis_Prep builds the moment
// maps/velocity profile that InitialAnalysis's ProjectedAnalysis consumes)
// and neither is large once the heavy lifting is delegated to the
// already-ported modules below.
//
// SCOPE DECISIONS (confirmed by reading the real call graph, not guessed):
//
// - InitialAnalysis_Prep: only the MomentMapSwitch==0 branch is ported
//   (ConstructProjectionsFromCube from the observed cube + SoFiA mask --
//   what every real bootstrap-realization / single-galaxy fit actually
//   uses). MomentMapSwitch==1/2 (load precomputed maps from files) and ==3
//   (secondary cube/mask path) aren't reachable from a DCP bootstrap
//   realization's own inputs.
//
// - GetNoise: initialAnalysis()'s own return value carries ONLY the RMS (via
//   estimateNoise) -- every bootstrap realization calls initialAnalysis() and
//   has no use for SN_Median/SN_Peak/SN_Avg/SN_Int (those feed ONLY
//   WriteBootstrappedFitOutputFile_Text's header block for the single
//   BEST-FIT model, Bootstrap_Outputs.py reading Model['RMS'][0] etc.), so
//   folding them into initialAnalysis() would be dead weight on the
//   per-realization hot path. computeSNDiagnostics() below ports that half
//   of Fortran's GetNoise (InitialAnalysis.f:244-273) as a separate function,
//   called explicitly only by the single initial/anchor JS fit
//   (RunInitialFitDCP), which does need them.
//
// - GetGalaxyShape: real Fortran's PFlags%ShapeSource branch (EstimateShape
//   or catalogue-ellipse-ratio) computes Incl/PA, but BOTH are then
//   unconditionally overwritten from CatItem%EllipseInc/EllipsePA (the SoFiA
//   catalogue) immediately after -- confirmed by reading InitialAnalysis.f
//   directly, no branch skips this override. So only the catalogue
//   conversion is ported; the dead branch (and EstimateShape.js's
//   estimateShape, which remains available/ported for other potential
//   callers) is not invoked from here.
//
// - EstimateFlatVDisp: Fortran hardcodes VDisp=10. (not data-derived at all).
//
// PRECISION
// ---------
// Fortran `real` -> Math.fround() per step, matching the -O0 build.
// =============================================================================

const f32 = Math.fround;
const Pi = f32(Math.PI);

const { copyDataCube } = require('../BootstrapSampler/FlipBootstrap.js');
const { maskCube } = require('../CompareCubes/MaskCube.js');
const { constructMomentMaps, makeVelProfile } = require('./GetMomentMaps.js');
const { estimateNoise } = require('./EstimateCubeNoise.js');
const { iterEstimateCenter } = require('./EstimateShape.js');
const { estimateVSysAndEdgesSimple } = require('./VelProfileAnalysis.js');
const { estimateProfiles } = require('./EstimateRadialProfiles.js');
const {
  setupTiltedRingStructures,
  setupTRFittingOptionsFromModelTR,
  setupTRParamLimitsFromModelTR,
  setInitialParamVector,
} = require('./ModellingInitializations.js');

// ---------------------------------------------------------------------------
// constructProjectionsFromCube
// Fortran: ConstructProjectionsFromCube(ObsDC, MaskDC)
//
// Masks a copy of the observed cube, then derives the moment maps and
// integrated velocity profile from it.
// ---------------------------------------------------------------------------
function constructProjectionsFromCube(obsDC, maskDC) {
  const maskedObservedDC = copyDataCube(obsDC);
  maskCube(maskedObservedDC, maskDC);

  const observedMaps = constructMomentMaps(maskedObservedDC);
  const observedVelocityProfile = makeVelProfile(maskedObservedDC);

  return { maskedObservedDC, observedMaps, observedVelocityProfile };
}

// ---------------------------------------------------------------------------
// getNoise
// Fortran: GetNoise(Noise, CatItem) -- RMS only, see module header.
// ---------------------------------------------------------------------------
function getNoise(observedDC) {
  return estimateNoise(observedDC);
}

// ---------------------------------------------------------------------------
// computeSNDiagnostics
// Fortran: the SN_Median/SN_Peak/SN_Avg/SN_Int half of GetNoise (InitialAnalysis.f
// :244-273), skipped when this file was first ported (see module header) since
// nothing in the per-bootstrap-realization path consumed them. Needed now for
// the single initial/anchor fit's output (RunInitialFitDCP's BestFitModel dict
// -- see Bootstrap_Outputs.py, which reads Model['RMS'][0] etc. off the ONE
// best-fit model, not per-realization). Deliberately NOT folded into
// initialAnalysis()'s own return value / call graph -- every bootstrap
// realization already calls initialAnalysis() and has no use for these, so
// keeping this a separate, explicitly-called function means the existing,
// verified bootstrap path is byte-for-byte unchanged.
//
// maskedObservedDC: result.maskedObservedDC from initialAnalysis() (observed
//   cube with everything outside the SoFiA mask zeroed -- Fortran's
//   MaskedObservedDC).
// maskDC: the 0/1 mask cube passed into initialAnalysis() (Fortran's
//   DataCubeMask) -- used ONLY for nCells (sum of mask=1 voxels), matching
//   Fortran's nCells=int(sum(DataCubeMask%Flux)) exactly -- NOT the same
//   count as "how many nonzero entries MakeMaskedFlatFluxArr finds" in the
//   (never hit in practice, since real flux is never exactly 0.0) edge case
//   where a masked-in voxel's observed flux happens to be exactly zero.
// noise: RMS, from getNoise()/initialAnalysis()'s `noise` field.
// beamAreaPixels: ObservedBeam%BeamAreaPixels (Beam2D.beamAreaPixels).
// ---------------------------------------------------------------------------
function computeSNDiagnostics(maskedObservedDC, maskDC, noise, beamAreaPixels) {
  const eps = f32(1.0e-20);
  const flux = maskedObservedDC.flux;

  // SPeak=maxval(...), SInt=sum(...) -- over the FULL cube array (masked-out
  // voxels are exactly 0.0, so this is equivalent to max/sum over the
  // nonzero-only set for SInt, but NOT for SPeak if every true value were
  // negative -- matching Fortran's maxval/sum over the whole array, not a
  // filtered one, preserves that edge case faithfully).
  let sPeak = flux[0], sInt = f32(0.0);
  for (let i = 0; i < flux.length; i++) {
    const v = flux[i];
    if (v > sPeak) sPeak = v;
    sInt = f32(sInt + f32(v));
  }

  // nCells=int(sum(DataCubeMask%Flux)) -- the 0/1 mask cube's own sum, NOT
  // a count of nonzero maskedObservedDC entries (see header comment).
  let maskSum = f32(0.0);
  const mflux = maskDC.flux;
  for (let i = 0; i < mflux.length; i++) maskSum = f32(maskSum + f32(mflux[i]));
  const nCells = Math.trunc(maskSum);

  // MakeMaskedFlatFluxArr: every |value|>eps entry of MaskedObservedDC%Flux,
  // in cube order (order doesn't matter -- sorted below).
  const flat = [];
  for (let i = 0; i < flux.length; i++) {
    if (Math.abs(flux[i]) > eps) flat.push(flux[i]);
  }
  flat.sort((a, b) => a - b);
  // MedIndx=nCells/2 (Fortran integer division, 1-indexed) -> 0-indexed
  // (nCells/2)-1. Ties are moot: duplicate float32 values return the same
  // number regardless of which duplicate's index a different sort algorithm
  // than Fortran's indexx would have placed there.
  const medIndx = Math.trunc(nCells / 2) - 1;
  const medianFlux = flat[medIndx];

  const snMedian = f32(medianFlux / noise);
  const snPeak = f32(sPeak / noise);
  const snAvg = f32(sInt / f32(nCells * noise));
  const snInt = f32(sInt / f32(Math.sqrt(f32(nCells * beamAreaPixels)) * noise));

  return { snMedian, snPeak, snAvg, snInt, nCells, rms: noise };
}

// ---------------------------------------------------------------------------
// getObjectCenter
// Fortran: GetObjectCenter(Center, CatItem)
//
// centerSource: 0 = iterative moment-map center estimate (the real analysis
//   path), 1 = observedDC's own RefLocation, 2 = observedMaps' PixelCent.
// ---------------------------------------------------------------------------
function getObjectCenter(observedMaps, centerSource, observedDC) {
  let center, centerFlag = 0;

  if (centerSource === 0) {
    const r = iterEstimateCenter(observedMaps);
    center = r.center;
    centerFlag = r.centerFlag;
  } else if (centerSource === 1) {
    center = new Float32Array([observedDC.dh.refLocation[0], observedDC.dh.refLocation[1]]);
  } else if (centerSource === 2) {
    center = new Float32Array([observedMaps.dh.pixelCent[0], observedMaps.dh.pixelCent[1]]);
  }

  return { center, centerFlag };
}

// ---------------------------------------------------------------------------
// getGalaxyShape
// Fortran: GetGalaxyShape(Incl, PA, Center, CatItem) -- catalogue conversion
// only, see module header for why the PFlags%ShapeSource branch is dead.
// ---------------------------------------------------------------------------
function getGalaxyShape(catalogueEllipseIncDeg, catalogueEllipsePADeg) {
  const incl = f32(f32(catalogueEllipseIncDeg) * Pi / f32(180.0));
  let pa = f32(f32(f32(catalogueEllipsePADeg) + f32(90.0)) * Pi / f32(180.0));
  if (pa < f32(0.0)) pa = f32(pa + f32(2.0) * Pi);
  if (pa > f32(2.0) * Pi) pa = f32(pa - f32(2.0) * Pi);
  return { incl, pa };
}

// ---------------------------------------------------------------------------
// estimateFlatVDisp
// Fortran: EstimateFlatVDisp(VDisp) -- hardcoded constant, not data-derived.
// ---------------------------------------------------------------------------
function estimateFlatVDisp() {
  return f32(10.0);
}

// ---------------------------------------------------------------------------
// convertFlatDiskProfilesToTR
// Fortran: ConvertFlatDiskProfilesToTR(nRings, Center, Incl, PA, VSys, VDisp,
//                                       EstimatedRadialProfiles)
//
// Builds a fresh ModelTiltedRing/TiltedRingFittingOptions pair at nRings and
// fills every ring from the (shared) shape/center/VSys/VDisp plus this
// ring's own [radius, sd, vRot] from estimatedRadialProfiles (the
// EstimateRadialProfiles.js RadialProfile representation).
// ---------------------------------------------------------------------------
function convertFlatDiskProfilesToTR(nRings, center, incl, pa, vSys, vDisp, estimatedRadialProfiles, nRingsPerBeam, beamMajorAxis, sizeLims) {
  const { modelTiltedRing, trFittingOptions, sizeFlag } =
    setupTiltedRingStructures(nRings, nRingsPerBeam, beamMajorAxis, sizeLims);

  for (let i = 0; i < nRings; i++) {
    const r = modelTiltedRing.r[i];
    const erp = estimatedRadialProfiles[i];
    r.rmid = erp[0];
    r.centPos[0] = center[0];
    r.centPos[1] = center[1];
    r.inclination = incl;
    r.positionAngle = pa;
    r.vSys = vSys;
    r.vRot = erp[2];
    r.vRad = f32(0.0);
    r.vDisp = vDisp;
    r.vvert = f32(0.0);
    r.dvdz = f32(0.0);
    r.sigma = erp[1];
    r.logSigma = f32(Math.log10(erp[1]));
    r.z0 = f32(0.0);
    r.zGradiantStart = f32(f32(5.0) * r.z0);
  }

  return { modelTiltedRing, trFittingOptions, sizeFlag };
}

// ---------------------------------------------------------------------------
// projectedAnalysis
// Fortran: ProjectedAnalysis(CatItem)
//
// opts: {
//   centerSource,                          // PFlags%CenterSource
//   catalogueEllipseIncDeg, catalogueEllipsePADeg,  // SoFiA catalogue, deg
//   nRingsPerBeam, nTargRings, radGridArcsec,       // ring grid (pinned for bootstrap fits)
//   sdSwitch,                              // PFlags%Linear_Log_SDSwitch
//   linSDLims, logSDLims, vRotLims, sizeLims,       // config limit pairs [lo,hi]
//   noiseSigmaLim, beamMajorAxis,
// }
// ---------------------------------------------------------------------------
function projectedAnalysis(observedDC, observedMaps, observedVelocityProfile, opts) {
  const { center } = getObjectCenter(observedMaps, opts.centerSource, observedDC);
  let { incl, pa } = getGalaxyShape(opts.catalogueEllipseIncDeg, opts.catalogueEllipsePADeg);

  let { vSys } = estimateVSysAndEdgesSimple(observedVelocityProfile);

  let minV = observedVelocityProfile.vel[0], maxV = observedVelocityProfile.vel[0];
  for (let i = 1; i < observedVelocityProfile.vel.length; i++) {
    const v = observedVelocityProfile.vel[i];
    if (v < minV) minV = v;
    if (v > maxV) maxV = v;
  }
  if (vSys < minV || vSys > maxV) {
    vSys = f32(f32(minV + maxV) / f32(2.0));
  }

  // MomentMapSwitch===0 branch: propagate the cube's own noise onto the maps.
  observedMaps.dh.uncertainty = observedDC.dh.uncertainty;

  const sdLimsForFit = opts.sdSwitch === 0 ? opts.linSDLims : opts.logSDLims;
  const { radialProfile, pa: paFinal, nRings } = estimateProfiles(
    observedMaps, opts.beamMajorAxis, center, pa, vSys,
    { nRingsPerBeam: opts.nRingsPerBeam, nTargRings: opts.nTargRings, radGridArcsec: opts.radGridArcsec, incl },
    sdLimsForFit, opts.vRotLims, opts.noiseSigmaLim, observedDC
  );
  pa = paFinal;

  const vDisp = estimateFlatVDisp();

  const { modelTiltedRing, trFittingOptions } = convertFlatDiskProfilesToTR(
    nRings, center, incl, pa, vSys, vDisp, radialProfile,
    opts.nRingsPerBeam, opts.beamMajorAxis, opts.sizeLims
  );

  // Which of the 13 ring properties are constant-across-rings / permanently
  // fixed (not varying, not optimized) comes from WRKP_Options.txt's fitting
  // flags -- Fortran-side configuration this function has no way to derive
  // on its own. trFittingOptions.constParams/fixedParams default to all-false
  // (tiltRingFittingOptions_Allocate), which means "every property free per
  // ring" -- silently wrong (13x too many free params: found via a direct
  // numerical cross-check against Fortran's real ProbeSwitch chi2, where this
  // produced nParams=65 instead of Fortran's real 16, and a garbage chi2).
  // opts.constParams/opts.fixedParams (13-element boolean arrays, same order
  // as TiltedRing.js's Ring fields) must be threaded in from the real
  // Fortran fixture dump (diskfit_fixture.json's own constParams/fixedParams
  // fields) -- required, not optional, since there is no safe default here
  // (the three previously-scattered hardcoded conventions in
  // ParameterToTiltedRingVector.js/FullModelComparison.js/TiltedRing.js
  // disagree with each other and with the real fixture).
  if (!opts.constParams || !opts.fixedParams) {
    throw new Error('initialAnalysis: opts.constParams/opts.fixedParams are required '
      + '(13-element boolean arrays from the real Fortran fixture dump) -- '
      + 'without them this silently computes the wrong number of free parameters.');
  }
  for (let i = 0; i <= 12; i++) {
    trFittingOptions.constParams[i] = !!opts.constParams[i];
    trFittingOptions.fixedParams[i] = !!opts.fixedParams[i];
  }

  setupTRFittingOptionsFromModelTR(opts.sdSwitch, modelTiltedRing, trFittingOptions);
  setupTRParamLimitsFromModelTR(
    opts.sdSwitch, f32(f32(2.0) * opts.beamMajorAxis), modelTiltedRing, trFittingOptions,
    {
      minChannel: minV, maxChannel: maxV,
      vRotLims: opts.vRotLims, linSDLims: opts.linSDLims, logSDLims: opts.logSDLims,
      pixelSize0: observedDC.dh.pixelSize[0], pixelSize1: observedDC.dh.pixelSize[1],
      channelSize: observedDC.dh.channelSize,
    }
  );
  const pvIni = setInitialParamVector(trFittingOptions);

  return { modelTiltedRing, trFittingOptions, pvIni, center, incl, pa, vSys, vDisp, nRings };
}

// ---------------------------------------------------------------------------
// initialAnalysis
// Fortran: InitialAnalysis(CatItem), preceded by InitialAnalysis_Prep()
// (matching PreGalaxyAnalysis.f's real call order: Prep runs BEFORE
// InitialAnalysis, since InitialAnalysis's ProjectedAnalysis needs the
// moment maps/velocity profile Prep builds).
// ---------------------------------------------------------------------------
function initialAnalysis(observedDC, maskDC, opts) {
  const { maskedObservedDC, observedMaps, observedVelocityProfile } =
    constructProjectionsFromCube(observedDC, maskDC);

  const noise = getNoise(observedDC);
  observedDC.dh.uncertainty = noise;

  const result = projectedAnalysis(observedDC, observedMaps, observedVelocityProfile, opts);

  return { ...result, noise, maskedObservedDC, observedMaps };
}

// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  constructProjectionsFromCube,
  getNoise,
  computeSNDiagnostics,
  getObjectCenter,
  getGalaxyShape,
  estimateFlatVDisp,
  convertFlatDiskProfilesToTR,
  projectedAnalysis,
  initialAnalysis,
};


// ---------------------------------------------------------------------------
// Self-test (node InitialAnalysis.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { DataCube, allocateDataCube } = require('../ObjectDefinitions/DataCube.js');

  console.log('=== getGalaxyShape ===');
  const { incl, pa } = getGalaxyShape(45.0, 30.0);
  console.log('  incl (expect 45deg in rad):', incl, f32(45 * Math.PI / 180));
  console.log('  pa   (expect (30+90)deg in rad):', pa, f32(120 * Math.PI / 180));

  console.log('\n=== estimateFlatVDisp ===');
  console.log('  VDisp (expect 10):', estimateFlatVDisp());

  console.log('\n=== full initialAnalysis on a synthetic rotating-disc cube ===');
  // A small cube with an actual spectral line per pixel (Gaussian in
  // channel), whose central velocity shifts linearly across x (a crude
  // rotating-disc velocity field) and whose peak flux declines with radius
  // (so outer rings fail modelability realistically, unlike a flat-top
  // disc with a spectrally-flat profile -- degenerate in a way no real
  // galaxy cube is, which is what a first version of this test used and
  // which hit an edge case in EstimateRadialProfiles' ring accounting; see
  // that file's module header for the "all rings modelable" bounds
  // question this surfaced).
  const n = 21, nCh = 30;
  const observedDC = new DataCube();
  observedDC.dh.nPixels[0] = n; observedDC.dh.nPixels[1] = n; observedDC.dh.nChannels = nCh;
  observedDC.dh.pixelSize[0] = f32(6.0); observedDC.dh.pixelSize[1] = f32(6.0);
  observedDC.dh.channelSize = f32(4.0);
  allocateDataCube(observedDC);
  const mid = Math.floor(n / 2);
  const Rmax = 9;
  const { flatIndxCalc } = require('../ObjectDefinitions/DataCube.js');
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const dx = i - mid, dy = j - mid;
      const r = Math.sqrt(dx * dx + dy * dy);
      if (r > Rmax) continue;
      const peak = Math.exp(-r / 4);       // radially-declining peak flux
      const vCent = mid + dx * 0.6;         // "rotation": center channel shifts with x
      for (let k = 0; k < nCh; k++) {
        const g = peak * Math.exp(-0.5 * ((k - vCent) / 2.5) ** 2);
        observedDC.flux[flatIndxCalc(i, j, k, observedDC.dh)] = f32(g);
      }
    }
  }

  const maskDC = new DataCube();
  maskDC.dh.nPixels[0] = n; maskDC.dh.nPixels[1] = n; maskDC.dh.nChannels = nCh;
  allocateDataCube(maskDC);
  maskDC.flux.fill(f32(1.0)); // trivial all-pass mask

  const opts = {
    centerSource: 0,
    catalogueEllipseIncDeg: 30.0,
    catalogueEllipsePADeg: 45.0,
    nRingsPerBeam: 2,
    nTargRings: -1,
    radGridArcsec: null,
    sdSwitch: 0,
    linSDLims: [f32(0.001), f32(10.0)],
    logSDLims: [f32(-3.0), f32(1.0)],
    vRotLims: [f32(0.0), f32(300.0)],
    sizeLims: [f32(1.0), f32(100.0)],
    noiseSigmaLim: f32(3.0),
    // Small relative to the disc (Rmax=9px): gives enough rings that the
    // outermost ones fall in the zero-flux region beyond the disc and fail
    // modelability naturally -- avoids the "all rings modelable" bounds
    // edge case documented in EstimateRadialProfiles.js's header (hit with
    // a larger beam during development of this test).
    beamMajorAxis: f32(3.0),
  };

  try {
    const result = initialAnalysis(observedDC, maskDC, opts);
    console.log('  noise (RMS from corner sub-boxes, all-inside-R disc, expect ~0):', result.noise);
    console.log('  center:', Array.from(result.center), `(expect ~[${mid},${mid}])`);
    console.log('  nRings:', result.nRings);
    console.log('  pvIni.nParams:', result.pvIni.nParams);
    console.log('  RESULT: pipeline ran end-to-end without throwing.');
  } catch (e) {
    console.log('  RESULT: threw ->', e.message);
    console.log('  (a throw here on synthetic/degenerate test data is not necessarily wrong --');
    console.log('   e.g. FillInMisingRings\' ring-0 bug can legitimately fire on a toy cube.)');
  }
}
