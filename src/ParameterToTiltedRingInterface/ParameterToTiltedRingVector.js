'use strict';

// =============================================================================
// ParameterToTiltedRingVector.js
// High-fidelity port of
//   src/ParameterToTiltedRingInterface/ParameterToTiltedRingVector.f
// (ParameterVectorToTiltedRingMod)
//
// PORTING NOTES
// -------------
// This module is the serializer/deserializer between the flat parameter
// vector used by the optimizer and the structured TiltedRing objects.
//
// Parameter index map (0:12):
//   0=CentPos[0]  1=CentPos[1]  2=Inclination  3=PositionAngle
//   4=VSys        5=VRot        6=VRad          7=VDisp
//   8=Vvert       9=dvdz        10=SigUse       11=z0
//   12=zGradiantStart
//
// constParams/fixedParams flag matrix (per parameter i):
//   const=T, fixed=T → not in PV; all rings get trfo.radialProfiles[0].field
//   const=T, fixed=F → 1 slot in PV; all rings get that one value
//   const=F, fixed=T → not in PV; each ring gets trfo.radialProfiles[r].field
//   const=F, fixed=F → nRings slots in PV; each ring gets its own value
//
// CRITICAL: SetParamLimsFromFittingOptions does NOT advance CurrParam.
// Only SetParamFromFittingOptions advances CurrParam. The four setters
// (lowerLims, upperLims, paramRange, cyclics, then values) share the same
// CurrParam and only the final values call advances it.
// =============================================================================

const f32 = Math.fround;


// ---------------------------------------------------------------------------
// getParamField
// Returns a getter/setter pair for parameter index i on a Ring object.
// Maps the 13 parameter indices to Ring field names.
// ---------------------------------------------------------------------------
function getParamField(i) {
  switch (i) {
    case  0: return { get: r => r.centPos[0],       set: (r,v) => { r.centPos[0]      = f32(v); } };
    case  1: return { get: r => r.centPos[1],       set: (r,v) => { r.centPos[1]      = f32(v); } };
    case  2: return { get: r => r.inclination,      set: (r,v) => { r.inclination     = f32(v); } };
    case  3: return { get: r => r.positionAngle,    set: (r,v) => { r.positionAngle   = f32(v); } };
    case  4: return { get: r => r.vSys,             set: (r,v) => { r.vSys            = f32(v); } };
    case  5: return { get: r => r.vRot,             set: (r,v) => { r.vRot            = f32(v); } };
    case  6: return { get: r => r.vRad,             set: (r,v) => { r.vRad            = f32(v); } };
    case  7: return { get: r => r.vDisp,            set: (r,v) => { r.vDisp           = f32(v); } };
    case  8: return { get: r => r.vvert,            set: (r,v) => { r.vvert           = f32(v); } };
    case  9: return { get: r => r.dvdz,             set: (r,v) => { r.dvdz            = f32(v); } };
    case 10: return { get: r => r.sigUse,           set: (r,v) => { r.sigUse          = f32(v); } };
    case 11: return { get: r => r.z0,               set: (r,v) => { r.z0              = f32(v); } };
    case 12: return { get: r => r.zGradiantStart,   set: (r,v) => { r.zGradiantStart  = f32(v); } };
  }
}


// ---------------------------------------------------------------------------
// setSpecificVector
// Fortran: SetSpecificVector(nParam, nRings, Const, Fixed, CurrParam,
//                            Param, AlternateSourceVector, TargVec)
//
// Sets ring fields (TargVec) from either the PV or the alternate source,
// based on const/fixed flags. Advances currParam when reading from PV.
//
// alternateSource: array of nRings values from trfo.radialProfiles
// target:          array of nRings ring objects to write to
// field:           { get, set } pair from getParamField
// Returns updated currParam.
// ---------------------------------------------------------------------------
function setSpecificVector(pv, nRings, isConst, isFixed, currParam,
                           alternateSource, target, field) {
  if (isConst) {
    if (isFixed) {
      // Broadcast trfo value[0] to all rings — not in PV
      const val = alternateSource[0];
      for (let r = 0; r < nRings; r++) field.set(target[r], val);
    } else {
      // One PV slot → broadcast to all rings, advance currParam
      const val = pv.param[currParam];
      for (let r = 0; r < nRings; r++) field.set(target[r], val);
      currParam++;
    }
  } else {
    if (isFixed) {
      // Per-ring trfo values → not in PV
      for (let r = 0; r < nRings; r++) field.set(target[r], alternateSource[r]);
    } else {
      // nRings PV slots → one per ring, advance currParam by nRings
      for (let r = 0; r < nRings; r++) {
        field.set(target[r], pv.param[currParam + r]);
      }
      currParam += nRings;
    }
  }
  return currParam;
}


// ---------------------------------------------------------------------------
// generalizedParamVectorToTiltedRing
// Fortran: GeneralizedParamVectorToTiltedRing(PV, TR, TRFO)
//
// Deserializes the flat PV%Param array into the TR ring structs,
// using TRFO constParams/fixedParams flags to determine layout.
// Called on every chi² evaluation.
// ---------------------------------------------------------------------------
function generalizedParamVectorToTiltedRing(pv, tr, trfo) {
  let currParam = 0;
  const nRings  = tr.nRings;

  for (let i = 0; i <= 12; i++) {
    const field          = getParamField(i);
    const isConst        = trfo.constParams[i];
    const isFixed        = trfo.fixedParams[i];
    const alternateSource = trfo.radialProfiles.map(r => field.get(r));

    currParam = setSpecificVector(
      pv, nRings, isConst, isFixed, currParam,
      alternateSource, tr.r, field
    );
  }
}


// ---------------------------------------------------------------------------
// setParamLimsFromFittingOptions
// Fortran: SetParamLimsFromFittingOptions(ParamNum, nParam, nRings,
//                                          Const, Fixed, CurrParam,
//                                          ParamLims, TargLims)
//
// Copies limit value from trfo.paramLowerLims/Upper/Range[paramNum]
// into pv array at CurrParam position(s).
// DOES NOT advance currParam — that is done only by setParamFromFittingOptions.
// ---------------------------------------------------------------------------
function setParamLimsFromFittingOptions(paramNum, nRings, isConst, isFixed,
                                         currParam, paramLims, targLims) {
  if (!isFixed) {
    if (isConst) {
      paramLims[currParam] = f32(targLims[paramNum]);
    } else {
      for (let r = 0; r < nRings; r++) {
        paramLims[currParam + r] = f32(targLims[paramNum]);
      }
    }
  }
}


// ---------------------------------------------------------------------------
// setParamCyclicsFromFittingOptions
// Fortran: SetParamCyclicsFromFittingOptions(...)
// Same pattern as setParamLimsFromFittingOptions but for integer cyclic array.
// ---------------------------------------------------------------------------
function setParamCyclicsFromFittingOptions(paramNum, nRings, isConst, isFixed,
                                            currParam, paramCyclics, targCyclics) {
  if (!isFixed) {
    if (isConst) {
      paramCyclics[currParam] = targCyclics[paramNum];
    } else {
      for (let r = 0; r < nRings; r++) {
        paramCyclics[currParam + r] = targCyclics[paramNum];
      }
    }
  }
}


// ---------------------------------------------------------------------------
// setParamFromFittingOptions
// Fortran: SetParamFromFittingOptions(nParam, nRings, Const, Fixed,
//                                      CurrParam, Param, TargVec)
//
// Copies trfo radialProfiles values into pv.param.
// THIS is the only function that advances currParam.
// Returns updated currParam.
// ---------------------------------------------------------------------------
function setParamFromFittingOptions(nRings, isConst, isFixed,
                                     currParam, param, targVec) {
  if (!isFixed) {
    if (isConst) {
      param[currParam] = f32(targVec[0]);
      currParam++;
    } else {
      for (let r = 0; r < nRings; r++) {
        param[currParam + r] = f32(targVec[r]);
      }
      currParam += nRings;
    }
  }
  return currParam;
}


// ---------------------------------------------------------------------------
// tiltedRingOptionsToPV
// Fortran: TiltedRingOptionsToPV(PV, TRFO)
//
// Serializes TRFO (radialProfiles + limits) into PV flat arrays.
// Called once before fitting to initialize PVIni.
//
// Note: all four setters (lowerLims, upperLims, paramRange, cyclics)
// share the same currParam value but only setParamFromFittingOptions
// advances it. The lims/cyclics setters read currParam but don't modify it.
// ---------------------------------------------------------------------------
function tiltedRingOptionsToPV(pv, trfo) {
  let currParam = 0;
  const nRings  = trfo.nRings;

  for (let i = 0; i <= 12; i++) {
    const field   = getParamField(i);
    const isConst = trfo.constParams[i];
    const isFixed = trfo.fixedParams[i];

    // Lower limits
    setParamLimsFromFittingOptions(i, nRings, isConst, isFixed,
      currParam, pv.paramLowerLims, trfo.paramLowerLims);

    // Upper limits
    setParamLimsFromFittingOptions(i, nRings, isConst, isFixed,
      currParam, pv.paramUpperLims, trfo.paramUpperLims);

    // Param ranges
    setParamLimsFromFittingOptions(i, nRings, isConst, isFixed,
      currParam, pv.paramRange, trfo.paramRange);

    // Cyclic switches
    setParamCyclicsFromFittingOptions(i, nRings, isConst, isFixed,
      currParam, pv.cyclicSwitch, trfo.cyclicSwitch);

    // Values — this one advances currParam
    const targVec = trfo.radialProfiles.map(r => field.get(r));
    currParam = setParamFromFittingOptions(
      nRings, isConst, isFixed, currParam, pv.param, targVec
    );
  }
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  generalizedParamVectorToTiltedRing,
  tiltedRingOptionsToPV,
  setSpecificVector,
  setParamLimsFromFittingOptions,
  setParamFromFittingOptions,
  setParamCyclicsFromFittingOptions
};


// ---------------------------------------------------------------------------
// Self-test (node ParameterToTiltedRingVector.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { TiltedRingModel, TiltedRingFittingOptions,
          tiltRing_Allocate, tiltRingFittingOptions_Allocate,
          logicalTiltedRingIndexing }
    = require('../ObjectDefinitions/TiltedRing.js');
  const { ParameterVector, allocateParamVector }
    = require('../ObjectDefinitions/ParameterVector.js');

  // Build a 2-ring TRFO with all params constant (simplest case)
  const trfo        = new TiltedRingFittingOptions();
  trfo.nRings       = 2;
  trfo.nRingsPerBeam = 2;
  trfo.nTargRings   = -1;

  // Make some params constant, some radial, one fixed
  for (let i = 0; i <= 12; i++) {
    trfo.constParams[i] = (i <= 4);   // 0-4 constant, 5-12 radial
    trfo.fixedParams[i] = (i === 6 || i === 8 || i === 9 || i === 11 || i === 12);
  }
  tiltRingFittingOptions_Allocate(trfo);
  logicalTiltedRingIndexing(trfo);

  // Set radial profile values
  for (let r = 0; r < 2; r++) {
    const rp = trfo.radialProfiles[r];
    rp.centPos[0]      = f32(32.0);
    rp.centPos[1]      = f32(32.0);
    rp.inclination     = f32(0.785);   // ~45 deg
    rp.positionAngle   = f32(0.524);   // ~30 deg
    rp.vSys            = f32(1000.0);
    rp.vRot            = f32(r === 0 ? 100.0 : 150.0);  // radial: different per ring
    rp.vRad            = f32(0.0);
    rp.vDisp           = f32(r === 0 ? 8.0 : 10.0);
    rp.vvert           = f32(0.0);
    rp.dvdz            = f32(0.0);
    rp.sigUse          = f32(r === 0 ? 0.01 : 0.02);
    rp.z0              = f32(0.0);
    rp.zGradiantStart  = f32(0.0);
  }

  // Set limits
  trfo.paramLowerLims[0] = f32(0.0);   trfo.paramUpperLims[0] = f32(64.0);
  trfo.paramLowerLims[1] = f32(0.0);   trfo.paramUpperLims[1] = f32(64.0);
  trfo.paramLowerLims[2] = f32(0.0);   trfo.paramUpperLims[2] = f32(1.57);
  trfo.paramLowerLims[3] = f32(0.0);   trfo.paramUpperLims[3] = f32(6.28);
  trfo.paramLowerLims[4] = f32(900.0); trfo.paramUpperLims[4] = f32(1100.0);
  trfo.paramLowerLims[5] = f32(0.0);   trfo.paramUpperLims[5] = f32(300.0);
  trfo.paramLowerLims[7] = f32(0.0);   trfo.paramUpperLims[7] = f32(20.0);
  trfo.paramLowerLims[10]= f32(0.0);   trfo.paramUpperLims[10]= f32(1.0);
  trfo.cyclicSwitch[3]   = 1;  // PA is cyclic

  console.log('=== tiltedRingOptionsToPV ===');
  console.log('nFittedParamsTotal:', trfo.nFittedParamsTotal);

  const pv    = new ParameterVector();
  pv.nParams  = trfo.nFittedParamsTotal;
  allocateParamVector(pv);
  tiltedRingOptionsToPV(pv, trfo);

  console.log('PV params:', Array.from(pv.param).map(v => v.toFixed(4)));
  console.log('PV lowerLims:', Array.from(pv.paramLowerLims).map(v => v.toFixed(2)));
  console.log('PV cyclicSwitch:', Array.from(pv.cyclicSwitch));

  // Now test roundtrip: PV → TR → check ring fields match
  console.log('\n=== generalizedParamVectorToTiltedRing roundtrip ===');
  const tr = new TiltedRingModel();
  tr.nRings = 2; tr.cmode = 1; tr.cloudBaseSurfDens = f32(1.0);
  tiltRing_Allocate(tr);

  generalizedParamVectorToTiltedRing(pv, tr, trfo);

  for (let r = 0; r < 2; r++) {
    const ring = tr.r[r];
    const prof = trfo.radialProfiles[r];
    console.log(`ring[${r}]:`);
    console.log(`  centPos:    [${ring.centPos[0].toFixed(2)}, ${ring.centPos[1].toFixed(2)}] (expect [32.00, 32.00])`);
    console.log(`  inclination: ${ring.inclination.toFixed(4)} (expect 0.7850)`);
    console.log(`  vRot:        ${ring.vRot.toFixed(2)} (expect ${prof.vRot.toFixed(2)})`);
    console.log(`  vDisp:       ${ring.vDisp.toFixed(2)} (expect ${prof.vDisp.toFixed(2)})`);
    console.log(`  sigUse:      ${ring.sigUse.toFixed(4)} (expect ${prof.sigUse.toFixed(4)})`);
  }
}
