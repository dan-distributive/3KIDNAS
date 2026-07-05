'use strict';

// =============================================================================
// TiltedRing.js
// High-fidelity port of src/ObjectDefinitions/TiltedRing.f (TiltedRingMod)
//
// PORTING NOTES
// -------------
// Fortran `real` → Float32Array / Math.fround().
// Fortran `logical` → JS boolean.
// Fortran `integer` → plain JS number (integer).
//
// Ring%P is an allocatable array of Particle — in JS this is a plain JS array
// of Particle objects, allocated by ring_ParticleAllocation().
//
// TiltedRingFittingOptions parameter index map (0:12):
//   0=Xcent  1=Ycent  2=Inclination  3=PA  4=VSys
//   5=VRot   6=VRad   7=VDisp        8=Vvert  9=dvdz
//   10=Sigma  11=z0   12=zGradiantStart
// =============================================================================

const f32 = Math.fround;
const { Particle } = require('./Particle.js');


// ---------------------------------------------------------------------------
// Ring
// Fortran TYPE Ring — a single tilted ring with its kinematic parameters
// and an array of model particles.
// ---------------------------------------------------------------------------
class Ring {
  constructor() {
    this.nParticles      = 0;
    this.rmid            = f32(0.0);
    this.rwidth          = f32(0.0);
    this.centPos         = new Float32Array(2);   // (0:1)
    this.inclination     = f32(0.0);
    this.positionAngle   = f32(0.0);
    this.vSys            = f32(0.0);
    this.vRot            = f32(0.0);
    this.vRad            = f32(0.0);
    this.vDisp           = f32(0.0);
    this.vvert           = f32(0.0);
    this.dvdz            = f32(0.0);
    this.sigma           = f32(0.0);
    this.logSigma        = f32(0.0);
    this.sigUse          = f32(0.0);
    this.z0              = f32(0.0);
    this.zGradiantStart  = f32(0.0);
    this.p               = null;   // Array of Particle, length nParticles
  }
}


// ---------------------------------------------------------------------------
// TiltedRingModel
// Fortran TYPE TiltedRingModel — array of rings + global model parameters.
// ---------------------------------------------------------------------------
class TiltedRingModel {
  constructor() {
    this.nRings            = 0;
    this.cmode             = 0;
    this.cloudBaseSurfDens = f32(0.0);
    this.r                 = null;   // Array of Ring, length nRings
    this.unitSwitchs       = new Int32Array(4);   // (0:3)
    this.sd_Switch         = 0;
  }
}


// ---------------------------------------------------------------------------
// TiltedRingFittingOptions
// Fortran TYPE TiltedRingFittingOptions — controls which parameters are
// free, fixed, or constant across rings.
// ---------------------------------------------------------------------------
class TiltedRingFittingOptions {
  constructor() {
    this.nRingsPerBeam          = 0;
    this.nRings                 = 0;
    this.nFittedParamsTotal     = 0;
    this.nTargRings             = 0;
    this.nFittedRadialParams    = 0;
    this.nFixedRadialParams     = 0;
    this.nFittedConstantParams  = 0;
    this.nFixedConstantParams   = 0;
    // Logical arrays (0:12) — false = fitted radial, true = constant
    this.constParams     = new Array(13).fill(false);
    this.fixedParams     = new Array(13).fill(false);
    // Parameter bounds and ranges (0:12)
    this.paramLowerLims  = new Float32Array(13);
    this.paramUpperLims  = new Float32Array(13);
    this.paramRange      = new Float32Array(13);
    this.cyclicSwitch    = new Int32Array(13);
    // Per-ring parameter profiles
    this.radialProfiles  = null;   // Array of Ring, length nRings
  }
}


// ---------------------------------------------------------------------------
// ring_ParticleAllocation
// Fortran: Ring_ParticleAllocation(R)
// Allocates R%P(0:nParticles-1) — a 0-indexed array of Particle objects.
// ---------------------------------------------------------------------------
function ring_ParticleAllocation(r) {
  r.p = new Array(r.nParticles);
  for (let i = 0; i < r.nParticles; i++) {
    r.p[i] = new Particle();
  }
}


// ---------------------------------------------------------------------------
// ring_ParticleDeAllocation
// Fortran: Ring_ParticleDeAllocation(R)
// ---------------------------------------------------------------------------
function ring_ParticleDeAllocation(r) {
  r.p = null;
}


// ---------------------------------------------------------------------------
// tiltRing_Allocate
// Fortran: TiltRing_Allocate(TR)
// Allocates TR%R(0:nRings-1).
// ---------------------------------------------------------------------------
function tiltRing_Allocate(tr) {
  tr.r = new Array(tr.nRings);
  for (let i = 0; i < tr.nRings; i++) {
    tr.r[i] = new Ring();
  }
}


// ---------------------------------------------------------------------------
// tiltRing_DeAllocate
// Fortran: TiltRing_DeAllocate(TR)
// Deallocates particles in each ring, then the ring array.
// ---------------------------------------------------------------------------
function tiltRing_DeAllocate(tr) {
  if (tr.r) {
    for (let i = 0; i < tr.nRings; i++) {
      ring_ParticleDeAllocation(tr.r[i]);
    }
  }
  tr.r = null;
}


// ---------------------------------------------------------------------------
// tiltRing_DeAllocateStruct
// Fortran: TiltRing_DeAllocateStruct(TR)
// Deallocates ring array only (no particle deallocation).
// ---------------------------------------------------------------------------
function tiltRing_DeAllocateStruct(tr) {
  tr.r = null;
}


// ---------------------------------------------------------------------------
// tiltRingFittingOptions_Allocate
// Fortran: TiltRingFittingOptions_Allocate(TRFO)
// Allocates radialProfiles(0:nRings-1) and initializes all kinematic fields
// to 0.
// ---------------------------------------------------------------------------
function tiltRingFittingOptions_Allocate(trfo) {
  trfo.radialProfiles = new Array(trfo.nRings);
  for (let i = 0; i < trfo.nRings; i++) {
    const r           = new Ring();
    r.centPos[0]      = f32(0.0);
    r.centPos[1]      = f32(0.0);
    r.inclination     = f32(0.0);
    r.positionAngle   = f32(0.0);
    r.vSys            = f32(0.0);
    r.vRot            = f32(0.0);
    r.vRad            = f32(0.0);
    r.vDisp           = f32(0.0);
    r.vvert           = f32(0.0);
    r.dvdz            = f32(0.0);
    r.sigma           = f32(0.0);
    r.z0              = f32(0.0);
    r.zGradiantStart  = f32(0.0);
    trfo.radialProfiles[i] = r;
  }
}


// ---------------------------------------------------------------------------
// tiltRingFittingOptions_DeAllocate
// Fortran: TiltRingFittingOptions_DeAllocate(TRFO)
// ---------------------------------------------------------------------------
function tiltRingFittingOptions_DeAllocate(trfo) {
  trfo.radialProfiles = null;
}


// ---------------------------------------------------------------------------
// logicalTiltedRingIndexing
// Fortran: LogicalTiltedRingIndexing(TRFO)
//
// Counts free/fixed/constant parameters from constParams and fixedParams
// flags, then computes nFittedParamsTotal.
//
// Parameter type matrix:
//   constParams=T, fixedParams=T → nFixedConstantParams++
//   constParams=T, fixedParams=F → nFittedConstantParams++  (1 slot in PV)
//   constParams=F, fixedParams=T → nFixedRadialParams++
//   constParams=F, fixedParams=F → nFittedRadialParams++    (nRings slots in PV)
//
// nFittedParamsTotal = nFittedConstantParams + nFittedRadialParams * nRings
// ---------------------------------------------------------------------------
function logicalTiltedRingIndexing(trfo) {
  trfo.nFittedConstantParams = 0;
  trfo.nFittedRadialParams   = 0;
  trfo.nFixedConstantParams  = 0;
  trfo.nFixedRadialParams    = 0;

  for (let i = 0; i <= 12; i++) {
    if (trfo.constParams[i]) {
      if (trfo.fixedParams[i]) {
        trfo.nFixedConstantParams++;
      } else {
        trfo.nFittedConstantParams++;
      }
    } else {
      if (trfo.fixedParams[i]) {
        trfo.nFixedRadialParams++;
      } else {
        trfo.nFittedRadialParams++;
      }
    }
  }

  trfo.nFittedParamsTotal = trfo.nFittedConstantParams
    + trfo.nFittedRadialParams * trfo.nRings;
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  Ring,
  TiltedRingModel,
  TiltedRingFittingOptions,
  ring_ParticleAllocation,
  ring_ParticleDeAllocation,
  tiltRing_Allocate,
  tiltRing_DeAllocate,
  tiltRing_DeAllocateStruct,
  tiltRingFittingOptions_Allocate,
  tiltRingFittingOptions_DeAllocate,
  logicalTiltedRingIndexing
};


// ---------------------------------------------------------------------------
// Self-test (node TiltedRing.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const {
    Ring, TiltedRingModel, TiltedRingFittingOptions,
    ring_ParticleAllocation, ring_ParticleDeAllocation,
    tiltRing_Allocate, tiltRing_DeAllocate,
    tiltRingFittingOptions_Allocate, logicalTiltedRingIndexing
  } = module.exports;

  // --- Ring particle allocation ---
  console.log('=== Ring particle allocation ===');
  const r = new Ring();
  r.nParticles = 5;
  ring_ParticleAllocation(r);
  console.log('nParticles:', r.nParticles, '(expect 5)');
  console.log('p.length:', r.p.length, '(expect 5)');
  console.log('p[0] is Particle:', r.p[0] instanceof Particle, '(expect true)');
  r.p[0].flux = f32(1.5);
  console.log('p[0].flux:', r.p[0].flux, '(expect 1.5)');
  ring_ParticleDeAllocation(r);
  console.log('p after dealloc:', r.p, '(expect null)');

  // --- TiltedRingModel allocation ---
  console.log('\n=== TiltedRingModel allocation ===');
  const tr = new TiltedRingModel();
  tr.nRings = 3;
  tiltRing_Allocate(tr);
  console.log('r.length:', tr.r.length, '(expect 3)');
  tr.r[0].rmid = f32(10.0);
  tr.r[2].vRot = f32(150.0);
  console.log('r[0].rmid:', tr.r[0].rmid, '(expect 10)');
  console.log('r[2].vRot:', tr.r[2].vRot, '(expect 150)');
  tiltRing_DeAllocate(tr);
  console.log('r after dealloc:', tr.r, '(expect null)');

  // --- LogicalTiltedRingIndexing ---
  console.log('\n=== logicalTiltedRingIndexing ===');
  const trfo = new TiltedRingFittingOptions();
  trfo.nRings = 4;
  tiltRingFittingOptions_Allocate(trfo);

  // All 13 params free and radial (default) → nFittedParamsTotal = 13*4 = 52
  logicalTiltedRingIndexing(trfo);
  console.log('all radial free: nFittedParamsTotal =', trfo.nFittedParamsTotal, '(expect 52)');
  console.log('nFittedRadialParams:', trfo.nFittedRadialParams, '(expect 13)');

  // Make params 0-4 constant (one slot each), rest radial
  for (let i = 0; i <= 4; i++) trfo.constParams[i] = true;
  logicalTiltedRingIndexing(trfo);
  // 5 const + 8 radial*4 = 5 + 32 = 37
  console.log('5 const, 8 radial: nFittedParamsTotal =', trfo.nFittedParamsTotal, '(expect 37)');

  // Fix param 5 (VRot) — should not appear in PV
  trfo.fixedParams[5] = true;
  logicalTiltedRingIndexing(trfo);
  // 5 const + 7 radial*4 = 5 + 28 = 33
  console.log('5 const, 1 fixed radial: nFittedParamsTotal =', trfo.nFittedParamsTotal, '(expect 33)');

  // Fix param 0 (const+fixed) — should not appear in PV
  trfo.fixedParams[0] = true;
  logicalTiltedRingIndexing(trfo);
  // 4 const + 7 radial*4 = 4 + 28 = 32
  console.log('4 const, 1 fixed const, 1 fixed radial: nFittedParamsTotal =', trfo.nFittedParamsTotal, '(expect 32)');
}
