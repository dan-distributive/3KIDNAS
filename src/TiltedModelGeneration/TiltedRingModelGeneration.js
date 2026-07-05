'use strict';

// =============================================================================
// TiltedRingModelGeneration.js
// High-fidelity port of src/TiltedRingModelGeneration/TiltedRingModelGeneration.f
// (TiltedRingGenerationMod)
//
// PORTING NOTES
// -------------
// Thin wrapper — all physics delegated to SingleRingGeneration.js.
// rng: makeRng() object from random.js, shared across all rings.
// =============================================================================

const { ring_ParticleAllocation } = require('../ObjectDefinitions/TiltedRing.js');
const {
  ring_CalcNumParticles,
  ring_ParticleGeneration
} = require('./SingleRingGeneration.js');


// ---------------------------------------------------------------------------
// buildTiltedRingModel
// Fortran: BuildTiltedRingModel(TR, idum)
//
// For each ring in TR:
//   1. Calculate number of particles
//   2. Allocate particle array
//   3. Generate particle positions, velocities, and fluxes
//
// TR must already be allocated (tiltRing_Allocate called) and all ring
// parameters set before calling this.
//
// rng: makeRng() object from random.js
// ---------------------------------------------------------------------------
function buildTiltedRingModel(tr, rng) {
  for (let i = 0; i < tr.nRings; i++) {
    ring_CalcNumParticles(tr.r[i], tr.cmode, tr.cloudBaseSurfDens);
    ring_ParticleAllocation(tr.r[i]);
    ring_ParticleGeneration(tr.r[i], rng);
  }
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = { buildTiltedRingModel };


// ---------------------------------------------------------------------------
// Self-test (node TiltedRingModelGeneration.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { makeRng }          = require('../StandardMath/random.js');
  const { TiltedRingModel, tiltRing_Allocate, tiltRing_DeAllocate }
                             = require('../ObjectDefinitions/TiltedRing.js');
  const f32 = Math.fround;
  const Pi  = f32(Math.PI);

  const tr               = new TiltedRingModel();
  tr.nRings              = 3;
  tr.cmode               = 1;
  tr.cloudBaseSurfDens   = f32(1.0);
  tiltRing_Allocate(tr);

  // Set up three rings with increasing radii
  for (let i = 0; i < 3; i++) {
    const r           = tr.r[i];
    r.rmid            = f32((i + 1) * 5.0);
    r.rwidth          = f32(2.0);
    r.sigma           = f32(0.01);
    r.inclination     = f32(45.0 * Math.PI / 180.0);
    r.positionAngle   = f32(30.0 * Math.PI / 180.0);
    r.centPos[0]      = f32(50.0);
    r.centPos[1]      = f32(50.0);
    r.vSys            = f32(1000.0);
    r.vRot            = f32(150.0);
    r.vRad            = f32(0.0);
    r.vDisp           = f32(8.0);
    r.vvert           = f32(0.0);
    r.dvdz            = f32(0.0);
    r.z0              = f32(0.0);
    r.zGradiantStart  = f32(0.0);
  }

  const rng = makeRng(-1);
  buildTiltedRingModel(tr, rng);

  console.log('=== buildTiltedRingModel ===');
  for (let i = 0; i < 3; i++) {
    const r = tr.r[i];
    let totalFlux = f32(0.0);
    for (let j = 0; j < r.nParticles; j++) totalFlux = f32(totalFlux + r.p[j].flux);
    const rmin = f32(r.rmid - r.rwidth / 2);
    const rmax = f32(r.rmid + r.rwidth / 2);
    const expectedFlux = f32(r.sigma * f32(Pi * f32(rmax * rmax - rmin * rmin)));
    console.log(`ring[${i}]: rmid=${r.rmid} nParticles=${r.nParticles} totalFlux=${totalFlux.toExponential(4)} expectedFlux=${expectedFlux.toExponential(4)} match=${Math.abs(totalFlux - expectedFlux) < 1e-5 ? 'OK' : 'FAIL'}`);
  }

  tiltRing_DeAllocate(tr);
  console.log('dealloc OK:', tr.r === null);
}
