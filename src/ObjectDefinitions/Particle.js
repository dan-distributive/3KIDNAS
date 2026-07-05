'use strict';

// =============================================================================
// Particle.js
// High-fidelity port of src/ObjectDefinitions/Particle.f (ParticleMod)
//
// PORTING NOTES
// -------------
// Fortran `real` → stored as Math.fround() values via Float32Array.
// All arrays are 0-indexed matching Fortran (0:2).
// No methods in the Fortran module — pure struct definition.
// =============================================================================

const f32 = Math.fround;


// ---------------------------------------------------------------------------
// Particle
// Fortran TYPE Particle — a single model particle in a tilted ring.
//
// Fields:
//   flux           — particle flux weight (Jy·pixel·channel)
//   pos[3]         — Cartesian position (x, y, z) in ring plane
//   vel[3]         — velocity vector (not used in standard fit)
//   angPos[3]      — cylindrical coordinates (R, theta, z)
//   projectedPos[3] — sky-plane position after inclination + PA rotation
//                     [0]=x_pixel, [1]=y_pixel (used for voxel lookup)
//   projectedVel[3] — projected velocity after inclination + PA rotation
//                     [2] = line-of-sight velocity (used for channel lookup)
// ---------------------------------------------------------------------------
class Particle {
  constructor() {
    this.flux          = f32(0.0);
    this.pos           = new Float32Array(3);   // (0:2)
    this.vel           = new Float32Array(3);   // (0:2)
    this.angPos        = new Float32Array(3);   // (0:2) R, theta, z
    this.projectedPos  = new Float32Array(3);   // (0:2)
    this.projectedVel  = new Float32Array(3);   // (0:2)
  }
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = { Particle };


// ---------------------------------------------------------------------------
// Self-test (node Particle.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { Particle } = module.exports;

  const p = new Particle();
  console.log('=== Particle defaults ===');
  console.log('flux:', p.flux, '(expect 0)');
  console.log('pos:', Array.from(p.pos), '(expect [0,0,0])');
  console.log('angPos:', Array.from(p.angPos), '(expect [0,0,0])');
  console.log('projectedPos:', Array.from(p.projectedPos), '(expect [0,0,0])');
  console.log('projectedVel:', Array.from(p.projectedVel), '(expect [0,0,0])');

  p.flux            = f32(0.005);
  p.angPos[0]       = f32(3.5);    // R
  p.angPos[1]       = f32(1.2);    // theta
  p.angPos[2]       = f32(0.1);    // z
  p.projectedPos[0] = f32(42.3);
  p.projectedPos[1] = f32(17.8);
  p.projectedVel[2] = f32(1250.0);

  console.log('\n=== after assignment ===');
  console.log('flux:', p.flux);
  console.log('angPos:', Array.from(p.angPos));
  console.log('projectedPos[0:1]:', p.projectedPos[0], p.projectedPos[1]);
  console.log('projectedVel[2]:', p.projectedVel[2]);
}
