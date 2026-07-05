'use strict';

// =============================================================================
// SingleRingGeneration.js
// High-fidelity port of src/TiltedRingModelGeneration/SingleRingGeneration.f
// (SingleRingGenerationMod)
//
// PORTING NOTES
// -------------
// Fortran `real` → Math.fround() at each arithmetic step.
// Trig functions (sin, cos, atanh, sqrt) operate on f32 inputs and their
// outputs are wrapped in f32 to match Fortran `real` precision.
//
// RNG: Fortran threads a single `idum` integer through ran2() and gasdev().
// In JS, pass the `rng` object from makeRng() (random.js) — rng.ran2() and
// rng.gasdev() share state exactly as Fortran's idum does.
//
// KNOWN ISSUE — delZ uninitialized in Ring_CalcParticle_VSys:
// The Fortran source assigns delZ in a commented-out line:
//   c      delZ=abs(R%P(PartID)%AngPos(2))-R%zGradiantStart
// but the branch that uses delZ:
//   if(abs(...AngPos(2)) .gt. R%zGradiantStart) VRotP=VRotP-R%dvdz*delZ
// references it uninitialized. In practice dvdz=0 in the standard fit so
// this branch has no effect. Ported exactly as-is with delZ=0 as the JS
// default (matching typical Fortran stack behaviour for local reals).
//
// Pi from CommonConsts: f32(Math.PI)
// =============================================================================

const f32  = Math.fround;
const Pi   = f32(Math.PI);

const { ring_ParticleAllocation } = require('../ObjectDefinitions/TiltedRing.js');


// ---------------------------------------------------------------------------
// ring_CalcNumParticles
// Fortran: Ring_CalcNumParticles(R, cmode, CloudSurfDens)
//
// Calculates the number of particles needed to represent the ring,
// proportional to ring area and surface density.
//
// nParticles = int(CloudSurfDens * Sigma^cmode * Pi*(Rh^2 - Rl^2)) + 1
// ---------------------------------------------------------------------------
function ring_CalcNumParticles(r, cmode, cloudSurfDens) {
  const rl          = f32(f32(r.rmid) - f32(f32(r.rwidth) / f32(2.0)));
  const rh          = f32(f32(r.rmid) + f32(f32(r.rwidth) / f32(2.0)));
  const pixelRing   = f32(Pi * f32(f32(rh * rh) - f32(rl * rl)));
  const densMulti   = f32(f32(cloudSurfDens) * f32(r.sigma ** cmode));
  r.nParticles      = Math.trunc(f32(densMulti * pixelRing)) + 1;
}


// ---------------------------------------------------------------------------
// ring_ParticleGeneration
// Fortran: Ring_ParticleGeneration(R, idum)
//
// Generates all particles for a ring:
//   1. Select random position (cylindrical coords + Cartesian)
//   2. Project to sky plane (inclination + PA rotation)
//   3. Shift to galaxy center
//   4. Calculate line-of-sight velocity
//   5. Assign flux weight
//
// rng: makeRng() object from random.js (shared ran2/gasdev state)
// ---------------------------------------------------------------------------
function ring_ParticleGeneration(r, rng) {
  const rmin = f32(f32(r.rmid) - f32(f32(r.rwidth) / f32(2.0)));
  const rmax = f32(f32(r.rmid) + f32(f32(r.rwidth) / f32(2.0)));
  const area = f32(Pi * f32(f32(rmax * rmax) - f32(rmin * rmin)));

  for (let i = 0; i < r.nParticles; i++) {
    ring_ParticlePosSelect(r, rng, rmin, rmax, i);
    particlePosProject(r.p[i], r.inclination, r.positionAngle);
    particlePos_NewCenter(r.p[i], r.centPos);
    ring_CalcParticle_VSys(r, i, rng);
    ring_CalcParticleFlux_Basic(r, i, area);
  }
}


// ---------------------------------------------------------------------------
// ring_ParticlePosSelect
// Fortran: Ring_ParticlePosSelect(R, idum, Rmin, Rmax, PartID)
//
// Selects a random position for particle PartID:
//   RR    = sqrt(ran2*(Rmax^2 - Rmin^2) + Rmin^2)  equal-area radius
//   Theta = ran2 * 2*Pi                              uniform angle
//   Z     = atanh(2*ran2 - 1) * z0                  sech^2 height
//
// Stores cylindrical (RR, Theta, Z) in angPos and Cartesian in pos.
// ---------------------------------------------------------------------------
function ring_ParticlePosSelect(r, rng, rmin, rmax, partID) {
  const p = r.p[partID];

  // Equal-area radius sampling
  const rr = f32(Math.sqrt(f32(
    f32(rng.ran2() * f32(f32(rmax * rmax) - f32(rmin * rmin)))
    + f32(rmin * rmin)
  )));

  // Uniform angle
  const theta = f32(f32(rng.ran2()) * f32(2.0 * Pi));

  // sech^2 height via atanh
  const z = f32(f32(Math.atanh(f32(f32(2.0) * f32(rng.ran2()) - f32(1.0))))
    * f32(r.z0));

  // Store cylindrical coordinates
  p.angPos[0] = rr;
  p.angPos[1] = theta;
  p.angPos[2] = z;

  // Store Cartesian coordinates
  p.pos[0] = f32(rr * f32(Math.cos(theta)));
  p.pos[1] = f32(rr * f32(Math.sin(theta)));
  p.pos[2] = z;
}


// ---------------------------------------------------------------------------
// particlePosProject
// Fortran: ParticlePosProject(P, Inclination, PositionAngle)
//
// Projects particle from ring plane to sky plane:
//   1. Incline: XTemp=Pos[0], YTemp=Pos[1]*cos(i) - Pos[2]*sin(i)
//   2. Rotate by PA:
//      ProjectedPos[0] = XTemp*cos(PA) - YTemp*sin(PA)
//      ProjectedPos[1] = XTemp*sin(PA) + YTemp*cos(PA)
// ---------------------------------------------------------------------------
function particlePosProject(p, inclination, positionAngle) {
  const xTemp = f32(p.pos[0]);
  const yTemp = f32(
    f32(f32(p.pos[1]) * f32(Math.cos(inclination)))
    - f32(f32(p.pos[2]) * f32(Math.sin(inclination)))
  );

  const cpa = f32(Math.cos(positionAngle));
  const spa = f32(Math.sin(positionAngle));

  p.projectedPos[0] = f32(f32(xTemp * cpa) - f32(yTemp * spa));
  p.projectedPos[1] = f32(f32(xTemp * spa) + f32(yTemp * cpa));
}


// ---------------------------------------------------------------------------
// particlePos_NewCenter
// Fortran: ParticlePos_NewCenter(P, NewCent)
//
// Shifts projected position by galaxy center offset.
// ---------------------------------------------------------------------------
function particlePos_NewCenter(p, newCent) {
  p.projectedPos[0] = f32(f32(p.projectedPos[0]) + f32(newCent[0]));
  p.projectedPos[1] = f32(f32(p.projectedPos[1]) + f32(newCent[1]));
}


// ---------------------------------------------------------------------------
// ring_CalcParticle_VSys
// Fortran: Ring_CalcParticle_VSys(R, PartID, idum)
//
// Calculates line-of-sight velocity for particle PartID:
//   vLOS = VSys
//        + VRot * cos(theta) * sin(incl)     (rotation)
//        + VRad * sin(theta) * sin(incl)     (radial)
//        + Vvert * cos(incl)                 (vertical bulk motion)
//        + gasdev() * VDisp                  (dispersion scatter)
//
// NOTE: The dvdz vertical gradient branch references `delZ` which is
// uninitialized in the Fortran source (assignment line is commented out).
// In the standard fit dvdz=0, so the branch is never entered. Ported
// with delZ=0 matching typical Fortran local variable initialisation.
// ---------------------------------------------------------------------------
function ring_CalcParticle_VSys(r, partID, rng) {
  const p = r.p[partID];

  let vRotP  = f32(r.vRot);
  const vRadP  = f32(r.vRad);
  const vDispP = f32(r.vDisp);

  // Vertical velocity gradient (dvdz branch)
  // delZ is uninitialized in Fortran — use 0 matching standard fit behaviour
  const delZ = f32(0.0);
  if (f32(Math.abs(p.angPos[2])) > f32(r.zGradiantStart)) {
    vRotP = f32(vRotP - f32(f32(r.dvdz) * delZ));
  }

  const cTheta = f32(Math.cos(p.angPos[1]));
  const sTheta = f32(Math.sin(p.angPos[1]));

  const vFromRotation = f32(f32(vRotP * cTheta) * f32(Math.sin(r.inclination)));
  const vFromRadial   = f32(f32(vRadP * sTheta) * f32(Math.sin(r.inclination)));
  const vFromVertical = f32(f32(r.vvert) * f32(Math.cos(r.inclination)));

  p.projectedVel[2] = f32(
    f32(f32(r.vSys) + vFromRotation)
    + f32(vFromRadial + vFromVertical)
  );

  // Gaussian dispersion scatter
  p.projectedVel[2] = f32(
    f32(p.projectedVel[2]) + f32(f32(rng.gasdev()) * vDispP)
  );
}


// ---------------------------------------------------------------------------
// ring_CalcParticleFlux_Basic
// Fortran: Ring_CalcParticleFlux_Basic(R, PartID, Area)
//
// Assigns flux to particle — equal share of total ring flux.
// Flux = Sigma * Area / nParticles
// ---------------------------------------------------------------------------
function ring_CalcParticleFlux_Basic(r, partID, area) {
  const sigmaP = f32(r.sigma);
  r.p[partID].flux = f32(f32(sigmaP * area) / f32(r.nParticles));
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  ring_CalcNumParticles,
  ring_ParticleGeneration,
  ring_ParticlePosSelect,
  particlePosProject,
  particlePos_NewCenter,
  ring_CalcParticle_VSys,
  ring_CalcParticleFlux_Basic
};


// ---------------------------------------------------------------------------
// Self-test (node SingleRingGeneration.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { makeRng }             = require('../StandardMath/random.js');
  const { Ring }                = require('../ObjectDefinitions/TiltedRing.js');
  const { ring_ParticleAllocation } = require('../ObjectDefinitions/TiltedRing.js');

  // Build a simple test ring
  const r           = new Ring();
  r.rmid            = f32(10.0);
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

  console.log('=== ring_CalcNumParticles ===');
  ring_CalcNumParticles(r, 1, f32(1.0));
  console.log('nParticles:', r.nParticles);

  console.log('\n=== ring_ParticleGeneration (seed -1) ===');
  const rng = makeRng(-1);
  ring_ParticleAllocation(r);
  ring_ParticleGeneration(r, rng);

  console.log('particle[0]:');
  console.log('  angPos:       ', Array.from(r.p[0].angPos).map(v => v.toFixed(6)));
  console.log('  projectedPos: ', Array.from(r.p[0].projectedPos).map(v => v.toFixed(6)));
  console.log('  projectedVel[2]:', r.p[0].projectedVel[2].toFixed(4));
  console.log('  flux:         ', r.p[0].flux.toExponential(6));

  // Flux conservation check
  let totalFlux = f32(0.0);
  for (let i = 0; i < r.nParticles; i++) totalFlux = f32(totalFlux + r.p[i].flux);
  const area = f32(Pi * f32(f32(f32(r.rmid + r.rwidth/2)**2) - f32(f32(r.rmid - r.rwidth/2)**2)));
  const expectedFlux = f32(r.sigma * area);
  console.log('\n=== flux conservation ===');
  console.log('total flux:   ', totalFlux.toExponential(6));
  console.log('sigma*area:   ', expectedFlux.toExponential(6));
  console.log('match:', Math.abs(totalFlux - expectedFlux) < 1e-5 ? 'OK' : 'FAIL');
}
