'use strict';

// =============================================================================
// FillDataCubeByTiltedRing.js
// High-fidelity port of src/TiltedRingToDataCube/FillDataCubeByTiltedRing.f
// (FillDataCubeWithTiltedRingMod)
//
// PORTING NOTES
// -------------
// Fortran `real` → Math.fround(). Integer truncation → Math.trunc().
//
// COORDINATE MAPPING — two different conventions per axis:
//   Spatial (x, y): projectedPos is already in pixel coordinates.
//     CellIndex = int(projectedPos + 0.5)  ← nearest-pixel rounding
//   Velocity:       projectedVel[2] is in km/s, must convert to channel.
//     CellIndex = int((vel - start) / channelSize + 0.5)
//
// Flux array: DC%Flux(i,j,k) stored as flat Float32Array via flatIndxCalc.
// DC%Flux=0. at start of FillDataCubeWithTiltedRing — resets every call.
// =============================================================================

const f32 = Math.fround;
const { flatIndxCalc } = require('../ObjectDefinitions/DataCube.js');


// ---------------------------------------------------------------------------
// findParticleCellLocation
// Fortran: FindParticleCellLocation(P, DC, CellIndex)
//
// Returns [i, j, k] cell indices for a particle.
// Spatial: nearest pixel via int(pos + 0.5)
// Velocity: int((vel - start) / channelSize + 0.5)
// ---------------------------------------------------------------------------
function findParticleCellLocation(p, dc) {
  const dh = dc.dh;
  const i  = Math.trunc(f32(p.projectedPos[0]) + f32(0.5));
  const j  = Math.trunc(f32(p.projectedPos[1]) + f32(0.5));
  const k  = Math.trunc(
    f32(f32(f32(p.projectedVel[2]) - f32(dh.start[2]))
    / f32(dh.channelSize)) + f32(0.5)
  );
  return [i, j, k];
}


// ---------------------------------------------------------------------------
// checkIfInCube
// Fortran: CheckIfInCube(CellIndex, DC, InBounds)
//
// Returns true if [i, j, k] is within the cube bounds.
// ---------------------------------------------------------------------------
function checkIfInCube(cellIndex, dc) {
  const dh = dc.dh;
  if (cellIndex[0] < 0 || cellIndex[0] >= dh.nPixels[0]) return false;
  if (cellIndex[1] < 0 || cellIndex[1] >= dh.nPixels[1]) return false;
  if (cellIndex[2] < 0 || cellIndex[2] >= dh.nChannels)  return false;
  return true;
}


// ---------------------------------------------------------------------------
// fillDataCubeWithTiltedRing
// Fortran: FillDataCubeWithTiltedRing(DC, TR)
//
// Resets DC flux to 0, then splats every particle from every ring
// into the appropriate voxel.
// ---------------------------------------------------------------------------
function fillDataCubeWithTiltedRing(dc, tr) {
  // Reset cube to zero — Fortran: DC%Flux=0.
  dc.flux.fill(0);

  for (let i = 0; i < tr.nRings; i++) {
    const ring = tr.r[i];
    for (let j = 0; j < ring.nParticles; j++) {
      const cell = findParticleCellLocation(ring.p[j], dc);
      if (checkIfInCube(cell, dc)) {
        const idx      = flatIndxCalc(cell[0], cell[1], cell[2], dc.dh);
        dc.flux[idx]   = f32(dc.flux[idx] + f32(ring.p[j].flux));
      }
    }
  }
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  findParticleCellLocation,
  checkIfInCube,
  fillDataCubeWithTiltedRing
};


// ---------------------------------------------------------------------------
// Self-test (node FillDataCubeByTiltedRing.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { makeRng }          = require('../StandardMath/random.js');
  const { DataCube, DataCubeHeader, allocateDataCube, flatIndxCalc }
                             = require('../ObjectDefinitions/DataCube.js');
  const { TiltedRingModel, tiltRing_Allocate, tiltRing_DeAllocate }
                             = require('../ObjectDefinitions/TiltedRing.js');
  const { buildTiltedRingModel }
                             = require('../TiltedModelGeneration/TiltedRingModelGeneration.js');
  const Pi = f32(Math.PI);

  // Build a small cube — 100x100 spatial, 200 channels
  const dc = new DataCube();
  const dh = dc.dh;
  dh.nPixels[0]     = 100;
  dh.nPixels[1]     = 100;
  dh.nChannels      = 200;
  dh.pixelSize[0]   = f32(-1.0);
  dh.pixelSize[1]   = f32(1.0);
  dh.channelSize    = f32(10.0);
  dh.refLocation[0] = f32(50.0);
  dh.refLocation[1] = f32(50.0);
  dh.refLocation[2] = f32(100.0);
  dh.refVal[0]      = f32(0.0);
  dh.refVal[1]      = f32(0.0);
  dh.refVal[2]      = f32(1000.0);
  allocateDataCube(dc);

  // Build a 3-ring model centred on the cube
  const tr             = new TiltedRingModel();
  tr.nRings            = 3;
  tr.cmode             = 1;
  tr.cloudBaseSurfDens = f32(1.0);
  tiltRing_Allocate(tr);

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
  fillDataCubeWithTiltedRing(dc, tr);

  // Sum total flux in cube — should match sum of particle fluxes
  let cubeFlux = f32(0.0);
  for (let v of dc.flux) cubeFlux = f32(cubeFlux + f32(v));

  let particleFlux = f32(0.0);
  for (let i = 0; i < tr.nRings; i++)
    for (let j = 0; j < tr.r[i].nParticles; j++)
      particleFlux = f32(particleFlux + f32(tr.r[i].p[j].flux));

  console.log('=== fillDataCubeWithTiltedRing ===');
  console.log('cube total flux:     ', cubeFlux.toExponential(6));
  console.log('particle total flux: ', particleFlux.toExponential(6));
  console.log('match:', Math.abs(cubeFlux - particleFlux) < 1e-4 ? 'OK' : 'FAIL');
  console.log('non-zero voxels:',
    Array.from(dc.flux).filter(v => v > 0).length);

  // Check findParticleCellLocation on a known particle
  const p0   = tr.r[0].p[0];
  const cell = findParticleCellLocation(p0, dc);
  console.log('\n=== findParticleCellLocation (ring 0, particle 0) ===');
  console.log('projectedPos:', Array.from(p0.projectedPos).map(v => v.toFixed(4)));
  console.log('projectedVel[2]:', p0.projectedVel[2].toFixed(4));
  console.log('cellIndex:', cell);
  console.log('inBounds:', checkIfInCube(cell, dc));

  tiltRing_DeAllocate(tr);
}
