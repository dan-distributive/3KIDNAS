module.declare(["./DataCube.js"], function (require, exports, module) {
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
const { flatIndxCalc } = require('./DataCube.js');

// ---------------------------------------------------------------------------
// roundForBinStability
// Same technique as SingleRingGeneration.js's roundForParticleStability /
// GenerateBootstrap.js's roundForInterpStability (masks off the low 12 bits
// of x's float32 mantissa, keeping the top 11 of 23 bits -- ~0.05% relative
// precision), applied here to a per-particle nearest-pixel rounding
// decision instead. Matches Fortran's RoundForBinStability
// (FillDataCubeByTiltedRing.f) bit-for-bit.
//
// Root cause: p.projectedPos/p.projectedVel (the product of a chain of
// float32 trig/rotation arithmetic during particle generation) can differ
// from Fortran by a handful of float32 ULPs -- the same already-known,
// unavoidable cross-platform difference documented at
// roundForParticleStability. Normally harmless, but
// findParticleCellLocation's int(pos+0.5) is a hard round-to-nearest-pixel
// decision, applied independently to EVERY particle (tens of thousands per
// ring): whenever one particle's position sits within that noise-distance
// of a half-integer boundary, Fortran and this port round it into
// different, adjacent cells. Since each particle's full flux is added (not
// interpolated) to whichever single cell it lands in, one straddling
// particle yanks a chunk of flux out of one cell and into its neighbor on
// only one platform -- confirmed directly (2026-08-17): with particle
// COUNTS already verified identical per ring between two independent
// Fortran/JS initial fits, and the model cube isolated as the only
// diverging input (observed cube confirmed bit-identical) to a downstream
// bootstrap-resampling pixel diff, this per-particle binning round is the
// remaining discrete decision capable of producing that divergence at this
// scale.
// ---------------------------------------------------------------------------
const _binF32buf = new Float32Array(1);
const _binU32buf = new Uint32Array(_binF32buf.buffer);
function roundForBinStability(x) {
  _binF32buf[0] = x;
  _binU32buf[0] = _binU32buf[0] & 0xFFFFF000;
  return _binF32buf[0];
}


// ---------------------------------------------------------------------------
// findParticleCellLocation
// Fortran: FindParticleCellLocation(P, DC, CellIndex)
//
// Writes [i, j, k] cell indices for a particle into `out` (also returned,
// for chaining) and returns it. Fortran's CellIndex is an INTENT(OUT) array
// the caller owns and reuses across every particle -- no per-call
// allocation. This used to `return [i, j, k]`, allocating a fresh array on
// every call; with hundreds of thousands of particles per ring across
// hundreds of evaluations per fit, that was a real, measured hot spot (see
// git history / session notes: this function's containing call,
// fillDataCubeWithTiltedRing, ran ~19.6x slower than Fortran's equivalent,
// vs. a ~2.8-3.7x baseline everywhere else that doesn't allocate per
// particle). `out` defaults to a fresh plain array (not a typed array --
// keeps identical semantics to the old `return [i, j, k]` for any
// out-of-range/NaN edge case Math.trunc could produce, rather than an
// Int32Array silently wrapping/zeroing those) so this still works as a
// plain two-argument call (e.g. the self-test below) when the caller
// doesn't have a scratch buffer to hand in.
//
// Spatial: nearest pixel via int(pos + 0.5)
// Velocity: int((vel - start) / channelSize + 0.5)
// ---------------------------------------------------------------------------
function findParticleCellLocation(p, dc, out = [0, 0, 0]) {
  const dh = dc.dh;
  // Fortran rounds the `+0.5` itself to real4 before int() truncates it
  // (int((...)+0.5)) -- every operation on a real4 rounds separately. Leaving
  // the addition in JS's native double precision before Math.trunc() skips
  // that rounding step, which can flip the truncated result for values
  // landing within ~1 float32 ULP of an integer/half-integer boundary.
  out[0] = Math.trunc(roundForBinStability(f32(f32(p.projectedPos[0]) + f32(0.5))));
  out[1] = Math.trunc(roundForBinStability(f32(f32(p.projectedPos[1]) + f32(0.5))));
  out[2] = Math.trunc(roundForBinStability(
    f32(f32(f32(f32(p.projectedVel[2]) - f32(dh.start[2]))
    / f32(dh.channelSize)) + f32(0.5))
  ));
  return out;
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

  // One scratch buffer reused for every particle -- see findParticleCellLocation's
  // header for why this matters (was a real, measured ~19.6x-vs-Fortran hot spot).
  const cell = [0, 0, 0];

  for (let i = 0; i < tr.nRings; i++) {
    const ring = tr.r[i];
    for (let j = 0; j < ring.nParticles; j++) {
      findParticleCellLocation(ring.p[j], dc, cell);
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
  const { makeRng }          = require('./random.js');
  const { DataCube, DataCubeHeader, allocateDataCube, flatIndxCalc }
                             = require('./DataCube.js');
  const { TiltedRingModel, tiltRing_Allocate, tiltRing_DeAllocate }
                             = require('./TiltedRing.js');
  const { buildTiltedRingModel }
                             = require('./TiltedRingModelGeneration.js');
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
  const fakeBeam = { beamMajorAxis: f32(3.0) };
  const noise = f32(0.01);
  buildTiltedRingModel(tr, rng, noise, dc, fakeBeam);
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

});
