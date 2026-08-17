module.declare(["./TiltedRing.js","./SingleRingGeneration.js"], function (require, exports, module) {
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
// dc/beam: the observed DataCube/Beam2D (js/src/ObjectDefinitions/{DataCube,
// Beam}.js) — feed calcAvgChanPerPix below, matching Fortran's BuildTiltedRingModel
// taking DC/BUse as of upstream commit 76ade48.
// =============================================================================

const f32 = Math.fround;

const { ring_ParticleAllocation } = require('./TiltedRing.js');
const {
  ring_CalcNumParticles,
  ring_ParticleGeneration,
  warmUpWasmTrig
} = require('./SingleRingGeneration.js');


// ---------------------------------------------------------------------------
// calcAvgChanPerPix
// Fortran: CalcAvgChanPerPix(ringIndx, TR, DC, BUse, AvgChanPerPix)
//
// How many spectral channels a ring's own velocity spread (rotation +
// dispersion) covers, relative to how many beam-widths wide the ring is --
// rings with a steep velocity gradient or wide dispersion smear their flux
// across more channels per spatial pixel, so need proportionally more
// particles to sample adequately. Floored at dDisp+1 so a ring with ~zero
// rotation still gets at least dispersion-width coverage.
// ---------------------------------------------------------------------------
function calcAvgChanPerPix(ring, dc, beam) {
  const chanSize = f32(Math.abs(dc.dh.channelSize));
  const dv        = f32(f32(2.0 * f32(ring.vRot)) / chanSize);
  const dDisp      = f32(f32(f32(2.0 * f32(Math.sqrt(2.0))) * f32(ring.vDisp)) / chanSize);
  const dr         = f32(2.0 * f32(ring.rmid));

  let avgChanPerPix = f32(f32(f32(beam.beamMajorAxis) / dr) * f32(dv + dDisp));
  if (avgChanPerPix < f32(dDisp + 1.0)) {
    avgChanPerPix = f32(dDisp + 1.0);
  }
  return avgChanPerPix;
}


// ---------------------------------------------------------------------------
// buildTiltedRingModel
// Fortran: BuildTiltedRingModel(TR, idum, Noise, DC, BUse)
//
// For each ring in TR:
//   1. Calculate the ring's average channels-per-pixel (calcAvgChanPerPix)
//   2. Calculate number of particles
//   3. Allocate particle array
//   4. Generate particle positions, velocities, and fluxes
//
// TR must already be allocated (tiltRing_Allocate called) and all ring
// parameters set before calling this.
//
// rng: makeRng() object from random.js
// noise: DC%DH%Uncertainty*abs(DC%DH%ChannelSize), computed by the caller
//   (matches FullModelComparison.f / FitOutput.f's NoiseSpec/SpecNoise).
// dc, beam: the observed DataCube/Beam2D, fed to calcAvgChanPerPix.
// ---------------------------------------------------------------------------
function buildTiltedRingModel(tr, rng, noise, dc, beam) {
  for (let i = 0; i < tr.nRings; i++) {
    const avgChanPerPix = calcAvgChanPerPix(tr.r[i], dc, beam);
    ring_CalcNumParticles(tr.r[i], tr.cmode, tr.cloudBaseSurfDens, noise, avgChanPerPix);
    ring_ParticleAllocation(tr.r[i]);
    ring_ParticleGeneration(tr.r[i], rng);
  }
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = { buildTiltedRingModel, calcAvgChanPerPix, warmUpWasmTrig };


// ---------------------------------------------------------------------------
// Self-test (node TiltedRingModelGeneration.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { makeRng }          = require('./random.js');
  const { TiltedRingModel, tiltRing_Allocate, tiltRing_DeAllocate }
                             = require('./TiltedRing.js');
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
  const fakeDc = { dh: { channelSize: f32(4.0) } };
  const fakeBeam = { beamMajorAxis: f32(3.0) };
  const noise = f32(0.01);
  buildTiltedRingModel(tr, rng, noise, fakeDc, fakeBeam);

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

});
