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
// Faithful fdlibm port of sin/cos/atanh -- Fortran now also calls fdlibm
// (instead of the system libm) for these, closing a confirmed nonzero
// per-call divergence. See ../StandardMath/fdlibm.js.
//
// Default to the JS port (correct everywhere, no build step). Optionally
// upgraded to the real compiled-wasm fdlibm via warmUpWasmTrig() below --
// verified bit-exact against this same JS port (0/200000 mismatches,
// incl. the specific FMA-fusion regression input; see that wasm build's
// own header for why a naive wasm build of fdlibm's kernel files does
// NOT get this for free) and measured ~1.2x faster per real-pipeline
// objective-function evaluation (real SoFiA-derived 5-ring galaxy fit,
// not a synthetic benchmark). Bindings are `let`, not `const`,
// specifically so warmUpWasmTrig() can swap them in place -- every call
// site below reads the same three names, so no call site changes when
// the backend is upgraded.
const jsFdlibm = require('../StandardMath/fdlibm.js');
let fdSin = jsFdlibm.fdSin;
let fdCos = jsFdlibm.fdCos;
let fdAtanh = jsFdlibm.fdAtanh;

let wasmTrigModulePromise = null;

/**
 * Warms up the real compiled-wasm fdlibm build (../StandardMath/
 * fdlibm-wasm.js, bundled in-tree here -- not yet published as a DCP
 * package, see that file's own header) and, once ready, swaps
 * fdSin/fdCos/fdAtanh (module-level, used by every call site below) to
 * its wasm-backed functions. Same warmUp-once-then-sync pattern as
 * ../ConvolveCube/CubeKernelConvolution.js's FFTW3 integration -- call
 * this once (e.g. alongside `await fftwWarmUpPromise` in the launcher)
 * before running any fit. Safe to never call: the JS fallback above is
 * already correct, just slower.
 */
async function warmUpWasmTrig() {
  if (!wasmTrigModulePromise) {
    wasmTrigModulePromise = (async () => {
      const fdlibmWasm = require('../StandardMath/fdlibm-wasm.js');
      await fdlibmWasm.warmUp();
      fdSin = fdlibmWasm.fdSinWasm;
      fdCos = fdlibmWasm.fdCosWasm;
      fdAtanh = fdlibmWasm.fdAtanhWasm;
    })();
  }
  await wasmTrigModulePromise;
}

// `process` doesn't exist in the DCP worker sandbox (browser-like Worker
// environment, not Node.js) -- must guard, not just check process.env
// directly, or this throws ReferenceError the moment either debug branch
// below is reached.
const TRACE_DEBUG = typeof process !== 'undefined' && process.env && process.env.TRACE_DEBUG === '1';


// ---------------------------------------------------------------------------
// ring_CalcNumParticles
// Fortran: Ring_CalcNumParticles(R, cmode, CloudSurfDens, Noise, AvgChannelsPerPix)
//
// Calculates the number of particles needed to represent the ring,
// proportional to ring area and surface density (now relative to the
// noise level), and scaled by AvgChannelsPerPix so rings whose velocity
// spread covers many spectral channels get proportionally more particles
// to sample them adequately (see calcAvgChanPerPix in
// TiltedRingModelGeneration.js). Matches upstream commit 76ade48
// ("Changed how the number of particles in each ring is calculated").
//
// nParticles = int(CloudSurfDens * (Sigma/Noise)^cmode * Pi*(Rh^2-Rl^2) * AvgChannelsPerPix) + 1
// ---------------------------------------------------------------------------
const _f32buf = new Float32Array(1);
const _u32buf = new Uint32Array(_f32buf.buffer);
function hexF32(x) {
  _f32buf[0] = x;
  return _u32buf[0].toString(16).toUpperCase().padStart(8, '0');
}

// Masks off the low 12 bits of x's float32 mantissa (keeping the top 11 of
// 23 bits -- ~0.05% relative precision). Matches Fortran's
// RoundForParticleStability (SingleRingGeneration.f) bit-for-bit -- see its
// header comment for why this exists: R%Sigma (and everything derived from
// it) can differ from this port's sigma by a handful of float32 ULPs, an
// already-known unavoidable cross-platform difference that's normally
// harmless everywhere else in this pipeline, but ring_CalcNumParticles's
// truncation to an integer particle count turns that smooth tiny
// difference into a whole-particle difference whenever the true value
// lands near an integer boundary -- which then permanently desyncs
// ran2()/gasdev()'s idum for the rest of the fit. ~900x safety margin over
// the observed ~5e-7 relative noise.
function roundForParticleStability(x) {
  _f32buf[0] = x;
  _u32buf[0] = _u32buf[0] & 0xFFFFF000;
  return _f32buf[0];
}

function ring_CalcNumParticles(r, cmode, cloudSurfDens, noise, avgChannelsPerPix) {
  const rl          = f32(f32(r.rmid) - f32(f32(r.rwidth) / f32(2.0)));
  const rh          = f32(f32(r.rmid) + f32(f32(r.rwidth) / f32(2.0)));
  const pixelRing   = f32(Pi * f32(f32(rh * rh) - f32(rl * rl)));
  // DELIBERATELY on the pre-upstream-commit-76ade48 formula for now (not
  // using noise/avgChannelsPerPix below) -- Nathan's own attached example
  // config used cdens=10, but his email text says the actual default is
  // cdens=100 (10 looks like a typo in the example, not the intended
  // default); revisit once that's confirmed and re-verify Fortran/JS
  // agreement before switching this back on. noise/avgChannelsPerPix stay
  // threaded through this function's signature (unused here) so
  // re-enabling is a one-line change -- see calcAvgChanPerPix in
  // TiltedRingModelGeneration.js, already wired and confirmed
  // bit-matching against Fortran.
  const densMulti   = f32(f32(cloudSurfDens) * f32(f32(r.sigma) ** cmode));
  r.nParticles      = Math.trunc(roundForParticleStability(f32(densMulti * pixelRing))) + 1;
  if (TRACE_DEBUG) {
    console.error('NPTRACE', r.rmid, r.sigma, densMulti, pixelRing, avgChannelsPerPix, r.nParticles);
    console.error('NPHEX', hexF32(r.sigma), hexF32(noise), hexF32(densMulti), hexF32(pixelRing), hexF32(avgChannelsPerPix));
  }
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

  // One-off diagnostic ("bisection paradox" investigation, Dan 2026-08-18):
  // see SingleRingGeneration.f's matching RINGGEOMTRACE comment for why.
  if (typeof process !== 'undefined' && process.env && process.env.TRACE_OVERRIDE_IDUM) {
    console.error('RINGGEOMTRACE', r.rmid.toExponential(17), r.rwidth.toExponential(17),
      r.sigma.toExponential(17), r.nParticles, area.toExponential(17));
  }

  // inclination/positionAngle are the same for every particle in this ring
  // -- fdCos/fdSin are pure, deterministic functions of these, so computing
  // them once here and passing the results into particlePosProject/
  // ring_CalcParticle_VSys (rather than each recomputing them per particle,
  // and ring_CalcParticle_VSys separately recomputing sin/cos(inclination)
  // that particlePosProject already computed) produces bit-identical
  // results -- same inputs, same pure function, just not redundantly
  // called up to ~4x per particle. Real measured cost on rings with tens
  // of thousands of particles, not a micro-optimization: see the perf
  // investigation that found this.
  const cosIncl = f32(fdCos(r.inclination));
  const sinIncl = f32(fdSin(r.inclination));
  const cosPA   = f32(fdCos(r.positionAngle));
  const sinPA   = f32(fdSin(r.positionAngle));

  for (let i = 0; i < r.nParticles; i++) {
    ring_ParticlePosSelect(r, rng, rmin, rmax, i);
    particlePosProject(r.p[i], cosIncl, sinIncl, cosPA, sinPA);
    particlePos_NewCenter(r.p[i], r.centPos);
    ring_CalcParticle_VSys(r, i, rng, sinIncl, cosIncl);
    ring_CalcParticleFlux_Basic(r, i, area);
    if (TRACE_DEBUG && i < 5) {
      const p = r.p[i];
      console.error('PARTTRACE', r.rmid, i, rng.state.ran2State.idum,
        p.pos[0], p.pos[1], p.pos[2], p.projectedPos[0], p.projectedPos[1], p.projectedVel[2]);
    }
    // One-off diagnostic (Fortran-vs-JS gasdev-desync bisection, Dan
    // 2026-08-18): checkpoint every 200th particle (plus the very last
    // one) across the FULL ring, not just the first 5 PARTTRACE covers --
    // to localize exactly where a desync (re)starts, if one does, deeper
    // into a ring's particle loop. Gated on TRACE_OVERRIDE_IDUM. Matches
    // Fortran's equivalent in this same function.
    if (typeof process !== 'undefined' && process.env && process.env.TRACE_OVERRIDE_IDUM
        && (i % 200 === 0 || i === r.nParticles - 1)) {
      console.error('BISECTTRACE', r.rmid.toFixed(6), i, rng.state.ran2State.idum,
        r.p[i].projectedVel[2].toFixed(8));
    }
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
  const z = f32(f32(fdAtanh(f32(f32(2.0) * f32(rng.ran2()) - f32(1.0))))
    * f32(r.z0));

  // Store cylindrical coordinates
  p.angPos[0] = rr;
  p.angPos[1] = theta;
  p.angPos[2] = z;

  // Store Cartesian coordinates
  p.pos[0] = f32(rr * f32(fdCos(theta)));
  p.pos[1] = f32(rr * f32(fdSin(theta)));
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
//
// Takes precomputed cos/sin(inclination)/cos/sin(positionAngle) -- see
// ring_ParticleGeneration's call site for why (ring-constant, hoisted
// out of the per-particle loop; bit-identical result either way).
// ---------------------------------------------------------------------------
function particlePosProject(p, cosIncl, sinIncl, cosPA, sinPA) {
  const xTemp = f32(p.pos[0]);
  const yTemp = f32(
    f32(f32(p.pos[1]) * cosIncl)
    - f32(f32(p.pos[2]) * sinIncl)
  );

  p.projectedPos[0] = f32(f32(xTemp * cosPA) - f32(yTemp * sinPA));
  p.projectedPos[1] = f32(f32(xTemp * sinPA) + f32(yTemp * cosPA));
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
//
// Takes precomputed sin/cos(inclination) -- see ring_ParticleGeneration's
// call site: ring-constant, and particlePosProject already needs the same
// two values, so passing them in avoids a second redundant fdSin/fdCos
// pair per particle. Bit-identical result either way (pure function, same
// input).
// ---------------------------------------------------------------------------
function ring_CalcParticle_VSys(r, partID, rng, sinIncl, cosIncl) {
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

  const cTheta = f32(fdCos(p.angPos[1]));
  const sTheta = f32(fdSin(p.angPos[1]));

  const vFromRotation = f32(f32(vRotP * cTheta) * sinIncl);
  const vFromRadial   = f32(f32(vRadP * sTheta) * sinIncl);
  const vFromVertical = f32(f32(r.vvert) * cosIncl);

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
  ring_CalcParticleFlux_Basic,
  warmUpWasmTrig
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
  // noise=1.0, avgChannelsPerPix=1.0: neutral stand-ins so this demo's
  // particle count stays directly comparable to the pre-upstream-sync
  // formula (sigma/noise=sigma, *avgChannelsPerPix is a no-op).
  ring_CalcNumParticles(r, 1, f32(1.0), f32(1.0), f32(1.0));
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
