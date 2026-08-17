module.declare(["./DataCube.js","./CubeDifference.js","./PhysCoordTransform.js","./GenerateBootstrap.js"], function (require, exports, module) {
'use strict';

// =============================================================================
// FlipBootstrap.js
// Port of src/BootstrapSampler/FlippingBootstrap.f (FlippingBootstrapMod)
//
// This is the resampling method the native `BootStrapSampler` binary
// ACTUALLY runs (src/ProgramMains/BootStrapGenerator.f calls
// GenFlipBootstrapSample(), not GenBootstrapSample() -- that call is
// commented out there). GenerateBootstrap.js in this same directory ports
// GenerateBootstrap.f's block-resample method instead, which is real,
// tested code, but not what native BootStrapSampler currently executes --
// left in place, unused by the launcher, rather than deleted.
//
// ALGORITHM
// ---------
// 1. residual = observed - model (CubeDifference.js, unchanged)
// 2. physical coord array built on the OBSERVED cube's header
//    (PhysCoordTransform.js, unchanged)
// 3. bootstrapCube = residual (full copy)
// 4. axisFlip(bootstrapCube, FlipType=1): per channel-block, flip a coin;
//    if heads, mirror-reflect that block's spatial map across the major
//    axis (theta -> 2*pi - theta, R and dV unchanged)
// 5. axisFlip(bootstrapCube, FlipType=2): same coin-flip loop, mirroring
//    through the disk centre instead (theta -> theta + pi, dV -> -dV) --
//    i.e. swapping which side is the receding/approaching side
// 6. bootstrapCube += model
//
// RNG
// ---
// Fortran's AxisFlip calls plain RANDOM_NUMBER, one draw per channel-or-
// block, in a loop that self-terminates once `count` (the channel pointer,
// which advances faster than the loop counter whenever a flip is taken)
// reaches nChannels -- see the exact count/goto replication in axisFlip()
// below. The TOTAL number of draws per call is therefore RNG-dependent
// (more flips taken -> fewer outer-loop iterations needed), but perfectly
// deterministic given a matched seed: this is not a rejection/retry loop
// (no draw is ever discarded), so bit-exact draw-for-draw matching against
// Fortran is achievable once both sides use the same ran2 (see
// src/StandardMath/random.js/.f), unlike GenerateBootstrap.f's rejection-
// sampling block resampler.
//
// `rng` here is the object returned by src/StandardMath/random.js's
// makeRng(idum) -- same convention bootstrap-fit-launcher.js already uses
// for round 2 (`rng.ran2()`). Passing plain Math.random() in its place
// (an object with a `ran2` method) works too for callers that don't care
// about seed matching.
// =============================================================================

const f32 = Math.fround;
const Pi  = f32(Math.PI);

const { DataCube, allocateDataCube, flatIndxCalc } =
  require('./DataCube.js');
const { constructDiffCube } = require('./CubeDifference.js');
const { buildPhysCoordsArray, getCubeCoords, coordGet } =
  require('./PhysCoordTransform.js');
const { getFluxAtPoint, simpleBoundCheckReal } = require('./GenerateBootstrap.js');


// ---------------------------------------------------------------------------
// copyDataCube
// Fortran derived-type assignment (`TempCube=Cube`, `BootstrapCube=DifferenceCube`)
// deep-copies header + flux + valid indices. No JS equivalent of that
// assignment operator, so this helper does it explicitly -- used at both
// struct-copy points in this file (mirrors CubeDifference.js's own header-copy
// block, which computes a subtraction instead of a plain copy).
// ---------------------------------------------------------------------------
function copyDataCube(src) {
  const dst = new DataCube();
  const dh  = src.dh;
  dst.dh.nPixels[0]     = dh.nPixels[0];
  dst.dh.nPixels[1]     = dh.nPixels[1];
  dst.dh.nChannels      = dh.nChannels;
  dst.dh.pixelSize[0]   = dh.pixelSize[0];
  dst.dh.pixelSize[1]   = dh.pixelSize[1];
  dst.dh.channelSize    = dh.channelSize;
  dst.dh.start[0]       = dh.start[0];
  dst.dh.start[1]       = dh.start[1];
  dst.dh.start[2]       = dh.start[2];
  dst.dh.refLocation[0] = dh.refLocation[0];
  dst.dh.refLocation[1] = dh.refLocation[1];
  dst.dh.refLocation[2] = dh.refLocation[2];
  dst.dh.refVal[0]      = dh.refVal[0];
  dst.dh.refVal[1]      = dh.refVal[1];
  dst.dh.refVal[2]      = dh.refVal[2];
  dst.dh.uncertainty    = dh.uncertainty;
  dst.dh.nValid         = dh.nValid;

  allocateDataCube(dst);
  dst.dh.start[0] = dh.start[0];
  dst.dh.start[1] = dh.start[1];
  dst.dh.start[2] = dh.start[2];

  dst.flux.set(src.flux);
  dst.flattendValidIndices.set(src.flattendValidIndices);
  return dst;
}


// ---------------------------------------------------------------------------
// flipChannelSpatial
// Fortran: FlipChannelSpatial_Spectral_Dir(Cube,TempCube,CoordArr,ChanID,
//            XC,YC,VSys,PA,Inc,FlipType)
//
// FlipType 1: mirror across the major axis -- theta -> 2*pi - theta,
//   R and dV unchanged.
// FlipType 2: mirror through the disk centre -- theta -> theta + pi
//   (wrapped into [0, 2*pi)), dV -> -dV (swaps receding/approaching side).
//
// Reads from `cube` (unmodified source), writes the flipped channel into
// `tempCube` (accumulator across both AxisFlip passes).
// ---------------------------------------------------------------------------
function flipChannelSpatial(cube, tempCube, coordArr, chanID, xc, yc, vSys, pa, inc, flipType) {
  const dh  = cube.dh;
  const nx  = dh.nPixels[0];
  const ny  = dh.nPixels[1];
  const nch = dh.nChannels;

  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      const r     = coordGet(coordArr, 0, i, j, chanID, nx, ny, nch);
      const theta = coordGet(coordArr, 1, i, j, chanID, nx, ny, nch);
      const dv    = coordGet(coordArr, 2, i, j, chanID, nx, ny, nch);

      let physCoordFlip;
      if (flipType === 1) {
        const thetaNew = f32(f32(2.0) * Pi - theta);
        physCoordFlip = [r, thetaNew, dv];
      } else {
        let thetaNew = f32(theta + Pi);
        if (thetaNew > f32(2.0) * Pi) thetaNew = f32(thetaNew - f32(2.0) * Pi);
        physCoordFlip = [r, thetaNew, f32(-dv)];
      }

      const cubeCoordFlip = getCubeCoords(xc, yc, vSys, pa, inc, physCoordFlip);

      let flux;
      if (simpleBoundCheckReal(cube, cubeCoordFlip[0], cubeCoordFlip[1], cubeCoordFlip[2])) {
        flux = getFluxAtPoint(cube, cubeCoordFlip);
      } else {
        flux = cube.flux[flatIndxCalc(i, j, chanID, dh)];
      }
      tempCube.flux[flatIndxCalc(i, j, chanID, dh)] = flux;
    }
  }
}


// ---------------------------------------------------------------------------
// axisFlip
// Fortran: AxisFlip(Cube,CoordArr,XC,YC,VSys,PA,Inc,FlipType,VBlockSize)
//
// Per-channel-block coin flip: for each block of vBlockInt channels, draw
// one RandVal; if >= 0.5, flip every channel in that block via
// flipChannelSpatial. `count` (the channel pointer) and `i` (the outer loop
// counter) advance at different rates -- replicated exactly, including the
// unconditional `count += 1` after the if-block (present in Fortran whether
// or not the flip was taken) and the early exit once count reaches
// nChannels (Fortran's `goto 100`).
// ---------------------------------------------------------------------------
function axisFlip(cube, coordArr, xc, yc, vSys, pa, inc, flipType, vBlockSize, rng) {
  const dh        = cube.dh;
  const nChannels = dh.nChannels;

  // Fortran: TempCube=Cube (deep copy, including flux -- unflipped channels
  // must retain Cube's original values in the final result).
  const tempCube = copyDataCube(cube);

  let vBlockInt = Math.trunc(vBlockSize);
  if (vBlockInt < 1) vBlockInt = 1;

  let count = 0;
  outer:
  for (let i = 0; i < nChannels; i++) {
    const randVal = rng.ran2();
    if (randVal >= 0.5) {
      for (let j = 1; j <= vBlockInt; j++) {
        flipChannelSpatial(cube, tempCube, coordArr, count, xc, yc, vSys, pa, inc, flipType);
        count += 1;
        if (count >= nChannels) break outer;
      }
    }
    count += 1;
    if (count >= nChannels) break outer;
  }

  // Fortran: Cube%Flux=TempCube%Flux
  cube.flux.set(tempCube.flux);
}


// ---------------------------------------------------------------------------
// genFlipBootstrapSample
// Fortran: GenFlipBootstrapSample() using BootstrapGlobals
//
// Arguments:
//   observedDC   — DataCube (observed)
//   modelDC      — DataCube (best-fit model, post-convolution)
//   bsCent       — { centX, centY, centV, pa, inc }
//   velBlockSize — real (channels per block)
//   rng          — { ran2: () => number }, e.g. makeRng(idum) from
//                  src/StandardMath/random.js
//
// Returns: DataCube (bootstrap realization)
// ---------------------------------------------------------------------------
function genFlipBootstrapSample(observedDC, modelDC, bsCent, velBlockSize, rng) {
  // Step 1: residual cube
  const diffCube = new DataCube();
  constructDiffCube(observedDC, modelDC, diffCube);

  // Step 2: physical coordinate array, built on the observed cube's header
  const coordArr = buildPhysCoordsArray(
    bsCent.centX, bsCent.centY, bsCent.centV,
    bsCent.pa, bsCent.inc,
    observedDC.dh
  );

  // Step 3: BootstrapCube = DifferenceCube (full copy)
  const bootstrapCube = copyDataCube(diffCube);

  // Step 4/5: flip across the major axis, then through the disk centre
  axisFlip(bootstrapCube, coordArr, bsCent.centX, bsCent.centY, bsCent.centV,
    bsCent.pa, bsCent.inc, 1, velBlockSize, rng);
  axisFlip(bootstrapCube, coordArr, bsCent.centX, bsCent.centY, bsCent.centV,
    bsCent.pa, bsCent.inc, 2, velBlockSize, rng);

  // Step 6: add model back
  const nCells = bootstrapCube.flux.length;
  for (let l = 0; l < nCells; l++) {
    bootstrapCube.flux[l] = f32(bootstrapCube.flux[l] + modelDC.flux[l]);
  }

  return bootstrapCube;
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  genFlipBootstrapSample,
  axisFlip,
  flipChannelSpatial,
  copyDataCube,
};

});
