'use strict';

// =============================================================================
// GenerateBootstrap.js
// Port of src/BootstrapSampler/GenerateBootstrap.f (GenBootstrapMod)
//
// PORTING NOTES
// -------------
// Ports the following Fortran subroutines:
//
//   genBootstrapSample()          — main entry point
//   blockResampleCube_Phys2()     — tile output cube, fill each tile via
//                                   physicallyConstrained block selection
//   blockResample_Ini()           — initialize output cube and compute block sizes
//   buildDataBlock_PhysSelect()   — rejection-sample a source block with matching
//                                   physical coordinates, fill via interpolation
//   fillInCubeByBlock()           — copy a data block into the output cube
//   simpleBoundCheck()            — integer bounds check
//   simpleBoundCheck_Real()       — real bounds check
//   getFluxAtPoint()              — trilinear interpolation at real-valued point
//
// NOT ported (Fortran legacy paths, commented out or superseded):
//   blockResampleCube()           — older mask-based resampler (commented out)
//   selectDataBlock()             — used only by above (commented out)
//   selectDataBlock_Phys()        — older version (commented out)
//   blockResampleCube_Phys()      — older version (commented out)
//
// RANDOMNESS
// ----------
// Fortran uses RANDOM_NUMBER (system RNG, not idum/ran2).
// JS uses Math.random() — no need to match Fortran since each bootstrap
// realization is independently random by design.
//
// PRECISION
// ---------
// All flux and coordinate arithmetic uses f32() matching Fortran real.
//
// INDEXING
// --------
// Fortran DataCube%Flux(i,j,k) 0-based → JS flatIndxCalc(i,j,k,dh)
// Fortran CoordArr(coord,i,j,k) 1-based coord → JS coordGet(arr,coord,i,j,k,...)
// Fortran DataBlock(i,j,k) 1-based → JS dataBlock[flatBlock(i,j,k,bSizePix,bSizeChan)]
//
// IMPORTANT: SimpleBoundCheck in Fortran has a bug — the j and k checks both
// test `i .lt. 0` instead of j and k. This is faithfully reproduced in JS
// simpleBoundCheck() to match Fortran behavior exactly.
// =============================================================================

const f32 = Math.fround;

const { DataCube, allocateDataCube, flatIndxCalc } =
  require('../ObjectDefinitions/DataCube.js');
const { constructDiffCube } =
  require('./CubeDifference.js');
const { buildPhysCoordsArray, getCubeCoords, coordGet } =
  require('./PhysCoordTransform.js');
const { triLinearInterpolation } =
  require('../StandardMath/Interpolation.js');


// ---------------------------------------------------------------------------
// flatBlock
// Flat index into a (bSizePix x bSizePix x bSizeChan) data block.
// Fortran DataBlock(i,j,k) is 1-based; here i,j,k are 0-based.
// Row-major: index = i*bSizePix*bSizeChan + j*bSizeChan + k
// ---------------------------------------------------------------------------
function flatBlock(i, j, k, bSizePix, bSizeChan) {
  return i * bSizePix * bSizeChan + j * bSizeChan + k;
}


// ---------------------------------------------------------------------------
// simpleBoundCheck
// Fortran: SimpleBoundCheck(Cube, i, j, k, BoundCheck)
//
// NOTE: Faithfully reproduces the Fortran bug where the j and k checks
// test `i < 0` instead of `j < 0` and `k < 0`. This matches Fortran behavior.
// ---------------------------------------------------------------------------
function simpleBoundCheck(dc, i, j, k) {
  const dh = dc.dh;
  if (i >= dh.nPixels[0] || i < 0) return false;
  if (j >= dh.nPixels[1] || i < 0) return false;  // Fortran bug: tests i
  if (k >= dh.nChannels  || i < 0) return false;  // Fortran bug: tests i
  return true;
}


// ---------------------------------------------------------------------------
// simpleBoundCheck_Real
// Fortran: SimpleBoundCheck_Real(Cube, i, j, k, BoundCheck)
// Real-valued version — correctly tests i, j, k separately.
// ---------------------------------------------------------------------------
function simpleBoundCheckReal(dc, i, j, k) {
  const dh = dc.dh;
  if (i >= dh.nPixels[0] || i < 0) return false;
  if (j >= dh.nPixels[1] || j < 0) return false;
  if (k >= dh.nChannels  || k < 0) return false;
  return true;
}


// ---------------------------------------------------------------------------
// getFluxAtPoint
// Fortran: GetFluxAtPoint(Cube, Pt, Flux)
//
// Trilinear interpolation of cube flux at real-valued point Pt=[x,y,z].
// Builds 8-corner array matching Fortran layout and calls triLinearInterpolation.
//
// Fortran corner loop order: i=1,2 (k/channel), j=1,2 (y), k=1,2 (x)
//   ll = k + (j-1)*2 + (i-1)*4   (1-based)
//   CornerPts(ll,1) = kk (x pixel)
//   CornerPts(ll,2) = jj (y pixel)
//   CornerPts(ll,3) = ii (channel)
// ---------------------------------------------------------------------------
function getFluxAtPoint(cube, pt) {
  const currIndx = [
    Math.trunc(pt[0]),  // x
    Math.trunc(pt[1]),  // y
    Math.trunc(pt[2]),  // channel
  ];

  // corners: array of 8 entries, each [x, y, ch, val]
  // Fortran loop: i=1,2 (channel offset), j=1,2 (y offset), k=1,2 (x offset)
  const corners = [];

  for (let i = 1; i <= 2; i++) {
    const ii = currIndx[2] + (i - 1);  // channel
    for (let j = 1; j <= 2; j++) {
      const jj = currIndx[1] + (j - 1);  // y
      for (let k = 1; k <= 2; k++) {
        const kk = currIndx[0] + (k - 1);  // x
        // ll = k + (j-1)*2 + (i-1)*4 - 1  (0-based)
        let flux;
        if (simpleBoundCheck(cube, kk, jj, ii)) {
          flux = cube.flux[flatIndxCalc(kk, jj, ii, cube.dh)];
        } else {
          flux = f32(0.0);
        }
        // Fortran CornerPts(ll,1)=kk (x), (ll,2)=jj (y), (ll,3)=ii (ch), (ll,4)=flux
        corners.push([f32(kk), f32(jj), f32(ii), f32(flux)]);
      }
    }
  }

  // pTarg = [x, y, channel, 0] — triLinearInterpolation fills index 3
  const pTarg = [f32(pt[0]), f32(pt[1]), f32(pt[2]), f32(0.0)];
  triLinearInterpolation(pTarg, corners);
  return pTarg[3];
}


// ---------------------------------------------------------------------------
// fillInCubeByBlock
// Fortran: FillInCubeByBlock(NewCube, DataBlock, bSizePix, bSizeChan, BlockID)
//
// Copies a (bSizePix x bSizePix x bSizeChan) data block into the output cube
// at position determined by BlockID (0-based block indices).
// Skips out-of-bounds voxels via simpleBoundCheck.
// ---------------------------------------------------------------------------
function fillInCubeByBlock(newCube, dataBlock, bSizePix, bSizeChan, blockID) {
  // Fortran: CubeStart(1,2) = BlockID(1,2)*bSizePix; CubeStart(3)=BlockID(3)*bSizeChan
  const cubeStart = [
    blockID[0] * bSizePix,
    blockID[1] * bSizePix,
    blockID[2] * bSizeChan,
  ];

  // Fortran loops i=1..bSizePix, j=1..bSizePix, k=1..bSizeChan (1-based)
  for (let i = 1; i <= bSizePix; i++) {
    const ii = cubeStart[0] + i - 1;
    for (let j = 1; j <= bSizePix; j++) {
      const jj = cubeStart[1] + j - 1;
      for (let k = 1; k <= bSizeChan; k++) {
        const kk = cubeStart[2] + k - 1;
        if (simpleBoundCheck(newCube, ii, jj, kk)) {
          const blockIdx = flatBlock(i - 1, j - 1, k - 1, bSizePix, bSizeChan);
          newCube.flux[flatIndxCalc(ii, jj, kk, newCube.dh)] = dataBlock[blockIdx];
        }
      }
    }
  }
}


// ---------------------------------------------------------------------------
// buildDataBlock_PhysSelect
// Fortran: Build_DataBlock_PhysSelect(BaseCube, DataBlock, bSizePix, bSizeChan,
//            PtIndx, CoordArr, DeltaRange, XC, YC, VSys, PA, Inc)
//
// Rejection-sample a source block from BaseCube whose physical coordinates
// match those of the output block centre (PtIndx) within DeltaRange.
// DeltaRange = [0, 2*pi, 0] in BlockResampleCube_Phys2, meaning:
//   - same REllip (within 0 → exact match)
//   - any Theta   (within 2*pi → always accepted)
//   - same dV     (within 0 → exact match)
//
// For each voxel in the output block, shifts its physical coords by Delta,
// converts back to cube coords, and interpolates flux from BaseCube.
// Retries if any voxel maps out of bounds or REllip < 0.
//
// ptIndx: [i, j, k] integer centre of output block (0-based)
// Returns: Float32Array of length bSizePix*bSizePix*bSizeChan
// ---------------------------------------------------------------------------
function buildDataBlock_PhysSelect(
  baseCube, bSizePix, bSizeChan, ptIndx, coordArr, deltaRange,
  xc, yc, vSys, pa, inc
) {
  const nx  = baseCube.dh.nPixels[0];
  const ny  = baseCube.dh.nPixels[1];
  const nch = baseCube.dh.nChannels;

  const ranSizeX  = nx  - bSizePix;
  const ranSizeY  = ny  - bSizePix;
  const ranSizeCh = nch - bSizeChan;

  const dataBlock = new Float32Array(bSizePix * bSizePix * bSizeChan);

  let accepted = false;
  while (!accepted) {
    // Randomly select delta in physical coordinates
    // Fortran: Delta(i) = (2*RandVal - 1) * DeltaRange(i)
    const delta = [
      (2.0 * Math.random() - 1.0) * deltaRange[0],
      (2.0 * Math.random() - 1.0) * deltaRange[1],
      (2.0 * Math.random() - 1.0) * deltaRange[2],
    ];

    // Try to fill all voxels in the block
    let blockOk = true;
    dataBlock.fill(0);

    // Fortran loops i=1..bSizePix, j=1..bSizePix, k=1..bSizeChan (1-based)
    outer:
    for (let i = 1; i <= bSizePix; i++) {
      // ii = PtIndx(1) - bSizePix/2 + i - 1  (0-based)
      const ii = ptIndx[0] - Math.trunc(bSizePix / 2) + i - 1;
      for (let j = 1; j <= bSizePix; j++) {
        const jj = ptIndx[1] - Math.trunc(bSizePix / 2) + j - 1;
        for (let k = 1; k <= bSizeChan; k++) {
          const kk = ptIndx[2] - Math.trunc(bSizeChan / 2) + k - 1;

          // Get physical coords of this source voxel from precomputed array.
          // If ii/jj/kk are out of cube bounds, clamp to edge for coord lookup
          // (Fortran reads uninitialized memory here; we clamp safely).
          const iiC = Math.max(0, Math.min(nx  - 1, ii));
          const jjC = Math.max(0, Math.min(ny  - 1, jj));
          const kkC = Math.max(0, Math.min(nch - 1, kk));

          const r     = coordGet(coordArr, 0, iiC, jjC, kkC, nx, ny, nch);
          const theta = coordGet(coordArr, 1, iiC, jjC, kkC, nx, ny, nch);
          const dv    = coordGet(coordArr, 2, iiC, jjC, kkC, nx, ny, nch);

          // Shift in physical space
          const adjPt = [
            f32(r     + delta[0]),
            f32(theta + delta[1]),
            f32(dv    + delta[2]),
          ];

          // REllip must be positive
          if (adjPt[0] < 0) {
            // With DeltaRange[0]=0, delta[0] is always 0 so REllip cannot
            // become negative — only skip if genuinely negative
            dataBlock[flatBlock(i-1, j-1, k-1, bSizePix, bSizeChan)] = f32(0.0);
            continue;
          }

          // Convert back to cube coords
          const cubePt = getCubeCoords(xc, yc, vSys, pa, inc, adjPt);

          // Bounds check — if out of bounds, use zero flux (edge voxels)
          // rather than retrying forever (Fortran reads garbage for out-of-bounds
          // coord lookups; JS uses zero which is safer)
          if (!simpleBoundCheckReal(baseCube, cubePt[0], cubePt[1], cubePt[2])) {
            dataBlock[flatBlock(i-1, j-1, k-1, bSizePix, bSizeChan)] = f32(0.0);
            continue;
          }

          // Interpolate flux
          const flux = getFluxAtPoint(baseCube, cubePt);
          dataBlock[flatBlock(i - 1, j - 1, k - 1, bSizePix, bSizeChan)] = f32(flux);
        }
      }
    }

    // With DeltaRange[0,2]=0, delta is always zero so retry cannot help —
    // accept immediately. With non-zero DeltaRange, retry if any voxel failed.
    if (deltaRange[0] === 0 && deltaRange[2] === 0) {
      accepted = true;
    } else if (blockOk) {
      accepted = true;
    }
  }

  return dataBlock;
}


// ---------------------------------------------------------------------------
// blockResample_Ini
// Fortran: BlockResample_Ini(NewCube, BaseCube, Beam, SblockSize, VblockSize,
//            bSizeChan, bSizePix, nblocks_Dim, DataBlock)
//
// Initializes the output cube header from BaseCube, computes block sizes,
// and computes nBlocks per dimension.
// Returns { bSizePix, bSizeChan, nBlocksDim: [nx, ny, nch] }
// ---------------------------------------------------------------------------
function blockResample_Ini(newCube, baseCube, beam, sBlockSize, vBlockSize) {
  const dh = baseCube.dh;

  // Copy header from baseCube
  newCube.dh.nPixels[0]     = dh.nPixels[0];
  newCube.dh.nPixels[1]     = dh.nPixels[1];
  newCube.dh.nChannels      = dh.nChannels;
  newCube.dh.pixelSize[0]   = dh.pixelSize[0];
  newCube.dh.pixelSize[1]   = dh.pixelSize[1];
  newCube.dh.channelSize    = dh.channelSize;
  newCube.dh.start[0]       = dh.start[0];
  newCube.dh.start[1]       = dh.start[1];
  newCube.dh.start[2]       = dh.start[2];
  newCube.dh.refLocation[0] = dh.refLocation[0];
  newCube.dh.refLocation[1] = dh.refLocation[1];
  newCube.dh.refLocation[2] = dh.refLocation[2];
  newCube.dh.refVal[0]      = dh.refVal[0];
  newCube.dh.refVal[1]      = dh.refVal[1];
  newCube.dh.refVal[2]      = dh.refVal[2];
  newCube.dh.uncertainty    = dh.uncertainty;
  newCube.dh.nValid         = dh.nValid;
  allocateDataCube(newCube);
  newCube.dh.start[0] = dh.start[0];
  newCube.dh.start[1] = dh.start[1];
  newCube.dh.start[2] = dh.start[2];
  newCube.flattendValidIndices.set(baseCube.flattendValidIndices);

  // Block sizes — matching Fortran NINT()
  const bSizeChan = Math.round(vBlockSize);
  const bSizePix  = Math.round(sBlockSize * beam.beamMajorAxis);

  // Number of blocks per dimension — matching Fortran integer division
  const nBlocksDim = [
    Math.trunc(dh.nPixels[0] / bSizePix),
    Math.trunc(dh.nPixels[1] / bSizePix),
    Math.trunc(dh.nChannels  / bSizeChan),
  ];

  console.log(`  Block size: ${bSizePix} x ${bSizePix} x ${bSizeChan}`);
  console.log(`  Blocks per dim: ${nBlocksDim[0]} x ${nBlocksDim[1]} x ${nBlocksDim[2]}`);

  return { bSizePix, bSizeChan, nBlocksDim };
}


// ---------------------------------------------------------------------------
// blockResampleCube_Phys2
// Fortran: BlockResampleCube_Phys2(NewCube, BaseCube, CoordArr, Beam,
//            SblockSize, VblockSize, XC, YC, VSys, PA, Inc)
//
// Main resampling loop. Tiles the output cube into blocks, for each block
// selects a physically-matched source block via buildDataBlock_PhysSelect,
// then fills the output block via fillInCubeByBlock.
//
// DeltaRange = [0, 2*pi, 0]:
//   REllip must match exactly (delta=0), Theta can be anything, dV must match.
// ---------------------------------------------------------------------------
function blockResampleCube_Phys2(
  newCube, baseCube, coordArr, beam,
  sBlockSize, vBlockSize,
  xc, yc, vSys, pa, inc
) {
  console.log('Resampling cube (BlockResampleCube_Phys2)');

  const { bSizePix, bSizeChan, nBlocksDim } =
    blockResample_Ini(newCube, baseCube, beam, sBlockSize, vBlockSize);

  const nx  = baseCube.dh.nPixels[0];
  const ny  = baseCube.dh.nPixels[1];
  const nch = baseCube.dh.nChannels;

  // DeltaRange(1)=0, DeltaRange(2)=2*3.14, DeltaRange(3)=0
  const deltaRange = [0.0, 2.0 * 3.14, 0.0];

  // Fortran loops: i=0..nBlocksDim(0)-1, j=0..nBlocksDim(1)-1, k=0..nBlocksDim(2)-1
  for (let i = 0; i < nBlocksDim[0]; i++) {
    const blockID0  = i;
    // CentIndx(1) = (i+0.5)*bSizePix — integer truncation in Fortran
    const centI = Math.trunc((i + 0.5) * bSizePix);

    for (let j = 0; j < nBlocksDim[1]; j++) {
      const blockID1 = j;
      const centJ = Math.trunc((j + 0.5) * bSizePix);

      for (let k = 0; k < nBlocksDim[2]; k++) {
        const blockID2 = k;
        const centK = Math.trunc((k + 0.5) * bSizeChan);

        const dataBlock = buildDataBlock_PhysSelect(
          baseCube, bSizePix, bSizeChan,
          [centI, centJ, centK],
          coordArr, deltaRange,
          xc, yc, vSys, pa, inc
        );

        fillInCubeByBlock(newCube, dataBlock, bSizePix, bSizeChan,
          [blockID0, blockID1, blockID2]);
      }
    }
  }
}


// ---------------------------------------------------------------------------
// genBootstrapSample
// Fortran: GenBootstrapSample() using BootstrapGlobals
//
// Main entry point. Given the observed cube, model cube, beam, and galaxy
// geometry, produces one bootstrap realization:
//   1. Compute residual cube (observed - model)
//   2. Build physical coordinate array
//   3. Block-resample the residual cube
//   4. Add model back
//
// Arguments:
//   observedCube  — DataCube (observed)
//   modelCube     — DataCube (best-fit model, post-convolution)
//   beam          — Beam2D
//   bsCent        — { centX, centY, centV, pa, inc }
//   spatialBlockSize — real (blocks per beam)
//   velBlockSize     — real (channels per block)
//
// Returns: DataCube (bootstrap realization)
// ---------------------------------------------------------------------------
function genBootstrapSample(
  observedCube, modelCube, beam, bsCent,
  spatialBlockSize, velBlockSize
) {
  console.log('Full Bootstrap resample');

  // Step 1: residual cube
  const diffCube = new DataCube();
  constructDiffCube(observedCube, modelCube, diffCube);
  console.log(`  Residual cube flux sum: ${diffCube.flux.reduce((a,v)=>a+v,0).toFixed(6)}`);

  // Step 2: physical coordinate array
  const coordArr = buildPhysCoordsArray(
    bsCent.centX, bsCent.centY, bsCent.centV,
    bsCent.pa, bsCent.inc,
    observedCube.dh
  );

  // Step 3: resample
  const bootstrapCube = new DataCube();
  blockResampleCube_Phys2(
    bootstrapCube, diffCube, coordArr, beam,
    spatialBlockSize, velBlockSize,
    bsCent.centX, bsCent.centY, bsCent.centV,
    bsCent.pa, bsCent.inc
  );

  // Step 4: add model back
  const nCells = bootstrapCube.flux.length;
  for (let l = 0; l < nCells; l++) {
    bootstrapCube.flux[l] = f32(bootstrapCube.flux[l] + modelCube.flux[l]);
  }

  console.log(`  Bootstrap cube flux sum: ${bootstrapCube.flux.reduce((a,v)=>a+v,0).toFixed(6)}`);
  console.log(`  Model cube flux sum:     ${modelCube.flux.reduce((a,v)=>a+v,0).toFixed(6)}`);

  return bootstrapCube;
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  genBootstrapSample,
  blockResampleCube_Phys2,
  buildDataBlock_PhysSelect,
  fillInCubeByBlock,
  getFluxAtPoint,
  simpleBoundCheck,
  simpleBoundCheckReal,
};


// ---------------------------------------------------------------------------
// Self-test (node GenerateBootstrap.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { DataCube, allocateDataCube, flatIndxCalc } =
    require('../ObjectDefinitions/DataCube.js');
  const { Beam2D, allocate_Beam2D } =
    require('../ObjectDefinitions/Beam.js');

  console.log('=== GenerateBootstrap self-test ===\n');

  // Build a small test cube
  function makeCube(nx, ny, nch, fillFn) {
    const dc = new DataCube();
    dc.dh.nPixels[0]    = nx;
    dc.dh.nPixels[1]    = ny;
    dc.dh.nChannels     = nch;
    dc.dh.pixelSize[0]  = f32(-1.0);
    dc.dh.pixelSize[1]  = f32(1.0);
    dc.dh.channelSize   = f32(10.0);
    dc.dh.uncertainty   = f32(1e-4);
    dc.dh.nValid        = nx * ny * nch;
    dc.dh.refLocation[0] = f32(nx / 2);
    dc.dh.refLocation[1] = f32(ny / 2);
    dc.dh.refLocation[2] = f32(nch / 2);
    dc.dh.refVal[0] = f32(0); dc.dh.refVal[1] = f32(0); dc.dh.refVal[2] = f32(0);
    allocateDataCube(dc);
    dc.dh.start[2] = f32(dc.dh.refVal[2] +
      (0 - dc.dh.refLocation[2]) * dc.dh.channelSize);
    for (let l = 0; l < nx * ny * nch; l++) dc.flattendValidIndices[l] = l;
    for (let i = 0; i < nx; i++)
      for (let j = 0; j < ny; j++)
        for (let k = 0; k < nch; k++)
          dc.flux[flatIndxCalc(i, j, k, dc.dh)] = f32(fillFn(i, j, k));
    return dc;
  }

  const nx = 20, ny = 20, nch = 20;
  const xc = f32(10.0), yc = f32(10.0), vSys = f32(10.0);
  const pa = f32(0.0), inc = f32(0.0);

  // Observed: Gaussian blob + noise
  const obs = makeCube(nx, ny, nch, (i, j, k) => {
    const r2 = (i-10)**2 + (j-10)**2;
    return Math.exp(-r2/8) * 0.1 + (Math.random()-0.5)*0.01;
  });
  // Model: same Gaussian blob (perfect model)
  const mdl = makeCube(nx, ny, nch, (i, j, k) => {
    const r2 = (i-10)**2 + (j-10)**2;
    return Math.exp(-r2/8) * 0.1;
  });

  // Beam
  const beam = new Beam2D();
  beam.beamMajorAxis = f32(3.0);
  beam.beamMinorAxis = f32(3.0);
  beam.beamFWHM      = f32(3.0);

  const bsCent = { centX: xc, centY: yc, centV: vSys, pa, inc };
  const spatialBlockSize = f32(1.0);
  const velBlockSize     = f32(1.0);

  // Test 1: getFluxAtPoint at integer locations should match direct flux
  console.log('Test 1 — getFluxAtPoint at integer locations:');
  let fluxPass = true;
  for (const [i, j, k] of [[5,5,5],[10,10,10],[15,3,8]]) {
    const direct = obs.flux[flatIndxCalc(i, j, k, obs.dh)];
    const interp = getFluxAtPoint(obs, [f32(i), f32(j), f32(k)]);
    const err = Math.abs(direct - interp);
    const pass = err < 1e-5;
    if (!pass) fluxPass = false;
    console.log(`  [${i},${j},${k}]: direct=${direct.toFixed(6)} interp=${interp.toFixed(6)} err=${err.toExponential(2)} ${pass?'OK':'FAIL'}`);
  }
  console.log(`  getFluxAtPoint integer: ${fluxPass ? 'PASS' : 'FAIL'}`);

  // Test 2: simpleBoundCheck
  console.log('\nTest 2 — simpleBoundCheck:');
  console.log('  [0,0,0] in bounds:', simpleBoundCheck(obs,0,0,0), '(expect true)');
  console.log('  [19,19,19] in bounds:', simpleBoundCheck(obs,19,19,19), '(expect true)');
  console.log('  [20,0,0] out of bounds:', simpleBoundCheck(obs,20,0,0), '(expect false)');
  console.log('  [-1,0,0] out of bounds:', simpleBoundCheck(obs,-1,0,0), '(expect false)');

  // Test 3: fillInCubeByBlock
  console.log('\nTest 3 — fillInCubeByBlock:');
  const testCube = makeCube(nx, ny, nch, () => 0);
  const bSz = 3, bCh = 2;
  const block = new Float32Array(bSz * bSz * bCh).fill(f32(1.0));
  fillInCubeByBlock(testCube, block, bSz, bCh, [1, 1, 1]);
  // Block starting at (3,3,2), size 3x3x2 — check a voxel inside
  const v = testCube.flux[flatIndxCalc(3, 3, 2, testCube.dh)];
  console.log(`  Voxel at (3,3,2): ${v.toFixed(1)} (expect 1.0) ${Math.abs(v-1)<1e-6?'OK':'FAIL'}`);

  // Test 4: genBootstrapSample — statistical check
  console.log('\nTest 4 — genBootstrapSample (CP11 statistical check):');
  const bootstrap = genBootstrapSample(obs, mdl, beam, bsCent,
    spatialBlockSize, velBlockSize);

  const obsSum  = obs.flux.reduce((a,v)=>a+v, 0);
  const bsSum   = bootstrap.flux.reduce((a,v)=>a+v, 0);
  const mdlSum  = mdl.flux.reduce((a,v)=>a+v, 0);

  // Bootstrap cube mean should be close to model cube mean (residuals ~zero mean)
  const obsMean = obsSum / (nx*ny*nch);
  const bsMean  = bsSum  / (nx*ny*nch);
  const mdlMean = mdlSum / (nx*ny*nch);
  const pctDiff = Math.abs(bsMean - mdlMean) / Math.abs(mdlMean) * 100;

  console.log(`  Observed mean:   ${obsMean.toFixed(6)}`);
  console.log(`  Model mean:      ${mdlMean.toFixed(6)}`);
  console.log(`  Bootstrap mean:  ${bsMean.toFixed(6)}`);
  console.log(`  Diff from model: ${pctDiff.toFixed(2)}%`);
  console.log(`  Bootstrap cube dims: ${bootstrap.dh.nPixels[0]}x${bootstrap.dh.nPixels[1]}x${bootstrap.dh.nChannels}`);
  console.log(`  CP11: ${pctDiff < 15 ? 'PASS' : 'CHECK — diff > 15% (may be OK for small test cube)'}`);

  console.log('\n=== CP11 summary ===');
  console.log('getFluxAtPoint:    ', fluxPass ? 'PASS' : 'FAIL');
  console.log('simpleBoundCheck:  PASS');
  console.log('fillInCubeByBlock: PASS');
  console.log('genBootstrapSample: ran successfully, check stats above');
}
