'use strict';

// =============================================================================
// CubeComparison.js
// High-fidelity port of src/CompareCubes/CubeComparison.f (CubeCompareMod)
//
// PORTING NOTES
// -------------
// Fortran `real` → Math.fround(). Accumulation in f32 matching Fortran.
//
// Only non-NaN valid voxels are compared, via Cube1.flattendValidIndices.
// These are pre-computed at cube load time (allocateDataCube sets all cells
// valid by default; the pipeline overwrites with actual valid indices).
//
// LikePoint is a function pointer in Fortran — defaults to chi2Calc.
// In JS we pass it as an explicit argument with chi2Calc as default,
// matching the pipeline's standard usage.
//
// ThreeDIndxCalc is used to convert flat indices to (i,j,k) — imported
// from DataCube.js.
// =============================================================================

const f32 = Math.fround;
const { threeDIndxCalc, flatIndxCalc } = require('../ObjectDefinitions/DataCube.js');
const { chi2Calc } = require('../CompareCubes/LikelihoodFunctions.js');


// ---------------------------------------------------------------------------
// cubeCompare
// Fortran: CubeCompare(Cube1, Cube2, chi2, sigma)
//
// Gathers flux values at all valid voxels (via flattendValidIndices),
// then computes chi² = sum((f1-f2)^2 / sigma^2) over those voxels.
//
// sigma: uniform noise level (scalar, same for all voxels)
// likePoint: function matching chi2Calc signature — defaults to chi2Calc.
//
// Returns chi2 as f32.
// ---------------------------------------------------------------------------
function cubeCompare(cube1, cube2, sigma, likePoint = chi2Calc) {
  const nValid = cube1.dh.nValid;
  const f1     = new Float32Array(nValid);
  const f2     = new Float32Array(nValid);
  const unc    = new Float32Array(nValid).fill(f32(sigma));

  for (let l = 0; l < nValid; l++) {
    const lUse      = cube1.flattendValidIndices[l];
    const { i, j, k } = threeDIndxCalc(lUse, cube1.dh);
    const idx       = flatIndxCalc(i, j, k, cube1.dh);
    f1[l]           = f32(cube1.flux[idx]);
    f2[l]           = f32(cube2.flux[idx]);
  }

  return likePoint(nValid, f1, f2, unc);
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = { cubeCompare };


// ---------------------------------------------------------------------------
// Self-test (node CubeComparison.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { DataCube, allocateDataCube, flatIndxCalc: fi }
    = require('../ObjectDefinitions/DataCube.js');

  // Build two identical small cubes then introduce a known difference
  function makeCube() {
    const dc = new DataCube();
    const dh = dc.dh;
    dh.nPixels[0] = 4; dh.nPixels[1] = 4; dh.nChannels = 4;
    dh.pixelSize[0] = -1; dh.pixelSize[1] = 1;
    dh.channelSize  = 10;
    dh.refLocation[0] = 2; dh.refLocation[1] = 2; dh.refLocation[2] = 2;
    dh.refVal[0] = 0; dh.refVal[1] = 0; dh.refVal[2] = 1000;
    allocateDataCube(dc);
    return dc;
  }

  const c1 = makeCube();
  const c2 = makeCube();

  // Fill c1 and c2 with known values
  for (let i = 0; i < 4; i++)
    for (let j = 0; j < 4; j++)
      for (let k = 0; k < 4; k++) {
        c1.flux[fi(i,j,k,c1.dh)] = f32(1.0);
        c2.flux[fi(i,j,k,c2.dh)] = f32(1.1);
      }

  const sigma = f32(0.1);

  console.log('=== cubeCompare ===');
  const chi2 = cubeCompare(c1, c2, sigma);
  // Expected: nValid=64, diff=0.1, sigma=0.1 → chi2 = 64*(0.1^2/0.1^2) = 64
  console.log('chi2:', chi2.toFixed(4), '(expect 64.0)');
  console.log('nValid:', c1.dh.nValid);

  // Test with identical cubes — chi2 should be 0
  const chi2Zero = cubeCompare(c1, c1, sigma);
  console.log('chi2 (identical cubes):', chi2Zero.toFixed(4), '(expect 0.0)');
}
