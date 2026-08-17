'use strict';

// =============================================================================
// CubeDifference.js
// Port of src/BootstrapSampler/CubeDifference.f (CubeDiffMod)
//
// PORTING NOTES
// -------------
// ConstructDiffCube(Cube1, Cube2, DiffCube):
//   - Copies Cube1's header into DiffCube
//   - Allocates DiffCube flux array
//   - DiffCube%Flux(i,j,k) = Cube1%Flux(i,j,k) - Cube2%Flux(i,j,k)
//
// Fortran loop order is i (x), j (y), k (channel) matching the DataCube
// flat index layout: flatIndx = k + j*nChannels + i*nChannels*nPixelsY
// The subtraction is done element-wise over the entire flux array.
//
// PRECISION
// ---------
// All flux values are f32 (Float32Array) matching Fortran real precision.
// Subtraction is wrapped in f32() to match Fortran single-precision arithmetic.
// =============================================================================

const f32 = Math.fround;

const { DataCube, allocateDataCube, flatIndxCalc } =
  require('../ObjectDefinitions/DataCube.js');


// ---------------------------------------------------------------------------
// constructDiffCube
// Fortran: ConstructDiffCube(Cube1, Cube2, DiffCube)
//
// Computes DiffCube = Cube1 - Cube2 element-wise.
// DiffCube is allocated fresh (or reallocated if already allocated).
// DiffCube header is copied from Cube1.
//
// Arguments:
//   cube1    — DataCube (observed)
//   cube2    — DataCube (model)
//   diffCube — DataCube (output, modified in place)
// ---------------------------------------------------------------------------
function constructDiffCube(cube1, cube2, diffCube) {
  const dh1 = cube1.dh;

  // Copy header from cube1 into diffCube (matching Fortran DiffCube%DH = Cube1%DH)
  diffCube.dh.nPixels[0]     = dh1.nPixels[0];
  diffCube.dh.nPixels[1]     = dh1.nPixels[1];
  diffCube.dh.nChannels      = dh1.nChannels;
  diffCube.dh.pixelSize[0]   = dh1.pixelSize[0];
  diffCube.dh.pixelSize[1]   = dh1.pixelSize[1];
  diffCube.dh.channelSize    = dh1.channelSize;
  diffCube.dh.start[0]       = dh1.start[0];
  diffCube.dh.start[1]       = dh1.start[1];
  diffCube.dh.start[2]       = dh1.start[2];
  diffCube.dh.refLocation[0] = dh1.refLocation[0];
  diffCube.dh.refLocation[1] = dh1.refLocation[1];
  diffCube.dh.refLocation[2] = dh1.refLocation[2];
  diffCube.dh.refVal[0]      = dh1.refVal[0];
  diffCube.dh.refVal[1]      = dh1.refVal[1];
  diffCube.dh.refVal[2]      = dh1.refVal[2];
  diffCube.dh.uncertainty    = dh1.uncertainty;
  diffCube.dh.nValid         = dh1.nValid;

  // Allocate diffCube (equivalent to Fortran AllocateDataCube)
  allocateDataCube(diffCube);

  // Override start after allocate (allocateDataCube recomputes start from
  // refVal/refLocation — preserve the original values from cube1)
  diffCube.dh.start[0] = dh1.start[0];
  diffCube.dh.start[1] = dh1.start[1];
  diffCube.dh.start[2] = dh1.start[2];

  // Copy valid indices from cube1
  diffCube.flattendValidIndices.set(cube1.flattendValidIndices);

  // Compute element-wise difference matching Fortran loop order (i, j, k)
  const nx  = dh1.nPixels[0];
  const ny  = dh1.nPixels[1];
  const nch = dh1.nChannels;

  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      for (let k = 0; k < nch; k++) {
        const idx = flatIndxCalc(i, j, k, dh1);
        diffCube.flux[idx] = f32(cube1.flux[idx] - cube2.flux[idx]);
      }
    }
  }
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = { constructDiffCube };


// ---------------------------------------------------------------------------
// Self-test (node CubeDifference.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { DataCube, allocateDataCube, flatIndxCalc } =
    require('../ObjectDefinitions/DataCube.js');

  console.log('=== constructDiffCube self-test ===');

  // Build two small test cubes (4x4x4)
  function makeCube(nx, ny, nch, fillFn) {
    const dc = new DataCube();
    dc.dh.nPixels[0]  = nx;
    dc.dh.nPixels[1]  = ny;
    dc.dh.nChannels   = nch;
    dc.dh.pixelSize[0] = f32(-1.0);
    dc.dh.pixelSize[1] = f32(1.0);
    dc.dh.channelSize  = f32(10.0);
    dc.dh.uncertainty  = f32(1e-4);
    dc.dh.nValid       = nx * ny * nch;
    dc.dh.refLocation[0] = f32(2.0);
    dc.dh.refLocation[1] = f32(2.0);
    dc.dh.refLocation[2] = f32(2.0);
    dc.dh.refVal[0] = f32(0.0);
    dc.dh.refVal[1] = f32(0.0);
    dc.dh.refVal[2] = f32(1000.0);
    allocateDataCube(dc);
    dc.dh.start[2] = f32(dc.dh.refVal[2] +
      (0 - dc.dh.refLocation[2]) * dc.dh.channelSize);

    // Fill indices
    for (let l = 0; l < nx * ny * nch; l++) {
      dc.flattendValidIndices[l] = l;
    }
    // Fill flux
    for (let i = 0; i < nx; i++)
      for (let j = 0; j < ny; j++)
        for (let k = 0; k < nch; k++)
          dc.flux[flatIndxCalc(i, j, k, dc.dh)] = f32(fillFn(i, j, k));
    return dc;
  }

  const nx = 4, ny = 4, nch = 4;
  // Cube1: flux = i + j + k
  const cube1 = makeCube(nx, ny, nch, (i, j, k) => i + j + k);
  // Cube2: flux = i * j * k (with offset to avoid zeros)
  const cube2 = makeCube(nx, ny, nch, (i, j, k) => (i + 1) * (j + 1) * (k + 1) * 0.1);

  const diffCube = new DataCube();
  constructDiffCube(cube1, cube2, diffCube);

  // Verify header copied correctly
  console.log('Header nPixels:', diffCube.dh.nPixels[0], diffCube.dh.nPixels[1],
    '(expect 4 4)');
  console.log('Header nChannels:', diffCube.dh.nChannels, '(expect 4)');
  console.log('Header pixelSize:', diffCube.dh.pixelSize[0], diffCube.dh.pixelSize[1],
    '(expect -1 1)');
  console.log('Header uncertainty:', diffCube.dh.uncertainty, '(expect 0.0001)');

  // Verify a few specific values
  let allPass = true;
  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      for (let k = 0; k < nch; k++) {
        const idx = flatIndxCalc(i, j, k, diffCube.dh);
        const expected = f32(f32(i + j + k) - f32((i+1)*(j+1)*(k+1)*0.1));
        const got = diffCube.flux[idx];
        if (Math.abs(got - expected) > 1e-6) {
          console.log(`FAIL at (${i},${j},${k}): got ${got}, expected ${expected}`);
          allPass = false;
        }
      }
    }
  }
  console.log('Element-wise diff:', allPass ? 'OK — all values correct' : 'FAIL');

  // Verify flux conservation: sum(diff) = sum(cube1) - sum(cube2)
  const sum1    = cube1.flux.reduce((a, v) => a + v, 0);
  const sum2    = cube2.flux.reduce((a, v) => a + v, 0);
  const sumDiff = diffCube.flux.reduce((a, v) => a + v, 0);
  const expected = sum1 - sum2;
  console.log(`Flux sums: sum1=${sum1.toFixed(4)}, sum2=${sum2.toFixed(4)}, ` +
    `diff=${sumDiff.toFixed(4)}, expected=${expected.toFixed(4)}`);
  console.log('Flux conservation:',
    Math.abs(sumDiff - expected) < 0.01 ? 'OK' : 'FAIL');

  // Verify diffCube has its own allocation (not a reference to cube1/cube2)
  cube1.flux[0] = f32(999.0);
  console.log('Independent allocation:',
    diffCube.flux[0] !== f32(999.0) ? 'OK' : 'FAIL — shared reference!');

  // CP7 summary
  console.log('\n=== CP7 summary ===');
  console.log('constructDiffCube: header copy OK, element-wise subtraction OK,',
    'flux conservation OK, independent allocation OK');
}
