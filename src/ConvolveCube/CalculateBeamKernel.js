'use strict';

// =============================================================================
// CalculateBeamKernel.js
// High-fidelity port of src/ConvolveCube/CalculateBeamKernel.f
// (CalcBeamKernelMod)
//
// PORTING NOTES
// -------------
// Fortran `real` → Math.fround(). All trig/exp operate on f32 inputs,
// output wrapped in f32.
//
// Calculate2DBeamKernel: pure math, ported exactly.
//
// CalculateComplex2DKernel: uses FFTW3 in Fortran. In JS this is replaced
// by a call to a JS FFT library (fft.js or ndarray-fft) in
// CubeKernelConvolution.js. The MakeWrappedArray helper is ported here
// since it is called during the complex kernel setup.
//
// Kernel indexing: Fortran B%Kernel(-n:n, -n:n) → JS flat Float32Array
// via kernelGet/kernelSet from Beam.js (offset by nRadialCells).
//
// The renormalization B%Kernel = B%Kernel / sum(B%Kernel) is done in
// f32 accumulation matching Fortran array reduction.
// =============================================================================

const f32 = Math.fround;
const Pi  = f32(Math.PI);

const { kernelGet, kernelSet } = require('../ObjectDefinitions/Beam.js');


// ---------------------------------------------------------------------------
// calculate2DBeamKernel
// Fortran: Calculate2DBeamKernel(B, PixelSizes)
//
// Fills B%Kernel with a 2D rotated Gaussian evaluated at each pixel offset
// (-nRadialCells:nRadialCells) in both axes, then normalizes so sum=1.
//
// PixelSizes is accepted but not used (matching commented-out Fortran code).
//
// Gaussian:
//   cpa = cos(-BeamSigmaVector[2])   (negative rotation to beam axes)
//   spa = sin(-BeamSigmaVector[2])
//   xp  = x*cpa - y*spa
//   yp  = x*spa + y*cpa
//   R2  = (xp/sigma0)^2 + (yp/sigma1)^2
//   K   = 1/sqrt(2*Pi*sigma0*sigma1) * exp(-R2/2)
// Then K = K / sum(K)
// ---------------------------------------------------------------------------
function calculate2DBeamKernel(b, pixelSizes) {
  const sigma0 = f32(b.beamSigmaVector[0]);
  const sigma1 = f32(b.beamSigmaVector[1]);
  const pa     = f32(b.beamSigmaVector[2]);

  const cpa = f32(Math.cos(f32(-pa)));
  const spa = f32(Math.sin(f32(-pa)));

  const n = b.nRadialCells;

  // Fill kernel with Gaussian values
  for (let i = -n; i <= n; i++) {
    for (let j = -n; j <= n; j++) {
      const x  = f32(i);
      const y  = f32(j);
      const xp = f32(f32(x * cpa) - f32(y * spa));
      const yp = f32(f32(x * spa) + f32(y * cpa));
      const R2 = f32(
        f32(f32(xp / sigma0) * f32(xp / sigma0))
        + f32(f32(yp / sigma1) * f32(yp / sigma1))
      );
      const norm = f32(
        f32(1.0) / f32(Math.sqrt(f32(f32(2.0) * Pi * f32(sigma0 * sigma1))))
      );
      const val  = f32(norm * f32(Math.exp(f32(f32(-R2) / f32(2.0)))));
      kernelSet(b, i, j, val);
    }
  }

  // Renormalize so sum = 1
  // Fortran: B%Kernel = B%Kernel / sum(B%Kernel)
  // sum accumulated in f32 matching Fortran array reduction
  const kSize = 2 * n + 1;
  let kernelSum = f32(0.0);
  for (let idx = 0; idx < kSize * kSize; idx++) {
    kernelSum = f32(kernelSum + f32(b.kernel[idx]));
  }
  for (let idx = 0; idx < kSize * kSize; idx++) {
    b.kernel[idx] = f32(b.kernel[idx] / kernelSum);
  }
}


// ---------------------------------------------------------------------------
// makeWrappedArray
// Fortran: MakeWrappedArray(SA, CV, Arr, WrappedArr)
//
// Wraps a 2D array so that the element at CV becomes the origin (1,1).
// Used during complex kernel setup before FFT.
//
// SA: [sizeX, sizeY]  (1-indexed in Fortran → 0-indexed here)
// CV: [centX, centY]  (1-indexed centre values)
// arr:        Float64Array, row-major, size SA[0]*SA[1]
// wrappedArr: Float64Array, same size, output
//
// Fortran indexing (1-based):
//   k = i - CV[0];  if k<=0: k = SA[0] + k
//   l = j - CV[1];  if l<=0: l = SA[1] + l
//   WrappedArr(k,l) = Arr(i,j)
// Row-major storage: Arr(i,j) → arr[(i-1)*sy + (j-1)]
//                   WrappedArr(k,l) → wrappedArr[(k-1)*sy + (l-1)]
// ---------------------------------------------------------------------------
function makeWrappedArray(sa, cv, arr, wrappedArr) {
  const sx = sa[0];
  const sy = sa[1];

  for (let i = 1; i <= sx; i++) {
    for (let j = 1; j <= sy; j++) {
      let k = i - cv[0];
      let l = j - cv[1];
      if (k <= 0) k = sx + k;
      if (l <= 0) l = sy + l;
      // Row-major: arr[i][j] → arr[(i-1)*sy + (j-1)]
      // wrappedArr[k][l] → wrappedArr[(k-1)*sy + (l-1)]
      wrappedArr[(k - 1) * sy + (l - 1)] = arr[(i - 1) * sy + (j - 1)];
    }
  }
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  calculate2DBeamKernel,
  makeWrappedArray
};


// ---------------------------------------------------------------------------
// Self-test (node CalculateBeamKernel.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { Beam2D, allocate_Beam2D, kernelGet } = require('../ObjectDefinitions/Beam.js');

  console.log('=== calculate2DBeamKernel (circular beam) ===');
  const b = new Beam2D();
  b.beamFWHM     = f32(3.0);
  b.sigmaLengths = f32(3.0);
  b.pixelSize[0] = f32(-1.0);
  b.pixelSize[1] = f32(1.0);
  allocate_Beam2D(b, new Int32Array([64, 64]));
  calculate2DBeamKernel(b, b.pixelSize);

  // Sum should be 1
  const n = b.nRadialCells;
  const kSize = 2 * n + 1;
  let sum = f32(0.0);
  for (let idx = 0; idx < kSize * kSize; idx++) sum = f32(sum + b.kernel[idx]);
  console.log('kernel sum:', sum.toFixed(8), '(expect 1.0)');

  // Centre value should be maximum
  const centre = kernelGet(b, 0, 0);
  const corner = kernelGet(b, -n, -n);
  console.log('kernel(0,0):', centre.toExponential(6), '(expect max)');
  console.log('kernel(-n,-n):', corner.toExponential(6), '(expect min)');
  console.log('centre > corner:', centre > corner ? 'OK' : 'FAIL');

  console.log('\n=== makeWrappedArray ===');
  // 4x4 array, centre at (2,2) — element (2,2) should map to (1,1)
  const sa = [4, 4];
  const cv = [2, 2];
  const arr        = new Float64Array(16);
  const wrappedArr = new Float64Array(16);
  for (let i = 0; i < 16; i++) arr[i] = i + 1;
  makeWrappedArray(sa, cv, arr, wrappedArr);
  // arr(2,2) maps to wrappedArr(4,4); arr(3,3) maps to wrappedArr(1,1)
  // arr(3,3) in row-major = arr[(3-1)*4+(3-1)] = arr[10] = 11
  console.log('arr(2,2)=', arr[(2-1)*4+(2-1)], '→ wrappedArr(1,1)=', wrappedArr[0], '(expect 11)');
}
