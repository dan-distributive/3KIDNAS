'use strict';

// =============================================================================
// Beam.js
// High-fidelity port of src/ObjectDefinitions/Beam.f (BeamMod)
//
// PORTING NOTES
// -------------
// Fortran `real` → Float32Array / Math.fround().
// Fortran `double complex` 2D array (ComplexKernel) → two Float64Array
//   (real and imaginary parts interleaved or separate — stored here as
//   a single Float64Array of length 2*size with [re, im, re, im, ...] layout,
//   matching the convention expected by most JS FFT libraries).
//
// CommonConsts Pi → Math.fround(Math.PI) for f32 arithmetic,
//   plain Math.PI for f64 arithmetic.
//
// Kernel is a 2D array indexed (-nRadialCells:nRadialCells) in Fortran.
// In JS we store it as a Float32Array of size (2*nRadialCells+1)^2 with
// helper functions kernelGet/kernelSet that offset the index by nRadialCells.
//
// PaddedSize and ComplexSize use 1-based indexing in Fortran (integer arrays
// declared as PaddedSize(2), ComplexSize(2)) — ported as 0-indexed [0],[1].
// Fortran PaddedSize(1) → JS paddedSize[0], PaddedSize(2) → paddedSize[1].
// =============================================================================

const f32 = Math.fround;
const Pi  = f32(Math.PI);   // Fortran CommonConsts Pi as real


// ---------------------------------------------------------------------------
// Beam2D
// Fortran TYPE Beam2D — 2D Gaussian beam descriptor + kernel arrays.
// ---------------------------------------------------------------------------
class Beam2D {
  constructor() {
    this.beamMajorAxis       = f32(0.0);
    this.beamMinorAxis       = f32(0.0);
    this.beamPositionAngle   = f32(0.0);
    this.beamFWHM            = f32(0.0);
    this.sigmaLengths        = f32(0.0);
    this.beamSigmaVector     = new Float32Array(3);   // (0:2)
    this.pixelSize           = new Float32Array(2);   // (0:1)
    this.nRadialCells        = 0;
    this.kernel              = null;   // Float32Array, (2n+1)^2, offset-indexed
    this.velocitySmoothSwitch = 0;
    this.velocitySmoothSigma = f32(0.0);
    this.beamUnitsSwitch     = new Int32Array(2);     // (0:1)
    this.beamAreaPixels      = f32(0.0);
    this.beamAreaUnits       = f32(0.0);
    this.paddedSize          = new Int32Array(2);     // (2) → [0],[1]
    this.complexSize         = new Int32Array(2);     // (2) → [0],[1]
    this.complexKernel       = null;   // Float64Array interleaved [re,im,...]
    this.complexKernelCreated = false;
  }
}


// ---------------------------------------------------------------------------
// allocate_Beam2D
// Fortran: Allocate_Beam2D(B, nCubePixels)
//
// nCubePixels: Int32Array or array of length 2 — [nPixX, nPixY]
//
// Sets sigma values from FWHM (if beamFWHM > 0) or from major/minor axes.
// Calculates beam area in pixels and in world units.
// Allocates kernel and complexKernel arrays.
//
// Fortran kernel is indexed (-nRadialCells:nRadialCells, -nRadialCells:nRadialCells).
// JS kernel is a flat Float32Array of size (2n+1)^2 — use kernelGet/kernelSet.
//
// Fortran ComplexKernel(ComplexSize(1), ComplexSize(2)) is double complex.
// JS complexKernel is Float64Array of length 2 * complexSize[0] * complexSize[1]
// with interleaved [re, im] pairs.
// ---------------------------------------------------------------------------
function allocate_Beam2D(b, nCubePixels) {
  const FWHM_TO_SIGMA = f32(2.355);

  if (b.beamFWHM > f32(0.0)) {
    // FWHM mode — circular beam
    b.beamSigmaVector[0]  = f32(b.beamFWHM / FWHM_TO_SIGMA);
    b.beamSigmaVector[1]  = f32(b.beamFWHM / FWHM_TO_SIGMA);
    b.beamSigmaVector[2]  = f32(0.0);
    b.beamMajorAxis       = b.beamFWHM;
    b.beamMinorAxis       = b.beamFWHM;
  } else {
    // Axis mode — elliptical beam, axes set prior to call
    b.beamSigmaVector[0]  = f32(b.beamMajorAxis / FWHM_TO_SIGMA);
    b.beamSigmaVector[1]  = f32(b.beamMinorAxis / FWHM_TO_SIGMA);
    b.beamSigmaVector[2]  = b.beamPositionAngle;
  }

  // Beam area in pixels
  // Fortran: 2.*Pi/2.355**2. * BeamMajorAxis * BeamMinorAxis
  b.beamAreaPixels = f32(
    f32(f32(2.0) * Pi) / f32(FWHM_TO_SIGMA * FWHM_TO_SIGMA)
    * f32(b.beamMajorAxis * b.beamMinorAxis)
  );

  // Beam area in world units
  b.beamAreaUnits = f32(
    b.beamAreaPixels
    * f32(Math.abs(b.pixelSize[0]) * Math.abs(b.pixelSize[1]))
  );

  // Number of radial cells — half-width of kernel in pixels
  b.nRadialCells = Math.trunc(
    Math.abs(f32(b.beamSigmaVector[0] * b.sigmaLengths))
  ) + 1;

  // Allocate real kernel: (2n+1) x (2n+1)
  const kSize   = 2 * b.nRadialCells + 1;
  b.kernel      = new Float32Array(kSize * kSize);

  b.complexKernelCreated = false;

  // PaddedSize = 2*nRadialCells + 1 + nCubePixels  (per axis)
  // Fortran: B%PaddedSize = 2*B%nRadialCells+1 + nCubePixels
  // Fortran arrays PaddedSize(2), ComplexSize(2) are 1-indexed
  // → stored here 0-indexed: [0]=axis1, [1]=axis2
  b.paddedSize[0] = 2 * b.nRadialCells + 1 + nCubePixels[0];
  b.paddedSize[1] = 2 * b.nRadialCells + 1 + nCubePixels[1];

  // ComplexSize(1) = PaddedSize(1)/2 + 1  (r2c FFT output)
  // ComplexSize(2) = PaddedSize(2)
  b.complexSize[0] = Math.trunc(b.paddedSize[0] / 2) + 1;
  b.complexSize[1] = b.paddedSize[1];

  // Allocate complex kernel: interleaved [re, im] Float64Array
  const cLen       = b.complexSize[0] * b.complexSize[1];
  b.complexKernel  = new Float64Array(2 * cLen);
}


// ---------------------------------------------------------------------------
// deAllocate_Beam2D
// Fortran: DeAllocate_Beam2D(B)
// ---------------------------------------------------------------------------
function deAllocate_Beam2D(b) {
  b.kernel               = null;
  b.complexKernel        = null;
  b.complexKernelCreated = false;
}


// ---------------------------------------------------------------------------
// Kernel index helpers
// Fortran kernel is indexed B%Kernel(-n:n, -n:n).
// In JS, offset by nRadialCells to map to 0-based flat array.
//
// kernelIdx(i, j, n) → flat index for kernel[i][j] where i,j ∈ [-n, n]
// ---------------------------------------------------------------------------
function kernelIdx(i, j, nRadialCells) {
  const kSize = 2 * nRadialCells + 1;
  return (i + nRadialCells) * kSize + (j + nRadialCells);
}

function kernelGet(b, i, j) {
  return b.kernel[kernelIdx(i, j, b.nRadialCells)];
}

function kernelSet(b, i, j, val) {
  b.kernel[kernelIdx(i, j, b.nRadialCells)] = f32(val);
}


// ---------------------------------------------------------------------------
// Complex kernel index helpers
// Fortran: ComplexKernel(row, col) 1-indexed, double complex.
// JS: interleaved Float64Array, 0-indexed.
//
// complexKernelGet(b, row, col) → { re, im }  (row, col 0-indexed)
// complexKernelSet(b, row, col, re, im)
// ---------------------------------------------------------------------------
function complexKernelIdx(row, col, b) {
  return 2 * (row * b.complexSize[1] + col);
}

function complexKernelGet(b, row, col) {
  const idx = complexKernelIdx(row, col, b);
  return { re: b.complexKernel[idx], im: b.complexKernel[idx + 1] };
}

function complexKernelSet(b, row, col, re, im) {
  const idx = complexKernelIdx(row, col, b);
  b.complexKernel[idx]     = re;
  b.complexKernel[idx + 1] = im;
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  Beam2D,
  allocate_Beam2D,
  deAllocate_Beam2D,
  kernelGet,
  kernelSet,
  kernelIdx,
  complexKernelGet,
  complexKernelSet
};


// ---------------------------------------------------------------------------
// Self-test (node Beam.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { Beam2D, allocate_Beam2D, kernelGet, kernelSet } = module.exports;

  console.log('=== FWHM mode ===');
  const b1 = new Beam2D();
  b1.beamFWHM      = f32(3.0);   // 3 pixel FWHM circular beam
  b1.sigmaLengths  = f32(3.0);   // kernel extends 3 sigma
  b1.pixelSize[0]  = f32(-1.5);
  b1.pixelSize[1]  = f32(1.5);
  allocate_Beam2D(b1, new Int32Array([64, 64]));

  const sigma = f32(3.0 / 2.355);
  console.log('beamSigmaVector[0]:', b1.beamSigmaVector[0].toFixed(6),
              '(expect', sigma.toFixed(6), ')');
  console.log('beamSigmaVector[2]:', b1.beamSigmaVector[2], '(expect 0)');
  console.log('beamMajorAxis:', b1.beamMajorAxis, '(expect 3)');
  console.log('beamAreaPixels:', b1.beamAreaPixels.toFixed(6));
  console.log('nRadialCells:', b1.nRadialCells);
  console.log('paddedSize:', Array.from(b1.paddedSize));
  console.log('complexSize:', Array.from(b1.complexSize));
  console.log('kernel length:', b1.kernel.length,
              '(expect', (2*b1.nRadialCells+1)**2, ')');

  console.log('\n=== axis mode ===');
  const b2 = new Beam2D();
  b2.beamFWHM          = f32(-1.0);   // negative → axis mode
  b2.beamMajorAxis     = f32(4.0);
  b2.beamMinorAxis     = f32(2.0);
  b2.beamPositionAngle = f32(0.5);
  b2.sigmaLengths      = f32(3.0);
  b2.pixelSize[0]      = f32(-1.0);
  b2.pixelSize[1]      = f32(1.0);
  allocate_Beam2D(b2, new Int32Array([32, 32]));

  console.log('beamSigmaVector[0]:', b2.beamSigmaVector[0].toFixed(6),
              '(expect', f32(4.0/2.355).toFixed(6), ')');
  console.log('beamSigmaVector[1]:', b2.beamSigmaVector[1].toFixed(6),
              '(expect', f32(2.0/2.355).toFixed(6), ')');
  console.log('beamSigmaVector[2]:', b2.beamSigmaVector[2], '(expect 0.5)');

  console.log('\n=== kernel indexing ===');
  kernelSet(b1, 0, 0, 1.0);
  kernelSet(b1, -b1.nRadialCells, -b1.nRadialCells, 99.0);
  kernelSet(b1,  b1.nRadialCells,  b1.nRadialCells, 77.0);
  console.log('kernel(0,0):', kernelGet(b1, 0, 0), '(expect 1)');
  console.log('kernel(-n,-n):', kernelGet(b1, -b1.nRadialCells, -b1.nRadialCells), '(expect 99)');
  console.log('kernel(n,n):', kernelGet(b1, b1.nRadialCells, b1.nRadialCells), '(expect 77)');
}
