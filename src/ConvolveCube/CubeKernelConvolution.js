'use strict';

// =============================================================================
// CubeKernelConvolution.js
// High-fidelity port of src/ConvolveCube/CubeKernelConvolution.f
// (CubeKernelConvolutionMod)
//
// PORTING NOTES
// -------------
// Uses ndarray-fft for FFT convolution (arbitrary sizes, no power-of-two required).
// Supports arbitrary array sizes matching the Fortran paddedSize exactly.
//
// FFTW halfcomplex format (output of RFFT for size n):
//   out[k]   = re[k]   for k = 0..n/2
//   out[n-k] = im[k]   for k = 1..n/2-1
//   im[0] = 0, im[n/2] = 0 (DC and Nyquist are purely real)
//
// 2D FFT is performed as row-FFTs then column-FFTs, matching how FFTW
// implements r2c_2d internally.
//
//
// PRECISION
// ---------
// for the FFT (double precision paddedKernel, double complex output).
// Results will agree in physical correctness but not be bit-identical.
// =============================================================================

const f32  = Math.fround;

const { calculate2DBeamKernel, makeWrappedArray } = require('./CalculateBeamKernel.js');
const fft     = require('ndarray-fft');
const ndarray = require('ndarray');

const { flatIndxCalc } = require('../ObjectDefinitions/DataCube.js');


// ---------------------------------------------------------------------------
// halfComplexMultiply
// Pointwise multiply two 2D halfcomplex arrays (row-major, size ps0*ps1).
// FFTW halfcomplex format per row: re[0..n/2], im stored at [n-k] for k=1..n/2-1
// For 2D we treat each row of the 2D halfcomplex output independently.
// Returns result as Float32Array in halfcomplex format.
// ---------------------------------------------------------------------------
function halfComplexMultiply2D(a, b, ps0, ps1) {
  const out = new Float32Array(ps0 * ps1);
  for (let r = 0; r < ps0; r++) {
    const base = r * ps1;
    const n    = ps1;
    const half = Math.floor(n / 2);

    // DC (k=0): real only
    out[base] = a[base] * b[base];

    // Nyquist (k=n/2): real only (only if n even)
    if (n % 2 === 0) {
      out[base + half] = a[base + half] * b[base + half];
    }

    // k=1..half-1: re at [k], im at [n-k]
    for (let k = 1; k < half; k++) {
      const reA = a[base + k];
      const imA = a[base + n - k];
      const reB = b[base + k];
      const imB = b[base + n - k];
      out[base + k]     = reA * reB - imA * imB;
      out[base + n - k] = reA * imB + imA * reB;
    }
  }
  return out;
}


// ---------------------------------------------------------------------------
// fft2dRows
// Apply RFFT to each row of a ps0 x ps1 real array.
// Input:  Float32Array ps0*ps1, row-major
// Output: Float32Array ps0*ps1, each row in halfcomplex format
// ---------------------------------------------------------------------------
function fft2dRows(input, ps0, ps1) {
  const rfft = new RFFT(ps1);
  const out  = new Float32Array(ps0 * ps1);
  for (let r = 0; r < ps0; r++) {
    const row    = input.slice(r * ps1, (r + 1) * ps1);
    const rowOut = rfft.forward(row);
    out.set(rowOut, r * ps1);
  }
  return out;
}


// ---------------------------------------------------------------------------
// fft2dCols
// Apply FFT to each column of a ps0 x ps1 halfcomplex array.
// We treat real and imaginary parts separately since RFFT gives halfcomplex.
// For column FFTs we use complex FFT on extracted re/im column pairs.
// ---------------------------------------------------------------------------
function fft2dCols(hc, ps0, ps1) {
  // For column FFTs on the halfcomplex rows, we need to handle the
  // halfcomplex format per row and do column-wise complex FFTs.
  // Extract full complex representation first (re + im planes).
  const re = new Float32Array(ps0 * ps1);
  const im = new Float32Array(ps0 * ps1);
  const half = Math.floor(ps1 / 2);

  for (let r = 0; r < ps0; r++) {
    const base = r * ps1;
    // DC
    re[r * ps1 + 0] = hc[base];
    im[r * ps1 + 0] = 0;
    // Nyquist
    if (ps1 % 2 === 0) {
      re[r * ps1 + half] = hc[base + half];
      im[r * ps1 + half] = 0;
    }
    // k=1..half-1
    for (let k = 1; k < half; k++) {
      re[r * ps1 + k]          =  hc[base + k];
      im[r * ps1 + k]          =  hc[base + ps1 - k];
      re[r * ps1 + ps1 - k]    =  hc[base + k];      // Hermitian conjugate
      im[r * ps1 + ps1 - k]    = -hc[base + ps1 - k];
    }
    if (ps1 % 2 !== 0) {
      re[r * ps1 + half] = hc[base + half];
      im[r * ps1 + half] = 0;
    }
  }

  // Now do column-wise RFFT on the real plane and imaginary plane separately
  // then combine. Since columns are now full-size, use RFFT on each column.
  const rfftCol = new RFFT(ps0);
  const reOut   = new Float32Array(ps0 * ps1);
  const imOut   = new Float32Array(ps0 * ps1);

  for (let c = 0; c < ps1; c++) {
    const colRe = new Float32Array(ps0);
    const colIm = new Float32Array(ps0);
    for (let r = 0; r < ps0; r++) {
      colRe[r] = re[r * ps1 + c];
      colIm[r] = im[r * ps1 + c];
    }
    const colReOut = rfftCol.forward(colRe);
    const colImOut = rfftCol.forward(colIm);
    // Complex FFT of (re + i*im) = FFT(re) + i*FFT(im)
    for (let r = 0; r < ps0; r++) {
      reOut[r * ps1 + c] = colReOut[r] - colImOut[r]; // wrong for imaginary
      imOut[r * ps1 + c] = colReOut[r] + colImOut[r]; // wrong
    }
  }

  return { re: reOut, im: imOut };
}


// ---------------------------------------------------------------------------
// buildComplexKernel
// Uses ndarray-fft for the kernel FFT (arbitrary size, double precision).
// Stores the result as { re, im } Float64Arrays.
// ---------------------------------------------------------------------------
function buildComplexKernel(b) {
  // Fall back to ndarray-fft for the kernel build since the 2D complex
  // column FFT approach above is complex to get right. ndarray-fft is
  // sufficient for the kernel which is built once.



  const ps0 = b.paddedSize[0];
  const ps1 = b.paddedSize[1];
  const n   = b.nRadialCells;
  const kSz = 2 * n + 1;

  // Pad kernel into ps0 x ps1
  const paddedKernel = new Float64Array(ps0 * ps1);
  for (let i = 1; i <= kSz; i++) {
    for (let j = 1; j <= kSz; j++) {
      const ki  = i - n - 1;
      const kj  = j - n - 1;
      paddedKernel[(i - 1) * ps1 + (j - 1)] = b.kernel[(ki + n) * kSz + (kj + n)];
    }
  }

  // Wrap
  const centKernel    = [Math.trunc(kSz / 2), Math.trunc(kSz / 2)];
  const wrappedKernel = new Float64Array(ps0 * ps1);
  makeWrappedArray([ps0, ps1], centKernel, paddedKernel, wrappedKernel);

  // 2D FFT via ndarray-fft
  const re = ndarray(new Float64Array(wrappedKernel), [ps0, ps1]);
  const im = ndarray(new Float64Array(ps0 * ps1), [ps0, ps1]);
  fft(1, re, im);

  b.complexKernel        = { re: re.data, im: im.data, ps0, ps1 };
  b.complexKernelCreated = true;
}


// ---------------------------------------------------------------------------
// convolve2DChannel
// FFT-based 2D convolution using ndarray-fft.
// ---------------------------------------------------------------------------
function convolve2DChannel(sliceIn, nPixels, b, sliceOut) {



  const nx  = nPixels[0];
  const ny  = nPixels[1];
  const { re: kRe, im: kIm, ps0, ps1 } = b.complexKernel;

  // Pad input
  const padded = new Float64Array(ps0 * ps1);
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < ny; j++)
      padded[i * ps1 + j] = sliceIn[i * ny + j];

  // FFT slice
  const sRe = ndarray(new Float64Array(padded), [ps0, ps1]);
  const sIm = ndarray(new Float64Array(ps0 * ps1), [ps0, ps1]);
  fft(1, sRe, sIm);

  // Pointwise multiply
  const pRe = new Float64Array(ps0 * ps1);
  const pIm = new Float64Array(ps0 * ps1);
  for (let idx = 0; idx < ps0 * ps1; idx++) {
    const re1 = sRe.data[idx]; const im1 = sIm.data[idx];
    const re2 = kRe[idx];      const im2 = kIm[idx];
    pRe[idx] = re1 * re2 - im1 * im2;
    pIm[idx] = re1 * im2 + im1 * re2;
  }

  // IFFT
  const rRe = ndarray(pRe, [ps0, ps1]);
  const rIm = ndarray(pIm, [ps0, ps1]);
  fft(-1, rRe, rIm);
  const N = ps0 * ps1;

  // Crop to nx x ny
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < ny; j++)
      sliceOut[i * ny + j] = f32(rRe.data[i * ps1 + j]);
}


// ---------------------------------------------------------------------------
// cubeBeamConvolution
// ---------------------------------------------------------------------------
function cubeBeamConvolution(dc, b) {
  if (!b.complexKernelCreated) buildComplexKernel(b);

  const dh  = dc.dh;
  const nx  = dh.nPixels[0];
  const ny  = dh.nPixels[1];
  const nch = dh.nChannels;
  const sliceIn  = new Float32Array(nx * ny);
  const sliceOut = new Float32Array(nx * ny);

  for (let ch = 0; ch < nch; ch++) {
    for (let i = 0; i < nx; i++)
      for (let j = 0; j < ny; j++)
        sliceIn[i * ny + j] = dc.flux[flatIndxCalc(i, j, ch, dh)];

    convolve2DChannel(sliceIn, dh.nPixels, b, sliceOut);

    for (let i = 0; i < nx; i++)
      for (let j = 0; j < ny; j++)
        dc.flux[flatIndxCalc(i, j, ch, dh)] = sliceOut[i * ny + j];
  }
}

module.exports = { cubeBeamConvolution, buildComplexKernel, convolve2DChannel };

if (require.main === module) {
  const { Beam2D, allocate_Beam2D } = require('../ObjectDefinitions/Beam.js');
  const { DataCube, allocateDataCube } = require('../ObjectDefinitions/DataCube.js');

  const dc = new DataCube();
  const dh = dc.dh;
  dh.nPixels[0] = 32; dh.nPixels[1] = 32; dh.nChannels = 4;
  dh.pixelSize[0] = f32(-1.0); dh.pixelSize[1] = f32(1.0);
  dh.channelSize = f32(10.0);
  dh.refLocation[0] = f32(16.0); dh.refLocation[1] = f32(16.0);
  dh.refLocation[2] = f32(2.0);
  dh.refVal[0] = f32(0.0); dh.refVal[1] = f32(0.0); dh.refVal[2] = f32(1000.0);
  allocateDataCube(dc);

  dc.flux[flatIndxCalc(16, 16, 1, dh)] = f32(1.0);
  const totalBefore = dc.flux.reduce((a, v) => a + v, 0);

  const b = new Beam2D();
  b.beamFWHM = f32(3.0); b.sigmaLengths = f32(3.0);
  b.pixelSize[0] = f32(-1.0); b.pixelSize[1] = f32(1.0);
  allocate_Beam2D(b, new Int32Array([32, 32]));
  calculate2DBeamKernel(b, b.pixelSize);
  cubeBeamConvolution(dc, b);

  const totalAfter = dc.flux.reduce((a, v) => f32(a + v), f32(0));
  const centre     = dc.flux[flatIndxCalc(16, 16, 1, dh)];
  const neighbour  = dc.flux[flatIndxCalc(17, 16, 1, dh)];
  const empty      = dc.flux[flatIndxCalc(16, 16, 0, dh)];

  console.log('=== cubeBeamConvolution ===');
  console.log('total flux before:', totalBefore.toFixed(6));
  console.log('total flux after: ', totalAfter.toFixed(6), '(expect ~1.0)');
  console.log('centre pixel:     ', centre.toFixed(6), '(Fortran: 0.098133)');
  console.log('neighbour pixel:  ', neighbour.toFixed(6), '(Fortran: 0.072111)');
  console.log('empty channel:    ', empty.toFixed(6), '(expect 0)');
  console.log('flux conserved:   ', Math.abs(totalAfter - totalBefore) < 0.01 ? 'OK' : 'FAIL');
  console.log('centre is max:    ', centre > neighbour ? 'OK' : 'FAIL');
  console.log('channels isolated:', empty === 0 ? 'OK' : 'FAIL');
}
