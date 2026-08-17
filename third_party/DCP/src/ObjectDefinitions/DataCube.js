'use strict';

// =============================================================================
// DataCube.js
// High-fidelity port of src/ObjectDefinitions/DataCube.f (DataCubeMod)
//
// PORTING NOTES
// -------------
// Fortran `real` is 32-bit single precision — Math.fround() used throughout.
// Fortran integer arithmetic is exact — standard JS integer ops used.
//
// Fortran TYPE → JS class:
//   DataCubeHeader → class DataCubeHeader
//   DataCube       → class DataCube
//
// Array indexing: all arrays are 0-indexed matching Fortran (0:n-1).
//
// Float32Array is used for all real-valued arrays to match Fortran `real`.
// Int32Array is used for integer arrays.
//
// AllocateDataCube initializes FlattendValidIndices to the full cube
// (all cells valid) matching the Fortran comment: "Always initialize the
// number of valid indices to the full cube in case not used elsewhere."
// The caller is responsible for populating FlattendValidIndices and
// setting nValid when masking is applied.
// =============================================================================

const f32 = Math.fround;


// ---------------------------------------------------------------------------
// DataCubeHeader
// Fortran TYPE DataCubeHeader — metadata for a 3D data cube.
// ---------------------------------------------------------------------------
class DataCubeHeader {
  constructor() {
    this.nPixels         = new Int32Array(2);     // (0:1) spatial dimensions
    this.nChannels       = 0;                      // velocity channels
    this.pixelSize       = new Float32Array(2);   // (0:1) arcsec/pixel
    this.channelSize     = f32(0.0);              // km/s per channel
    this.pixelCent       = new Float32Array(2);   // (0:1) centre pixel
    this.channelCent     = f32(0.0);
    this.start           = new Float32Array(3);   // (0:2) cube origin x,y,vel
    this.axisType        = ['', '', ''];           // (0:2) character(8)
    this.units           = ['', '', ''];           // (0:2) character(8)
    this.dimensionUnitSwitch = new Int32Array(4); // (0:3)
    this.uncertaintyUnitSwitch = 0;
    this.uncertainty     = f32(0.0);
    this.fUnit           = '';                     // character(10)
    this.fType           = '';                     // character(10)
    this.epoch           = f32(0.0);
    this.pixelCenterIndx = new Int32Array(2);     // (0:1)
    this.channelCenterIndx = 0;
    this.refLocation     = new Float32Array(3);   // (0:2)
    this.refVal          = new Float32Array(3);   // (0:2)
    this.startFreq       = f32(0.0);
    this.restFreq        = f32(0.0);
    this.delFreq         = f32(0.0);
    this.nValid          = 0;                      // non-null cell count
    this.maskSwitch      = false;
    this.sn_Peak         = f32(0.0);
    this.sn_Avg          = f32(0.0);
    this.sn_Median       = f32(0.0);
    this.sn_Int          = f32(0.0);
  }
}


// ---------------------------------------------------------------------------
// DataCube
// Fortran TYPE DataCube — header + flux array + valid index list.
// ---------------------------------------------------------------------------
class DataCube {
  constructor() {
    this.dh                   = new DataCubeHeader();
    this.pixels               = null;   // Float32Array [2][mPix]
    this.channels             = null;   // Float32Array [nChannels]
    this.flux                 = null;   // Float32Array [nPixX][nPixY][nChan]
    this.flattendValidIndices = null;   // Int32Array [nCells]
  }
}


// ---------------------------------------------------------------------------
// allocateDataCube
// Fortran: AllocateDataCube(DC)
//
// Allocates and initializes all arrays. Sets pixel/channel coordinate arrays
// from RefLocation, RefVal, PixelSize, ChannelSize. Initializes flux to 0.
// Sets FlattendValidIndices to 0..nCells-1 (all valid).
//
// Modifies dc in place.
// ---------------------------------------------------------------------------
function allocateDataCube(dc) {
  const dh   = dc.dh;
  const mPix = Math.max(dh.nPixels[0], dh.nPixels[1]);

  // Allocate coordinate and flux arrays
  // Fortran: ALLOCATE(DC%Pixels(0:1, 0:mPix-1))
  dc.pixels   = new Float32Array(2 * mPix);          // [axis][pixel]
  dc.channels = new Float32Array(dh.nChannels);
  dc.flux     = new Float32Array(
    dh.nPixels[0] * dh.nPixels[1] * dh.nChannels
  );

  // All cells valid by default
  const nCells = dh.nPixels[0] * dh.nPixels[1] * dh.nChannels;
  dh.nValid    = nCells;
  dc.flattendValidIndices = new Int32Array(nCells);
  for (let l = 0; l < nCells; l++) dc.flattendValidIndices[l] = l;

  // Flux initialized to 0 (Float32Array default)

  // Spatial axes (j=0,1)
  for (let j = 0; j <= 1; j++) {
    const delta0    = f32(f32(0.0) - f32(dh.refLocation[j]));
    dh.start[j]     = f32(f32(dh.refVal[j]) + f32(delta0 * f32(dh.pixelSize[j])));

    for (let i = 0; i < mPix; i++) {
      const delta    = f32(f32(i) - f32(dh.refLocation[j]));
      dc.pixels[j * mPix + i] = f32(f32(dh.refVal[j]) + f32(delta * f32(dh.pixelSize[j])));
    }
  }

  // Velocity axis (j=2)
  const delta0v  = f32(f32(0.0) - f32(dh.refLocation[2]));
  dh.start[2]    = f32(f32(dh.refVal[2]) + f32(delta0v * f32(dh.channelSize)));

  for (let i = 0; i < dh.nChannels; i++) {
    const delta    = f32(f32(i) - f32(dh.refLocation[2]));
    dc.channels[i] = f32(f32(dh.refVal[2]) + f32(delta * f32(dh.channelSize)));
  }

  // Centre indices (integer division matching Fortran)
  dh.pixelCenterIndx[0]  = Math.trunc(dh.nPixels[0] / 2);
  dh.pixelCenterIndx[1]  = Math.trunc(dh.nPixels[1] / 2);
  dh.channelCenterIndx   = Math.trunc(dh.nChannels / 2);
}


// ---------------------------------------------------------------------------
// deAllocateDataCube
// Fortran: DeAllocateDataCube(DC)
// In JS, just null out the arrays (GC handles the rest).
// ---------------------------------------------------------------------------
function deAllocateDataCube(dc) {
  dc.pixels               = null;
  dc.channels             = null;
  dc.flux                 = null;
  dc.flattendValidIndices = null;
}


// ---------------------------------------------------------------------------
// flatIndxCalc
// Fortran: FlatIndxCalc(i, j, k, DH, l)
// Converts 3D (x, y, channel) indices to a flat array index.
// Loop order: x → y → channel (k is innermost).
//
// l = k + j*nChannels + i*nChannels*nPixels[1]
// ---------------------------------------------------------------------------
function flatIndxCalc(i, j, k, dh) {
  return k + j * dh.nChannels + i * dh.nChannels * dh.nPixels[1];
}


// ---------------------------------------------------------------------------
// threeDIndxCalc
// Fortran: ThreeDIndxCalc(l, DH, i, j, k)
// Converts a flat array index to 3D (x, y, channel) indices.
// Inverse of flatIndxCalc.
//
// i = floor(l / (nChannels * nPixels[1]))
// j = floor(remainder / nChannels)
// k = remainder
//
// Returns { i, j, k }.
// ---------------------------------------------------------------------------
function threeDIndxCalc(l, dh) {
  let ltemp = l;
  const i   = Math.trunc(ltemp / (dh.nChannels * dh.nPixels[1]));
  ltemp    -= i * dh.nChannels * dh.nPixels[1];
  const j   = Math.trunc(ltemp / dh.nChannels);
  ltemp    -= j * dh.nChannels;
  const k   = ltemp;
  return { i, j, k };
}


// ---------------------------------------------------------------------------
// Flux array accessors
// The Fortran flux array is DC%Flux(0:nPixX-1, 0:nPixY-1, 0:nChan-1).
// In JS it's a flat Float32Array — use these helpers for clean access.
// ---------------------------------------------------------------------------
function fluxGet(dc, i, j, k) {
  return dc.flux[flatIndxCalc(i, j, k, dc.dh)];
}

function fluxSet(dc, i, j, k, val) {
  dc.flux[flatIndxCalc(i, j, k, dc.dh)] = f32(val);
}

function fluxAdd(dc, i, j, k, val) {
  const idx     = flatIndxCalc(i, j, k, dc.dh);
  dc.flux[idx]  = f32(dc.flux[idx] + f32(val));
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  DataCubeHeader,
  DataCube,
  allocateDataCube,
  deAllocateDataCube,
  flatIndxCalc,
  threeDIndxCalc,
  fluxGet,
  fluxSet,
  fluxAdd
};


// ---------------------------------------------------------------------------
// Self-test (node DataCube.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { DataCube, DataCubeHeader, allocateDataCube,
          flatIndxCalc, threeDIndxCalc, fluxGet, fluxSet } = module.exports;

  // Build a small 4x4x8 cube
  const dc = new DataCube();
  const dh = dc.dh;
  dh.nPixels[0]    = 4;
  dh.nPixels[1]    = 4;
  dh.nChannels     = 8;
  dh.pixelSize[0]  = -1.0;
  dh.pixelSize[1]  =  1.0;
  dh.channelSize   =  10.0;
  dh.refLocation[0] = 2.0;
  dh.refLocation[1] = 2.0;
  dh.refLocation[2] = 4.0;
  dh.refVal[0]     = 100.0;
  dh.refVal[1]     = 200.0;
  dh.refVal[2]     = 1000.0;

  allocateDataCube(dc);

  console.log('=== allocateDataCube ===');
  console.log('nValid:', dh.nValid, '(expect 128)');
  console.log('start[0]:', dh.start[0].toFixed(4), '(expect 102.0)');
  console.log('start[1]:', dh.start[1].toFixed(4), '(expect 198.0)');
  console.log('start[2]:', dh.start[2].toFixed(4), '(expect 960.0)');
  console.log('channels[0]:', dc.channels[0].toFixed(4), '(expect 960.0)');
  console.log('channels[4]:', dc.channels[4].toFixed(4), '(expect 1000.0)');
  console.log('pixelCenterIndx:', dh.pixelCenterIndx[0], dh.pixelCenterIndx[1], '(expect 2 2)');

  console.log('\n=== flatIndxCalc / threeDIndxCalc roundtrip ===');
  for (const [i,j,k] of [[0,0,0],[3,3,7],[1,2,5],[2,0,3]]) {
    const l   = flatIndxCalc(i, j, k, dh);
    const ijk = threeDIndxCalc(l, dh);
    const ok  = ijk.i===i && ijk.j===j && ijk.k===k;
    console.log(`(${i},${j},${k}) -> l=${l} -> (${ijk.i},${ijk.j},${ijk.k}) ${ok?'OK':'FAIL'}`);
  }

  console.log('\n=== fluxGet / fluxSet ===');
  fluxSet(dc, 1, 2, 3, 42.5);
  console.log('fluxGet(1,2,3):', fluxGet(dc, 1, 2, 3).toFixed(4), '(expect 42.5)');
  console.log('fluxGet(0,0,0):', fluxGet(dc, 0, 0, 0).toFixed(4), '(expect 0.0)');
}
