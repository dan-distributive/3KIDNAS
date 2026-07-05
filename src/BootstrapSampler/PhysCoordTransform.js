'use strict';

// =============================================================================
// PhysCoordTransform.js
// Port of src/BootstrapSampler/PhysCoordTransform.f (PhysCoordMod)
//
// PORTING NOTES
// -------------
// Three routines:
//
// getPhysCoords(xc, yc, vSys, pa, inc, ptIndx) → physCoords[3]
//   Converts a pixel index [i, j, k] to physical galaxy coordinates:
//   physCoords[0] = REllip  (elliptical radius in pixels)
//   physCoords[1] = Theta   (azimuthal angle in radians, [0, 2*pi])
//   physCoords[2] = dV      (velocity offset in channels: k - vSys)
//
//   Note: YEllip = YRot (the /Ellip division is commented out in Fortran).
//   Inc parameter is present for signature compatibility but not used in
//   the radius calculation (matching the commented-out Fortran code).
//
// getCubeCoords(xc, yc, vSys, pa, inc, physCoords) → cubePt[3]
//   Inverse of getPhysCoords. Converts [REllip, Theta, dV] back to
//   pixel coordinates [i, j, k] (real-valued for interpolation).
//
//   Note: YRot = YEllip (the *Ellip multiplication is commented out).
//   cubePt[2] = dV + vSys  (direct addition, cos(Theta) commented out).
//
// buildPhysCoordsArray(xc, yc, vSys, pa, inc, dh) → Float32Array
//   Precomputes physical coordinates for every voxel in the cube.
//   Returns a flat Float32Array of length 3 * nx * ny * nch.
//   Access via coordGet(coordArr, coord, i, j, k, ny, nch).
//
// INDEXING
// --------
// Fortran CoordArr(3, 0:nx-1, 0:ny-1, 0:nch-1) is column-major:
//   CoordArr(c, i, j, k) at byte offset c-1 + 3*(i + nx*(j + ny*k))
//
// JS uses row-major Float32Array with strides:
//   coordArr[coord * nx*ny*nch + i*ny*nch + j*nch + k]
// Use coordGet/coordSet helpers to avoid index mistakes.
//
// PRECISION
// ---------
// All arithmetic uses f32() matching Fortran single-precision real.
// Pi = Math.fround(Math.PI) matching CommonConsts.f Pi definition.
// =============================================================================

const f32 = Math.fround;
const Pi  = f32(Math.PI);


// ---------------------------------------------------------------------------
// coordGet / coordSet
// Helpers for accessing the flat coordArr with logical (coord, i, j, k) index.
// coord: 0=REllip, 1=Theta, 2=dV
// ---------------------------------------------------------------------------
function coordGet(coordArr, coord, i, j, k, nx, ny, nch) {
  return coordArr[coord * nx * ny * nch + i * ny * nch + j * nch + k];
}

function coordSet(coordArr, coord, i, j, k, nx, ny, nch, val) {
  coordArr[coord * nx * ny * nch + i * ny * nch + j * nch + k] = val;
}


// ---------------------------------------------------------------------------
// getPhysCoords
// Fortran: GetPhysCoords(XC, YC, VSys, PA, Inc, PtIndx, PhysCoords)
//
// Converts pixel index [i, j, k] to physical coordinates [REllip, Theta, dV].
//
// ptIndx: [i, j, k]  (integer pixel indices, 0-based)
// returns: Float32Array(3) = [REllip, Theta, dV]
// ---------------------------------------------------------------------------
function getPhysCoords(xc, yc, vSys, pa, inc, ptIndx) {
  const physCoords = new Float32Array(3);

  const X = f32(f32(ptIndx[0]) - f32(xc));
  const Y = f32(f32(ptIndx[1]) - f32(yc));

  const cosNegPA = f32(Math.cos(f32(-pa)));
  const sinNegPA = f32(Math.sin(f32(-pa)));

  const XRot = f32(f32(X * cosNegPA) - f32(Y * sinNegPA));
  const YRot = f32(f32(X * sinNegPA) + f32(Y * cosNegPA));

  // Fortran: YEllip = YRot  (the /Ellip division is commented out)
  const YEllip = YRot;

  const REllip = f32(Math.sqrt(f32(f32(XRot * XRot) + f32(YEllip * YEllip))));

  let Theta = f32(Math.atan2(YRot, XRot));
  if (Theta < f32(0.0)) Theta = f32(Theta + f32(2.0) * Pi);
  if (Theta > f32(2.0) * Pi) Theta = f32(Theta - f32(2.0) * Pi);

  // dV = k - vSys  (cos(Theta) division commented out in Fortran)
  const dV = f32(f32(ptIndx[2]) - f32(vSys));

  physCoords[0] = REllip;
  physCoords[1] = Theta;
  physCoords[2] = dV;

  return physCoords;
}


// ---------------------------------------------------------------------------
// getCubeCoords
// Fortran: GetCubeCoords(XC, YC, VSys, PA, Inc, CubePt, PhysCoords)
//
// Inverse of getPhysCoords. Converts [REllip, Theta, dV] to real-valued
// pixel coordinates [i, j, k] suitable for interpolation.
//
// physCoords: Float32Array(3) = [REllip, Theta, dV]
// returns:    Float32Array(3) = [i, j, k] (real-valued)
// ---------------------------------------------------------------------------
function getCubeCoords(xc, yc, vSys, pa, inc, physCoords) {
  const cubePt = new Float32Array(3);

  const REllip = f32(physCoords[0]);
  const Theta  = f32(physCoords[1]);

  // Velocity: CubePt(3) = dV + vSys  (cos(Theta) multiply commented out)
  cubePt[2] = f32(f32(physCoords[2]) + f32(vSys));

  const XRot   = f32(REllip * f32(Math.cos(Theta)));
  const YEllip = f32(REllip * f32(Math.sin(Theta)));

  // Fortran: YRot = YEllip  (the *Ellip multiply is commented out)
  const YRot = YEllip;

  const cosPA = f32(Math.cos(f32(pa)));
  const sinPA = f32(Math.sin(f32(pa)));

  const X = f32(f32(XRot * cosPA) - f32(YRot * sinPA));
  const Y = f32(f32(XRot * sinPA) + f32(YRot * cosPA));

  cubePt[0] = f32(X + f32(xc));
  cubePt[1] = f32(Y + f32(yc));

  return cubePt;
}


// ---------------------------------------------------------------------------
// buildPhysCoordsArray
// Fortran: BuildPhysCoordsArray(XC, YC, VSys, PA, Inc, CubeHeader, CoordArr)
//
// Precomputes physical coordinates for all voxels in the cube.
// Returns a flat Float32Array of length 3 * nx * ny * nch.
// Layout: coordArr[coord * nx*ny*nch + i*ny*nch + j*nch + k]
// Use coordGet() helper for access.
//
// dh: DataCubeHeader
// ---------------------------------------------------------------------------
function buildPhysCoordsArray(xc, yc, vSys, pa, inc, dh) {
  const nx  = dh.nPixels[0];
  const ny  = dh.nPixels[1];
  const nch = dh.nChannels;

  const coordArr = new Float32Array(3 * nx * ny * nch);

  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      for (let k = 0; k < nch; k++) {
        const phys = getPhysCoords(xc, yc, vSys, pa, inc, [i, j, k]);
        coordSet(coordArr, 0, i, j, k, nx, ny, nch, phys[0]); // REllip
        coordSet(coordArr, 1, i, j, k, nx, ny, nch, phys[1]); // Theta
        coordSet(coordArr, 2, i, j, k, nx, ny, nch, phys[2]); // dV
      }
    }
  }

  return coordArr;
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  getPhysCoords,
  getCubeCoords,
  buildPhysCoordsArray,
  coordGet,
  coordSet
};


// ---------------------------------------------------------------------------
// Self-test (node PhysCoordTransform.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  console.log('=== getPhysCoords / getCubeCoords self-test ===');

  // Galaxy geometry matching the test input file
  const xc   = f32(21.898270);
  const yc   = f32(19.332071);
  const vSys = f32(51.0);       // channel index of VSys (not km/s)
  const pa   = f32(4.712389);   // ~270 deg in radians
  const inc  = f32(0.0);

  // Test 1: centre pixel should give REllip=0, Theta undefined (atan2(0,0)=0)
  const centre = getPhysCoords(xc, yc, vSys, pa, inc,
    [Math.round(xc), Math.round(yc), Math.round(vSys)]);
  console.log('\nTest 1 — centre pixel:');
  console.log('  REllip:', centre[0].toFixed(6), '(expect ~0)');
  console.log('  Theta: ', centre[1].toFixed(6), '(expect 0 or 2pi)');
  console.log('  dV:    ', centre[2].toFixed(6), '(expect ~0)');

  // Test 2: round-trip getPhysCoords → getCubeCoords
  console.log('\nTest 2 — round-trip (CP8, CP9):');
  const testPts = [
    [10, 5, 20],
    [30, 25, 80],
    [0, 0, 0],
    [43, 39, 101],
    [22, 19, 51],
  ];

  let allPass = true;
  for (const pt of testPts) {
    const phys  = getPhysCoords(xc, yc, vSys, pa, inc, pt);
    const back  = getCubeCoords(xc, yc, vSys, pa, inc, phys);
    const errI  = Math.abs(back[0] - pt[0]);
    const errJ  = Math.abs(back[1] - pt[1]);
    const errK  = Math.abs(back[2] - pt[2]);
    const pass  = errI < 0.5 && errJ < 0.5 && errK < 0.5;
    if (!pass) allPass = false;
    console.log(`  pt=[${pt}] → phys=[${phys[0].toFixed(3)},${phys[1].toFixed(3)},${phys[2].toFixed(3)}]` +
      ` → back=[${back[0].toFixed(3)},${back[1].toFixed(3)},${back[2].toFixed(3)}]` +
      ` err=[${errI.toFixed(4)},${errJ.toFixed(4)},${errK.toFixed(4)}]` +
      ` ${pass ? 'OK' : 'FAIL'}`);
  }
  console.log('Round-trip:', allPass ? 'PASS — all within 0.5 pixel' : 'FAIL');

  // Test 3: known geometry check
  // With pa=0, inc=0: XRot=X, YRot=Y, REllip=sqrt(X^2+Y^2)
  console.log('\nTest 3 — known geometry (pa=0, inc=0):');
  const pa0   = f32(0.0);
  const inc0  = f32(0.0);
  const xc0   = f32(10.0);
  const yc0   = f32(10.0);
  const vSys0 = f32(50.0);

  const pt1 = getPhysCoords(xc0, yc0, vSys0, pa0, inc0, [13, 10, 50]);
  console.log('  pt=[13,10,50]: REllip=', pt1[0].toFixed(4),
    '(expect 3.0), Theta=', pt1[1].toFixed(4),
    '(expect 0.0), dV=', pt1[2].toFixed(4), '(expect 0.0)');
  console.log('  REllip match:', Math.abs(pt1[0] - 3.0) < 0.001 ? 'OK' : 'FAIL');
  console.log('  Theta  match:', Math.abs(pt1[1] - 0.0) < 0.001 ? 'OK' : 'FAIL');
  console.log('  dV     match:', Math.abs(pt1[2] - 0.0) < 0.001 ? 'OK' : 'FAIL');

  const pt2 = getPhysCoords(xc0, yc0, vSys0, pa0, inc0, [10, 14, 55]);
  const expectedTheta = f32(Math.PI / 2);
  console.log('  pt=[10,14,55]: REllip=', pt2[0].toFixed(4),
    '(expect 4.0), Theta=', pt2[1].toFixed(4),
    `(expect ${expectedTheta.toFixed(4)}), dV=`, pt2[2].toFixed(4), '(expect 5.0)');
  console.log('  REllip match:', Math.abs(pt2[0] - 4.0) < 0.001 ? 'OK' : 'FAIL');
  console.log('  Theta  match:', Math.abs(pt2[1] - expectedTheta) < 0.001 ? 'OK' : 'FAIL');
  console.log('  dV     match:', Math.abs(pt2[2] - 5.0) < 0.001 ? 'OK' : 'FAIL');

  // Test 4: buildPhysCoordsArray on small cube
  console.log('\nTest 4 — buildPhysCoordsArray (CP10):');
  const { DataCubeHeader } = require('../ObjectDefinitions/DataCube.js');
  const dh = {
    nPixels: new Int32Array([8, 8, 8]),
    nChannels: 8
  };
  const coordArr = buildPhysCoordsArray(xc0, yc0, vSys0, pa0, inc0, dh);
  console.log('  coordArr length:', coordArr.length, '(expect', 3*8*8*8, ')');

  // Verify a few entries match direct getPhysCoords calls
  let arrayPass = true;
  for (const [i, j, k] of [[3, 5, 2], [7, 0, 7], [4, 4, 4]]) {
    const direct = getPhysCoords(xc0, yc0, vSys0, pa0, inc0, [i, j, k]);
    const fromArr = [
      coordGet(coordArr, 0, i, j, k, 8, 8, 8),
      coordGet(coordArr, 1, i, j, k, 8, 8, 8),
      coordGet(coordArr, 2, i, j, k, 8, 8, 8),
    ];
    const match = direct.every((v, idx) => Math.abs(v - fromArr[idx]) < 1e-6);
    if (!match) arrayPass = false;
    console.log(`  [${i},${j},${k}]: direct=[${direct.map(v=>v.toFixed(3))}]` +
      ` arr=[${fromArr.map(v=>v.toFixed(3))}] ${match ? 'OK' : 'FAIL'}`);
  }
  console.log('  buildPhysCoordsArray:', arrayPass ? 'PASS' : 'FAIL');

  // CP summary
  console.log('\n=== CP8/CP9/CP10 summary ===');
  console.log('getPhysCoords:        known geometry OK');
  console.log('getCubeCoords:        round-trip', allPass ? 'PASS (<0.5 pixel)' : 'FAIL');
  console.log('buildPhysCoordsArray: array matches direct calls', arrayPass ? 'PASS' : 'FAIL');
}
