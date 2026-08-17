'use strict';

// =============================================================================
// Interpolation.js
// High-fidelity port of src/StandardMath/Interpolation.f (InterpolateMod)
//
// PORTING NOTES
// -------------
// Fortran `real` is 32-bit single precision. All arithmetic uses Math.fround()
// to match Fortran rounding at each step.
//
// Fortran arrays are 1-indexed. All array arguments here use 0-indexed JS
// arrays with the same layout — callers must adjust accordingly.
//
// Point convention (matches Fortran):
//   P1, P2 = [x, y]  (2-element arrays, 0-indexed)
//   Corners (bilinear) = 4x3 array: corners[i] = [x, y, val]
//   Corners (trilinear) = 8x4 array: corners[i] = [x, y, z, val]
//   PTarg (bilinear) = [x, y, val]  (val is output, index 2)
//   PTarg (trilinear) = [x, y, z, val]  (val is output, index 3)
// =============================================================================

const f32 = Math.fround;


// ---------------------------------------------------------------------------
// simpleInterpolateX
// Linear interpolation between two points to find the X value at a target Y.
// Fortran: SimpleInterpolateX(P1, P2, YTarg, XIntVal)
//
// P1 = [x1, y1], P2 = [x2, y2]
// slope = (y2 - y1) / (x2 - x1)
// XIntVal = (YTarg - y1) / slope + x1
//
// Returns XIntVal as f32.
// ---------------------------------------------------------------------------
function simpleInterpolateX(p1, p2, yTarg) {
  const slope = f32(f32(p2[1] - p1[1]) / f32(p2[0] - p1[0]));
  const delY  = f32(f32(yTarg) - f32(p1[1]));
  return f32(f32(delY / slope) + f32(p1[0]));
}


// ---------------------------------------------------------------------------
// simpleInterpolateY
// Linear interpolation between two points to find the Y value at a target X.
// Fortran: SimpleInterpolateY(P1, P2, XTarg, YIntVal)
//
// P1 = [x1, y1], P2 = [x2, y2]
// slope = (y2 - y1) / (x2 - x1)
// YIntVal = (XTarg - x1) * slope + y1
//
// Returns YIntVal as f32.
// ---------------------------------------------------------------------------
function simpleInterpolateY(p1, p2, xTarg) {
  const slope = f32(f32(p2[1] - p1[1]) / f32(p2[0] - p1[0]));
  const delX  = f32(f32(xTarg) - f32(p1[0]));
  return f32(f32(delX * slope) + f32(p1[1]));
}


// ---------------------------------------------------------------------------
// biLinearInterpolation
// Bilinear interpolation on a square grid with 4 corners.
//
// Fortran layout (1-indexed):
//   Corners(4,3): corners[i][j] where i=corner(1-4), j=x/y/val(1-3)
//   Corner ordering:
//     C1 = corners[0] (bottom-left)   C2 = corners[1] (bottom-right)
//     C3 = corners[2] (top-left)      C4 = corners[3] (top-right)
//
// JS layout (0-indexed):
//   corners[i] = [x, y, val]   i = 0..3
//   pTarg = [x, y, val]        pTarg[2] is output
//
// Algorithm:
//   1. Interpolate along bottom edge (C1-C2) at target X → point A
//   2. Interpolate along top edge (C3-C4) at target X → point B
//   3. Interpolate between A and B at target Y → output val
//
// Modifies pTarg[2] in place and returns it.
// ---------------------------------------------------------------------------
function biLinearInterpolation(pTarg, corners) {
  // CInterpolate[i] = [x, y, val] for the two edge-interpolated points
  // Fortran: CInterpolate(1:2, 1) = PTarg(1)  → x is fixed to pTarg[0]
  const cInterpolate = [
    [f32(pTarg[0]), 0.0, 0.0],
    [f32(pTarg[0]), 0.0, 0.0]
  ];

  const cTempX = f32(pTarg[0]);

  // Loop over the two edges (bottom: i=0, top: i=1)
  // Fortran: do i=1,2  →  j=(i-1)*2+1  gives j=1,3 (corner pair starts)
  for (let i = 0; i < 2; i++) {
    const j = i * 2;   // corner pair start index: 0 (bottom) or 2 (top)

    // P1Temp = [corners[j].x, corners[j].val]
    // P2Temp = [corners[j+1].x, corners[j+1].val]
    // Fortran: kk=(k-1)*2+1 gives kk=1,3 → indices x(0) and val(2)
    const p1Temp = [f32(corners[j][0]),     f32(corners[j][2])];
    const p2Temp = [f32(corners[j + 1][0]), f32(corners[j + 1][2])];

    // Interpolate along the edge to target X
    const interpVal = simpleInterpolateY(p1Temp, p2Temp, cTempX);

    // Store: y coordinate of this edge, interpolated val
    // Fortran: CInterpolate(i,2)=Corners(j,2)  → edge y coordinate
    //          CInterpolate(i,3)=CTemp(2)       → interpolated val
    cInterpolate[i][1] = f32(corners[j][1]);   // y of this edge
    cInterpolate[i][2] = interpVal;
  }

  // Interpolate between the two edge points at target Y
  // P1Temp = [cInterpolate[0].y, cInterpolate[0].val]
  // P2Temp = [cInterpolate[1].y, cInterpolate[1].val]
  const p1Final = [cInterpolate[0][1], cInterpolate[0][2]];
  const p2Final = [cInterpolate[1][1], cInterpolate[1][2]];
  const yTarg   = f32(pTarg[1]);

  pTarg[2] = simpleInterpolateY(p1Final, p2Final, yTarg);
  return pTarg;
}


// ---------------------------------------------------------------------------
// triLinearInterpolation
// Trilinear interpolation on a cubic grid with 8 corners.
//
// Fortran layout (1-indexed):
//   Corners(8,4): corners[i][j] where i=corner(1-8), j=x/y/z/val(1-4)
//   Corner ordering (bottom face first, then top face):
//     C1=[0], C2=[1], C3=[2], C4=[3]  (bottom z face)
//     C5=[4], C6=[5], C7=[6], C8=[7]  (top z face)
//
// JS layout (0-indexed):
//   corners[i] = [x, y, z, val]   i = 0..7
//   pTarg = [x, y, z, val]        pTarg[3] is output
//
// Algorithm:
//   1. Bilinear interpolation on bottom face (corners 0-3) at (x,y) → valBottom
//   2. Bilinear interpolation on top face (corners 4-7) at (x,y) → valTop
//   3. Linear interpolation between valBottom and valTop at target Z → output
//
// Modifies pTarg[3] in place and returns it.
// ---------------------------------------------------------------------------
function triLinearInterpolation(pTarg, corners) {
  // surfacePoints[i] = [x, y, val] for the two face-interpolated points
  // Fortran: SurfacePoints(1:2,1)=PTarg(1), SurfacePoints(1:2,2)=PTarg(2)
  const surfacePoints = [
    [f32(pTarg[0]), f32(pTarg[1]), 0.0],
    [f32(pTarg[0]), f32(pTarg[1]), 0.0]
  ];

  // Loop through the two Z faces
  // Fortran: jj tracks global corner index (0-based here)
  let jj = 0;
  for (let i = 0; i < 2; i++) {
    // Build 4-corner temp array for this face
    // Fortran: CornersTemp(ii,1:2)=Corners(jj,1:2), CornersTemp(ii,3)=Corners(jj,4)
    const cornersTemp = [];
    for (let j = 0; j < 2; j++) {
      for (let k = 0; k < 2; k++) {
        cornersTemp.push([
          f32(corners[jj][0]),   // x
          f32(corners[jj][1]),   // y
          f32(corners[jj][3])    // val (index 3, not z)
        ]);
        jj++;
      }
    }
    // Bilinear interpolation on this face
    biLinearInterpolation(surfacePoints[i], cornersTemp);
  }

  // Linear interpolation between the two face results at target Z
  // Fortran: PTemp(i,1)=Corners((i-1)*4+1, 3)  → z coordinate of face
  //          PTemp(i,2)=SurfacePoints(i,3)       → interpolated val
  const p1Final = [f32(corners[0][2]), surfacePoints[0][2]];   // [z_bottom, val_bottom]
  const p2Final = [f32(corners[4][2]), surfacePoints[1][2]];   // [z_top,    val_top]
  const zTarg   = f32(pTarg[2]);

  pTarg[3] = simpleInterpolateY(p1Final, p2Final, zTarg);
  return pTarg;
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  simpleInterpolateX,
  simpleInterpolateY,
  biLinearInterpolation,
  triLinearInterpolation
};


// ---------------------------------------------------------------------------
// Self-test (node Interpolation.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  console.log('=== simpleInterpolateX ===');
  // P1=(0,0), P2=(1,1), YTarg=0.5 → XIntVal=0.5
  console.log(simpleInterpolateX([0,0],[1,1], 0.5).toFixed(10));

  console.log('\n=== simpleInterpolateY ===');
  // P1=(0,0), P2=(2,4), XTarg=1 → slope=2, YIntVal=2
  console.log(simpleInterpolateY([0,0],[2,4], 1.0).toFixed(10));

  console.log('\n=== biLinearInterpolation ===');
  // Unit square corners, value = x+y
  // C1=(0,0,0), C2=(1,0,1), C3=(0,1,1), C4=(1,1,2)
  // Target (0.5, 0.5) → expected val = 1.0
  const corners4 = [
    [0,0,0], [1,0,1],
    [0,1,1], [1,1,2]
  ];
  const pTarg2 = [0.5, 0.5, 0.0];
  biLinearInterpolation(pTarg2, corners4);
  console.log(pTarg2[2].toFixed(10));

  console.log('\n=== triLinearInterpolation ===');
  // Unit cube corners, value = x+y+z
  const corners8 = [
    [0,0,0,0],[1,0,0,1],[0,1,0,1],[1,1,0,2],
    [0,0,1,1],[1,0,1,2],[0,1,1,2],[1,1,1,3]
  ];
  const pTarg3 = [0.5, 0.5, 0.5, 0.0];
  triLinearInterpolation(pTarg3, corners8);
  console.log(pTarg3[3].toFixed(10));
}
