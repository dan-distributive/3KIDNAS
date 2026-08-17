'use strict';

// =============================================================================
// MaskCube.js
// Port of src/CompareCubes/MaskCube.f (MaskCubeMod) -- MaskCube only.
// MakeMaskCube (build a 0/1 mask from a flux tolerance) isn't part of the
// InitialAnalysis call chain (SoFiA already supplies the mask) and isn't
// ported here.
// =============================================================================

// ---------------------------------------------------------------------------
// maskCube
// Fortran: MaskCube(Cube, Mask)
//
// Elementwise Cube.Flux *= Mask.Flux (mask assumed 1 where there's signal,
// 0 elsewhere). Mutates cube.flux in place. Both cubes must share the same
// dimensions/flatIndxCalc layout -- a pointwise product needs no particular
// iteration order, so this just walks the flat array directly.
// ---------------------------------------------------------------------------
function maskCube(cube, mask) {
  const f32 = Math.fround;
  const flux = cube.flux;
  const maskFlux = mask.flux;
  for (let l = 0; l < flux.length; l++) {
    flux[l] = f32(flux[l] * maskFlux[l]);
  }
}

module.exports = { maskCube };


// ---------------------------------------------------------------------------
// Self-test (node MaskCube.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { DataCube, allocateDataCube } = require('../ObjectDefinitions/DataCube.js');

  console.log('=== maskCube ===');
  const cube = new DataCube();
  cube.dh.nPixels[0] = 2; cube.dh.nPixels[1] = 2; cube.dh.nChannels = 2;
  allocateDataCube(cube);
  cube.flux.fill(5.0);

  const mask = new DataCube();
  mask.dh.nPixels[0] = 2; mask.dh.nPixels[1] = 2; mask.dh.nChannels = 2;
  allocateDataCube(mask);
  mask.flux.set([1, 0, 1, 0, 0, 1, 0, 1]);

  maskCube(cube, mask);
  console.log('  masked flux:', Array.from(cube.flux));
  const expected = [5, 0, 5, 0, 0, 5, 0, 5];
  const ok = cube.flux.every((v, i) => v === expected[i]);
  console.log('  ', ok ? 'OK' : 'FAIL');
}
