'use strict';

// =============================================================================
// EstimateCubeNoise.js
// Port of src/PreAnalysis/EstimateCubeNoise.f (EstimateCubeNoiseMod)
//
// Estimates the cube's noise RMS from the mean flux^2 over the cube's 8
// corner sub-boxes: the outer frac=0.12 slab of each axis, both sides, all
// 2^3 combinations -- assumes the source is confined away from the edges,
// so these corners are pure noise.
//
// PRECISION
// ---------
// Fortran `real` -> Math.fround() at each arithmetic step, matching the
// -O0 build (src/makeflags: no FMA/reassociation) -- see random.js's header
// for the project-wide rationale. `Cube%Flux(i,j,k)**2.` (real ** real
// literal 2.0) constant-folds to a plain multiply under gfortran -O0 (same
// reasoning already confirmed for X**2. in PhysCoordTransform.f/.js), so
// this ports as v*v, not Math.pow.
//
// No RNG usage anywhere in this file (confirmed against the Fortran source).
// =============================================================================

const f32 = Math.fround;
const { flatIndxCalc } = require('../ObjectDefinitions/DataCube.js');

// `process` isn't guaranteed to exist in a real DCP worker sandbox --
// confirmed directly (a bare `process.env.TRACE_DEBUG` check threw
// "process is not defined" on a real worker, while working fine locally
// under real Node). This file is a proper require()'d module (bundled by
// webpack via job.requires()), not the directly-shipped work function
// itself, so a module-scope helper here is fine -- unlike
// bootstrap-realization-launcher.js's own runBootstrapRealization, which
// needs this declared inside itself instead (see that file's comment).
function isTraceDebug() {
  return typeof process !== 'undefined' && process.env && process.env.TRACE_DEBUG === '1';
}

// ---------------------------------------------------------------------------
// setLims
// Fortran: SetLims(Cube, frac, DimSwitch, SideSwitch, Lims)
//
// dimSize: size of the dimension (nPixels[0]/nPixels[1]/nChannels, an
//          integer -> converted to real once, matching Fortran's
//          UpperLim=Cube%DH%nPixels(0) etc. implicit int->real conversion).
// sideSwitch: 1 = low side  [0, trunc(frac*upperLim)]
//             2 = high side [trunc((1-frac)*(upperLim-1)), upperLim-1]
//
// Fortran assigns these two REAL expressions directly into INTEGER Lims(1)/
// Lims(2) -- truncating-towards-zero on assignment, same as int(); ported
// here as explicit Math.trunc() at the same point.
// ---------------------------------------------------------------------------
function setLims(dimSize, frac, sideSwitch) {
  const upperLim = f32(dimSize);
  if (sideSwitch === 1) {
    return [0, Math.trunc(f32(frac * upperLim))];
  }
  const oneMinusFrac    = f32(f32(1.0) - frac);
  const upperLimMinus1  = f32(upperLim - f32(1.0));
  const lo = Math.trunc(f32(oneMinusFrac * upperLimMinus1));
  const hi = Math.trunc(upperLim) - 1;
  return [lo, hi];
}

// ---------------------------------------------------------------------------
// getSubBox
// Fortran: GetSubBox(Cube, CountLims, count, sum)
//
// Accumulates sum(flux^2) and count over one corner sub-box (inclusive index
// ranges per axis). acc is mutated in place, matching Fortran's INTENT(INOUT)
// count/sum threaded across all 8 corner-box calls in estimateNoise below.
// ---------------------------------------------------------------------------
function getSubBox(flux, dh, lims, acc) {
  const [i0, i1] = lims[0];
  const [j0, j1] = lims[1];
  const [k0, k1] = lims[2];
  for (let i = i0; i <= i1; i++) {
    for (let j = j0; j <= j1; j++) {
      for (let k = k0; k <= k1; k++) {
        const v = flux[flatIndxCalc(i, j, k, dh)];
        acc.sum = f32(acc.sum + f32(v * v));
        acc.count += 1;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// estimateNoise
// Fortran: EstimateNoise(Cube, RMS)
//
// RMS = sqrt(sum(flux^2) / count) over all 8 corner sub-boxes combined (sum
// and count accumulate across all 8, not reset per box -- matching Fortran's
// single count/sum threaded through the whole triple loop).
// ---------------------------------------------------------------------------
function estimateNoise(cube) {
  const dh   = cube.dh;
  const frac = f32(0.12);
  const acc  = { sum: f32(0.0), count: 0 };

  for (let ii = 1; ii <= 2; ii++) {
    const limsX = setLims(dh.nPixels[0], frac, ii);
    for (let jj = 1; jj <= 2; jj++) {
      const limsY = setLims(dh.nPixels[1], frac, jj);
      for (let kk = 1; kk <= 2; kk++) {
        const limsK = setLims(dh.nChannels, frac, kk);
        getSubBox(cube.flux, dh, [limsX, limsY, limsK], acc);
      }
    }
  }

  if (isTraceDebug()) {
    console.error('TRACE Noise square count', acc.count, acc.sum, f32(Math.sqrt(f32(acc.sum / f32(acc.count)))));
  }

  // Fortran: RMS=sqrt(sum/count) -- count is INTEGER, promoted to real for
  // the mixed-mode division before the sqrt.
  return f32(Math.sqrt(f32(acc.sum / f32(acc.count))));
}

// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
// Bundle-hash cache-buster: dcp-client's Job._publishLocalModules names its
// package from a hash of the POST-WEBPACK-MINIFIED bundle, not raw source
// -- a comment-only edit here got stripped by minification and produced an
// identical hash/name, still colliding with a NAMECONFLICT from a
// differently-owned prior publish of the same bundle. This exported
// constant can't be minified away (it's a real binding), so it changes the
// bundle's actual content. Bump the value if a NAMECONFLICT recurs.
const BUNDLE_CACHE_BUSTER = 1;

module.exports = { estimateNoise, setLims, getSubBox, BUNDLE_CACHE_BUSTER };


// ---------------------------------------------------------------------------
// Self-test (node EstimateCubeNoise.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { DataCube, allocateDataCube } = require('../ObjectDefinitions/DataCube.js');

  console.log('=== setLims ===');
  // nPixels=44, frac=0.12: low=[0, trunc(0.12*44)]=[0,5]; high=[trunc(0.88*43),43]=[37,43]
  const low44  = setLims(44, f32(0.12), 1);
  const high44 = setLims(44, f32(0.12), 2);
  console.log('  dim=44 low :', low44,  '(expect [0, 5])');
  console.log('  dim=44 high:', high44, '(expect [37, 43])');

  console.log('\n=== estimateNoise on a constant cube ===');
  // RMS of a cube filled with a single constant C must equal |C| exactly,
  // regardless of which 8 corners get sampled.
  const cube = new DataCube();
  cube.dh.nPixels[0] = 44; cube.dh.nPixels[1] = 40; cube.dh.nChannels = 20;
  allocateDataCube(cube);
  const C = f32(2.5);
  cube.flux.fill(C);
  const rms = estimateNoise(cube);
  console.log('  RMS:', rms, '(expect', C, ')', rms === C ? 'OK' : 'FAIL');

  console.log('\n=== estimateNoise ignores the interior ===');
  // Corrupting only the cube's interior (well inside every corner box's
  // range) must not change the RMS -- proves the corner-only sampling.
  const cube2 = new DataCube();
  cube2.dh.nPixels[0] = 44; cube2.dh.nPixels[1] = 40; cube2.dh.nChannels = 20;
  allocateDataCube(cube2);
  cube2.flux.fill(C);
  const { flatIndxCalc: fic } = require('../ObjectDefinitions/DataCube.js');
  cube2.flux[fic(22, 20, 10, cube2.dh)] = f32(9999.0);
  const rms2 = estimateNoise(cube2);
  console.log('  RMS:', rms2, '(expect', C, ')', rms2 === C ? 'OK' : 'FAIL');
}
