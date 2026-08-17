'use strict';

// =============================================================================
// FullCircTrig.js
// Port of src/StandardMath/FullCircTrig.f (FullCircTrig module)
//
// Hand-rolled full-circle arctangent: atan(Y/X) plus a quadrant correction,
// NOT Fortran's ATAN2 intrinsic and NOT JS's Math.atan2 -- a different
// (if numerically similar) code path, ported literally rather than swapped
// for the "obviously equivalent" builtin, since the whole point of this
// project is matching Fortran's actual arithmetic, not just its intent.
//
// FIXED FORTRAN BUG (2026, flagged by Dan, reported to Nathan) -- Fortran's
// plain ATAN() intrinsic here was not routed through this project's
// vendored fdlibm the way SIN/COS/ATANH elsewhere in this project are.
// Unlike sin/cos (where V8's Math.sin/Math.cos already happen to be
// bit-exact fdlibm derivatives at f32 precision -- confirmed empirically,
// zero mismatches across 200k random samples), V8's Math.atan is NOT: real
// disagreements against gfortran's native atan() were found for ordinary
// inputs (e.g. atan(-0.3), atan(0.618)), and this is exactly the function
// FullCircATan uses for essentially every angle computed anywhere in this
// pipeline (particle azimuth, shape/Fourier-moment sums, ring geometry).
// Fixed by porting fdlibm's actual atan algorithm (s_atan.c) -- see fdAtan
// in fdlibm.js, verified bit-exact against the real compiled fdlibm_atan.c
// (same FMA fusion as this project's actual C build) across 500k random
// values, 0 mismatches.
// =============================================================================

const f32 = Math.fround;
const Pi  = f32(Math.PI);
const { fdAtan } = require('./fdlibm.js');

// ---------------------------------------------------------------------------
// fullCircATan
// Fortran: FullCircATan(X, Y, Theta)
//
// Theta = atan(Y/X), then quadrant-corrected into [0, 2*Pi):
//   X>=0, Y>=0: no change
//   X<0        : + Pi   (both Y>=0 and Y<0 cases, per the Fortran source)
//   X>=0, Y<0  : + 2*Pi
//
// X=0,Y=0 -> Y/X = 0/0 = NaN -> Theta stays NaN (X>=0 && Y>=0 matches, and
// that branch is a no-op) -- matches Fortran's IEEE754 comparison semantics
// exactly (NaN comparisons are false in both languages, so this is the one
// combination where the "no correction" branch is chosen and NaN survives
// unmodified). Not special-cased -- this is what Fortran actually returns.
// ---------------------------------------------------------------------------
function fullCircATan(x, y) {
  let theta = f32(fdAtan(f32(y / x)));

  if (x >= f32(0.0) && y >= f32(0.0)) {
    // no correction
  } else if (x < f32(0.0) && y >= f32(0.0)) {
    theta = f32(theta + Pi);
  } else if (x < f32(0.0) && y < f32(0.0)) {
    theta = f32(theta + Pi);
  } else if (x >= f32(0.0) && y < f32(0.0)) {
    theta = f32(theta + f32(2.0) * Pi);
  }

  return theta;
}

module.exports = { fullCircATan };


// ---------------------------------------------------------------------------
// Self-test (node FullCircTrig.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const cases = [
    [1, 1, 'Q1'], [-1, 1, 'Q2'], [-1, -1, 'Q3'], [1, -1, 'Q4'],
    [1, 0, '+X axis'], [0, 1, '+Y axis'], [-1, 0, '-X axis'], [0, -1, '-Y axis'],
  ];
  console.log('=== fullCircATan quadrants ===');
  for (const [x, y, label] of cases) {
    const theta = fullCircATan(f32(x), f32(y));
    console.log(`  (${x},${y}) ${label}: theta=${theta.toFixed(6)} (${(theta * 180 / Math.PI).toFixed(2)} deg)`);
  }
  console.log('\n=== degenerate (0,0) -> NaN ===');
  const nanCase = fullCircATan(f32(0), f32(0));
  console.log('  theta:', nanCase, Number.isNaN(nanCase) ? 'OK (NaN, matches Fortran)' : 'FAIL');
}
