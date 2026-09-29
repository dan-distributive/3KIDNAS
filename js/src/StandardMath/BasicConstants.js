'use strict';

// =============================================================================
// BasicConstants.js
// Port of src/StandardMath/BasicConstants.f (CommonConsts module)
//
// All Fortran `real,parameter` compile-time constants. JyAS_To_MsolPC's
// exact expression (including operation order) matches what
// bootstrap-fit-launcher.js already computes inline for WRKP output
// conversion -- centralized here for reuse by new code (ModellingInitializations.js)
// rather than duplicated again.
// =============================================================================

const f32 = Math.fround;

const Pi           = f32(4.0 * Math.atan(1.0));
const HIRestFreq    = f32(1.42040575179e9);
const Lightspeed    = f32(2.99792458e5);
const Degree_To_AS  = f32(3600.0);
const Radian_To_AS  = f32(206265.0);
// BUG FIX (Dan, 2026): combining this whole expression inside one f32()
// call computes it at full double precision and rounds only ONCE, which
// does NOT match gfortran's compile-time folding of the equivalent real*4
// PARAMETER expression (verified directly: a standalone Fortran program
// printing this constant's hex bit pattern gives 0x38D11999 at
// -O0/-ffp-contract=off, gfortran 15.2.0 arm64; the single-double-then-
// round form here gave 0x38D11998, ONE ULP off -- traced all the way back
// from a real Sigma-parameter-range divergence at optimizer call 3 on
// WALLABY_J100336-262923). js/src/StandardMath/CommonConsts.js already has
// this right, by nesting f32() around each sub-operation in the same
// left-to-right/parenthesized order Fortran evaluates -- reused that exact
// form here instead of introducing a second, differently-computed copy.
const JyAS_To_MsolPC = f32(
  f32(1.24756e20) / f32(
    f32(6.0574e5) * f32(
      f32(1.823e18) * f32(
        f32(2.0) * Pi / f32(Math.log(256.0))
      )
    )
  )
);
const H0            = f32(70.0);

module.exports = {
  Pi, HIRestFreq, Lightspeed, Degree_To_AS, Radian_To_AS, JyAS_To_MsolPC, H0,
};


// ---------------------------------------------------------------------------
// Self-test (node BasicConstants.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  console.log('Pi:', Pi, '(expect', Math.PI, ')');
  console.log('JyAS_To_MsolPC:', JyAS_To_MsolPC);
}
