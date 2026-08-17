module.declare([], function (require, exports, module) {
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
const JyAS_To_MsolPC = f32(1.24756e20 / (6.0574e5 * 1.823e18 * (2.0 * Pi / Math.log(256.0))));
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

});
