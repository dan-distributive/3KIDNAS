'use strict';

// =============================================================================
// CommonConsts.js
// High-fidelity port of src/StandardMath/BasicConstants.f (CommonConsts module)
//
// PORTING NOTES
// -------------
// All constants are Fortran `real` (32-bit single precision) parameters
// evaluated at compile time. Values are hardcoded here as their exact
// IEEE 754 f32 representations, computed via Python struct.pack/unpack.
//
// Fortran: Pi = 4.*atan(1.)
//   → f32(4.0 * atan(1.0)) = f32(3.14159265...) = 3.141592741012573e+00
//
// Fortran: JyAS_To_MsolPC = 1.24756e+20/(6.0574E5*1.823E18*(2.*Pi/log(256.)))
//   → evaluated entirely in f32 at compile time
//   → 9.970664541469887e-05
// =============================================================================

const f32 = Math.fround;

// Pi = 4.*atan(1.) as Fortran real
const Pi            = f32(Math.PI);                  // 3.141592741012573e+00

// HI 21cm rest frequency in Hz
const HIRestFreq    = f32(1.42040575179e+09);        // 1.420405760000000e+09

// Speed of light in km/s
const lightspeed    = f32(2.99792458e5);             // 2.997924687500000e+05

// Degree to arcsecond conversion
const Degree_To_AS  = f32(3600.0);                   // 3.600000000000000e+03

// Radian to arcsecond conversion
const Radian_To_AS  = f32(206265.0);                 // 2.062650000000000e+05

// Jy/arcsec² to M☉/pc² conversion factor
// Fortran: 1.24756e+20 / (6.0574E5 * 1.823E18 * (2.*Pi/log(256.)))
const JyAS_To_MsolPC = f32(
  f32(1.24756e+20) / f32(
    f32(6.0574e5) * f32(
      f32(1.823e18) * f32(
        f32(2.0) * Pi / f32(Math.log(256.0))
      )
    )
  )
);                                                   // 9.970664541469887e-05

// Hubble constant in km/s/Mpc
const H0            = f32(70.0);                     // 7.000000000000000e+01


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  Pi,
  HIRestFreq,
  lightspeed,
  Degree_To_AS,
  Radian_To_AS,
  JyAS_To_MsolPC,
  H0
};


// ---------------------------------------------------------------------------
// Self-test (node CommonConsts.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const c = module.exports;
  console.log('Pi:             ', c.Pi.toFixed(15));
  console.log('HIRestFreq:     ', c.HIRestFreq.toFixed(2));
  console.log('lightspeed:     ', c.lightspeed.toFixed(4));
  console.log('Degree_To_AS:   ', c.Degree_To_AS);
  console.log('Radian_To_AS:   ', c.Radian_To_AS);
  console.log('JyAS_To_MsolPC: ', c.JyAS_To_MsolPC.toExponential(10));
  console.log('H0:             ', c.H0);
}
