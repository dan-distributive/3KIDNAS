module.declare(["./BasicConstants.js"], function (require, exports, module) {
'use strict';

// =============================================================================
// UnitConversions.js
// Port of src/UnitConversions/UnitConversions.f (InputUnitConversionsMod) --
// msolPc2ToJyAS2/jyAS2ToMsolPc2 only, the two routines this project's
// InitialAnalysis-chain port currently needs. RedshiftCalc (also in this
// Fortran file) belongs to the WALLABY-frequency-cubelet path this port
// doesn't reach (see GetMomentMaps.js's header) and isn't ported here.
// =============================================================================

const f32 = Math.fround;
const { JyAS_To_MsolPC } = require('./BasicConstants.js');

// ---------------------------------------------------------------------------
// jyAS2ToMsolPc2
// Fortran: JyAS2_To_MSolPc2(SD_J, SD_M)
// ---------------------------------------------------------------------------
function jyAS2ToMsolPc2(sdJ) {
  return f32(sdJ / JyAS_To_MsolPC);
}

// ---------------------------------------------------------------------------
// msolPc2ToJyAS2
// Fortran: MSolPc2_To_JyAS2(SD_M, SD_J)
// ---------------------------------------------------------------------------
function msolPc2ToJyAS2(sdM) {
  return f32(sdM * JyAS_To_MsolPC);
}

module.exports = { jyAS2ToMsolPc2, msolPc2ToJyAS2 };


// ---------------------------------------------------------------------------
// Self-test (node UnitConversions.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const sdM = f32(2.0);
  const sdJ = msolPc2ToJyAS2(sdM);
  const back = jyAS2ToMsolPc2(sdJ);
  console.log('msolPc2ToJyAS2(2.0):', sdJ);
  console.log('round-trip back    :', back, '(expect ~2.0)', Math.abs(back - 2.0) < 1e-4 ? 'OK' : 'FAIL');
}

});
