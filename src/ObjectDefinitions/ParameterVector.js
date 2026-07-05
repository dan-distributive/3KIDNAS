'use strict';

// =============================================================================
// ParameterVector.js
// High-fidelity port of src/ObjectDefinitions/ParameterVector.f
// (ParameterVectorMod)
//
// PORTING NOTES
// -------------
// Fortran `real` arrays → Float32Array (32-bit single precision).
// Fortran `integer` arrays → Int32Array.
// Fortran scalars `real` → stored as Math.fround() values.
//
// AllocateParamVector initializes PV%Param = -1. (all elements).
// All other arrays are left at their Float32Array/Int32Array defaults (0).
//
// The Fortran TYPE has no constructor — the caller sets nParams before
// calling AllocateParamVector. In JS, allocate() takes nParams directly.
// =============================================================================

const f32 = Math.fround;


// ---------------------------------------------------------------------------
// ParameterVector
// Fortran TYPE ParameterVector — flat parameter array + metadata for optimizer.
//
// Fields:
//   nParams        — number of free parameters
//   param          — current parameter values (initialized to -1)
//   bestParam      — best-fit parameter values
//   paramErr       — parameter uncertainties
//   paramLowerLims — lower bounds per parameter
//   paramUpperLims — upper bounds per parameter
//   paramRange     — search range per parameter (for simplex initialization)
//   cyclicSwitch   — 1 if parameter is cyclic (e.g. PA), 0 otherwise
//   currLike       — current likelihood / chi²
//   bestLike       — best likelihood / chi² found so far
// ---------------------------------------------------------------------------
class ParameterVector {
  constructor() {
    this.nParams        = 0;
    this.param          = null;   // Float32Array(nParams)
    this.bestParam      = null;   // Float32Array(nParams)
    this.paramErr       = null;   // Float32Array(nParams)
    this.paramLowerLims = null;   // Float32Array(nParams)
    this.paramUpperLims = null;   // Float32Array(nParams)
    this.paramRange     = null;   // Float32Array(nParams)
    this.cyclicSwitch   = null;   // Int32Array(nParams)
    this.currLike       = f32(0.0);
    this.bestLike       = f32(0.0);
  }
}


// ---------------------------------------------------------------------------
// allocateParamVector
// Fortran: AllocateParamVector(PV)
// Allocates all arrays to length nParams.
// Initializes param = -1. (matching Fortran: PV%Param=-1.)
// All other arrays default to 0.
// ---------------------------------------------------------------------------
function allocateParamVector(pv) {
  const n         = pv.nParams;
  pv.param        = new Float32Array(n).fill(f32(-1.0));
  pv.bestParam    = new Float32Array(n);
  pv.paramErr     = new Float32Array(n);
  pv.paramLowerLims = new Float32Array(n);
  pv.paramUpperLims = new Float32Array(n);
  pv.paramRange   = new Float32Array(n);
  pv.cyclicSwitch = new Int32Array(n);
}


// ---------------------------------------------------------------------------
// deAllocateParamVector
// Fortran: DeAllocateParamVector(PV)
// In JS, null out arrays (GC handles memory).
// ---------------------------------------------------------------------------
function deAllocateParamVector(pv) {
  pv.param        = null;
  pv.bestParam    = null;
  pv.paramErr     = null;
  pv.paramLowerLims = null;
  pv.paramUpperLims = null;
  pv.paramRange   = null;
  pv.cyclicSwitch = null;
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  ParameterVector,
  allocateParamVector,
  deAllocateParamVector
};


// ---------------------------------------------------------------------------
// Self-test (node ParameterVector.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { ParameterVector, allocateParamVector, deAllocateParamVector }
    = module.exports;

  const pv    = new ParameterVector();
  pv.nParams  = 5;
  allocateParamVector(pv);

  console.log('=== allocateParamVector ===');
  console.log('nParams:', pv.nParams, '(expect 5)');
  console.log('param:', Array.from(pv.param), '(expect all -1)');
  console.log('bestParam:', Array.from(pv.bestParam), '(expect all 0)');
  console.log('cyclicSwitch:', Array.from(pv.cyclicSwitch), '(expect all 0)');
  console.log('currLike:', pv.currLike, '(expect 0)');
  console.log('bestLike:', pv.bestLike, '(expect 0)');

  // Assign some values
  pv.param[0]          = f32(1.5);
  pv.paramLowerLims[0] = f32(0.0);
  pv.paramUpperLims[0] = f32(3.0);
  pv.paramRange[0]     = f32(0.5);
  pv.cyclicSwitch[0]   = 1;
  pv.bestLike          = f32(42.0);

  console.log('\n=== after assignment ===');
  console.log('param[0]:', pv.param[0], '(expect 1.5)');
  console.log('cyclicSwitch[0]:', pv.cyclicSwitch[0], '(expect 1)');
  console.log('bestLike:', pv.bestLike, '(expect 42)');

  deAllocateParamVector(pv);
  console.log('\n=== after dealloc ===');
  console.log('param:', pv.param, '(expect null)');
}
