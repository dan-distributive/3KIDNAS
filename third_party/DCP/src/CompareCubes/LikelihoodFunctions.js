'use strict';

// =============================================================================
// likelihood.js
// High-fidelity port of src/CompareCubes/LikelihoodFunctions.f (LikelihoodMod)
//
// PORTING NOTES
// -------------
// Fortran `real` is 32-bit single precision. All arithmetic uses Math.fround()
// to match Fortran rounding behaviour at each step.
//
// The Fortran module declares a function pointer `LikePoint` that gets
// assigned elsewhere in the pipeline (to Chi2Calc in the standard fit).
// In JS this is represented as a module-level variable `likePoint` that
// can be reassigned before use.
//
// Array indexing: Fortran arrays are 0-indexed here (Model(0:nElements-1))
// which maps directly to JS arrays — no offset needed.
//
// FUNCTION POINTER INTERFACE
// --------------------------
// All three functions share the same signature:
//   fn(chi2Out, nElements, model, obs, uncertainties)
// where chi2Out is returned as a value (Fortran INTENT(OUT) scalar).
// In JS we return chi2 directly rather than using an output parameter.
// =============================================================================

const f32 = Math.fround;


// ---------------------------------------------------------------------------
// chi2Calc
// Standard chi-squared statistic.
// Fortran: Chi2Calc(chi2, nElements, Model, Obs, Uncertainties)
//
// chi2 = sum_i [ (Model(i) - Obs(i))^2 / Uncertainties(i)^2 ]
//
// Arrays are 0-indexed, length nElements (Fortran: 0:nElements-1).
// Returns chi2 as a 32-bit float matching Fortran `real` accumulation.
// ---------------------------------------------------------------------------
function chi2Calc(nElements, model, obs, uncertainties) {
  let chi2 = f32(0.0);
  for (let i = 0; i < nElements; i++) {
    const diff  = f32(f32(model[i]) - f32(obs[i]));
    const sigma = f32(uncertainties[i]);
    chi2 = f32(chi2 + f32(f32(diff * diff) / f32(sigma * sigma)));
  }
  return chi2;
}


// ---------------------------------------------------------------------------
// logChi2Calc
// Log10 of the standard chi-squared statistic.
// Fortran: LogChi2Calc — calls Chi2Calc then takes log10.
// ---------------------------------------------------------------------------
function logChi2Calc(nElements, model, obs, uncertainties) {
  const chi2 = chi2Calc(nElements, model, obs, uncertainties);
  return f32(Math.log10(chi2));
}


// ---------------------------------------------------------------------------
// chi2Calc_logElements
// Chi-squared in log10 space — compares log10(Model) vs log10(Obs).
// Used when surface densities span many orders of magnitude.
//
// Small = 1e-10 floor applied to both Model and Obs before log.
// chi2 = sum_i [ (log10(M_i) - log10(O_i))^2 / Uncertainties(i)^2 ]
// ---------------------------------------------------------------------------
function chi2Calc_logElements(nElements, model, obs, uncertainties) {
  const Small = f32(1.0e-10);
  let chi2 = f32(0.0);
  for (let i = 0; i < nElements; i++) {
    let M = f32(model[i]);
    let O = f32(obs[i]);
    if (M <= Small) M = Small;
    if (O <= Small) O = Small;
    const logDiff = f32(f32(Math.log10(M)) - f32(Math.log10(O)));
    const sigma   = f32(uncertainties[i]);
    chi2 = f32(chi2 + f32(f32(logDiff * logDiff) / f32(sigma * sigma)));
  }
  return chi2;
}


// ---------------------------------------------------------------------------
// likePoint — function pointer (matches Fortran LikePoint procedure pointer)
// Assign before use. Default: chi2Calc (standard fit).
//
// Usage:
//   likePoint = chi2Calc;                  // standard chi^2
//   likePoint = logChi2Calc;               // log chi^2
//   likePoint = chi2Calc_logElements;      // log-space chi^2
//
//   const chi2 = likePoint(nElements, model, obs, uncertainties);
// ---------------------------------------------------------------------------
let likePoint = chi2Calc;


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  chi2Calc,
  logChi2Calc,
  chi2Calc_logElements,
  get likePoint()  { return likePoint; },
  set likePoint(fn) { likePoint = fn; }
};


// ---------------------------------------------------------------------------
// Self-test (node likelihood.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const n = 5;
  const model        = new Float32Array([1.0, 2.0, 3.0, 4.0, 5.0]);
  const obs          = new Float32Array([1.1, 1.9, 3.2, 3.8, 5.1]);
  const uncertainties = new Float32Array([0.1, 0.1, 0.1, 0.1, 0.1]);

  console.log('=== chi2Calc ===');
  console.log(chi2Calc(n, model, obs, uncertainties).toFixed(10));

  console.log('\n=== logChi2Calc ===');
  console.log(logChi2Calc(n, model, obs, uncertainties).toFixed(10));

  console.log('\n=== chi2Calc_logElements ===');
  console.log(chi2Calc_logElements(n, model, obs, uncertainties).toFixed(10));

  console.log('\n=== likePoint (default = chi2Calc) ===');
  console.log(likePoint(n, model, obs, uncertainties).toFixed(10));
}
