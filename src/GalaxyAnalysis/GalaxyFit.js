'use strict';

// =============================================================================
// GalaxyFit.js
// High-fidelity port of src/GalaxyAnalysis/GalaxyFit.f (GalaxyFitMod)
//
// PORTING NOTES
// -------------
// Fortran `real` → Math.fround(). All param arrays are Float32Array.
//
// The Fortran uses PipelineGlobals for all shared state. In JS everything
// is passed explicitly via a `state` context object — making the optimizer
// a pure function suitable for DCP worker use.
//
// amoeba (Nelder-Mead) is ported directly from Numerical Recipes Fortran
// rather than using an npm library, for exact algorithm fidelity.
//
// state object (same as FullModelComparison.js plus additional fields):
//   pvIni             — ParameterVector (initial params, read-only)
//   pvModel           — ParameterVector (working copy, modified)
//   pvFirstFit        — ParameterVector (stores pass 1 result)
//   modelTiltedRing   — TiltedRingModel
//   modelDC           — DataCube (model cube)
//   observedDC        — DataCube (observed, read-only)
//   observedBeam      — Beam2D
//   trFittingOptions  — TiltedRingFittingOptions
//   rng               — makeRng() object
//   linearLogSDSwitch — 0 or 1
//   paramToTiltedRing — generalizedParamVectorToTiltedRing
//   ftol              — convergence tolerance (default 0.005)
//   iniGuessWidth     — simplex perturbation scale (default 1.0)
// =============================================================================

const f32 = Math.fround;

const { allocateDataCube }              = require('../ObjectDefinitions/DataCube.js');
const { allocateParamVector, ParameterVector } = require('../ObjectDefinitions/ParameterVector.js');
const { calculate2DBeamKernel }         = require('../ConvolveCube/CalculateBeamKernel.js');
const { tiltedRingModelComparison }     = require('../CompareCubes/FullModelComparison.js');


// ---------------------------------------------------------------------------
// amoeba
// Numerical Recipes Nelder-Mead downhill simplex minimizer.
// Ported directly from the Fortran source for algorithm fidelity.
//
// p:       Float32Array (nParams+1) x nParams — simplex vertices (row-major)
//          p[k][i] = vertex k, parameter i  (1-indexed in Fortran → 0-indexed here)
// y:       Float32Array (nParams+1) — function values at each vertex
// ftol:    convergence tolerance on fractional range of y values
// objFn:   function(paramsArray) → chi2 scalar
//
// Modifies p and y in place. Returns number of iterations used.
//
// Fortran uses 1-indexed arrays: p(nParams+1, nParams), y(nParams+1)
// JS uses 0-indexed: p[k][i] for vertex k, param i
// ---------------------------------------------------------------------------
function amoeba(p, y, nParams, ftol, objFn) {
  const ALPHA = 1.0;    // reflection coefficient
  const BETA  = 0.5;    // contraction coefficient
  const GAMMA = 2.0;    // expansion coefficient
  const ITMAX = 5000;

  const ndim = nParams;
  const mpts = ndim + 1;
  let nfunk  = 0;

  const psum = new Float32Array(ndim);

  // Compute initial psum (sum of all vertex coords per dimension)
  function getPsum() {
    for (let j = 0; j < ndim; j++) {
      let sum = f32(0.0);
      for (let i = 0; i < mpts; i++) sum = f32(sum + f32(p[i][j]));
      psum[j] = sum;
    }
  }

  getPsum();

  // Helper: try a new point by reflecting/expanding/contracting the simplex
  function amotry(ihi, fac) {
    const fac1 = f32(f32(1.0 - fac) / f32(ndim));
    const fac2 = f32(f32(fac1) - fac);
    const ptry = new Float32Array(ndim);
    for (let j = 0; j < ndim; j++) {
      ptry[j] = f32(f32(psum[j]) * fac1 - f32(p[ihi][j]) * fac2);
    }
    const ytry = f32(objFn(ptry));
    nfunk++;
    if (ytry < y[ihi]) {
      y[ihi] = ytry;
      for (let j = 0; j < ndim; j++) {
        psum[j] = f32(f32(psum[j]) + f32(ptry[j]) - f32(p[ihi][j]));
        p[ihi][j] = ptry[j];
      }
    }
    return ytry;
  }

  for (let iter = 0; iter < ITMAX; iter++) {
    // DCP Worker Progress call
    progress(iter/ITMAX);

    // Find highest, second highest, lowest
    let ilo = 0;
    let ihi, inhi;
    if (y[0] > y[1]) { ihi = 0; inhi = 1; }
    else             { ihi = 1; inhi = 0; }
    for (let i = 0; i < mpts; i++) {
      if (y[i] <= y[ilo]) ilo = i;
      if (y[i] >  y[ihi]) { inhi = ihi; ihi = i; }
      else if (y[i] > y[inhi] && i !== ihi) inhi = i;
    }

    // Convergence check: fractional range of y values
    const rtol = f32(f32(2.0) * f32(Math.abs(y[ihi] - y[ilo])) /
      f32(Math.abs(y[ihi]) + Math.abs(y[ilo]) + f32(1.0e-10)));
    if (rtol < f32(ftol)) break;

    if (nfunk >= ITMAX) break;

    // Reflect
    let ytry = amotry(ihi, -ALPHA);
    nfunk++;

    if (ytry <= y[ilo]) {
      // Expansion
      amotry(ihi, GAMMA);
      nfunk++;
    } else if (ytry >= y[inhi]) {
      // Contraction
      const ysave = y[ihi];
      ytry = amotry(ihi, BETA);
      nfunk++;
      if (ytry >= ysave) {
        // Contract around lowest point
        for (let i = 0; i < mpts; i++) {
          if (i !== ilo) {
            for (let j = 0; j < ndim; j++) {
              psum[j] = f32(f32(0.5) * f32(p[i][j] + p[ilo][j]));
              p[i][j] = psum[j];
            }
            y[i] = f32(objFn(p[i]));
            nfunk++;
          }
        }
        getPsum();
      }
    }
  }

  return nfunk;
}


// ---------------------------------------------------------------------------
// makeParamGuessArray
// Fortran: MakeParamGuessArray(PredictedPV, ParamGuesses, ndim, idum,
//                               lambda, StrictEstimate)
//
// Generates nParams+1 starting points for the simplex:
//   - Row 0: current best estimate (PVModel.param)
//   - Rows 1..nParams: random perturbations within lambda * paramRange
//
// Handles cyclic parameters (wrapping) and bounds rejection.
// StrictEstimate=1: reject out-of-bounds always
// StrictEstimate=0: accept out-of-bounds after 200 tries per parameter
//
// Fortran layout: ParamGuesses(nParams+1, nParams) — (row, col), 1-indexed
// JS layout: paramGuesses[k][i] — (vertex, param), 0-indexed
// ---------------------------------------------------------------------------
function makeParamGuessArray(pv, rng, lambda, strictEstimate) {
  const ndim    = pv.nParams;
  // paramGuesses[k][i]: k=0..ndim, i=0..ndim-1
  const guesses = Array.from({ length: ndim + 1 }, () => new Float32Array(ndim));

  // First row: current best estimate
  for (let i = 0; i < ndim; i++) {
    guesses[0][i] = f32(pv.param[i]);
  }

  // Remaining rows: random perturbations
  for (let k = 1; k <= ndim; k++) {
    for (let i = 0; i < ndim; i++) {
      let parCounter        = 0;
      let acceptBeyondLimits = 0;
      let val;

      do {
        const lambdaPar = f32(f32(f32(2.0) * f32(rng.ran2()) - f32(1.0))
          * f32(lambda * f32(pv.paramRange[i])));
        val = f32(f32(pv.param[i]) + lambdaPar);
        parCounter++;

        if (strictEstimate === 0 && parCounter >= 200) {
          acceptBeyondLimits = 1;
        }

        // Cyclic wrapping (e.g. PA)
        if (pv.cyclicSwitch[i] === 1) {
          while (val < f32(pv.paramLowerLims[i])) {
            val = f32(val + f32(pv.paramUpperLims[i]));
          }
          while (val > f32(pv.paramUpperLims[i])) {
            val = f32(val - f32(pv.paramUpperLims[i]));
          }
        }

        if (acceptBeyondLimits === 1) break;
      } while (
        val < f32(pv.paramLowerLims[i]) ||
        val > f32(pv.paramUpperLims[i])
      );

      guesses[k][i] = val;
    }
  }

  return guesses;
}


// ---------------------------------------------------------------------------
// downhillSimplexRun
// Fortran: DownhillSimplexRun(nParams, paramGuesses, chiArray)
//
// Evaluates chi² at each simplex vertex, then runs amoeba.
// Updates pvModel with the best result.
// ---------------------------------------------------------------------------
function downhillSimplexRun(paramGuesses, chiArray, state) {
  const pv  = state.pvModel;
  const n   = pv.nParams;

  // Evaluate chi² at each starting vertex
  for (let i = 0; i <= n; i++) {
    for (let j = 0; j < n; j++) {
      pv.param[j] = f32(paramGuesses[i][j]);
    }
    chiArray[i] = f32(tiltedRingModelComparison(paramGuesses[i], state));
  }

  // Run Nelder-Mead
  amoeba(paramGuesses, chiArray, n, state.ftol,
    params => f32(tiltedRingModelComparison(params, state))
  );

  // Store best result in pvModel
  pv.bestLike = f32(chiArray[0]);
  for (let j = 0; j < n; j++) {
    pv.param[j] = f32(paramGuesses[0][j]);
  }
}


// ---------------------------------------------------------------------------
// galaxyFit_Simple
// Fortran: GalaxyFit_Simple(CatItem)
//
// Two-pass Nelder-Mead optimizer:
//   Pass 1: wide search (iniGuessWidth=1.0, ftol=0.005, strictEstimate=1)
//   Pass 2: refined search (iniGuessWidth=0.5, ftol/5, strictEstimate=0)
//
// Returns pvModel (best-fit parameter vector after both passes).
//
// state fields used:
//   pvIni, pvModel, pvFirstFit, observedDC, observedBeam,
//   modelDC, modelTiltedRing, trFittingOptions, rng,
//   linearLogSDSwitch, paramToTiltedRing
// ---------------------------------------------------------------------------
function galaxyFit_Simple(state) {
  const { pvIni, observedDC, observedBeam } = state;

  // Step 1: calculate beam kernel (if not already done)
  calculate2DBeamKernel(observedBeam, observedDC.dh.pixelSize);

  // Step 2: allocate model cube matching observed
  state.modelDC.dh = Object.assign(
    Object.create(Object.getPrototypeOf(observedDC.dh)),
    observedDC.dh
  );
  state.modelDC.dh.nPixels     = new Int32Array(observedDC.dh.nPixels);
  state.modelDC.dh.pixelSize   = new Float32Array(observedDC.dh.pixelSize);
  state.modelDC.dh.refLocation = new Float32Array(observedDC.dh.refLocation);
  state.modelDC.dh.refVal      = new Float32Array(observedDC.dh.refVal);
  allocateDataCube(state.modelDC);

  // Step 3: copy pvIni → pvModel
  const pvModel      = state.pvModel;
  pvModel.nParams    = pvIni.nParams;
  allocateParamVector(pvModel);
  pvModel.param.set(pvIni.param);
  pvModel.paramLowerLims.set(pvIni.paramLowerLims);
  pvModel.paramUpperLims.set(pvIni.paramUpperLims);
  pvModel.cyclicSwitch.set(pvIni.cyclicSwitch);
  pvModel.paramRange.set(pvIni.paramRange);

  const n = pvModel.nParams;

  // Step 4: evaluate initial chi²
  const chi2Ini = f32(tiltedRingModelComparison(Array.from(pvModel.param), state));
  console.log('Initial model fit:', chi2Ini);

  // ---- Pass 1: wide search ----
  state.ftol         = f32(0.005);
  state.iniGuessWidth = f32(1.0);

  let paramGuesses = makeParamGuessArray(pvModel, state.rng, state.iniGuessWidth, 1);
  let chiArray     = new Float32Array(n + 1);

  downhillSimplexRun(paramGuesses, chiArray, state);

  // Store pass 1 result
  const pvFirstFit    = state.pvFirstFit;
  pvFirstFit.nParams  = pvIni.nParams;
  allocateParamVector(pvFirstFit);
  for (let j = 0; j < n; j++) pvFirstFit.param[j] = f32(paramGuesses[0][j]);
  pvFirstFit.bestLike = f32(chiArray[0]);

  // ---- Pass 2: refined search ----
  state.iniGuessWidth = f32(0.5);
  state.ftol          = f32(state.ftol / 5.0);

  paramGuesses = makeParamGuessArray(pvModel, state.rng, state.iniGuessWidth, 0);
  chiArray     = new Float32Array(n + 1);

  downhillSimplexRun(paramGuesses, chiArray, state);

  // Final best-fit is in pvModel (updated by downhillSimplexRun)
  return pvModel;
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  galaxyFit_Simple,
  downhillSimplexRun,
  makeParamGuessArray,
  amoeba
};


// ---------------------------------------------------------------------------
// Self-test (node GalaxyFit.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const f32 = Math.fround;

  // Simple 2D quadratic test — no astronomy needed
  // Minimize f(x,y) = (x-3)^2 + (y-2)^2, minimum at [3,2]
  console.log('=== amoeba self-test (2D quadratic) ===');
  const n = 2;
  const p = [
    new Float32Array([0.0, 0.0]),
    new Float32Array([1.0, 0.0]),
    new Float32Array([0.0, 1.0])
  ];
  const y = new Float32Array([
    (0-3)**2 + (0-2)**2,
    (1-3)**2 + (0-2)**2,
    (0-3)**2 + (1-2)**2
  ]);

  const nfunk = amoeba(p, y, n, 1e-6,
    params => f32(f32(params[0]-3)**2 + f32(params[1]-2)**2)
  );

  console.log('minimum at:', Array.from(p[0]).map(v=>v.toFixed(4)), '(expect [3,2])');
  console.log('f(min):', y[0].toFixed(8), '(expect ~0)');
  console.log('nfunk:', nfunk);
  console.log('converged:', Math.abs(p[0][0]-3) < 0.01 && Math.abs(p[0][1]-2) < 0.01 ? 'OK' : 'FAIL');

  // Test makeParamGuessArray
  console.log('\n=== makeParamGuessArray ===');
  const { makeRng } = require('../StandardMath/random.js');
  const { ParameterVector, allocateParamVector } = require('../ObjectDefinitions/ParameterVector.js');
  const pv = new ParameterVector();
  pv.nParams = 3;
  allocateParamVector(pv);
  pv.param[0] = f32(1.0); pv.paramLowerLims[0] = f32(0.0); pv.paramUpperLims[0] = f32(2.0); pv.paramRange[0] = f32(0.5);
  pv.param[1] = f32(5.0); pv.paramLowerLims[1] = f32(0.0); pv.paramUpperLims[1] = f32(10.0); pv.paramRange[1] = f32(2.0);
  pv.param[2] = f32(3.0); pv.paramLowerLims[2] = f32(0.0); pv.paramUpperLims[2] = f32(6.28); pv.paramRange[2] = f32(1.0);
  pv.cyclicSwitch[2] = 1; // cyclic

  const rng     = makeRng(-1);
  const guesses = makeParamGuessArray(pv, rng, f32(1.0), 1);
  console.log('nGuesses:', guesses.length, '(expect 4 = nParams+1)');
  console.log('guess[0] (initial):', Array.from(guesses[0]).map(v=>v.toFixed(4)), '(expect [1,5,3])');
  let allInBounds = true;
  for (let k = 1; k <= 3; k++) {
    for (let i = 0; i < 3; i++) {
      if (guesses[k][i] < pv.paramLowerLims[i] || guesses[k][i] > pv.paramUpperLims[i]) {
        if (pv.cyclicSwitch[i] !== 1) { allInBounds = false; }
      }
    }
  }
  console.log('all guesses in bounds:', allInBounds ? 'OK' : 'FAIL');
}
