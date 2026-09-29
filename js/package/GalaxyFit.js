module.declare(["./DataCube.js","./ParameterVector.js","./CalculateBeamKernel.js","./FullModelComparison.js"], function (require, exports, module) {
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

const TRACE_DEBUG = typeof process !== 'undefined' && process.env && process.env.TRACE_DEBUG === '1';

const { allocateDataCube }              = require('./DataCube.js');
const { allocateParamVector, ParameterVector } = require('./ParameterVector.js');
const { calculate2DBeamKernel }         = require('./CalculateBeamKernel.js');
const { tiltedRingModelComparison }     = require('./FullModelComparison.js');


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
function amoeba(p, y, nParams, ftol, objFn, onProgress) {
  const ITMAX = 5000;
  const TINY  = f32(1.0e-10);
  const ndim  = nParams;
  const mpts  = ndim + 1;

  const psum = new Float32Array(ndim);
  let iter = 0;
  let noConvergence = false;

  // Stall detector (Dan, 2026-09-16; DEFAULTED OFF 2026-09-28): a practical
  // compute-saving early exit, NOT a numerical-parity fix, and NOT present
  // in Fortran's DownhillSimplex.f -- Fortran has no equivalent and always
  // grinds to ITMAX=5000 regardless of stalling. Root cause of the stalls
  // this was built to catch: Monte Carlo sampling noise in the
  // particle-placement objective function sitting at or above ftol's
  // threshold for some resampled datasets -- present in principle on EITHER
  // platform (Fortran's own rtol trajectory for one such case wobbled
  // non-monotonically in the 0.001-0.003 range before happening to dip
  // under ftol=0.001 by chance at iter=39; it isn't immune, it's just been
  // lucky across the runs checked so far). Confirmed 2026-09-28: at
  // cloudDensity=20, disabling this detector and grinding to ITMAX=5000
  // does NOT produce convergence either (rtol genuinely never dips below
  // ftol at that noise level) -- so it was never masking a bug, but
  // defaulting it on made JS's behavior diverge from Fortran's (JS reports
  // noConvergence sooner than Fortran would ever report it, since Fortran
  // has no early exit at all). Defaulted OFF so JS's amoeba matches
  // Fortran's iteration behavior exactly; opt back in with
  // JS_ENABLE_STALL_DETECTOR=1 if the compute savings are wanted for a
  // production/batch context where matching Fortran's iteration count
  // doesn't matter.
  // STALL_WINDOW/STALL_REL_EPS chosen conservatively (window wide relative
  // to a typical shrink-step's iter+=ndim cost; epsilon an order of
  // magnitude tighter than ftol) so this never fires during genuine slow
  // convergence -- only once real progress has demonstrably stopped.
  const stallDetectorEnabled = typeof process !== 'undefined' && process.env
    && process.env.JS_ENABLE_STALL_DETECTOR;
  const STALL_WINDOW = Math.max(300, 20 * ndim);
  const STALL_REL_EPS = f32(ftol) / 10;
  let stallBestY = Infinity;
  let stallBestIter = 0;

  function computePsum() {                    // Fortran label 1
    for (let j = 0; j < ndim; j++) {
      let s = f32(0.0);
      for (let i = 0; i < mpts; i++) s = f32(s + f32(p[i][j]));
      psum[j] = s;
    }
  }

  function amotry(ihi, fac) {                 // returns ytry; mutates in place
    const fac1 = f32(f32(1.0 - fac) / f32(ndim));
    const fac2 = f32(fac1 - fac);
    const ptry = new Float32Array(ndim);
    // Fortran: ptry(j)=psum(j)*fac1-p(ihi,j)*fac2 -- p/psum/fac1/fac2/ptry are
    // all REAL(4), so each individual multiply and the subtract is its own
    // float32-rounded hardware op, not "compute in double, round once at the
    // end". Each intermediate needs its own f32() wrap to match.
    for (let j = 0; j < ndim; j++)
      ptry[j] = f32(f32(f32(psum[j]) * fac1) - f32(f32(p[ihi][j]) * fac2));
    const ytry = f32(objFn(ptry));
    if (ytry < y[ihi]) {
      y[ihi] = ytry;
      for (let j = 0; j < ndim; j++) {
        // Fortran: psum(j)=psum(j)-p(ihi,j)+ptry(j) -- same per-op rounding.
        psum[j]   = f32(f32(f32(psum[j]) - f32(p[ihi][j])) + f32(ptry[j]));
        p[ihi][j] = ptry[j];
      }
    }
    return ytry;                              // NR: no nfunk bookkeeping here
  }

  computePsum();                              // label 1

  for (;;) {                                  // label 2 loop
    if (onProgress) onProgress(Math.min(iter / ITMAX, 1));

    let ilo = 0, ihi, inhi;
    if (y[0] > y[1]) { ihi = 0; inhi = 1; } else { ihi = 1; inhi = 0; }
    for (let i = 0; i < mpts; i++) {
      if (y[i] <= y[ilo]) ilo = i;
      if (y[i] >  y[ihi]) { inhi = ihi; ihi = i; }
      else if (y[i] > y[inhi] && i !== ihi) inhi = i;
    }

    // Fortran: rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY) --
    // y/TINY are REAL(4), so the subtract, the *2., the two-step sum, and
    // the final divide are each their own float32-rounded op (abs() itself
    // needs no extra rounding -- a sign flip on an already-exact value).
    const absDiff = Math.abs(f32(y[ihi] - y[ilo]));
    const numer   = f32(f32(2.0) * absDiff);
    const denom   = f32(f32(Math.abs(y[ihi]) + Math.abs(y[ilo])) + TINY);
    const rtol    = f32(numer / denom);
    if (typeof process !== 'undefined' && process.env && process.env.AMOEBA_DEBUG) {
      console.log('Current tolerance', iter, rtol, y[ihi], y[ilo]);
    }
    if (typeof process !== 'undefined' && process.env && process.env.TRACE_DUMP_PRECONV) {
      require('fs').appendFileSync('AmoebaTraceJS.txt',
        `iter=${iter} rtol=${rtol.toExponential(19)} yhi=${y[ihi].toExponential(19)} ylo=${y[ilo].toExponential(19)}\n`);
    }
    if (rtol < f32(ftol)) {                   // converged: swap best into slot 0
      let s = y[0]; y[0] = y[ilo]; y[ilo] = s;
      for (let j = 0; j < ndim; j++) { const t = p[0][j]; p[0][j] = p[ilo][j]; p[ilo][j] = t; }
      break;
    }
    if (iter >= ITMAX) { noConvergence = true; break; }

    // Stall detector -- see its declaration above for why this exists and
    // why it's safe. Relative improvement measured against |stallBestY| so
    // it scales with chi2's own magnitude, matching rtol's own normalization.
    const relImprovement = Number.isFinite(stallBestY)
      ? f32(f32(stallBestY - y[ilo]) / f32(Math.abs(stallBestY) + TINY))
      : Infinity;                              // first pass: always seed
    if (y[ilo] < stallBestY && relImprovement > STALL_REL_EPS) {
      stallBestY = y[ilo];
      stallBestIter = iter;
    } else if (stallDetectorEnabled && iter - stallBestIter > STALL_WINDOW) {
      let s = y[0]; y[0] = y[ilo]; y[ilo] = s;
      for (let j = 0; j < ndim; j++) { const t = p[0][j]; p[0][j] = p[ilo][j]; p[ilo][j] = t; }
      noConvergence = true;
      break;
    }

    iter += 2;
    let ytry = amotry(ihi, -1.0);
    if (ytry <= y[ilo]) {
      amotry(ihi, 2.0);                       // expansion
    } else if (ytry >= y[inhi]) {
      const ysave = y[ihi];
      ytry = amotry(ihi, 0.5);                // contraction
      if (ytry >= ysave) {                    // shrink toward ilo
        for (let i = 0; i < mpts; i++) {
          if (i !== ilo) {
            for (let j = 0; j < ndim; j++) {
              psum[j] = f32(f32(0.5) * f32(p[i][j] + p[ilo][j]));
              p[i][j] = psum[j];
            }
            // Historical bug fix, NOT a current JS-vs-Fortran divergence:
            // DownhillSimplex.f's shrink step used to call funk(psum,ytry)
            // here without ever storing ytry into y(i), leaving y() stale/
            // inconsistent with the just-shrunk p() positions. Once
            // triggered this could freeze rtol permanently (ihi/ilo/rtol
            // recomputed from stale y() next pass), burning iterations with
            // zero real progress until ITMAX -- confirmed live via
            // AMOEBA_DEBUG when this JS port was first written: y[ihi]/
            // y[ilo] frozen identical across dozens of iterations while
            // iter climbed by exactly 17 (2 + ndim) each pass. Canonical
            // Numerical Recipes amoeba DOES store the recomputed value here.
            // DownhillSimplex.f:98-103 has SINCE been fixed to match (its
            // own "BUG FIX" comment there) -- both implementations agree
            // now, so this is no longer a source of fortran/js divergence.
            y[i] = f32(objFn(psum));
          }
        }
        iter += ndim;
        computePsum();                        // goto 1: rebuild psum
      }
    } else {
      iter -= 1;
    }
  }

  return { iter, noConvergence };
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
        // Fortran: lambdaPar=lambda*ParamRange(j) [round]; lambdaPar=
        // (2*ran2(idum)-1.)*lambdaPar [round the 2*ran2, round the -1.,
        // round the final multiply] -- REAL(4) throughout, each op its own
        // hardware rounding, not one double-precision expression rounded once.
        const rangeTerm = f32(lambda * f32(pv.paramRange[i]));
        const randTerm  = f32(f32(f32(2.0) * f32(rng.ran2())) - f32(1.0));
        const lambdaPar = f32(randTerm * rangeTerm);
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
// Fortran: DownhillSimplexRun(nParams, paramGuesses, chiArray, onProgress)
//
// Evaluates chi² at each simplex vertex, then runs amoeba.
// Updates pvModel with the best result.
// ---------------------------------------------------------------------------
function downhillSimplexRun(paramGuesses, chiArray, state, onProgress) {
  const pv = state.pvModel, n = pv.nParams;
  for (let i = 0; i <= n; i++) {
    for (let j = 0; j < n; j++) pv.param[j] = f32(paramGuesses[i][j]);
    chiArray[i] = f32(tiltedRingModelComparison(paramGuesses[i], state));
  }
  if (typeof process !== 'undefined' && process.env && process.env.TRACE_DUMP_PRECONV) {
    const lines = [];
    for (let i = 0; i <= n; i++) lines.push(`vertex${i} ${chiArray[i].toExponential(19)}`);
    require('fs').appendFileSync('SimplexTraceJS.txt', lines.join('\n') + '\n---\n');
  }
  const { iter, noConvergence } = amoeba(
    paramGuesses, chiArray, n, state.ftol,
    params => f32(tiltedRingModelComparison(params, state)), onProgress);
  pv.bestLike = f32(chiArray[0]);
  for (let j = 0; j < n; j++) pv.param[j] = f32(paramGuesses[0][j]);
  return { iter, noConvergence };
}


// ---------------------------------------------------------------------------
// galaxyFit_Simple
// Fortran: GalaxyFit_Simple(CatItem)
//
// Two-pass Nelder-Mead optimizer:
//   Pass 1: wide search (iniGuessWidth=1.0, ftol=0.005, strictEstimate=1)
//   Pass 2: refined search (iniGuessWidth=0.25, ftol/5, strictEstimate=0)
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

  // Progress reporting: strictly non-decreasing + throttled to 0.1% steps.
  // Swallows the pass-1→pass-2 reset (pass 2 restarts at 0, so it stays
  // silent until it climbs back above pass 1's peak) and avoids firing
  // on every tight iteration.
  let progHigh = -1;                        // force the first call through
  const report = (frac) => {
    if (frac >= progHigh + 0.001) {         // only if advanced ≥ 0.1%
      progHigh = frac;
      if (typeof progress === 'function') progress(frac);
    }
  };

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

  // Diagnostic (Dan probe, 2026-09-27): force pvModel.param and the RNG's
  // idum to exact externally-supplied values before the single unperturbed
  // chi2Ini evaluation below. Proven result (WALLABY_J100336-262923,
  // cdens=20, seed=42, bootstrap realization 0): forcing call 1's inputs
  // to Fortran's exact bits, then letting the rest of the fit run
  // NATURALLY (unforced), reproduced all 127 of Fortran's evaluations
  // bit-for-bit (all 13 params, hex-for-hex) to convergence. This confirms
  // the per-evaluation math (particle generation, convolution, chi2,
  // amoeba, makeParamGuessArray's perturbation) has zero remaining bugs --
  // the ENTIRE naturally-observed divergence traces to this one upstream
  // cause: pvIni itself differing by 1-30 ULP in 3 of 13 params, inherited
  // from the resampled cube's already-irreducible ~1-ULP float32 noise
  // floor (see BuildPhysCoordsArray/GetPhysCoords history). No effect
  // unless JS_FORCE_PVINI_HEX_PATH is set. Scoped to
  // JS_OVERRIDE_REALIZATION_INDEX (same convention as
  // JS_SIMPLEX_OVERRIDE_PATH below) so it doesn't corrupt the anchor fit
  // or other realizations.
  const forceRealizationIndex = process.env.JS_OVERRIDE_REALIZATION_INDEX != null
    ? parseInt(process.env.JS_OVERRIDE_REALIZATION_INDEX, 10) : null;
  const forceScopeOk = forceRealizationIndex == null
    || state.realizationIndex === forceRealizationIndex;
  if (process.env.JS_FORCE_PVINI_HEX_PATH && forceScopeOk) {
    const hexLines = require('fs').readFileSync(process.env.JS_FORCE_PVINI_HEX_PATH, 'utf8')
      .trim().split('\n');
    const fBuf = new ArrayBuffer(4);
    const fU32 = new Uint32Array(fBuf);
    const fF32 = new Float32Array(fBuf);
    for (let j = 0; j < n && j < hexLines.length; j++) {
      fU32[0] = parseInt(hexLines[j], 16);
      pvModel.param[j] = fF32[0];
    }
    const fMsg = `JS_FORCE_PVINI applied: ${hexLines.length} values from ${process.env.JS_FORCE_PVINI_HEX_PATH}\n`;
    if (process.env.TRACE_DEBUG_FINALVEC_FILE) {
      require('fs').appendFileSync(process.env.TRACE_DEBUG_FINALVEC_FILE + '.debug', fMsg);
    } else {
      console.error(fMsg.trim());
    }
  }
  if (process.env.JS_FORCE_IDUM != null && forceScopeOk) {
    state.rng.state.ran2State.idum = parseInt(process.env.JS_FORCE_IDUM, 10);
    const iMsg = `JS_FORCE_IDUM applied: ${state.rng.state.ran2State.idum}\n`;
    if (process.env.TRACE_DEBUG_FINALVEC_FILE) {
      require('fs').appendFileSync(process.env.TRACE_DEBUG_FINALVEC_FILE + '.debug', iMsg);
    } else {
      console.error(iMsg.trim());
    }
  }

  // Step 4: evaluate initial chi²
  const chi2Ini = f32(tiltedRingModelComparison(Array.from(pvModel.param), state));
  console.log('Initial model fit:', chi2Ini);
  if (typeof process !== 'undefined' && process.env && process.env.TRACE_DUMP_PRECONV) {
    require('fs').appendFileSync('Chi2IniTraceJS.txt', `chi2Ini ${chi2Ini.toExponential(19)}\n`);
  }

  // ---- Pass 1: wide search ----
  state.ftol         = f32(0.005);
  state.iniGuessWidth = f32(1.0);

  let paramGuesses = makeParamGuessArray(pvModel, state.rng, state.iniGuessWidth, 1);
  let chiArray     = new Float32Array(n + 1);

  downhillSimplexRun(paramGuesses, chiArray, state, report);

  // Store pass 1 result
  const pvFirstFit    = state.pvFirstFit;
  pvFirstFit.nParams  = pvIni.nParams;
  allocateParamVector(pvFirstFit);
  for (let j = 0; j < n; j++) pvFirstFit.param[j] = f32(paramGuesses[0][j]);
  pvFirstFit.bestLike = f32(chiArray[0]);

  // ---- Pass 2: refined search ----
  state.iniGuessWidth = f32(0.25);
  state.ftol          = f32(state.ftol / 5.0);

  paramGuesses = makeParamGuessArray(pvModel, state.rng, state.iniGuessWidth, 0);
  chiArray     = new Float32Array(n + 1);

  // One-off diagnostic (Dan, 2026): force-feed Fortran's own exact pass-2
  // starting simplex (dumped by GalaxyFit.f's matching
  // FORTRAN_SIMPLEX_DUMP_PATH block -- one hex-encoded float32 per line,
  // row-major over [vertex][param], (n+1)*n lines total) instead of this
  // platform's own makeParamGuessArray() output. Controlled experiment:
  // if amoeba then takes the SAME iteration path/count to the SAME
  // answer Fortran got from this exact simplex, amoeba itself is fine and
  // any real-world divergence traces to the two platforms starting from
  // very slightly different points upstream (already-characterized ~1
  // float32 ULP noise). If it STILL diverges given an IDENTICAL starting
  // simplex, that's a genuine bug in amoeba or the objective function
  // evaluation itself, not just Nelder-Mead's known sensitivity to tiny
  // perturbations. No effect unless JS_SIMPLEX_OVERRIDE_PATH is set.
  // Scoped to JS_OVERRIDE_REALIZATION_INDEX (Dan, 2026-09-16 fix): without
  // this, the override applied to EVERY galaxyFit_Simple call -- the
  // initial fit AND all 5 bootstrap realizations, not just the one being
  // investigated -- corrupting the initial fit ("No best fit model made")
  // and cascading into every realization derived from it.
  const overrideRealizationIndex = process.env.JS_OVERRIDE_REALIZATION_INDEX != null
    ? parseInt(process.env.JS_OVERRIDE_REALIZATION_INDEX, 10) : null;
  const overrideScopeOk = overrideRealizationIndex == null
    || state.realizationIndex === overrideRealizationIndex;
  if (process.env.JS_SIMPLEX_OVERRIDE_PATH && overrideScopeOk) {
    const hexLines = require('fs').readFileSync(process.env.JS_SIMPLEX_OVERRIDE_PATH, 'utf8')
      .trim().split('\n');
    const ovBuf = new ArrayBuffer(4);
    const ovU32 = new Uint32Array(ovBuf);
    const ovF32 = new Float32Array(ovBuf);
    let hexIdx = 0;
    for (let i = 0; i <= n; i++) {
      for (let j = 0; j < n; j++) {
        ovU32[0] = parseInt(hexLines[hexIdx++], 16);
        paramGuesses[i][j] = ovF32[0];
      }
    }
    // Pass-2-only gate for FullModelComparison.js's idum-sequence override
    // (Dan, 2026-09-16 fix): realizationIndex alone isn't enough to scope
    // that override to pass 2's vertex loop -- pass 1 shares the same
    // realizationIndex and its own natural (unforced) evalCount range can
    // overlap the target range, corrupting pass 1's trajectory and shifting
    // where pass 2 actually starts. This flag makes the idum override
    // strictly pass-2-only regardless of evalCount.
    state._simplexOverrideActive = true;
    const overrideMsg = `JS_SIMPLEX_OVERRIDE applied: ${hexIdx} values from ${process.env.JS_SIMPLEX_OVERRIDE_PATH}\n`;
    // console.error from inside a worker_threads Worker races
    // worker.terminate() -- see FINALVEC's identical comment/fix below.
    if (process.env.TRACE_DEBUG_FINALVEC_FILE) {
      require('fs').appendFileSync(process.env.TRACE_DEBUG_FINALVEC_FILE, overrideMsg);
    } else {
      console.error(overrideMsg.trim());
    }
  }

  const { iter: pass2Iter, noConvergence } = downhillSimplexRun(paramGuesses, chiArray, state, report);
  if (process.env.JS_SIMPLEX_OVERRIDE_PATH || TRACE_DEBUG) {
    const iterMsg = `PASS2_ITER ${pass2Iter} noConvergence=${noConvergence} finalChi2=${chiArray[0]}\n`;
    if (process.env.TRACE_DEBUG_FINALVEC_FILE) {
      require('fs').appendFileSync(process.env.TRACE_DEBUG_FINALVEC_FILE, iterMsg);
    } else {
      console.error(iterMsg.trim());
    }
  }

  // One-off diagnostic (Dan, 2026): dump the ACTUAL final returned vector
  // (pvModel.param right after the SECOND/refined-pass
  // downhillSimplexRun) -- not inferred from the per-evaluation TRACE/
  // FULLVEC history, which doesn't distinguish "the vertex the simplex
  // settled on" from "whatever the last/lowest-chi2 individual evaluation
  // happened to probe". Mirrors GalaxyFit.f's matching FINALVEC dump.
  if (TRACE_DEBUG) {
    const hexBuf = new ArrayBuffer(4);
    const hexF32v = new Float32Array(hexBuf);
    const hexU32v = new Uint32Array(hexBuf);
    const lines = ['FINALVEC'];
    for (let i = 0; i < pvModel.nParams; i++) {
      hexF32v[0] = pvModel.param[i];
      lines.push(`FINALVECPARAM ${i + 1} ${hexU32v[0].toString(16).toUpperCase().padStart(8, '0')} ${pvModel.param[i]}`);
    }
    // console.error from inside a worker_threads Worker races
    // worker.terminate() (called by the pool right after it receives this
    // realization's result message) -- the buffered stderr pipe relay to
    // the parent can be cut off before it flushes, silently dropping this
    // whole block for every realization but the last one the pool ever
    // runs. Synchronous file write sidesteps that race entirely.
    if (process.env.TRACE_DEBUG_FINALVEC_FILE) {
      require('fs').appendFileSync(process.env.TRACE_DEBUG_FINALVEC_FILE, lines.join('\n') + '\n');
    } else {
      console.error(lines.join('\n'));
    }
  }

  // Final best-fit is in pvModel (updated by downhillSimplexRun)
  return { pvModel, noConvergence };
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

  const { iter, noConvergence } = amoeba(p, y, n, 1e-6,
    params => f32(f32(params[0]-3)**2 + f32(params[1]-2)**2));
  console.log('iter:', iter, 'converged:', !noConvergence);

  console.log('minimum at:', Array.from(p[0]).map(v=>v.toFixed(4)), '(expect [3,2])');
  console.log('f(min):', y[0].toFixed(8), '(expect ~0)');
  console.log('converged:', Math.abs(p[0][0]-3) < 0.01 && Math.abs(p[0][1]-2) < 0.01 ? 'OK' : 'FAIL');

  // Test makeParamGuessArray
  console.log('\n=== makeParamGuessArray ===');
  const { makeRng } = require('./random.js');
  const { ParameterVector, allocateParamVector } = require('./ParameterVector.js');
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

});
