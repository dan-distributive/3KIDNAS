'use strict';
// =============================================================================
// defaultFittingOptions.js
//
// Plain-data equivalent of the ~15 fixed pipeline-configuration values that,
// on the Python/Fortran side, live inside two legacy Fortran-format text
// templates and get extracted by replaying FittingOptionsInputs.f's exact
// line-by-line read order (FitDriverScripts/RunInitialFitDCP.py:
// ParseFittingOptionsForInitialFit / RunBootstrapsDCP.py:
// ParseFittingOptionsExtras). These values are NOT derived from any cube --
// they're fixed pipeline defaults, identical for every galaxy run today.
//
// Transcribed directly from both template files (read 2026-08-16), not
// assumed:
//   Inputs/SingleGalaxyTestFittingOptions_Base.txt (lines 1-58)
//   Inputs/SingleFitInput_Base.in (line 12: idum)
//
// A from-scratch JS/browser pipeline has no Fortran-format file to parse, so
// this object replaces that parsing step entirely -- it IS the config, not a
// parser for one.
//
// constParams/fixedParams order (13 ring properties, matches
// FittingOptionsInputs.f / DumpFittingFixture exactly): XCenter, YCenter,
// Inclination, PositionAngle, SystemicVelocity, RotationVelocity,
// RadialVelocity, VelocityDispersion, VerticalVelocity, dvdz,
// SurfaceDensity, VerticalHeight, VerticalGradientStart.
//
// Loadable as CommonJS (Node CLI) or a plain <script> global (browser) --
// see the export guard at the bottom. No `require`/`fs` anywhere in this
// file, so it's safe in both environments unconditionally.
// =============================================================================

const RING_PARAM_LABELS = [
  'X Center', 'Y Center', 'Inclination', 'Position Angle', 'Systemic Velocity',
  'Rotation Velocity', 'Radial Velocity', 'Velocity Dispersion',
  'Vertical Velocity', 'dvdz', 'Surface Density', 'Vertical Height',
  'Vertical Gradient Start',
];

// [constant-across-rings, fixed-at-initial-value] per ring property,
// transcribed verbatim from SingleGalaxyTestFittingOptions_Base.txt lines
// 34-58. Only Rotation Velocity and Surface Density are actually free
// per-ring fit parameters under these defaults -- everything else fits as a
// single global value (constant) or never varies from its initial guess
// (fixed).
const DEFAULT_CONST_PARAMS = [true, true, true, true, true, false, true, true, true, true, false, true, true];
const DEFAULT_FIXED_PARAMS = [false, false, false, false, false, false, true, true, true, true, false, true, true];

function defaultFittingOptions() {
  // Returns a fresh deep copy every call -- callers (the advanced-options UI,
  // a payload builder) mutate their own copy freely without corrupting the
  // shared default.
  return {
    likelihoodSwitch: 1,        // 1 = chi^2
    centerSource: 0,            // 0 = estimate center (no live SoFiA needed for the initial fit)
    sdSwitch: 0,                 // 0 = linear surface density, 1 = logarithmic
    // cmode=0, cloudBaseSurfDens=400 -- the pre-upstream-commit-76ade48
    // pairing, matching Inputs/SingleGalaxyTestFittingOptions_Base.txt
    // exactly (previously 500 here vs 400 there -- an old, undocumented
    // asymmetric fudge that happened to compensate for some other,
    // unrelated discrepancy. Removed: confirmed it was silently causing
    // Fortran and the JS port to generate different particle *counts*
    // per ring -- hence different ran2()/gasdev() draw counts -- from the
    // very first optimizer evaluation, every run. With matched cdens on
    // both platforms instead, plus SingleRingGeneration.js's
    // roundForParticleStability rounding, Fortran and the JS port ran in
    // bit-exact idum lockstep for an entire fit to convergence, chi2
    // agreeing to ~8e-8 relative).
    //
    // Deliberately NOT yet switched to upstream commit 76ade48's new
    // formula (Noise-normalization + AvgChannelsPerPix multiplier,
    // cmode=1/cloudBaseSurfDens=10) -- that cdens=10 comes from Nathan
    // Deg's own attached example config, but his email text says the
    // actual new default is cdens=100; 10 looks like a typo in the
    // example rather than the intended value. Revisit once that's
    // confirmed, and re-verify Fortran/JS agreement at whatever the real
    // default turns out to be before switching this back on -- see
    // SingleRingGeneration.f/.js's own comments at the same spot.
    cmode: 0,
    cloudBaseSurfDens: 400.0,
    sigmaLengths: 2.5,          // convolution kernel half-width, in beam sigmas
    noiseSigmaLim: 1.0,
    nTargRings: -1,              // -1 = derive ring count from nRingsPerBeam, not a fixed count
    radGridArcsec: null,         // only meaningful when nTargRings > 0 (a manual ring-grid override)
    nRingsPerBeam: 2,
    constParams: DEFAULT_CONST_PARAMS.slice(),
    fixedParams: DEFAULT_FIXED_PARAMS.slice(),
    vRotLims: [0.0, 400.0],       // PipelineGlobals.f hardcoded constant, not file-read
    sizeLims: [2.0, 10.0],        // PipelineGlobals.f hardcoded constant, not file-read
    fitIdum: -3,                  // static per-galaxy optimizer seed (SingleFitInput_Base.in line 12)
    ftol: 0.005,                  // galaxyFit_Simple's own pass-1 ftol -- hardcoded in Python too, never file-read
  };
}

// Static asset, sibling to this file's ultimate consumers -- see
// buildBootstrapPayload.js for how it's loaded (fetch() in the browser,
// fs.readFileSync in Node). Only the bootstrap loop needs it (the initial
// fit runs no live SoFiA at all).
const SOFIA_PAR_TEMPLATE_RELATIVE_PATH = 'sofia-template-par-file.par';

const defaultFittingOptionsApi = { defaultFittingOptions, RING_PARAM_LABELS, SOFIA_PAR_TEMPLATE_RELATIVE_PATH };
if (typeof module !== 'undefined' && module.exports) {
  module.exports = defaultFittingOptionsApi;
} else if (typeof globalThis !== 'undefined') {
  globalThis.DefaultFittingOptions = defaultFittingOptionsApi;
}
