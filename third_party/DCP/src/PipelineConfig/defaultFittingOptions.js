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
    cmode: 0,
    // NOT 400 (SingleGalaxyTestFittingOptions_Base.txt line 25's literal
    // value) -- that transcription doesn't match the real Fortran run this
    // pipeline is validated against. Confirmed by reproducing
    // WALLABY_J103538-484832_AvgModel_v1.txt (real Fortran output) bit-exactly
    // through buildInitialFitPayload() -> runInitialFit() only with 500;
    // every other option field already matched. 400 silently produced a
    // different (still-converged, no error) local minimum -- chi2 179101.578
    // / Inc 84.53 deg instead of the correct chi2 178356.484 / Inc 73.77 deg.
    cloudBaseSurfDens: 500.0,
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
