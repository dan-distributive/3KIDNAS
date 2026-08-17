'use strict';
// =============================================================================
// buildFitPayloads.js
//
// Pure, dependency-free (no wasm, no cfitsio, no fs) equivalent of
// FitDriverScripts/RunInitialFitDCP.BuildInitialFitPayload and
// RunBootstrapsDCP.BuildRealizationPayload -- builds the exact same payload
// shapes bootstrap-realization-launcher.js's runInitialFit/
// runBootstrapRealization expect, from plain JS values instead of a Python
// process reading Fortran-format text templates off disk.
//
// Isomorphic on purpose: every function here takes plain data in (raw FITS
// bytes as Uint8Array, numbers, the fitting-options object from
// defaultFittingOptions.js) and returns plain data out (a JSON-serializable
// payload object). No `require`, no `fs`, no wasm module of its own -- the
// caller (a Node CLI script, or galaxy-fit.html in a browser) is responsible
// for getting bytes off disk/a file input and for eventually handing the
// returned payload to `compute.for`. This is what lets the SAME file be
// `require()`'d by Node and loaded via a plain <script> tag in a browser.
//
// Field-for-field provenance is documented inline against the two Python
// files it replaces -- every field here was read directly off
// RunInitialFitDCP.py / RunBootstrapsDCP.py / bootstrap-realization-
// launcher.js's runInitialFit report shape, not assumed.
// =============================================================================

// atob/btoa, not Buffer -- Buffer isn't guaranteed to exist in a browser or
// in a DCP worker sandbox (see bootstrap-realization-launcher.js's own
// b64ToBytes/bytesToB64 for the identical reasoning). Node has had global
// atob/btoa since v16, so this needs no environment branch.
function bytesToB64(bytes) {
  let bin = '';
  const chunkSize = 8192;
  for (let i = 0; i < bytes.length; i += chunkSize) {
    bin += String.fromCharCode.apply(null, bytes.subarray(i, i + chunkSize));
  }
  return btoa(bin);
}

// -----------------------------------------------------------------------
// buildInitialFitPayload
//
// Mirrors RunInitialFitDCP.BuildInitialFitPayload exactly. Unlike the
// Python version, there's no Fortran-format options file to parse --
// `options` is already the plain object defaultFittingOptions.js exports
// (optionally edited by the advanced-settings panel), so this function is
// pure field mapping, no parsing.
//
// @param {Uint8Array} cubeBytes  Raw (un-brightness-converted) observed cube FITS bytes.
// @param {Uint8Array} maskBytes  Raw pre-existing catalogue mask FITS bytes.
// @param {number} paEstDeg
// @param {number} incEstDeg
// @param {object} options  A defaultFittingOptions()-shaped object.
// @returns {object} payload, ready for compute.for([0], runInitialFit, [payload]).
// -----------------------------------------------------------------------
// buildGeometryEstimatePayload
//
// Payload for the estimateGeometry work function -- a one-off SoFiA run on
// the raw observed cube (no mask, no resample, no fit) that produces the
// same PA/Inc starting-point estimate a survey's own source-finding
// catalogue would normally supply ahead of time (see estimateGeometry's own
// header comment in bootstrap-realization-launcher.js). Only needs the cube
// bytes and the SoFiA .par template -- no options object, since this step
// doesn't touch any fitting config at all.
// -----------------------------------------------------------------------
function buildGeometryEstimatePayload({ cubeBytes, sofiaParTemplateText }) {
  if (!(cubeBytes instanceof Uint8Array)) throw new Error('buildGeometryEstimatePayload: cubeBytes must be a Uint8Array');
  if (!sofiaParTemplateText) throw new Error('buildGeometryEstimatePayload: sofiaParTemplateText is required');

  return {
    observedCubeRawFitsB64: bytesToB64(cubeBytes),
    sofiaParTemplate: sofiaParTemplateText,
  };
}

// -----------------------------------------------------------------------
function buildInitialFitPayload({ cubeBytes, maskBytes, paEstDeg, incEstDeg, options }) {
  if (!(cubeBytes instanceof Uint8Array)) throw new Error('buildInitialFitPayload: cubeBytes must be a Uint8Array');
  if (!(maskBytes instanceof Uint8Array)) throw new Error('buildInitialFitPayload: maskBytes must be a Uint8Array');
  if (!Number.isFinite(paEstDeg)) throw new Error('buildInitialFitPayload: paEstDeg must be a number');
  if (!Number.isFinite(incEstDeg)) throw new Error('buildInitialFitPayload: incEstDeg must be a number');

  return {
    observedCubeRawFitsB64: bytesToB64(cubeBytes),
    observedMaskRawFitsB64: bytesToB64(maskBytes),
    paEstDeg,
    incEstDeg,
    centerSource: options.centerSource,
    nRingsPerBeam: options.nRingsPerBeam,
    nTargRings: options.nTargRings,
    radGridArcsec: options.radGridArcsec,
    sdSwitch: options.sdSwitch,
    constParams: options.constParams,
    fixedParams: options.fixedParams,
    vRotLims: options.vRotLims,
    sizeLims: options.sizeLims,
    noiseSigmaLim: options.noiseSigmaLim,
    sigmaLengths: options.sigmaLengths,
    cmode: options.cmode,
    cloudBaseSurfDens: options.cloudBaseSurfDens,
    fitIdum: options.fitIdum,
    ftol: options.ftol,
    likelihoodSwitch: options.likelihoodSwitch,
  };
}

// -----------------------------------------------------------------------
// computeBsCent
//
// Mirrors RunBootstrapsDCP.ComputeBsCent exactly (same formula, same
// units-in/units-out), reading its three raw spectral-header inputs off
// initialFitResult.rawSpectralHeader (see that field's own comment in
// bootstrap-realization-launcher.js's runInitialFit for why it's surfaced
// there rather than re-read from FITS here) instead of GalaxyDict
// ['CubeHeader']. Geometry (XCENTER/YCENTER/POSITIONANGLE/INCLINATION/VSYS)
// comes straight off initialFitResult's own report fields -- identical
// values to what Python's GalaxyDict['BestFitModel'] carries for these,
// since AdaptInitialFitResult passes them through unchanged (only RA/DEC
// go through a WCS transform, which this doesn't need).
// -----------------------------------------------------------------------
function computeBsCent(initialFitResult) {
  const { crval3, crpix3, cdelt3 } = initialFitResult.rawSpectralHeader;
  const dV = cdelt3 / 1000;
  const deltaV = initialFitResult.VSYS[0] - crval3 / 1000;
  const centV = deltaV / dV + crpix3;
  return {
    centX: initialFitResult.XCENTER[0],
    centY: initialFitResult.YCENTER[0],
    centV,
    pa: (initialFitResult.POSITIONANGLE[0] + 90) * Math.PI / 180,
    inc: initialFitResult.INCLINATION[0] * Math.PI / 180,
  };
}

// -----------------------------------------------------------------------
// computeSDLims
//
// Mirrors RunBootstrapsDCP.ComputeSDLims exactly (plain double precision,
// same as Python's own numpy float64 -- not f32-rounded, matching that
// function's own header comment that this only ever feeds optimizer
// parameter BOUNDS, not the objective function itself).
// -----------------------------------------------------------------------
function computeSDLims(pixelSizeX, pixelSizeY, channelSize) {
  const JyAS_To_MsolPC = 1.24756e20 / (6.0574e5 * 1.823e18 * (2.0 * Math.PI / Math.log(256.0)));
  const sdConv1 = 1.0 / Math.abs(pixelSizeX * pixelSizeY);
  const LIN_SD_LIMS_MSOL = [0.1, 50.0];
  const linSDLims = LIN_SD_LIMS_MSOL.map((vMsol) => {
    const sdJyAS2 = vMsol * JyAS_To_MsolPC;
    const sd = sdJyAS2 / sdConv1;
    return sd / Math.abs(channelSize);
  });
  const logSDLims = linSDLims.map((sd) => Math.log10(sd));
  return { linSDLims, logSDLims };
}

// -----------------------------------------------------------------------
// buildBootstrapPayload
//
// Mirrors RunBootstrapsDCP.BuildRealizationPayload, but reads everything it
// needs from initialFitResult (the plain object runInitialFit's report
// already is, in memory) instead of round-tripping diskfit_fixture.json /
// model_cube_bestfit.json through disk -- the fixture IS
// initialFitResult.fixtureJson already; the model cube IS
// initialFitResult.modelCubeFitsB64 already (already base64, no
// decode/re-encode needed).
//
// @param {object} initialFitResult  The report object runInitialFit() returned.
// @param {Uint8Array} cubeBytes  The SAME raw observed-cube bytes passed to buildInitialFitPayload.
// @param {string} sofiaParTemplateText  Contents of sofia-template-par-file.par.
// @param {number} [bootstrapSeed]  0 (default) = today's time-based-seed behaviour.
// @param {object} options  Same defaultFittingOptions()-shaped object used for the initial fit.
// @returns {object} payload, ready for compute.for(inputSet, runBootstrapRealization, [payload]).
// -----------------------------------------------------------------------
function buildBootstrapPayload({ initialFitResult, cubeBytes, sofiaParTemplateText, bootstrapSeed, options }) {
  if (!initialFitResult || !initialFitResult.fixtureJson) {
    throw new Error('buildBootstrapPayload: initialFitResult.fixtureJson is required -- run the initial fit first');
  }
  if (!(cubeBytes instanceof Uint8Array)) throw new Error('buildBootstrapPayload: cubeBytes must be a Uint8Array');
  if (!sofiaParTemplateText) throw new Error('buildBootstrapPayload: sofiaParTemplateText is required');

  const fixture = initialFitResult.fixtureJson;
  // Pinned ring grid: the bootstrap loop always uses the initial fit's own
  // best-fit radii (RunBootstrapsDCP.py: best_fit_radii = list(BestFitModel['R'])),
  // NOT whatever the generic options object says -- initialFitResult.R is
  // that exact array (arcsec, same units Python's BestFitModel['R'] carries).
  const radGridArcsec = Array.from(initialFitResult.R);

  const { validFlux, flattendValidIndices, ...observedDCHeader } = fixture.observedDC;
  void validFlux; void flattendValidIndices; // stripped intentionally -- see RunBootstrapsDCP.py's own comment on why (unread downstream, ~9.4MB -> ~1.4MB payload)

  const { linSDLims, logSDLims } = computeSDLims(
    fixture.observedDC.pixelSizeX, fixture.observedDC.pixelSizeY, fixture.observedDC.channelSize
  );

  return {
    observedDC: observedDCHeader,
    observedCubeRawFitsB64: bytesToB64(cubeBytes),
    modelCubeRawFitsB64: initialFitResult.modelCubeFitsB64, // already base64 -- no decode/re-encode
    beamMajorAxis: fixture.observedBeam.beamMajorAxis,
    beamMinorAxis: fixture.observedBeam.beamMinorAxis,
    bsCent: computeBsCent(initialFitResult),
    velBlockSize: 1,
    bootstrapSeed: bootstrapSeed || 0,
    sofiaParTemplate: sofiaParTemplateText,

    observedBeam: fixture.observedBeam,
    centerSource: options.centerSource,
    nRingsPerBeam: fixture.nRingsPerBeam,
    nTargRings: radGridArcsec.length,
    radGridArcsec,
    sdSwitch: fixture.linearLogSDSwitch,
    constParams: fixture.constParams,
    fixedParams: fixture.fixedParams,
    linSDLims,
    logSDLims,
    vRotLims: options.vRotLims,
    sizeLims: options.sizeLims,
    noiseSigmaLim: options.noiseSigmaLim,
    cmode: fixture.cmode,
    cloudBaseSurfDens: fixture.cloudBaseSurfDens,
    fitIdum: fixture.idum,
    ftol: fixture.ftol,
    likelihoodSwitch: fixture.likelihoodSwitch,
  };
}

const buildFitPayloadsApi = { buildGeometryEstimatePayload, buildInitialFitPayload, buildBootstrapPayload, computeBsCent, computeSDLims, bytesToB64 };
if (typeof module !== 'undefined' && module.exports) {
  module.exports = buildFitPayloadsApi;
} else if (typeof globalThis !== 'undefined') {
  globalThis.BuildFitPayloads = buildFitPayloadsApi;
}
