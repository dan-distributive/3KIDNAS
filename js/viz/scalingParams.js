'use strict';

// =============================================================================
// scalingParams.js
// Port of ExtractScalingParams.py (RHI/VHI scaling relations) and
// GeometryCorrection.py:GetGlobalPositionAngle (WCS-corrected global PA).
// See ~/.claude/plans/breezy-launching-nova.md Phase 4.
//
// Runs entirely client-side, on data index.html already holds after a
// bootstrap run (initialFitResult as the Python's BestFitModel, and the
// bootstrapResults array as BootstrapModels) -- no new server-side/DCP-
// dispatched computation needed for this phase.
//
// Deliberate simplification vs. the Python source: CalcRHI/CalcVHI there
// first check GalaxyDict['ExtendedSDProfile']['ProfileAcceptFlag'] -- a
// SoFiA-derived, wider-than-the-fitted-rings surface density profile this
// JS port doesn't have an equivalent of. Here the fit's own R_SD/
// SURFDENS_FACEON ring array is used directly as "the profile" and that
// gate is treated as always-true; RHI/VHI can end up less well-constrained
// for a galaxy whose HI extends past the last fitted ring than the Python
// pipeline's version would be, but the extraction algorithm itself
// (GetSD_Intecept/FindProfileIntersection/GetProfilePoint, and the
// bootstrap sqrt(std^2+diff^2) uncertainty combination) is ported exactly.
// =============================================================================

// FindProfileIntersection -- brentq(interp1d(X,YNorm), X[indx], X[indx+1])
// on a piecewise-LINEAR interpolant restricted to one bracket IS just the
// linear-interpolation root between those two points; no iterative solve
// needed for an exact match.
function findProfileIntersection(X, Y, lim) {
  const yNorm = Y.map((y) => y - lim);
  let idx = -1;
  for (let i = 0; i < yNorm.length - 1; i++) {
    if (Math.sign(yNorm[i]) !== Math.sign(yNorm[i + 1])) idx = i; // last (outermost) crossing wins
  }
  if (idx === -1) return { xInt: NaN, indx: -1 };
  const x0 = X[idx], x1 = X[idx + 1], y0 = yNorm[idx], y1 = yNorm[idx + 1];
  const xInt = x0 + (0 - y0) * ((x1 - x0) / (y1 - y0));
  return { xInt, indx: idx };
}

function getSDIntercept(R, SD, sdLim) {
  const finite = SD.filter(Number.isFinite);
  const minSD = Math.min(...finite), maxSD = Math.max(...finite);
  if (minSD > sdLim) {
    return { rhi: R[R.length - 1], found: false };
  }
  if (maxSD < sdLim) {
    let maxIdx = 0;
    for (let i = 1; i < SD.length; i++) if (SD[i] > SD[maxIdx]) maxIdx = i;
    return { rhi: R[maxIdx], found: false };
  }
  const { xInt, indx } = findProfileIntersection(R, SD, sdLim);
  if (indx === -1) return { rhi: R[R.length - 1], found: false };
  return { rhi: xInt, found: true };
}

function extractRHI(model, sdLim) {
  if (model.FITAchieved === false) return { rhi: NaN, found: false };
  const R = model.R_SD, SD = model.SURFDENS_FACEON;
  if (!R || R.length <= 1) return { rhi: R ? R[0] : NaN, found: false };
  return getSDIntercept(R, SD, sdLim);
}

// GetProfilePoint -- linear interpolation of Y(X) at X=xTarg, clamped-bracket
// at the array ends (matches Python's own out-of-range i=0 / i=len-2 clamp,
// not an extrapolation guard bolted on separately).
function getProfilePoint(X, Y, xTarg) {
  let i = 0;
  if (xTarg <= X[0]) {
    i = 0;
  } else if (xTarg > X[X.length - 1]) {
    i = X.length - 2;
  } else {
    for (let j = 0; j < X.length - 1; j++) {
      if (xTarg > X[j] && xTarg <= X[j + 1]) { i = j; break; }
    }
  }
  const x1 = X[i], x2 = X[i + 1], y1 = Y[i], y2 = Y[i + 1];
  const m = (y2 - y1) / (x2 - x1);
  return y1 + m * (xTarg - x1);
}

function getVHIFromProf(R, VProf, rhi) {
  if (!VProf || VProf.length <= 1) return VProf ? VProf[0] : NaN;
  return getProfilePoint(R, VProf, rhi);
}

function meanStd(arr) {
  const mean = arr.reduce((a, b) => a + b, 0) / arr.length;
  const variance = arr.reduce((a, b) => a + (b - mean) ** 2, 0) / arr.length; // ddof=0, matches np.std default
  return { mean, std: Math.sqrt(variance) };
}

// ---------------------------------------------------------------------------
// computeScalingParams(bestFit, bootstrapResults, sdLim=1.0)
// Returns RHI/VHI (best-fit value + bootstrap uncertainty) and RHI in kpc
// (H0=70 km/s/Mpc flat Hubble-flow distance, matching ExtractScalingParams.
// py's DistEst -- no peculiar-velocity or redshift correction, same as the
// Python source).
// ---------------------------------------------------------------------------
function computeScalingParams(bestFit, bootstrapResults, sdLim = 1.0) {
  const bestRHI = extractRHI(bestFit, sdLim);
  const bsRHI = bootstrapResults.map((r) => extractRHI(r, sdLim));

  let RHI_flag, RHI_err;
  if (!bestRHI.found) {
    RHI_flag = -1; RHI_err = NaN;
  } else {
    const nRequired = Math.trunc(0.6 * bootstrapResults.length);
    const found = bsRHI.filter((r) => r.found).map((r) => r.rhi);
    if (found.length < nRequired) {
      RHI_flag = -1; RHI_err = NaN;
    } else {
      const { mean, std } = meanStd(found);
      const diff = bestRHI.rhi - mean;
      RHI_err = Math.sqrt(std * std + diff * diff);
      RHI_flag = 0;
    }
  }

  const H0 = 70; // km/s/Mpc
  const distMpc = bestFit.VSYS[0] / H0;
  const ARCSEC_PER_RAD = 206265;
  const toKpc = (arcsec) => (arcsec / ARCSEC_PER_RAD) * distMpc * 1000;
  const RHI_kpc = RHI_flag === -1 ? NaN : toKpc(bestRHI.rhi);
  const RHI_kpc_err = RHI_flag === -1 ? NaN : toKpc(RHI_err);

  let VHI_flag, VHI_err, VHI = NaN;
  if (RHI_flag === -1) {
    VHI_flag = -1; VHI_err = NaN;
  } else {
    VHI = getVHIFromProf(bestFit.R, bestFit.VROT, bestRHI.rhi);
    const bsVHI = [];
    for (let i = 0; i < bootstrapResults.length; i++) {
      if (bsRHI[i].found) bsVHI.push(getVHIFromProf(bootstrapResults[i].R, bootstrapResults[i].VROT, bsRHI[i].rhi));
    }
    if (bsVHI.length === 0) {
      VHI_flag = -1; VHI_err = NaN;
    } else {
      const { mean, std } = meanStd(bsVHI);
      const diff = VHI - mean;
      VHI_err = Math.sqrt(std * std + diff * diff);
      VHI_flag = 0;
    }
  }

  return {
    RHI: bestRHI.rhi, RHI_err, RHI_flag,
    RHI_kpc, RHI_kpc_err,
    VHI, VHI_err, VHI_flag,
    distMpc,
  };
}

// ---------------------------------------------------------------------------
// GetGlobalPositionAngle -- WCS-corrected PA (pixel-frame PA -> sky-frame
// PA). Uses a linear tangent-plane pixel->sky approximation rather than
// astropy's full SIN-projection WCS: exact for this data (no CROTA2/PC/CD
// rotation terms in these cubes' headers -- confirmed directly), and the
// SIN projection's own spherical curvature is negligible at the ~4 arcmin
// field of view these cubes cover (deliberately flagged, not silently
// dropped -- see the plan's Phase 4 section for why a full WCS library
// wasn't pulled in for this one calculation).
// ---------------------------------------------------------------------------
function pixelToSkyDeg(xPix, yPix, wcsHeader) {
  const crval1Deg = wcsHeader.refValArcsec[0] / 3600;
  const crval2Deg = wcsHeader.refValArcsec[1] / 3600;
  const cdelt1Deg = wcsHeader.pixelSizeArcsec[0] / 3600;
  const cdelt2Deg = wcsHeader.pixelSizeArcsec[1] / 3600;
  const dx = xPix - wcsHeader.refPix0Indexed[0];
  const dy = yPix - wcsHeader.refPix0Indexed[1];
  const cosDec = Math.cos((crval2Deg * Math.PI) / 180);
  return { ra: crval1Deg + (dx * cdelt1Deg) / cosDec, dec: crval2Deg + dy * cdelt2Deg };
}

function computeGlobalPA(bestFit, wcsHeader, beamMajorAxisPix) {
  const xPix = bestFit.XCENTER[0], yPix = bestFit.YCENTER[0];
  const cent = pixelToSkyDeg(xPix, yPix, wcsHeader);
  const angURad = ((bestFit.POSITIONANGLE[0] + 90) * Math.PI) / 180;
  const xNew = xPix + beamMajorAxisPix * Math.cos(angURad);
  const yNew = yPix + beamMajorAxisPix * Math.sin(angURad);
  const newSky = pixelToSkyDeg(xNew, yNew, wcsHeader);
  const deltRA = newSky.ra - cent.ra;
  const deltDEC = newSky.dec - cent.dec;
  let newPA = (Math.atan2(deltDEC, -deltRA) * 180) / Math.PI - 90;
  if (newPA < 0) newPA += 360;
  else if (newPA > 360) newPA -= 360;
  return newPA;
}

// Per-model RHI/VHI -- exposed so a per-bootstrap-realization scatter plot
// (cornerPlot.js) can show each realization's OWN RHI/VHI, the same values
// computeScalingParams' bootstrap loop above computes internally, without
// duplicating the extraction logic.
function extractVHIForModel(model, sdLim = 1.0) {
  const rhiRes = extractRHI(model, sdLim);
  if (!rhiRes.found) return { rhi: rhiRes.rhi, vhi: NaN, found: false };
  return { rhi: rhiRes.rhi, vhi: getVHIFromProf(model.R, model.VROT, rhiRes.rhi), found: true };
}

if (typeof module !== 'undefined') {
  module.exports = {
    computeScalingParams, computeGlobalPA, pixelToSkyDeg,
    getSDIntercept, getProfilePoint, extractRHI, extractVHIForModel,
  };
}
if (typeof globalThis !== 'undefined') {
  globalThis.ScalingParams = {
    computeScalingParams, computeGlobalPA, pixelToSkyDeg,
    getSDIntercept, getProfilePoint, extractRHI, extractVHIForModel,
  };
}
