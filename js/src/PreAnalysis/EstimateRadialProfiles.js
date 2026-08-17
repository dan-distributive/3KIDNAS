'use strict';

// =============================================================================
// EstimateRadialProfiles.js
// Port of src/PreAnalysis/EstimateRadialProfiles.f (EstimateProfileMod)
//
// DATA REPRESENTATION
// --------------------
// A "RadialProfile"-shaped Fortran array (0:2, 0:nRings-1) is represented
// here as an Array of length nRings, each element a Float32Array(3) =
// [radius, surfDens, vRot] -- ring i is profile[i]. A "wide" profile
// (0:2, -nRings+1:nRings), used only for SliceMap's output (EstimatedProfile
// in the Fortran source), is returned as { nRings, get(idx) } where idx may
// be negative; get(idx) returns the Float32Array(3) for that ring index.
//
// LIKELY FORTRAN BUG -- "all rings modelable" overflows the wide profile's
// allocated bounds:
//   When FittingOptions%nTargRings==-1 (auto mode), nRings is reassigned to
//   `sum(ModelerableRingSwitch)+1` (an intentional "+1 extra ring trial" --
//   see the source comment right above it). ModelerableRingSwitch has
//   length nRingsMax, so if EVERY candidate ring is modelable, this
//   reassigned nRings becomes nRingsMax+1 -- one past what EstimatedProfile
//   (the wide sliceMap output, allocated at -nRingsMax+1:nRingsMax) actually
//   holds. The later leading/trailing-velocity loops then index
//   EstimatedProfile one step outside its allocated range. Reproduced here
//   directly (hit while exercising this port against a small, clean
//   synthetic disc where every ring happened to pass modelability -- see
//   InitialAnalysis.js's self-test history). Not specially handled -- this
//   throws (Array access returns undefined, TypeError on the next property
//   read) rather than silently clamping, so it surfaces the same way an
//   out-of-bounds access would in Fortran (compiled with -fbounds-check per
//   src/makeflags) instead of being masked. Flagged for Nathan alongside the
//   other issues in this file; not yet confirmed against a real Fortran run
//   with real (always-noisier) galaxy data, where 100% ring modelability
//   may simply never occur in practice.
//
// FIXED FORTRAN BUG (2026, flagged by Dan, reported to Nathan) --
// FillInMisingRings used to unconditionally abort when ring 0 needed
// filling in, even after successfully finding a fallback value:
//   The `if(i.eq.0) then <fill-in loop> ; print "no ring..."; stop; else ...`
//   structure had NO conditional guard between the fill-in loop and the
//   stop -- confirmed by matching every if/do/enddo/else/endif in the source
//   (2 nested `if`s, 1 `else`, 2 `endif`s: the print+stop sat inside the
//   i==0 branch, BEFORE its `else`, with nothing skipping them even when the
//   loop's `exit` fired). So whenever ring 0 was non-modelable, the real
//   Fortran binary ALWAYS stopped the program, whether or not a valid
//   fallback ring was found -- a missing "did we find one?" guard around
//   the abort. Fixed on both sides -- native Fortran now only aborts if no
//   fallback ring was found (src/PreAnalysis/EstimateRadialProfiles.f), and
//   this port matches.
//
// DEAD/DISCARDED VALUES (not ported -- confirmed zero effect on output):
//   - EstimateProfiles: the very first "Leading Vel Check" print block
//     (source lines 178-184) prints a debug expression and an UNINITIALIZED
//     LeadingVelProfile (its real assignment is commented out) -- pure
//     stack-garbage debug noise with no computational effect. The REAL
//     Leading/TrailingVelProfile computation (source lines 186-207, used for
//     the PA leading/trailing flip) is ported below.
//   - SliceMap: `l=nint(yPrime)` is computed and never read again. AvgX(i)
//     is accumulated and averaged solely for a debug print, never read by
//     anything that affects AvgQuantities (the actual output). Neither is
//     ported.
//   - BuildIniProfile: InclU (a would-be inclination clamp) is computed but
//     never referenced anywhere else in the subroutine -- Incl itself is
//     what's actually used in cos(Incl)/sin(Incl) below. Not ported.
//   - BeamCorrectProfile: dR is computed once and never used again. Not
//     ported.
//
// PRECISION
// ---------
// Fortran `real` -> Math.fround() per step, matching the -O0 build.
// `X**2.` (compile-time-literal real exponent) constant-folds to a plain
// multiply, as established elsewhere in this project.
// Fortran NINT (round half AWAY FROM ZERO) is NOT the same as JS's
// Math.round (round half towards +Infinity) -- SliceMap needs a real NINT,
// implemented locally as fortranNint() below.
//
// FIXED FORTRAN BUG (2026, flagged by Dan, reported to Nathan) -- every
// cos()/sin() call in this file (sliceMap's rotation, buildIniProfile's
// SD/VRot de-projection, checkProfilePointModelability's velocity-limit
// checks) used the native Math.cos()/Math.sin() intrinsic on this side and
// Fortran's native cos()/sin() intrinsic on the other -- unlike every other
// trig call in this pipeline, which already goes through the shared fdlibm
// port precisely to guarantee bit-identical results. Traced directly: a
// real 4-ULP disagreement in buildIniProfile's cos(Incl), for one specific
// bootstrap realization's inclination, propagated unchanged through the
// initial parameter vector for ~130 chained optimizer evaluations before
// tipping a near-tied Nelder-Mead simplex comparison and forking that
// realization's entire fit trajectory. Switched to fdSin/fdCos below.
// =============================================================================

const f32 = Math.fround;
const { simpleInterpolateY } = require('../StandardMath/Interpolation.js');
const { fdSin, fdCos } = require('../StandardMath/fdlibm.js');

// `process` isn't guaranteed to exist in a real DCP worker sandbox --
// confirmed directly against a real worker. See EstimateCubeNoise.js's
// identical helper for the full explanation.
function isTraceDebug() {
  return typeof process !== 'undefined' && process.env && process.env.TRACE_DEBUG === '1';
}

// ---------------------------------------------------------------------------
// fortranNint
// Fortran NINT(x): nearest integer, ties rounding AWAY from zero.
// (Math.round rounds ties towards +Infinity -- wrong for negative x.5 cases.)
// ---------------------------------------------------------------------------
function fortranNint(x) {
  return x >= 0 ? Math.floor(x + 0.5) : Math.ceil(x - 0.5);
}

// ---------------------------------------------------------------------------
// makeProfile / makeWideProfile
// Small local helpers for the two array shapes used throughout this file.
// ---------------------------------------------------------------------------
function makeProfile(nRings) {
  const arr = new Array(nRings);
  for (let i = 0; i < nRings; i++) arr[i] = new Float32Array(3);
  return arr;
}

function makeWideProfile(nRings) {
  const nCols = 2 * nRings;
  const arr = new Array(nCols);
  for (let i = 0; i < nCols; i++) arr[i] = new Float32Array(3);
  const idxOf = (idx) => idx + nRings - 1;
  return {
    nRings,
    get(idx) { return arr[idxOf(idx)]; },
    _arr: arr,
    _idxOf: idxOf,
  };
}

// ---------------------------------------------------------------------------
// sliceMap
// Fortran: SliceMap(Maps, Beam, Center, incl, Phi, Size, nRings,
//                    AvgQuantities, nRingsPerBeam)
//
// `incl` and `Size` are unused dummy parameters in the Fortran body (checked
// directly) -- dropped from this port's signature.
//
// Slices the moment maps into nRingsPerBeam-per-beam bins along the major
// axis (Phi), averaging moment-0/1 within +-0.5 beam of the slice line.
// Returns a wide profile (see module header) with columns
// [radius(arcsec-equivalent-in-beam-units), avg moment0/chanwidth, avg moment1].
// ---------------------------------------------------------------------------
function sliceMap(maps, beamMajorAxis, center, phi, nRings, nRingsPerBeam) {
  const dh = maps.dh;
  const nx = dh.nPixels[0];
  const ny = dh.nPixels[1];
  const sliceLimit = f32(0.5);
  const beamPix = f32(Math.abs(beamMajorAxis));
  const cx = center[0], cy = center[1];

  const wide = makeWideProfile(nRings);
  const nBin = new Int32Array(wide._arr.length);

  const cosNegPhi = f32(fdCos(f32(-phi)));
  const sinNegPhi = f32(fdSin(f32(-phi)));
  const xScale = f32(beamPix / f32(nRingsPerBeam));

  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      const x = f32(f32(i) - cx);
      const y = f32(f32(j) - cy);
      let xPrime = f32(f32(x * cosNegPhi) - f32(y * sinNegPhi));
      let yPrime = f32(f32(x * sinNegPhi) + f32(y * cosNegPhi));
      xPrime = f32(xPrime / xScale);
      yPrime = f32(yPrime / beamPix);
      const k = fortranNint(f32(xPrime + f32(0.5)));

      if (Math.abs(yPrime) <= sliceLimit) {
        if (k >= -nRings + 1 && k <= nRings) {
          const base = i * ny * 3 + j * 3;
          const flux0 = maps.flux[base];
          const flux1 = maps.flux[base + 1];
          if (flux0 >= f32(0.0) && !Number.isNaN(flux1)) {
            const col = wide._idxOf(k);
            wide._arr[col][1] = f32(wide._arr[col][1] + flux0);
            wide._arr[col][2] = f32(wide._arr[col][2] + flux1);
            nBin[col] += 1;
          }
        }
      }
    }
  }

  const chanAbs = f32(Math.abs(dh.channelSize));
  for (let idx = -nRings + 1; idx <= nRings; idx++) {
    const col = wide._idxOf(idx);
    const n = f32(nBin[col]);
    const cell = wide._arr[col];
    cell[0] = f32(f32(beamPix * f32(f32(idx) - f32(0.5))) / f32(nRingsPerBeam));
    cell[1] = f32(cell[1] / n);
    cell[1] = f32(cell[1] / chanAbs);
    cell[2] = f32(cell[2] / n);
  }

  return wide;
}

// ---------------------------------------------------------------------------
// buildIniProfile
// Fortran: BuildIniProfile(nRings, VSys, Incl, ModelerableRingSwitch,
//                           EstimatedProfile, RadialProfile, ProjectedProfile)
//
// Pairs each "leading" ring i (1..nRings) with its "trailing" mirror
// j=-i+1 from the wide sliceMap profile, averages them, and de-projects by
// inclination into an initial radius/surfDens/vRot guess per ring.
// ---------------------------------------------------------------------------
function buildIniProfile(nRings, vSys, incl, estimatedProfile) {
  const modelerableRingSwitch = new Int32Array(nRings);
  const radialProfile = makeProfile(nRings);
  const projectedProfile = makeProfile(nRings);

  const cosIncl = f32(fdCos(incl));
  const sinIncl = f32(fdSin(incl));

  for (let i = 1; i <= nRings; i++) {
    modelerableRingSwitch[i - 1] = 1;
    const j = -i + 1;
    const epI = estimatedProfile.get(i);
    const epJ = estimatedProfile.get(j);

    let t1 = f32(Math.abs(f32(epI[2] - vSys)));
    let t2 = f32(Math.abs(f32(epJ[2] - vSys)));
    if (epI[2] <= f32(0.0)) t1 = t2;
    if (epJ[2] <= f32(0.0)) t2 = t1;

    const pp = projectedProfile[i - 1];
    pp[2] = f32(f32(t1 + t2) / f32(2.0));
    pp[1] = f32(f32(epI[1] + epJ[1]) / f32(2.0));

    const rp = radialProfile[i - 1];
    rp[0] = f32(Math.abs(epJ[0]));
    rp[1] = f32(pp[1] * cosIncl);
    rp[2] = f32(pp[2] / sinIncl);
  }

  return { modelerableRingSwitch, radialProfile, projectedProfile };
}

// ---------------------------------------------------------------------------
// checkProfilePointModelability
// Fortran: CheckProfilePointModelability(SDLims, VLims, Noise_SDLim,
//              ProfilePt, ProjectedProfilePt, ModelSwitch, Incl, VSys, DC)
//
// Mutates profilePt (clamped to VLims / high-SDLims) in place, returns the
// (possibly-cleared-to-0) modelSwitch. The low-SDLims check is commented out
// in the Fortran source -- not ported (matches that dead branch faithfully
// by omission).
// ---------------------------------------------------------------------------
function checkProfilePointModelability(sdLims, vLims, noiseSDLim, profilePt, projectedProfilePt, modelSwitchIn, incl, vSys, channels) {
  let modelSwitch = modelSwitchIn;

  if (profilePt[1] <= f32(0.0)) modelSwitch = 0;
  if (profilePt[1] !== profilePt[1]) modelSwitch = 0; // NaN check
  if (projectedProfilePt[1] <= noiseSDLim) modelSwitch = 0;

  if (profilePt[1] >= sdLims[1]) profilePt[1] = sdLims[1];
  if (profilePt[2] <= vLims[0]) profilePt[2] = vLims[0];
  if (profilePt[2] >= vLims[1]) profilePt[2] = vLims[1];

  const sinIncl = f32(fdSin(incl));
  const vProjLow  = f32(vSys - f32(profilePt[2] * sinIncl));
  const vProjHigh = f32(vSys + f32(profilePt[2] * sinIncl));

  let minV = channels[0], maxV = channels[0];
  for (let i = 1; i < channels.length; i++) {
    if (channels[i] < minV) minV = channels[i];
    if (channels[i] > maxV) maxV = channels[i];
  }

  if (vProjLow < minV)  profilePt[2] = f32(f32(vSys - minV) / sinIncl);
  if (vProjHigh > maxV) profilePt[2] = f32(f32(maxV - vSys) / sinIncl);

  return modelSwitch;
}

// ---------------------------------------------------------------------------
// assignProfile
// Fortran: AssignProfile(nRings, nRingsMax, FittingOptions, RadialProfile,
//              TempRadialProfile, ModelerableRingSwitch,
//              RadialProfileModelSwitch, Maps)
//
// nTargRings === -1: auto mode, keep modelable rings only.
// nTargRings >= 1:   pinned-grid mode (radGridArcsec supplies the target
//                    radii, converted to pixels here), keep rings whose
//                    temp-profile radius matches the grid within 0.01 pixel.
// ---------------------------------------------------------------------------
function assignProfile(nRings, nRingsMax, nTargRings, radGridArcsec, pixelSizeX, tempRadialProfile, modelerableRingSwitch) {
  const radialProfile = makeProfile(nRings);
  const radialProfileModelSwitch = new Int32Array(nRings);

  let radGridPixels = null;
  if (nTargRings >= 1) {
    const pixSize = f32(Math.abs(pixelSizeX));
    radGridPixels = radGridArcsec.map((v) => f32(v / pixSize));
  }

  let j = 0;
  for (let i = 0; i < nRingsMax; i++) {
    if (nTargRings === -1) {
      if (modelerableRingSwitch[i] === 1) {
        radialProfile[j].set(tempRadialProfile[i]);
        radialProfileModelSwitch[j] = modelerableRingSwitch[i];
        j += 1;
      }
    } else {
      if (Math.abs(f32(tempRadialProfile[i][0] - radGridPixels[j])) < 0.01) {
        radialProfile[j].set(tempRadialProfile[i]);
        radialProfileModelSwitch[j] = modelerableRingSwitch[i];
        j += 1;
      }
    }
    if (j === nRings) break;
  }

  return { radialProfile, radialProfileModelSwitch };
}

// ---------------------------------------------------------------------------
// fillInMisingRings
// Fortran: FillInMisingRings(nRings, RadialProfile, ModelerableRingSwitch)
//
// See module header for the fixed unconditional-abort bug at ring 0.
// ---------------------------------------------------------------------------
function fillInMisingRings(nRings, radialProfile, modelerableRingSwitch) {
  for (let i = 0; i < nRings; i++) {
    if (modelerableRingSwitch[i] === 0) {
      if (i === 0) {
        let foundFillIn = false;
        for (let k = i; k < nRings; k++) {
          if (modelerableRingSwitch[k] === 1) {
            radialProfile[i][1] = radialProfile[k][1];
            radialProfile[i][2] = radialProfile[k][2];
            foundFillIn = true;
            break;
          }
        }
        // FIXED (2026, flagged by Dan, reported to Nathan): Fortran used to
        // STOP here unconditionally, even after the loop above found a
        // fallback value -- see module header note. Now only aborts if no
        // fallback ring was actually found, matching the fixed Fortran.
        if (!foundFillIn) {
          throw new Error('FillInMisingRings: no ring with an acceptable surface density found');
        }
      } else {
        for (let k = i; k >= 0; k--) {
          if (modelerableRingSwitch[k] === 1) {
            radialProfile[i][1] = radialProfile[k][1];
            radialProfile[i][2] = radialProfile[k][2];
            break;
          }
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// beamCorrectProfile
// Fortran: BeamCorrectProfile(nRings, RadialProfile, Beam, SDLims, VLims)
//
// Deconvolves each ring's radius by (progressively finer) beam fractions
// until the corrected radius-squared is non-negative, then re-interpolates
// SD/VRot at that corrected radius from the (uncorrected) profile.
//
// rIndx persists ACROSS the whole i loop (only ever increments) -- NOT reset
// per ring. Confirmed by reading the source: RIndx is initialized once
// before the loop, matching an assumption that RadialProfile(0,:) and the
// corrected radii are both non-decreasing in i, so search can resume where
// the previous ring's search left off. Faithfully preserved as loop-scoped
// state here, not re-initialized inside the per-ring iteration.
// ---------------------------------------------------------------------------
function beamCorrectProfile(nRings, radialProfile, beamMajorAxis, sdLims, vLims) {
  if (nRings < 2) return;

  const sdCorr = new Float32Array(nRings);
  const rcCorr = new Float32Array(nRings);
  let rIndx = 0;

  for (let i = 0; i < nRings; i++) {
    let k = 2;
    let rCorr;
    for (;;) {
      const bmOverK = f32(beamMajorAxis / f32(k));
      rCorr = f32(f32(radialProfile[i][0] * radialProfile[i][0]) - f32(bmOverK * bmOverK));
      if (rCorr >= f32(0.0)) break;
      k += 1;
    }
    rCorr = f32(Math.sqrt(rCorr));

    while (rCorr > radialProfile[rIndx][0]) {
      rIndx += 1;
    }
    if (rIndx === nRings - 1) rIndx -= 1;

    const sdL = [radialProfile[rIndx][0], radialProfile[rIndx][1]];
    const sdH = [radialProfile[rIndx + 1][0], radialProfile[rIndx + 1][1]];
    sdCorr[i] = simpleInterpolateY(sdL, sdH, rCorr);

    const rcL = [radialProfile[rIndx][0], radialProfile[rIndx][2]];
    const rcH = [radialProfile[rIndx + 1][0], radialProfile[rIndx + 1][2]];
    rcCorr[i] = simpleInterpolateY(rcL, rcH, rCorr);

    if (rcCorr[i] > vLims[1]) rcCorr[i] = vLims[1];
    else if (rcCorr[i] < vLims[0]) rcCorr[i] = vLims[0];

    if (sdCorr[i] > sdLims[1]) sdCorr[i] = sdLims[1];
    else if (sdCorr[i] < sdLims[0]) sdCorr[i] = sdLims[0];
  }

  for (let i = 0; i < nRings; i++) {
    radialProfile[i][1] = sdCorr[i];
    radialProfile[i][2] = rcCorr[i];
  }
}

// ---------------------------------------------------------------------------
// estimateProfiles
// Fortran: EstimateProfiles(Maps, Beam, Center, incl, PA, Size,
//              EstimatedProfile, RadialProfile, VSys, FittingOptions, nRings,
//              SDLims, VLims, NoiseFrac, DC)
//
// Main driver. opts: { nRingsPerBeam, nTargRings, radGridArcsec (only used
// when nTargRings>=1) }. Returns { radialProfile, radialProfileModelSwitch,
// paOut, nRings } -- nRings is the FINAL (post-reassignment) ring count,
// which can be much smaller than the geometry-derived nRingsMax used
// internally for slicing.
// ---------------------------------------------------------------------------
function estimateProfiles(maps, beamMajorAxis, center, pa, vSys, opts, sdLims, vLims, noiseFrac, dc) {
  if (isTraceDebug()) {
    console.error('TRACE Estimating shape check', f32(pa * f32(180.0 / Math.PI)), f32(opts.incl * f32(180.0 / Math.PI)), vSys, center[0], center[1]);
  }
  const dh = maps.dh;
  const nxf = f32(dh.nPixels[0]);
  const nyf = f32(dh.nPixels[1]);

  let rTest = f32(f32(nxf * nxf) + f32(nyf * nyf));
  rTest = f32(Math.sqrt(rTest));
  rTest = f32(rTest / f32(2.0));
  rTest = f32(rTest / f32(Math.abs(beamMajorAxis)));
  let nRings = Math.trunc(f32(rTest * f32(opts.nRingsPerBeam))) + 1;
  const nRingsMax = nRings;

  const estimatedProfile = sliceMap(maps, beamMajorAxis, center, pa, nRings, opts.nRingsPerBeam);
  const { modelerableRingSwitch, radialProfile: tempRadialProfile, projectedProfile } =
    buildIniProfile(nRings, vSys, opts.incl, estimatedProfile);

  const noiseSDLim = f32(f32(noiseFrac * dh.uncertainty) * f32(Math.abs(dh.channelSize)));
  for (let i = 1; i <= nRings; i++) {
    modelerableRingSwitch[i - 1] = checkProfilePointModelability(
      sdLims, vLims, noiseSDLim,
      tempRadialProfile[i - 1],
      projectedProfile[i - 1],
      modelerableRingSwitch[i - 1],
      opts.incl, vSys, dc.channels
    );
  }

  if (isTraceDebug()) {
    console.error('TRACE Total potentially modelable rings', nRings);
  }
  let paOut = pa;
  if (opts.nTargRings === -1) {
    let sum = 0;
    for (let i = 0; i < modelerableRingSwitch.length; i++) sum += modelerableRingSwitch[i];
    nRings = sum + 1;
  } else {
    nRings = opts.nTargRings;
  }
  if (isTraceDebug()) {
    console.error('TRACE Total modelable number of rings', nRings, nRingsMax);
  }
  if (nRings === 0) {
    throw new Error('EstimateProfiles: no modelable rings have been found');
  }

  const { radialProfile, radialProfileModelSwitch } = assignProfile(
    nRings, nRingsMax, opts.nTargRings, opts.radGridArcsec,
    dh.pixelSize[0], tempRadialProfile, modelerableRingSwitch
  );

  fillInMisingRings(nRings, radialProfile, radialProfileModelSwitch);

  radialProfile[nRings - 1][0] = f32(
    radialProfile[nRings - 2][0] + f32(beamMajorAxis / f32(opts.nRingsPerBeam))
  );
  radialProfile[nRings - 1][1] = radialProfile[nRings - 2][1];
  radialProfile[nRings - 1][2] = radialProfile[nRings - 2][2];
  radialProfileModelSwitch[nRings - 1] = 1;

  if (isTraceDebug()) {
    for (let i = 0; i < nRings; i++) {
      console.error('TRACE Post Fill-in', radialProfile[i][0], radialProfile[i][1], radialProfile[i][2]);
    }
  }

  beamCorrectProfile(nRings, radialProfile, beamMajorAxis, sdLims, vLims);

  if (isTraceDebug()) {
    console.error('TRACE PA Ini', paOut);
  }

  // Leading/trailing velocity comparison -> PA 180-degree flip decision.
  // Uses the WIDE estimatedProfile (from sliceMap, sized by nRingsMax),
  // indexed by the NEW (post-reassignment) nRings -- a subset of its range.
  let trailingVelProfile = f32(0.0);
  for (let i = -nRings + 1; i <= -1; i++) {
    const v = estimatedProfile.get(i - 1)[2];
    if (v === v) { // NaN check
      trailingVelProfile = f32(trailingVelProfile + f32(v - vSys));
    }
  }
  let leadingVelProfile = f32(0.0);
  for (let i = 1; i <= nRings - 1; i++) {
    const v = estimatedProfile.get(i - 1)[2];
    if (v === v) {
      leadingVelProfile = f32(leadingVelProfile + f32(v - vSys));
    }
  }

  const Pi = f32(Math.PI);
  if (leadingVelProfile <= trailingVelProfile) {
    paOut = f32(paOut + Pi);
    if (paOut >= f32(2.0) * Pi) paOut = f32(paOut - f32(2.0) * Pi);
  }

  if (isTraceDebug()) {
    console.error('TRACE Final PA Check', paOut);
  }

  return { radialProfile, radialProfileModelSwitch, pa: paOut, nRings };
}

// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  fortranNint,
  sliceMap,
  buildIniProfile,
  checkProfilePointModelability,
  assignProfile,
  fillInMisingRings,
  beamCorrectProfile,
  estimateProfiles,
};


// ---------------------------------------------------------------------------
// Self-test (node EstimateRadialProfiles.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  console.log('=== fortranNint (ties away from zero, unlike Math.round) ===');
  const cases = [[2.5, 3], [-2.5, -3], [0.4, 0], [-0.4, 0], [1.5, 2], [-1.5, -2]];
  for (const [x, expect] of cases) {
    const r = fortranNint(x);
    console.log(`  nint(${x}) = ${r} (expect ${expect})`, r === expect ? 'OK' : 'FAIL');
  }

  console.log('\n=== assignProfile: pinned grid ===');
  const nRingsMax = 5;
  const temp = makeProfile(nRingsMax);
  for (let i = 0; i < nRingsMax; i++) { temp[i][0] = i * 6; temp[i][1] = 1 + i; temp[i][2] = 10 + i; }
  const modelable = new Int32Array([1, 1, 1, 1, 1]);
  // Grid in arcsec, pixelSize=6 arcsec/pixel -> pixel radii [0,6,18] (skip ring 1 => 6/6=1... let's target rings 0,2 explicitly)
  const radGridArcsec = [0, 12]; // pixel radii [0, 2] -- won't match ring1(pix=6) but should match ring0(pix=0) and... let's just check basic behavior
  const { radialProfile, radialProfileModelSwitch } = assignProfile(2, nRingsMax, 2, radGridArcsec, 6, temp, modelable);
  console.log('  radialProfile:', radialProfile.map(r => Array.from(r)));
  console.log('  modelSwitch  :', Array.from(radialProfileModelSwitch));

  console.log('\n=== fillInMisingRings: interior fill ===');
  const rp = makeProfile(3);
  rp[0].set([0, 5, 50]); rp[1].set([6, 0, 0]); rp[2].set([12, 7, 70]);
  const sw = new Int32Array([1, 0, 1]);
  fillInMisingRings(3, rp, sw);
  console.log('  ring1 after fill (expect SD/VRot from ring0=[5,50] or ring2, ported rule uses leftward-first-found=ring0):',
    Array.from(rp[1]));

  console.log('\n=== fillInMisingRings: ring 0 fallback found (bug fixed, no throw) ===');
  const rp2 = makeProfile(2);
  rp2[0].set([0, 0, 0]); rp2[1].set([6, 9, 90]);
  const sw2 = new Int32Array([0, 1]);
  try {
    fillInMisingRings(2, rp2, sw2);
    console.log('  OK: did not throw, ring0 filled from ring1 (expect [9,90]):', Array.from(rp2[0]).slice(1));
  } catch (e) {
    console.log('  FAIL: threw unexpectedly ->', e.message);
  }

  console.log('\n=== fillInMisingRings: ring 0, no fallback anywhere (still throws) ===');
  const rp3 = makeProfile(2);
  rp3[0].set([0, 0, 0]); rp3[1].set([6, 0, 0]);
  const sw3 = new Int32Array([0, 0]);
  try {
    fillInMisingRings(2, rp3, sw3);
    console.log('  FAIL: expected a throw (no ring anywhere is modelable)');
  } catch (e) {
    console.log('  OK: threw as expected ->', e.message);
  }
}
