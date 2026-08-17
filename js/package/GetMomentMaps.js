module.declare(["./DataCube.js"], function (require, exports, module) {
'use strict';

// =============================================================================
// GetMomentMaps.js
// Port of src/PreAnalysis/GetMomentMaps.f (GetMomentMapsMod) +
// the MakeMomentMaps/CalculateGeneralMoment/MakeVelProfile subset of
// src/MomentMaps/CalculateMomentMaps.f (MakeMomentMapsMod) that
// InitialAnalysis_PreProcessing.f's ConstructProjectionsFromCube actually
// calls (MakeMomentMaps + CalculateGeneralMoment via ConstructMomentMaps,
// MakeVelProfile directly for ObservedVelocityProfile).
//
// CORRECTION: an earlier pass of this port claimed MakeVelProfile was
// unused -- wrong, ConstructProjectionsFromCube calls it directly (confirmed
// by reading InitialAnalysis_PreProcessing.f, not just GetMomentMaps.f's own
// callers). Added below.
//
// NOT ported here (confirmed unused by this call chain -- grepped
// InitialAnalysis.f / InitialAnalysis_PreProcessing.f, neither calls them):
//   - GetMomentMaps.f's FreqencyMomMapsConvert (WALLABY frequency-cubelet ->
//     km/s conversion). DataCube.js already carries restFreq/delFreq fields
//     for whenever this IS needed, but InitialAnalysis_Prep's real call
//     chain (-> ConstructMomentMaps) never reaches it.
//   - CalculateMomentMaps.f's RadioMomentMaps/CalcPVDiagrams (different
//     callers, not part of this chain).
//
// FIXED FORTRAN BUG (2026, flagged by Dan, reported to Nathan) -- last
// channel used to be silently dropped from every moment:
//   MakeMomentMaps called CalculateGeneralMoment(DC%DH%nChannels-1, ...),
//   passing n = nChannels-1 as the "number of elements" argument. Inside
//   CalculateGeneralMoment, the dummy arrays were declared Arr(0:n-1)/
//   WArr(0:n-1) using THAT n, so the `do i=0,n-1` loop only visited channel
//   indices 0..nChannels-2 -- the LAST channel (index nChannels-1) was never
//   summed into moment 0, 1, or 2, for every pixel, every moment map, every
//   InitialAnalysis run. An unintentional off-by-one (n was meant to be
//   nChannels). Fixed on both sides -- native Fortran now passes nChannels
//   (src/MomentMaps/CalculateMomentMaps.f), and this port now sums every
//   channel to match.
//
// PRECISION
// ---------
// Arr(i)**real(WMom): WMom is a RUNTIME integer (loop variable 0/1/2), not a
// compile-time literal like the X**2. cases already resolved elsewhere in
// this project (PhysCoordTransform.f, EstimateCubeNoise.f) -- gfortran
// cannot constant-fold this to a multiply, so it's a genuine runtime
// real**real power call. Verified empirically (gfortran -O0, standalone
// scratch test): x**real(1) and x**real(2) for x=3.14159265 (f32) match
// Math.pow(x,1)/Math.pow(x,2) bit-for-bit (0x40490FDB, 0x411DE9E7) -- also
// matching plain x*x for the WMom=2 case at this test value. Ported as
// Math.pow rather than a hand-rolled x*x special-case, since Math.pow
// already agrees and this stays correct if CalculateGeneralMoment is ever
// called with a WMom this project hasn't exercised yet.
// =============================================================================

const f32 = Math.fround;
const { DataCube, allocateDataCube, flatIndxCalc } = require('./DataCube.js');

// ---------------------------------------------------------------------------
// calculateGeneralMoment
// Fortran: CalculateGeneralMoment(n, Arr, WArr, WMom, CalcMoment)
//
// channels/flux: full-length arrays (length nChannels). Sums over every
// channel (see the module header's fixed-bug note -- this used to truncate
// to nChannels-1, matching a Fortran off-by-one that's now fixed on both
// sides).
//
// WMom=0: CalcMoment = sum(flux)                              (total flux)
// WMom>0: CalcMoment = sum(channel^WMom * flux) / sum(flux)    (weighted moment)
// ---------------------------------------------------------------------------
function calculateGeneralMoment(channels, flux, wMom) {
  const n = channels.length; // Fortran: n = DC%DH%nChannels (fixed, was nChannels-1)
  let calcMoment = f32(0.0);
  if (wMom === 0) {
    for (let i = 0; i < n; i++) {
      calcMoment = f32(calcMoment + f32(flux[i]));
    }
  } else {
    let wsum = f32(0.0);
    for (let i = 0; i < n; i++) {
      wsum = f32(wsum + f32(flux[i]));
      const term = f32(f32(Math.pow(channels[i], wMom)) * f32(flux[i]));
      calcMoment = f32(calcMoment + term);
    }
    calcMoment = f32(calcMoment / wsum);
  }
  return calcMoment;
}

// ---------------------------------------------------------------------------
// makeMomentMaps
// Fortran: MakeMomentMaps(DC, Maps)  (src/MomentMaps/CalculateMomentMaps.f)
//
// dc: DataCube (masked observed cube). maps: DataCube pre-allocated with
// nChannels=3 -- moment 0/1/2 stored as "channels" 0/1/2 of maps.flux.
// ---------------------------------------------------------------------------
function makeMomentMaps(dc, maps) {
  const dh  = dc.dh;
  const nx  = dh.nPixels[0];
  const ny  = dh.nPixels[1];
  const mdh = maps.dh;

  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      const chanStart = flatIndxCalc(i, j, 0, dh);
      const flux = dc.flux.subarray(chanStart, chanStart + dh.nChannels);

      const moments = new Float32Array(3);
      for (let moment = 0; moment <= 2; moment++) {
        let m = calculateGeneralMoment(dc.channels, flux, moment);
        if (moment === 0 && m < f32(0.0)) m = f32(0.0);
        moments[moment] = m;
      }

      // Moment 0 adjusted by channel width; moment 2 converted from a raw
      // flux-weighted <v^2> into variance (subtract <v>^2), then dispersion
      // (sqrt). moments[1]*moments[1] here IS a compile-time-literal-exponent
      // case (Fortran's `**2.`), so it constant-folds to a plain multiply --
      // unlike calculateGeneralMoment's runtime-WMom pow above.
      moments[0] = f32(moments[0] * f32(Math.abs(dh.channelSize)));
      moments[2] = f32(moments[2] - f32(moments[1] * moments[1]));
      moments[2] = f32(Math.sqrt(moments[2]));

      const outBase = flatIndxCalc(i, j, 0, mdh);
      maps.flux[outBase]     = moments[0];
      maps.flux[outBase + 1] = moments[1];
      maps.flux[outBase + 2] = moments[2];
    }
  }
}

// ---------------------------------------------------------------------------
// constructMomentMaps
// Fortran: ConstructMomentMaps(MaskedObsDC, ObsMaps)
//
// Allocates a 3-channel output cube (moment 0/1/2) with the same spatial
// header as the input (full header copy, then nChannels overridden to 3 --
// matching Fortran's `ObsMaps%DH=MaskedObsDC%DH; ObsMaps%DH%nChannels=3`),
// then fills it via makeMomentMaps. Header-copy field list matches
// FlipBootstrap.js's copyDataCube (this project's established convention
// for a Fortran derived-type header assignment).
// ---------------------------------------------------------------------------
function constructMomentMaps(maskedObsDC) {
  const obsMaps = new DataCube();
  const dh  = obsMaps.dh;
  const src = maskedObsDC.dh;

  dh.nPixels[0]     = src.nPixels[0];
  dh.nPixels[1]     = src.nPixels[1];
  dh.nChannels      = 3;
  dh.pixelSize[0]   = src.pixelSize[0];
  dh.pixelSize[1]   = src.pixelSize[1];
  dh.channelSize    = src.channelSize;
  dh.refLocation[0] = src.refLocation[0];
  dh.refLocation[1] = src.refLocation[1];
  dh.refLocation[2] = src.refLocation[2];
  dh.refVal[0]      = src.refVal[0];
  dh.refVal[1]      = src.refVal[1];
  dh.refVal[2]      = src.refVal[2];
  dh.uncertainty    = src.uncertainty;

  allocateDataCube(obsMaps);
  dh.start[0] = src.start[0];
  dh.start[1] = src.start[1];
  dh.start[2] = src.start[2];

  makeMomentMaps(maskedObsDC, obsMaps);
  return obsMaps;
}

// ---------------------------------------------------------------------------
// makeVelProfile
// Fortran: MakeVelProfile(DC, VelProf)  (src/MomentMaps/CalculateMomentMaps.f)
//
// Integrated 1D flux-vs-channel profile: for each channel, sums flux over
// every spatial pixel. Returns { vel, flux } in VelProfileAnalysis.js's
// working-array shape (vel = channel velocity axis, flux = spatial sum).
// ---------------------------------------------------------------------------
function makeVelProfile(dc) {
  const dh  = dc.dh;
  const nx  = dh.nPixels[0];
  const ny  = dh.nPixels[1];
  const nCh = dh.nChannels;

  const vel  = new Float32Array(nCh);
  const flux = new Float32Array(nCh);

  for (let i = 0; i < nCh; i++) {
    vel[i] = dc.channels[i];
    let sum = f32(0.0);
    for (let j = 0; j < nx; j++) {
      for (let k = 0; k < ny; k++) {
        sum = f32(sum + f32(dc.flux[flatIndxCalc(j, k, i, dh)]));
      }
    }
    flux[i] = sum;
  }

  return { vel, flux };
}

// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = { constructMomentMaps, makeMomentMaps, calculateGeneralMoment, makeVelProfile };


// ---------------------------------------------------------------------------
// Self-test (node GetMomentMaps.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { DataCube: DC, allocateDataCube: alloc } = require('./DataCube.js');

  console.log('=== calculateGeneralMoment: all channels now included (bug fixed) ===');
  // 4 channels, all flux=1: moment 0 should be 4 (every channel), not 3.
  const channels = new Float32Array([10, 20, 30, 40]);
  const flux     = new Float32Array([1, 1, 1, 1]);
  const m0 = calculateGeneralMoment(channels, flux, 0);
  console.log('  moment0 (expect 4):', m0, m0 === 4 ? 'OK' : 'FAIL');
  // moment1 = sum(chan*flux)/sum(flux) over all 4 channels = (10+20+30+40)/4 = 25
  const m1 = calculateGeneralMoment(channels, flux, 1);
  console.log('  moment1 (expect 25)      :', m1, m1 === 25 ? 'OK' : 'FAIL');

  console.log('\n=== makeMomentMaps / constructMomentMaps: uniform cube ===');
  // A spectrally-flat cube (constant flux F across all channels) at every
  // pixel, with a fully symmetric-about-zero 5-channel set so moment1/moment2
  // have a clean closed form and the test would fail loudly if any channel
  // (including the last) were ever dropped again.
  const dc = new DC();
  dc.dh.nPixels[0] = 3; dc.dh.nPixels[1] = 2; dc.dh.nChannels = 5;
  dc.dh.channelSize = f32(2.0);
  alloc(dc);
  // Symmetric channel values around 0, all 5 channels: [-4,-2,0,2,4] --
  // mean = 0, so moment2 (dispersion) = sqrt(mean(chan^2)) =
  // sqrt((16+4+0+4+16)/5) = sqrt(8).
  dc.channels.set([-4, -2, 0, 2, 4]);
  const F = f32(4.0);
  dc.flux.fill(F);

  const maps = constructMomentMaps(dc);
  console.log('  maps nChannels (expect 3):', maps.dh.nChannels);
  const base = flatIndxCalc(0, 0, 0, maps.dh);
  const expectedM0 = f32(F * 5 * f32(Math.abs(dc.dh.channelSize))); // all 5 channels summed, * |channelSize|
  const expectedM2 = f32(Math.sqrt(8));
  console.log('  moment0[0,0] (expect', expectedM0, '):', maps.flux[base],
    maps.flux[base] === expectedM0 ? 'OK' : 'FAIL');
  console.log('  moment1[0,0] (expect 0)          :', maps.flux[base + 1],
    maps.flux[base + 1] === 0 ? 'OK' : 'FAIL');
  console.log('  moment2[0,0] (expect', expectedM2.toFixed(6), '):', maps.flux[base + 2],
    Math.abs(maps.flux[base + 2] - expectedM2) < 1e-6 ? 'OK' : 'FAIL');

  console.log('\n=== makeVelProfile ===');
  // Reuse dc from the block above: 3x2 spatial, 5 channels, flux=4 everywhere.
  const { vel, flux: velFlux } = makeVelProfile(dc);
  console.log('  vel  :', Array.from(vel), '(expect the channel array itself)');
  console.log('  flux :', Array.from(velFlux), '(expect 4*6=24 at every channel, 3x2 spatial pixels)');
  const velOk = velFlux.every((v) => v === 24);
  console.log('  ', velOk ? 'OK' : 'FAIL');
}

});
