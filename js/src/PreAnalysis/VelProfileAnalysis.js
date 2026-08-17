'use strict';

// =============================================================================
// VelProfileAnalysis.js
// Port of src/PreAnalysis/VelProfileAnalysis.f (VelProfileAnalysisMod)
//
// A Fortran "VelProfile(0:1, 0:nChannels-1)" (row 0 = velocity/channel axis,
// row 1 = summed flux) is represented here as a plain object
// { vel: Float32Array(nCh), flux: Float32Array(nCh) } rather than replicating
// the 2D column-major layout -- this array is a local working structure, not
// something read from/written to a shared DataCube-like layout elsewhere, so
// there's no cross-file indexing convention to preserve.
//
// CALLER NOTE: InitialAnalysis.f's real ProjectedAnalysis flow calls
// EstimateVSysAndEdges_Simple (the last, simplest subroutine below), NOT the
// full peak-finding EstimateVSysAndEdges pipeline (confirmed by reading
// InitialAnalysis.f directly). The whole file is ported anyway (BoxCarAvg /
// FluxLims / FindPeakFlux / FindEdges / FirstDeriv / LinearInterp_YTarg /
// EstimateVSysAndEdges), since EstimateRadialProfiles.f (next file in this
// port) may call the more complex pipeline internally for per-ring work --
// to be confirmed when that file is read.
//
// DEAD/UNUSED VALUES:
//   - EstimateVSysAndEdges: VGuess=VelProfile(1,VGuessIndx) is computed then
//     never read again anywhere in the subroutine. Not ported.
//   - EstimateVSysAndEdges_Simple: FracLimit is a dummy argument the Fortran
//     body never references. Dropped from this port's signature entirely.
//
// PRECISION
// ---------
// Fortran `real` -> Math.fround() per step, matching the -O0 build.
// Fortran's `sum()` intrinsic (FluxLims) ported as a plain sequential
// left-to-right accumulation, consistent with how every other running sum
// in this project is ported (no reduction-tree/vectorized reordering).
// =============================================================================

const f32 = Math.fround;

// ---------------------------------------------------------------------------
// boxCarAvg
// Fortran: BoxCarAvg(n, Arr, SmoothArr, nBoxCar)
//
// Simple boxcar smoothing; averages only over in-bounds neighbors near the
// array edges (nUse shrinks there), not zero-padded.
// ---------------------------------------------------------------------------
function boxCarAvg(arr, nBoxCar) {
  const n = arr.length;
  const smoothArr = new Float32Array(n);
  const j = Math.trunc(nBoxCar / 2);

  for (let i = 0; i < n; i++) {
    let nUse = 0;
    let sum = f32(0.0);
    for (let k = -j; k <= j; k++) {
      if (i + k >= 0 && i + k <= n - 1) {
        nUse += 1;
        sum = f32(sum + f32(arr[i + k]));
      }
    }
    smoothArr[i] = f32(sum / f32(nUse));
  }
  return smoothArr;
}

// ---------------------------------------------------------------------------
// fluxLims
// Fortran: FluxLims(n, SmoothFlux, FracLim, LimGuess)
//
// First/last index where SmoothFlux exceeds FracLim*mean(SmoothFlux).
// Returns Int32Array([lo, hi]); either stays -1 if no index qualifies
// (matching Fortran's LimGuess=-1 initialization with no subsequent set).
// ---------------------------------------------------------------------------
function fluxLims(smoothFlux, fracLim) {
  const n = smoothFlux.length;
  let sum = f32(0.0);
  for (let i = 0; i < n; i++) sum = f32(sum + f32(smoothFlux[i]));
  const avgFlux = f32(sum / f32(n));
  const fluxLim = f32(avgFlux * fracLim);

  const limGuess = new Int32Array([-1, -1]);
  for (let i = 0; i < n; i++) {
    if (smoothFlux[i] > fluxLim) { limGuess[0] = i; break; }
  }
  for (let i = n - 1; i >= 0; i--) {
    if (smoothFlux[i] > fluxLim) { limGuess[1] = i; break; }
  }
  return limGuess;
}

// ---------------------------------------------------------------------------
// firstDeriv
// Fortran: FirstDeriv(n, X, Y, YP)
//
// Simple forward finite difference between consecutive points; length n-1.
// ---------------------------------------------------------------------------
function firstDeriv(x, y) {
  const n = x.length;
  const yp = new Float32Array(n - 1);
  for (let i = 0; i < n - 1; i++) {
    const dX = f32(x[i + 1] - x[i]);
    const dY = f32(y[i + 1] - y[i]);
    yp[i] = f32(dY / dX);
  }
  return yp;
}

// ---------------------------------------------------------------------------
// linearInterpYTarg
// Fortran: LinearInterp_YTarg(X, Y, YTarg, XFound)
//
// Linear interpolation for the X where Y=YTarg, given a 2-point bracket
// (x0,y0)-(x1,y1).
// ---------------------------------------------------------------------------
function linearInterpYTarg(x0, x1, y0, y1, yTarg) {
  const dX   = f32(x1 - x0);
  const dY   = f32(y1 - y0);
  const delY = f32(yTarg - y0);
  const delX = f32(f32(delY / dY) * dX);
  return f32(delX + x0);
}

// ---------------------------------------------------------------------------
// findPeakFlux
// Fortran: FindPeakFlux(n, nPeaks, VelProfile, EdgeGuessIndx, VGuessIndx,
//                        PeakIndx, PeakVel, PeakFlux)
//
// Scans outward from vGuessIndx in mirrored pairs (i increasing, j
// decreasing). nPeaks=1: both sides compete for the single peak slot 0.
// Otherwise: i feeds slot 1, j feeds slot 0 (two independent peaks, one per
// side of the profile).
// ---------------------------------------------------------------------------
function findPeakFlux(velProfile, edgeGuessIndx, vGuessIndx, nPeaks) {
  const { vel, flux } = velProfile;
  const peakFlux = new Float32Array([-1, -1]);
  const peakIndx = new Int32Array(2);
  const peakVel  = new Float32Array(2);

  for (let i = vGuessIndx; i <= edgeGuessIndx[1] + 1; i++) {
    const j = vGuessIndx - (i - vGuessIndx);
    if (nPeaks === 1) {
      if (flux[i] > peakFlux[0]) { peakFlux[0] = flux[i]; peakIndx[0] = i; peakVel[0] = vel[i]; }
      if (flux[j] > peakFlux[0]) { peakFlux[0] = flux[j]; peakIndx[0] = j; peakVel[0] = vel[j]; }
    } else {
      if (flux[i] > peakFlux[1]) { peakFlux[1] = flux[i]; peakIndx[1] = i; peakVel[1] = vel[i]; }
      if (flux[j] > peakFlux[0]) { peakFlux[0] = flux[j]; peakIndx[0] = j; peakVel[0] = vel[j]; }
    }
  }
  return { peakIndx, peakVel, peakFlux };
}

// ---------------------------------------------------------------------------
// findEdges
// Fortran: FindEdges(n, nPeaks, VelProfile, PeakIndx, PeakFlux, FracLimit,
//                     EdgeIndx, Edges)
//
// From each peak, scans toward the profile's own edge (j=0: backward from
// peak 0 to channel 0; j=1: forward from peak k to the last channel) for the
// first point at/below FracLimit*peakFlux, then linearly interpolates the
// exact crossing velocity. edgeIndx[j]/edges[j] stay undefined if no
// crossing is found (matching Fortran's uninitialized-if-never-set OUT args).
// ---------------------------------------------------------------------------
function findEdges(velProfile, peakIndx, peakFlux, fracLimit, nPeaks) {
  const { vel, flux } = velProfile;
  const n = vel.length;
  const k = (nPeaks === 1) ? 0 : 1;
  const edgeIndx = new Int32Array(2);
  const edges = new Float32Array(2).fill(NaN);

  for (let j = 0; j <= 1; j++) {
    let dir, start, fin, fluxTarg;
    if (j === 0) {
      dir = -1; start = peakIndx[0]; fin = 0;
      fluxTarg = f32(peakFlux[0] * fracLimit);
    } else {
      dir = 1; start = peakIndx[k]; fin = n - 1;
      fluxTarg = f32(peakFlux[k] * fracLimit);
    }

    for (let i = start; dir === 1 ? i <= fin : i >= fin; i += dir) {
      if (flux[i] <= fluxTarg) {
        edgeIndx[j] = i;
        const l = (j === 0) ? i : i - 1;
        edges[j] = linearInterpYTarg(vel[l], vel[l + 1], flux[l], flux[l + 1], fluxTarg);
        break;
      }
    }
  }
  return { edgeIndx, edges };
}

// ---------------------------------------------------------------------------
// estimateVSysAndEdges
// Fortran: EstimateVSysAndEdges(nChannels, VelProfile, FracLimit, VSys, Edges)
//
// Full peak-finding pipeline: smooth -> derivative sign changes (# peaks) ->
// locate peaks -> locate edges at FracLimit*peakFlux -> VSys = edge midpoint.
// ---------------------------------------------------------------------------
function estimateVSysAndEdges(velProfile, fracLimit) {
  const { vel, flux } = velProfile;
  const nBoxCar = 5;
  const avgLim  = f32(1.0);

  const smoothFlux = boxCarAvg(flux, nBoxCar);
  const smoothFluxDeriv = firstDeriv(vel, smoothFlux);
  const edgeGuessIndx = fluxLims(smoothFlux, avgLim);
  const vGuessIndx = Math.trunc((edgeGuessIndx[1] + edgeGuessIndx[0]) / 2);

  let nPeaks = 0;
  for (let i = edgeGuessIndx[0] - 1; i <= edgeGuessIndx[1]; i++) {
    if (f32(smoothFluxDeriv[i] * smoothFluxDeriv[i + 1]) < f32(0.0)) nPeaks += 1;
  }

  const { peakIndx, peakFlux } = findPeakFlux(velProfile, edgeGuessIndx, vGuessIndx, nPeaks);
  const { edges } = findEdges(velProfile, peakIndx, peakFlux, fracLimit, nPeaks);
  const vSys = f32(f32(edges[0] + edges[1]) / f32(2.0));
  return { vSys, edges };
}

// ---------------------------------------------------------------------------
// estimateVSysAndEdgesSimple
// Fortran: EstimateVSysAndEdges_Simple(nChannels, VelProfile, FracLimit,
//                                       VSys, Edges)
//
// The version InitialAnalysis.f's real ProjectedAnalysis flow actually
// calls. First/last channel with flux >= 0; VSys = their midpoint.
// ---------------------------------------------------------------------------
function estimateVSysAndEdgesSimple(velProfile) {
  const { vel, flux } = velProfile;
  const n = vel.length;
  let edge0, edge1;

  for (let i = 0; i < n; i++) {
    if (flux[i] >= f32(0.0)) { edge0 = vel[i]; break; }
  }
  for (let i = n - 1; i >= 0; i--) {
    if (flux[i] >= f32(0.0)) { edge1 = vel[i]; break; }
  }

  const vSys = f32(f32(edge0 + edge1) / f32(2.0));
  return { vSys, edges: new Float32Array([edge0, edge1]) };
}

// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  boxCarAvg,
  fluxLims,
  firstDeriv,
  linearInterpYTarg,
  findPeakFlux,
  findEdges,
  estimateVSysAndEdges,
  estimateVSysAndEdgesSimple,
};


// ---------------------------------------------------------------------------
// Self-test (node VelProfileAnalysis.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  console.log('=== boxCarAvg ===');
  const arr = new Float32Array([0, 0, 10, 0, 0]);
  const smooth = boxCarAvg(arr, 3); // j=1: window [i-1,i+1], edge-clipped
  console.log('  in :', Array.from(arr));
  console.log('  out:', Array.from(smooth), '(expect center-weighted smear of the spike)');

  console.log('\n=== firstDeriv / linearInterpYTarg ===');
  const x = new Float32Array([0, 1, 2, 3]);
  const y = new Float32Array([0, 2, 4, 6]);
  console.log('  deriv (expect all 2):', Array.from(firstDeriv(x, y)));
  console.log('  interp y=3 in [1,2]->[2,4] (expect x=1.5):', linearInterpYTarg(1, 2, 2, 4, 3));

  console.log('\n=== estimateVSysAndEdgesSimple: synthetic double-horn profile ===');
  // Flux >=0 only in channels 10..40 of a 50-channel profile; VSys should be
  // the midpoint of those two edge velocities.
  const n = 50;
  const vel  = new Float32Array(n);
  const flux = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    vel[i] = i * 2; // arbitrary velocity axis
    flux[i] = (i >= 10 && i <= 40) ? 5 - Math.abs(i - 25) * 0.1 : -1;
  }
  const { vSys, edges } = estimateVSysAndEdgesSimple({ vel, flux });
  console.log('  edges:', Array.from(edges), '(expect [20, 80])');
  console.log('  VSys :', vSys, '(expect 50)', vSys === 50 ? 'OK' : 'FAIL');
}
