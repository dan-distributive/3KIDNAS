'use strict';

// =============================================================================
// EstimateShape.js
// Port of src/PreAnalysis/EstimateShape.f (EstimateShapeMod)
//
// NOT ported here (confirmed dead code, not just unused-by-this-chain):
//   - EstimateSize: grepped its name across all of src/ -- defined but never
//     called anywhere in the Fortran codebase.
//   - CalculateFourerA2Moments: EstimateShape.f's own call site for it is
//     commented out AND the receiving variables (A2, Phi2) aren't even
//     declared in that caller -- it could not compile if uncommented as-is.
//
// DEAD/DISCARDED VALUES (kept faithful to what's actually computed and
// actually used, not to every intermediate Fortran assigns):
//   - calculateFlatMoments: the loop accumulates FTot=sum(flux) per iterated
//     pixel, but immediately after the loop Fortran reassigns
//     `FTot=Moments(0)+Moments(1)`, fully discarding the accumulated value
//     before it's ever read. Not ported (zero effect on output).
//   - estimateShape: computes Q,U -> an intermediate Phi via FullCircATan,
//     scales/wraps it, and this Phi is then unconditionally overwritten by
//     Phi1 (from calculateFourerA1VelMoments) before return. Q and U still
//     feed D -> ellip -> incl (a real, live dependency chain) and are kept;
//     the intermediate Phi computation itself is skipped since it has no
//     effect on either returned value (incl doesn't depend on it, and Phi
//     gets clobbered).
//
// PRECISION
// ---------
// Fortran `real` -> Math.fround() per step, matching the -O0 build.
// `X**2.` (compile-time-literal real exponent) and `X**2` (integer exponent)
// both constant-fold to a plain multiply under gfortran -- confirmed
// separately for each case elsewhere in this project (PhysCoordTransform.f,
// EstimateCubeNoise.f for the real-literal case; integer-exponent `**N` is
// always a repeated-squaring intrinsic in Fortran regardless of optimization
// level, never a generic real**real pow call) -- so all `**2.`/`**2` in this
// file port as v*v, not Math.pow. Contrast with GetMomentMaps.js's
// Arr(i)**real(WMom), where the exponent is a RUNTIME value and a genuine
// pow() call is unavoidable.
// =============================================================================

const f32 = Math.fround;
const { fullCircATan } = require('../StandardMath/FullCircTrig.js');
// FIXED FORTRAN BUG (2026, flagged by Dan, reported to Nathan) --
// calculateFourerA1VelMoments used to call the native Math.cos()/Math.sin()
// intrinsic (Fortran side: native cos()/sin()) inside its per-pixel Fourier
// sum -- the dominant entry point for a real Fortran/JS shape/inclination
// divergence, since it accumulates a tiny native-libm mismatch over every
// pixel in the map rather than a single call. Switched to the shared
// fdlibm port used elsewhere in this pipeline for exactly this reason.
const { fdSin, fdCos } = require('../StandardMath/fdlibm.js');

// ---------------------------------------------------------------------------
// estimateCenter
// Fortran: EstimateCenter(Maps, Center)
//
// Flux-weighted (moment-0 map) centre of brightness, in pixel-index units.
// ---------------------------------------------------------------------------
function estimateCenter(maps) {
  const dh = maps.dh;
  const nx = dh.nPixels[0];
  const ny = dh.nPixels[1];
  let cx = f32(0.0), cy = f32(0.0), fSum = f32(0.0);

  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      const flux0 = maps.flux[i * ny * 3 + j * 3 + 0];
      cx   = f32(cx + f32(f32(i) * flux0));
      cy   = f32(cy + f32(f32(j) * flux0));
      fSum = f32(fSum + flux0);
    }
  }
  return new Float32Array([f32(cx / fSum), f32(cy / fSum)]);
}

// ---------------------------------------------------------------------------
// rLimitedEstimateCenter
// Fortran: RLimited_EstimateCenter(Maps, Center, RLim, CentNew)
//
// Same flux-weighted centre, restricted to pixels within RLim of the
// current center estimate.
// ---------------------------------------------------------------------------
function rLimitedEstimateCenter(maps, center, rLim) {
  const dh = maps.dh;
  const nx = dh.nPixels[0];
  const ny = dh.nPixels[1];
  let cx = f32(0.0), cy = f32(0.0), fSum = f32(0.0);

  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      const x = f32(i), y = f32(j);
      const xx = f32(x - center[0]);
      const yy = f32(y - center[1]);
      const r  = f32(Math.sqrt(f32(f32(xx * xx) + f32(yy * yy))));
      if (r <= rLim) {
        const flux0 = maps.flux[i * ny * 3 + j * 3 + 0];
        cx   = f32(cx + f32(x * flux0));
        cy   = f32(cy + f32(y * flux0));
        fSum = f32(fSum + flux0);
      }
    }
  }
  return new Float32Array([f32(cx / fSum), f32(cy / fSum)]);
}

// ---------------------------------------------------------------------------
// iterEstimateCenter
// Fortran: Iter_EstimateCenter(Maps, Center, CenterFlag)
//
// Iteratively refines the flux-weighted centre with a shrinking radius cut,
// up to 10 iterations, converging once the shift between iterations is
// <=0.5 pixel (and at least 3 iterations have run). Returns { center,
// centerFlag } -- centerFlag=1 means it hit the iteration cap without
// converging (matching Fortran's own fallback/warning path).
// ---------------------------------------------------------------------------
function iterEstimateCenter(maps) {
  const dh = maps.dh;
  let centerFlag = 0;
  let center = estimateCenter(maps);

  const nx = f32(dh.nPixels[0]);
  const ny = f32(dh.nPixels[1]);
  let rMax = f32(Math.sqrt(f32(f32(nx * nx) + f32(ny * ny))));
  rMax = f32(rMax / f32(2.0));
  let rLim = f32(f32(0.5) * rMax);

  let centTemp = rLimitedEstimateCenter(maps, center, rLim);

  const imax = 10;
  for (let i = 1; i <= imax; i++) {
    const dx = f32(center[0] - centTemp[0]);
    const dy = f32(center[1] - centTemp[1]);
    const delR = f32(Math.sqrt(f32(f32(dx * dx) + f32(dy * dy))));

    center = centTemp;

    if (delR > f32(0.5) || i <= 3) {
      centTemp = rLimitedEstimateCenter(maps, center, rLim);
      rLim = f32(rMax / f32(i + 2));
    } else {
      return { center, centerFlag };
    }
  }

  // Iteration cap reached without converging -- Fortran's own fallback.
  center = centTemp;
  centerFlag = 1;
  return { center, centerFlag };
}

// ---------------------------------------------------------------------------
// calculateFlatMoments
// Fortran: CalculateFlatMoments(Maps, Moments, RMax, Center)
//
// Second spatial moments (xx, yy, xy) of the moment-0 map about `center`,
// restricted to R <= RLim = 0.25*sqrt(nx^2+ny^2) (half of half the diagonal),
// normalized by (Moments(0)+Moments(1)) -- NOT by total flux (see module
// header: the loop's own FTot accumulation is discarded before use).
// Returns { moments: Float32Array(3), rMax } (rMax capped at RLim by
// construction, since it's only updated inside the R<=RLim branch).
// ---------------------------------------------------------------------------
function calculateFlatMoments(maps, center) {
  const dh = maps.dh;
  const nx = dh.nPixels[0];
  const ny = dh.nPixels[1];
  const moments = new Float32Array(3);
  let rMax = f32(0.0);

  const nxf = f32(nx), nyf = f32(ny);
  let rLim = f32(Math.sqrt(f32(f32(nxf * nxf) + f32(nyf * nyf))));
  rLim = f32(rLim / f32(2.0));
  rLim = f32(f32(0.5) * rLim);
  const rLimSq = f32(rLim * rLim);

  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      const x = f32(f32(i) - center[0]);
      const y = f32(f32(j) - center[1]);
      const r2 = f32(f32(x * x) + f32(y * y));
      const flux0 = maps.flux[i * ny * 3 + j * 3 + 0];

      if (flux0 !== f32(0.0)) {
        if (r2 !== f32(0.0) && r2 <= rLimSq) {
          moments[0] = f32(moments[0] + f32(f32(x * x) * flux0));
          moments[1] = f32(moments[1] + f32(f32(y * y) * flux0));
          moments[2] = f32(moments[2] + f32(f32(x * y) * flux0));
          const r = f32(Math.sqrt(r2));
          if (r >= rMax) rMax = r;
        }
      }
    }
  }

  const fTot = f32(moments[0] + moments[1]);
  moments[0] = f32(moments[0] / fTot);
  moments[1] = f32(moments[1] / fTot);
  moments[2] = f32(moments[2] / fTot);

  return { moments, rMax };
}

// ---------------------------------------------------------------------------
// calculateFourerA1VelMoments
// Fortran: CalculateFourerA1VelMoments(Maps, A1, Phi, RMax, Center)
//
// First-harmonic Fourier decomposition of the moment-1 (velocity) map in
// azimuth, restricted to the same R<=RLim cut as calculateFlatMoments.
// Returns { a1, phi }.
// ---------------------------------------------------------------------------
function calculateFourerA1VelMoments(maps, center) {
  const dh = maps.dh;
  const nx = dh.nPixels[0];
  const ny = dh.nPixels[1];

  const nxf = f32(nx), nyf = f32(ny);
  let rLim = f32(Math.sqrt(f32(f32(nxf * nxf) + f32(nyf * nyf))));
  rLim = f32(rLim / f32(2.0));
  rLim = f32(f32(0.5) * rLim);
  const rLimSq = f32(rLim * rLim);

  let fTot = f32(0.0), aa1 = f32(0.0), bb1 = f32(0.0);

  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      const x = f32(f32(i) - center[0]);
      const y = f32(f32(j) - center[1]);
      const r2 = f32(f32(x * x) + f32(y * y));
      const theta = fullCircATan(x, y);
      const base = i * ny * 3 + j * 3;
      const flux0 = maps.flux[base];
      const flux1 = maps.flux[base + 1];

      if (flux0 !== f32(0.0)) {
        if (r2 !== f32(0.0) && r2 <= rLimSq) {
          aa1  = f32(aa1 + f32(flux1 * f32(fdCos(theta))));
          bb1  = f32(bb1 + f32(flux1 * f32(fdSin(theta))));
          fTot = f32(fTot + flux0);
        }
      }
    }
  }

  const a1  = f32(Math.sqrt(f32(f32(f32(aa1 * aa1) + f32(bb1 * bb1)) / f32(fTot * fTot))));
  const phi = fullCircATan(aa1, bb1);
  return { a1, phi };
}

// ---------------------------------------------------------------------------
// estimateShape
// Fortran: EstimateShape(Maps, Center, incl, Phi)
//
// incl: from the flat (0th-moment) shape's Q/U ellipticity.
// Phi:  from the 1st-harmonic velocity-map Fourier decomposition (the flat
//       moments' own Q/U-derived Phi is computed then discarded in Fortran
//       -- see module header -- so it's skipped here too).
// ---------------------------------------------------------------------------
function estimateShape(maps, center) {
  const { moments } = calculateFlatMoments(maps, center);
  const q = f32(moments[0] - moments[1]);
  const u = f32(f32(2.0) * moments[2]);
  const d = f32(Math.sqrt(f32(f32(q * q) + f32(u * u))));
  const ellip = f32(f32(f32(1.0) - d) / f32(f32(1.0) + d));
  const incl = f32(Math.acos(ellip));

  const { phi: phi1 } = calculateFourerA1VelMoments(maps, center);

  return { incl, phi: phi1 };
}

// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
  estimateCenter,
  rLimitedEstimateCenter,
  iterEstimateCenter,
  calculateFlatMoments,
  calculateFourerA1VelMoments,
  estimateShape,
};


// ---------------------------------------------------------------------------
// Self-test (node EstimateShape.js)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const { DataCube, allocateDataCube } = require('../ObjectDefinitions/DataCube.js');

  // Symmetric, centered "moment map" cube: a flat-topped square blob of
  // constant moment-0 flux centered exactly on the cube, zero elsewhere.
  // A perfectly symmetric blob should estimate its center at the geometric
  // center of the blob and inclination near 0 (circular/face-on: Q,U -> 0
  // -> D -> 0 -> ellip -> 1 -> incl -> 0).
  const n = 21; // odd, so there's an exact center pixel
  const maps = new DataCube();
  maps.dh.nPixels[0] = n; maps.dh.nPixels[1] = n; maps.dh.nChannels = 3;
  allocateDataCube(maps);
  const mid = Math.floor(n / 2);
  const R = 5;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const dx = i - mid, dy = j - mid;
      const inside = (dx * dx + dy * dy) <= R * R;
      const base = i * n * 3 + j * 3;
      maps.flux[base]     = inside ? 1.0 : 0.0; // moment 0
      maps.flux[base + 1] = 0.0;                // moment 1 (velocity) -- flat, no rotation signal
      maps.flux[base + 2] = 0.0;
    }
  }

  console.log('=== estimateCenter on a symmetric disc ===');
  const c = estimateCenter(maps);
  console.log('  center:', Array.from(c), `(expect [${mid}, ${mid}])`);

  console.log('\n=== iterEstimateCenter ===');
  const { center: cIter, centerFlag } = iterEstimateCenter(maps);
  console.log('  center:', Array.from(cIter), `(expect ~[${mid}, ${mid}])`, 'centerFlag:', centerFlag);

  console.log('\n=== estimateShape on a circularly symmetric, non-rotating disc ===');
  const { incl, phi } = estimateShape(maps, cIter);
  console.log('  incl (rad):', incl, '(expect ~0, symmetric disc)');
  console.log('  phi  (rad):', phi, '(velocity map is flat/zero -> aa1=bb1=0 -> phi undefined/NaN-ish, informational only)');
}
