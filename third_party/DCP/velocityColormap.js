'use strict';

// =============================================================================
// velocityColormap.js
// Port of the custom LCh-based "double velocity" diverging colormap
// MomentMapPlotFncs.py's SetCMap uses for the moment-1 (velocity field)
// panel: third_party/CosmosCanvas/velmap.py's create_cmap_doubleVelocity ->
// create_cmap_velocity -> third_party/colourspace/maps.py's
// make_cmap_segmented -> third_party/colourspace/convert.py's LCh->Lab->
// XYZ->sRGB chain.
//
// convert.py supports 3 backends (custom formulas / colorspacious /
// colour-science) but only ever calls set_convertor('custom') (line 60,
// module load, the only call anywhere in this project) -- so 'custom' is
// the actual pipeline in production, not a fallback; ported that one,
// verbatim (D65 illuminant, sRGB primaries, easyrgb.com-derived formulas).
//
// DELIBERATE SIMPLIFICATION (flagged, not silently dropped -- see
// ~/.claude/plans/breezy-launching-nova.md Phase 1): the Python side clips
// chroma to the actual sRGB gamut boundary per (L,H) pair before conversion
// (gamut.py's Cmax_for_LH, an iterative gamut-boundary search) so no RGB
// channel ever needs post-hoc clamping. This port skips that precomputed
// gamut-boundary table and just clamps the final RGB to [0,1] instead --
// the colormap's own Cval_max=35 is a moderate chroma, unlikely to blow out
// the gamut badly for these control points, but on a wide/extreme velocity
// range the two approaches could diverge slightly at the extremes. Revisit
// with a real gamut-boundary port if that's ever visually significant.
// =============================================================================

const D65 = { Xn: 95.047, Yn: 100.000, Zn: 108.883 };

function fReverse(x) {
  return x > 6 / 29 ? x ** 3 : 3 * (6 / 29) ** 2 * (x - 4 / 29);
}

function gammaForward(cLin) {
  return cLin <= 0.0031308 ? 12.92 * cLin : 1.055 * Math.pow(cLin, 1 / 2.4) - 0.055;
}

// LCh(ab) -> Lab (cylindrical -> cartesian; H in degrees)
function lchToLab(L, C, H) {
  const hRad = (H * 2 * Math.PI) / 360;
  return [L, C * Math.cos(hRad), C * Math.sin(hRad)];
}

// Lab -> XYZ (D65)
function labToXyz(L, a, b) {
  const X = D65.Xn * fReverse((L + 16) / 116 + a / 500);
  const Y = D65.Yn * fReverse((L + 16) / 116);
  const Z = D65.Zn * fReverse((L + 16) / 116 - b / 200);
  return [X, Y, Z];
}

// XYZ -> sRGB (D65), gamma-companded, in [0,1] before clamping
function xyzToRgb(X, Y, Z) {
  const R = gammaForward(3.2406 * (X / 100) - 1.5372 * (Y / 100) - 0.4986 * (Z / 100));
  const G = gammaForward(-0.9689 * (X / 100) + 1.8758 * (Y / 100) + 0.0415 * (Z / 100));
  const B = gammaForward(0.0557 * (X / 100) - 0.2040 * (Y / 100) + 1.0570 * (Z / 100));
  return [R, G, B];
}

function clamp01(x) {
  return Math.max(0, Math.min(1, x));
}

function lchToRgb255(L, C, H) {
  const [, a, b] = lchToLab(L, C, H);
  const [X, Y, Z] = labToXyz(L, a, b);
  const [r, g, bb] = xyzToRgb(X, Y, Z);
  return [Math.round(clamp01(r) * 255), Math.round(clamp01(g) * 255), Math.round(clamp01(bb) * 255)];
}

// ---------------------------------------------------------------------------
// buildDoubleVelocityColormap(minV, maxV, div, Cval_max)
//
// Reproduces velmap.py's create_cmap_doubleVelocity(minV, maxV, div,
// Cval_max) -> create_cmap_velocity(...)'s 7 (L,C,H) control points, at
// create_cmap_velocity's own DEFAULT values for every parameter
// create_cmap_doubleVelocity doesn't itself override (width=0.05,
// Lval_max=90, Lpoint_1=0.33, Lval_min=10, Cval_2=0.4). Hue control points
// are interpolated as PLAIN linear values, not circularly (matching
// make_cmap_segmented exactly -- it has no hue-wraparound special case,
// even though the H path here does jump 230->135 directly rather than the
// "short way around").
//
// Returns { colorAt(value) -> 'rgb(r,g,b)' string, or null for NaN/non-finite }.
// ---------------------------------------------------------------------------
function buildDoubleVelocityColormap(minV, maxV, div, Cval_max) {
  if (Cval_max === undefined) Cval_max = 35;
  const width = 0.05, Lpoint_1 = 0.33, Lval_max = 90, Lval_min = 10;
  const Lval_1 = 61, Lval_mid = 55, Lval_3 = 40, Lval_4 = 30;
  const Cval_1 = Cval_max, Cval_2 = 0.4, Cval_mid = 0.0;
  const Hval_L = 190, Hval_1 = 210, Hval_2 = 230, Hval_3 = 40, Hval_4 = 30, Hval_R = 10;
  const Hval_mid = (Hval_2 + Hval_3) / 2;

  const d0 = (div - minV) / (maxV - minV);
  const p2 = d0 - width / 2;
  const p3 = d0 + width / 2;
  // Lval_2 == None branch (velmap.py always passes Lval_2=None for the
  // double-velocity map): interpolated between Lval_mid and Lval_1.
  const Lval_2 = ((p2 - Lpoint_1) * Lval_mid + (d0 - p2) * Lval_1) / (d0 - Lpoint_1);

  const xs = [0, Lpoint_1, p2, d0, p3, 1 - Lpoint_1, 1];
  const Ls = [Lval_max, Lval_1, Lval_2, Lval_mid, Lval_3, Lval_4, Lval_min];
  const Cs = [Cval_max, Cval_1, Cval_2, Cval_mid, Cval_2, Cval_1, Cval_max];
  const Hs = [Hval_L, Hval_1, Hval_2, Hval_mid, Hval_3, Hval_4, Hval_R];

  function sampleAtT(t) {
    t = Math.max(0, Math.min(1, t));
    let seg = 0;
    while (seg < xs.length - 2 && t > xs[seg + 1]) seg++;
    const x0 = xs[seg], x1 = xs[seg + 1];
    const frac = x1 > x0 ? (t - x0) / (x1 - x0) : 0;
    const L = Ls[seg] + frac * (Ls[seg + 1] - Ls[seg]);
    const C = Cs[seg] + frac * (Cs[seg + 1] - Cs[seg]);
    const H = Hs[seg] + frac * (Hs[seg + 1] - Hs[seg]);
    return lchToRgb255(L, C, H);
  }

  // Small cache -- canvas rendering samples this once per pixel, and the
  // LCh->Lab->XYZ->sRGB chain isn't free; quantizing to 512 steps is far
  // finer than visually distinguishable and keeps repeated moment-map
  // renders (e.g. re-rendering after a window resize) cheap.
  const STEPS = 512;
  const cache = new Array(STEPS + 1);

  function colorAt(value) {
    if (!Number.isFinite(value)) return null;
    const t = (value - minV) / (maxV - minV || 1);
    const idx = Math.max(0, Math.min(STEPS, Math.round(Math.max(0, Math.min(1, t)) * STEPS)));
    if (!cache[idx]) cache[idx] = sampleAtT(idx / STEPS);
    const [r, g, b] = cache[idx];
    return `rgb(${r},${g},${b})`;
  }

  return { colorAt, d0, controlPoints: { xs, Ls, Cs, Hs } };
}

if (typeof module !== 'undefined') {
  module.exports = { buildDoubleVelocityColormap, lchToRgb255 };
}
if (typeof globalThis !== 'undefined') {
  globalThis.VelocityColormap = { buildDoubleVelocityColormap, lchToRgb255 };
}
