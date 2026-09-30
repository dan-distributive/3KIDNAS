'use strict';

// =============================================================================
// momentMapCanvas.js
// Browser-only canvas rendering of the moment-0/moment-1 maps runInitialFit's
// report now carries (report.momentMaps -- see bootstrap-realization-
// launcher.js's own comment at the point it's built, and
// ~/.claude/plans/breezy-launching-nova.md Phase 1).
//
// Deliberately the simpler "arcsec-extent" axes, not full WCS gridlines --
// see the plan file for why this is a legitimate simplification (the Python
// side has both a WCS-projected path, MakeMomPanel, and a simpler
// arcsec-extent path, MomentPlot; this ports the latter). Faithful to
// MomentMapPlotFncs.py's actual visual choices where it matters (grayscale
// Mom0 with vmin=0/vmax=1.05*max, beam drawn as a circle from BMAJ only,
// model contour at 1% of model peak, magenta center marker), not a
// from-scratch redesign.
//
// Loaded as a plain <script> tag (this project's convention -- no bundler,
// no npm). No dependencies.
// =============================================================================

// ---------------------------------------------------------------------------
// Marching squares -- single-level contour extraction on a regular grid.
// Returns an array of line segments [[x0,y0,x1,y1], ...] in GRID (pixel
// index, fractional) coordinates -- the caller maps those to canvas pixels.
// `nx`,`ny` are grid dimensions; `getVal(i,j)` reads the field (NaN treated
// as "below level", the usual marching-squares convention for masked data).
// ---------------------------------------------------------------------------
function marchingSquaresSegments(nx, ny, getVal, level) {
  const segs = [];
  function v(i, j) {
    const x = getVal(i, j);
    return Number.isFinite(x) ? x : -Infinity;
  }
  // Linear-interpolate the crossing point along a cell edge.
  function interp(a, b, va, vb) {
    const t = (level - va) / (vb - va);
    return a + t * (b - a);
  }
  for (let i = 0; i < nx - 1; i++) {
    for (let j = 0; j < ny - 1; j++) {
      const v00 = v(i, j), v10 = v(i + 1, j), v01 = v(i, j + 1), v11 = v(i + 1, j + 1);
      let idx = 0;
      if (v00 > level) idx |= 1;
      if (v10 > level) idx |= 2;
      if (v11 > level) idx |= 4;
      if (v01 > level) idx |= 8;
      if (idx === 0 || idx === 15) continue;

      const eBottom = () => [interp(i, i + 1, v00, v10), j];
      const eRight  = () => [i + 1, interp(j, j + 1, v10, v11)];
      const eTop    = () => [interp(i, i + 1, v01, v11), j + 1];
      const eLeft   = () => [i, interp(j, j + 1, v00, v01)];

      const push = (p1, p2) => segs.push([p1[0], p1[1], p2[0], p2[1]]);
      switch (idx) {
        case 1: case 14: push(eLeft(), eBottom()); break;
        case 2: case 13: push(eBottom(), eRight()); break;
        case 3: case 12: push(eLeft(), eRight()); break;
        case 4: case 11: push(eRight(), eTop()); break;
        case 6: case 9:  push(eBottom(), eTop()); break;
        case 7: case 8:  push(eLeft(), eTop()); break;
        case 5: push(eLeft(), eBottom()); push(eRight(), eTop()); break; // saddle
        case 10: push(eBottom(), eRight()); push(eLeft(), eTop()); break; // saddle
        default: break;
      }
    }
  }
  return segs;
}

// ---------------------------------------------------------------------------
// Grayscale colormap: vmin -> black, vmax -> white, clamped.
// ---------------------------------------------------------------------------
function grayColor(value, vmin, vmax) {
  if (!Number.isFinite(value)) return null; // transparent
  const t = Math.max(0, Math.min(1, (value - vmin) / (vmax - vmin || 1)));
  const g = Math.round(t * 255);
  return `rgb(${g},${g},${g})`;
}

// ---------------------------------------------------------------------------
// Draws one moment-0 panel: grayscale image, beam circle, magenta center
// marker, tan dashed model-Mom0 contour.
//
// @param canvas   an HTMLCanvasElement (its width/height attrs get set here
//                  to match the cube's own pixel aspect ratio -- caller just
//                  needs to give it a CSS box to sit in, same convention as
//                  drawChart()).
// @param momentMaps  report.momentMaps (see bootstrap-realization-launcher.js)
// @param centerPix   [x,y] pixel coords of the model center (report.XCENTER[0], report.YCENTER[0])
// ---------------------------------------------------------------------------
function drawMom0Panel(canvas, momentMaps, centerPix) {
  const { nPixX: nx, nPixY: ny, beamMajPix, observed, model } = momentMaps;
  const scale = 8; // canvas pixels per cube pixel
  canvas.width = nx * scale;
  canvas.height = ny * scale;
  const ctx = canvas.getContext('2d');
  ctx.fillStyle = '#000';
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  // FITS/astronomy convention: pixel row j=0 is the bottom of the image:
  // flip vertically when mapping cube-pixel (i,j) -> canvas (px,py).
  const toCanvasX = (i) => i * scale;
  const toCanvasY = (j) => (ny - 1 - j) * scale;

  const vmax = 1.05 * Math.max(0, ...observed.mom0.filter(Number.isFinite));
  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      const val = observed.mom0[i + j * nx];
      const color = grayColor(val, 0, vmax);
      if (color === null) continue;
      ctx.fillStyle = color;
      ctx.fillRect(toCanvasX(i), toCanvasY(j), scale, scale);
    }
  }

  // Beam circle -- MakeMomPanel draws this as a circle (BMAJ/CDELT1 only,
  // not a true ellipse) at a fixed corner-ish pixel (7,7); replicated as-is.
  const beamCenterI = 7, beamCenterJ = 7;
  ctx.strokeStyle = '#00e5ff';
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.arc(toCanvasX(beamCenterI) + scale / 2, toCanvasY(beamCenterJ) + scale / 2, (beamMajPix / 2) * scale, 0, 2 * Math.PI);
  ctx.stroke();

  // Model Mom0 outline contour -- single level at 1% of the model's own peak.
  const modelMax = Math.max(0, ...model.mom0.filter(Number.isFinite));
  const level = modelMax / 100;
  const segs = marchingSquaresSegments(nx, ny, (i, j) => model.mom0[i + j * nx], level);
  ctx.strokeStyle = '#deb887'; // burlywood/tan
  ctx.lineWidth = 2;
  ctx.setLineDash([4, 3]);
  ctx.beginPath();
  for (const [x0, y0, x1, y1] of segs) {
    ctx.moveTo(toCanvasX(x0) + scale / 2, toCanvasY(y0) + scale / 2);
    ctx.lineTo(toCanvasX(x1) + scale / 2, toCanvasY(y1) + scale / 2);
  }
  ctx.stroke();
  ctx.setLineDash([]);

  // Center marker -- magenta X.
  const [cx, cy] = centerPix;
  const mx = toCanvasX(cx) + scale / 2, my = toCanvasY(cy) + scale / 2;
  const markerSize = 8;
  ctx.strokeStyle = '#ff00ff';
  ctx.lineWidth = 3;
  ctx.beginPath();
  ctx.moveTo(mx - markerSize, my - markerSize); ctx.lineTo(mx + markerSize, my + markerSize);
  ctx.moveTo(mx - markerSize, my + markerSize); ctx.lineTo(mx + markerSize, my - markerSize);
  ctx.stroke();
}

// ---------------------------------------------------------------------------
// Draws the moment-1 (velocity field) panel: double-diverging LCh colormap
// (velocityColormap.js), black background (unmapped/NaN pixels stay black,
// matching MakeMomPanel's facecolor='black'), magenta center marker.
// Velocity contours + PA arrow are a separate, later increment (see
// ~/.claude/plans/breezy-launching-nova.md Phase 1) -- verify the base
// colormap image first.
//
// Colormap range is MODEL-derived, not data-derived (SetCMap's own choice):
// minV/maxV = VSYS +/- max(VROT)*sin(inc), NOT the observed mom1 array's own
// min/max -- get this backwards and the colorbar looks wrong for any
// off-center/noisy source (see MomentMapPlotFncs.js's header note on this
// same point, momentMapCanvas.js's own drawMom0Panel by contrast IS
// data-ranged, matching SetCMap's Mom0 case exactly).
//
// @param report  the full runInitialFit report (needs VSYS/VROT/INCLINATION
//                 in addition to momentMaps, unlike drawMom0Panel).
// ---------------------------------------------------------------------------
function drawMom1Panel(canvas, report, centerPix) {
  const momentMaps = report.momentMaps;
  const { nPixX: nx, nPixY: ny, observed } = momentMaps;
  const scale = 8;
  canvas.width = nx * scale;
  canvas.height = ny * scale;
  const ctx = canvas.getContext('2d');
  ctx.fillStyle = '#000';
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  const vsys = report.VSYS[0];
  const dV = Math.max(...report.VROT) * Math.sin((report.INCLINATION[0] * Math.PI) / 180);
  const cmap = VelocityColormap.buildDoubleVelocityColormap(vsys - dV, vsys + dV, vsys, 35);

  const toCanvasX = (i) => i * scale;
  const toCanvasY = (j) => (ny - 1 - j) * scale;

  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      const val = observed.mom1[i + j * nx];
      const color = cmap.colorAt(val);
      if (color === null) continue;
      ctx.fillStyle = color;
      ctx.fillRect(toCanvasX(i), toCanvasY(j), scale, scale);
    }
  }

  // Velocity contours -- from the MODEL cube's own Mom1 map (momentMaps.model.mom1),
  // not the observed one (AddVelContoursToMomentPlot's own choice). 11 levels,
  // VSYS + i*dV/7 for i=-5..5 (note: /7, not /5 -- NOT evenly spaced across
  // the full +-dV colormap range, an intentional asymmetry in the source,
  // replicated as-is), dashed tan, plus ONE extra thicker SOLID contour at
  // exactly VSYS on top (i=0 already draws VSYS once as part of the loop;
  // this is a deliberate second, heavier redraw at the same level, not a
  // duplicate to dedupe).
  // Flux-weighted mean velocity is numerically unstable wherever flux is
  // near zero (dividing a near-zero numerator by an even-smaller
  // denominator) -- the unmasked model cube's Mom1 has huge, meaningless
  // spikes (+-10^4 km/s) in its far field where "real" flux is essentially
  // zero but not exactly zero. NOT an issue for the colormap image (those
  // pixels just clamp to the colormap's min/max color, invisible), but
  // directly breaks marching-squares contour extraction, which reads raw
  // values everywhere -- feeding it produces a spider-web of spurious
  // crossings across the whole canvas instead of contours localized to the
  // galaxy. Fix: treat the same 1%-of-peak threshold already used for the
  // Mom0 contour as "no real flux here", and exclude those pixels from the
  // Mom1 contour source (NaN, so marching squares just skips them).
  const model = momentMaps.model;
  const modelMom0Max = Math.max(0, ...model.mom0.filter(Number.isFinite));
  const modelFluxFloor = modelMom0Max / 100;
  const modelMom1At = (i, j) => {
    const idx = i + j * nx;
    return model.mom0[idx] >= modelFluxFloor ? model.mom1[idx] : NaN;
  };
  ctx.strokeStyle = '#deb887';
  for (let k = -5; k <= 5; k++) {
    const level = vsys + (k * dV * 2) / 7; // VWidth = 2*VSinI = 2*dV
    const segs = marchingSquaresSegments(nx, ny, modelMom1At, level);
    ctx.lineWidth = 1.5;
    ctx.setLineDash([4, 3]);
    ctx.beginPath();
    for (const [x0, y0, x1, y1] of segs) {
      ctx.moveTo(toCanvasX(x0) + scale / 2, toCanvasY(y0) + scale / 2);
      ctx.lineTo(toCanvasX(x1) + scale / 2, toCanvasY(y1) + scale / 2);
    }
    ctx.stroke();
  }
  ctx.setLineDash([]);
  {
    const segs = marchingSquaresSegments(nx, ny, modelMom1At, vsys);
    ctx.lineWidth = 3.75;
    ctx.beginPath();
    for (const [x0, y0, x1, y1] of segs) {
      ctx.moveTo(toCanvasX(x0) + scale / 2, toCanvasY(y0) + scale / 2);
      ctx.lineTo(toCanvasX(x1) + scale / 2, toCanvasY(y1) + scale / 2);
    }
    ctx.stroke();
  }

  // PA arrow -- white, dotted, from the model center, length 1.1x the
  // outermost ring radius (arcsec -> pixels via pixelSizeX). AngleUse flips
  // the PA because RA increases leftward on the sky (AddArrowToMomMap).
  const [cx, cy] = centerPix;
  const arrowLenPix = (1.1 * Math.max(...report.R)) / Math.abs(momentMaps.pixelSizeX);
  const angleUse = 360 - report.POSITIONANGLE[0];
  const dXPix = -arrowLenPix * Math.cos(((angleUse + 90) * Math.PI) / 180);
  const dYPix = arrowLenPix * Math.sin(((angleUse + 90) * Math.PI) / 180);
  const arrowStartX = toCanvasX(cx) + scale / 2, arrowStartY = toCanvasY(cy) + scale / 2;
  const arrowEndX = toCanvasX(cx + dXPix) + scale / 2, arrowEndY = toCanvasY(cy + dYPix) + scale / 2;
  ctx.strokeStyle = '#fff';
  ctx.fillStyle = '#fff';
  ctx.lineWidth = 1.5;
  ctx.setLineDash([2, 3]);
  ctx.beginPath();
  ctx.moveTo(arrowStartX, arrowStartY);
  ctx.lineTo(arrowEndX, arrowEndY);
  ctx.stroke();
  ctx.setLineDash([]);
  const arrowAngle = Math.atan2(arrowEndY - arrowStartY, arrowEndX - arrowStartX);
  const headLen = 12;
  ctx.beginPath();
  ctx.moveTo(arrowEndX, arrowEndY);
  ctx.lineTo(arrowEndX - headLen * Math.cos(arrowAngle - Math.PI / 6), arrowEndY - headLen * Math.sin(arrowAngle - Math.PI / 6));
  ctx.lineTo(arrowEndX - headLen * Math.cos(arrowAngle + Math.PI / 6), arrowEndY - headLen * Math.sin(arrowAngle + Math.PI / 6));
  ctx.closePath();
  ctx.fill();

  // Center marker -- magenta X, drawn last/on top.
  const mx = toCanvasX(cx) + scale / 2, my = toCanvasY(cy) + scale / 2;
  const markerSize = 8;
  ctx.strokeStyle = '#ff00ff';
  ctx.lineWidth = 3;
  ctx.beginPath();
  ctx.moveTo(mx - markerSize, my - markerSize); ctx.lineTo(mx + markerSize, my + markerSize);
  ctx.moveTo(mx - markerSize, my + markerSize); ctx.lineTo(mx + markerSize, my - markerSize);
  ctx.stroke();
}

if (typeof module !== 'undefined') {
  module.exports = { drawMom0Panel, drawMom1Panel, marchingSquaresSegments, grayColor };
}
if (typeof globalThis !== 'undefined') {
  globalThis.MomentMapCanvas = { drawMom0Panel, drawMom1Panel, marchingSquaresSegments, grayColor };
}
