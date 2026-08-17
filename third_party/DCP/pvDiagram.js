'use strict';

// =============================================================================
// pvDiagram.js
// Renders the 4 position-velocity panels PVPlotFncs.py/BootstrapModelPlot.py
// produce (major/minor axis x data/residual), from report.pvDiagrams
// (bootstrap-realization-launcher.js -- see that file's own PV-construction
// comment for the ConstructModelBasedPVDiagram port this data comes from).
// See ~/.claude/plans/breezy-launching-nova.md Phase 3.
//
// X axis: BasePVPlot's own centering -- XMid = nSpatialPix/2 (float, NOT
// an integer floor -- matches the SAME half-index ConstructModelBasedPVDiagram
// itself used when placing flux into the array, so X=0 lines up with the
// true rotation-frame origin), X[k] = (k+1-XMid)*pixelSizeX for 0-indexed k.
// =============================================================================

// Simple linear blue-white-red diverging colormap for the residual panels --
// NOT matplotlib's exact 'coolwarm' (a real LCH/CIE-based colormap of its
// own); this is a deliberately simpler stand-in for a secondary/diagnostic
// panel, unlike Mom1's carefully-ported velocity colormap. Flagged, not
// silently approximated.
function coolwarmColor(value, limit) {
  if (!Number.isFinite(value) || limit <= 0) return null;
  const t = Math.max(-1, Math.min(1, value / limit));
  let r, g, b;
  if (t < 0) {
    const u = 1 + t; // 0 (most negative) .. 1 (zero)
    r = Math.round(60 + u * (255 - 60));
    g = Math.round(80 + u * (255 - 80));
    b = 255;
  } else {
    const u = t; // 0 (zero) .. 1 (most positive)
    r = 255;
    g = Math.round(255 - u * (255 - 60));
    b = Math.round(255 - u * (255 - 60));
  }
  return `rgb(${r},${g},${b})`;
}

function grayColorSigned(value, vmin, vmax) {
  if (!Number.isFinite(value)) return null;
  const t = Math.max(0, Math.min(1, (value - vmin) / (vmax - vmin || 1)));
  const g = Math.round(t * 255);
  return `rgb(${g},${g},${g})`;
}

// Evenly-spaced "nice" tick values (multiples of 1/2/5 x 10^n) covering
// [lo, hi] -- just enough to make the axes readable, not a full engine.
function niceTicks(lo, hi, targetCount) {
  const range = hi - lo || 1;
  const rawStep = range / targetCount;
  const mag = 10 ** Math.floor(Math.log10(rawStep));
  const norm = rawStep / mag;
  const step = (norm >= 5 ? 10 : norm >= 2 ? 5 : norm >= 1 ? 2 : 1) * mag;
  const start = Math.ceil(lo / step) * step;
  const ticks = [];
  for (let v = start; v <= hi + 1e-9; v += step) ticks.push(Math.round(v / step) * step);
  return ticks;
}

// ---------------------------------------------------------------------------
// drawPVPanel: one raster + contour + crosshair panel.
// pvFlat: flat [nSpatialPix*nChan] array (row-major, spatial-major -- matches
//   bootstrap-realization-launcher.js's pv[k*nChan+m] layout).
// modelFlat: same shape, used for the contour overlay.
// ---------------------------------------------------------------------------
function drawPVPanel(ctx, x0, y0, w, h, pvData, opts) {
  const { nSpatialPix, nChan, pixelSizeX, channelVelsKmS, modelFlat, noise, nLevels,
          contourColor, diverging, vsys, rotCurve, drawYTicks, drawXTicks, tickColor } = opts;
  const half = nSpatialPix / 2;
  const arcsecX = (k) => (k + 1 - half) * pixelSizeX;
  // channelVelsKmS (from observedDC.channels) is already km/s -- unlike the
  // Python source's CubeVels (raw m/s straight off CRVAL3/CDELT3, needing
  // BasePVPlot's own /1000.), no further conversion needed here.
  const velY = (m) => channelVelsKmS[m];

  const xMin = arcsecX(0), xMax = arcsecX(nSpatialPix - 1);
  const yMin = Math.min(...channelVelsKmS), yMax = Math.max(...channelVelsKmS);
  const toPx = (xArcsec) => x0 + ((xArcsec - xMin) / (xMax - xMin || 1)) * w;
  const toPy = (vKmS) => y0 + h - ((vKmS - yMin) / (yMax - yMin || 1)) * h;

  let vmin, vmax;
  if (diverging) {
    const limit = Math.max(1e-12, ...pvData.map(Math.abs).filter(Number.isFinite));
    vmin = -limit; vmax = limit;
  } else {
    const finite = pvData.filter(Number.isFinite);
    vmin = Math.min(...finite); vmax = Math.max(...finite);
  }

  // Raster (pcolormesh equivalent -- one filled rect per PV cell).
  const cellW = w / nSpatialPix, cellH = h / nChan;
  for (let k = 0; k < nSpatialPix; k++) {
    const px = toPx(arcsecX(k) - pixelSizeX / 2);
    for (let m = 0; m < nChan; m++) {
      const val = pvData[k * nChan + m];
      const color = diverging ? coolwarmColor(val, vmax) : grayColorSigned(val, vmin, vmax);
      if (color === null) continue;
      const py = toPy(velY(m) + (yMax - yMin) / nChan / 2);
      ctx.fillStyle = color;
      ctx.fillRect(px, py, cellW + 1, cellH + 1);
    }
  }

  // Model contours -- noise-based levels [1,3,5]*sigma (only nLevels used:
  // 3 for data panels, 1 for residual panels), line style dotted/dashed/
  // solid by level index (AddPVContoursToPlot).
  if (modelFlat && noise > 0) {
    const dashByLevel = [[2, 3], [6, 3], []]; // dotted, dashed, solid
    for (let lvlIdx = 0; lvlIdx < nLevels; lvlIdx++) {
      const level = (lvlIdx * 2 + 1) * noise; // 1,3,5 * noise
      const segs = MomentMapCanvas.marchingSquaresSegments(
        nSpatialPix, nChan, (k, m) => modelFlat[k * nChan + m], level
      );
      ctx.strokeStyle = contourColor;
      ctx.lineWidth = lvlIdx === 2 ? 2 : 1.5;
      ctx.setLineDash(dashByLevel[lvlIdx] || []);
      ctx.beginPath();
      for (const [kx0, my0, kx1, my1] of segs) {
        ctx.moveTo(toPx(arcsecX(kx0)), toPy(velY(my0)));
        ctx.lineTo(toPx(arcsecX(kx1)), toPy(velY(my1)));
      }
      ctx.stroke();
      ctx.setLineDash([]);
    }
  }

  // Crosshair at V=Vsys, X=0.
  ctx.strokeStyle = '#000';
  ctx.lineWidth = 1;
  ctx.setLineDash([2, 2]);
  ctx.beginPath();
  ctx.moveTo(toPx(0), y0); ctx.lineTo(toPx(0), y0 + h);
  ctx.moveTo(x0, toPy(vsys)); ctx.lineTo(x0 + w, toPy(vsys));
  ctx.stroke();
  ctx.setLineDash([]);

  // Rotation curve overlay (major axis panel only) -- two branches,
  // +R/-R, V = Vsys +- VROT*sin(inc) (AddSingleRC_to_PVPlot/DrawSingleRC_onPVPanel).
  if (rotCurve) {
    const { R, VROT, incRad } = rotCurve;
    ctx.strokeStyle = '#1f6fa8';
    ctx.fillStyle = '#1f6fa8';
    ctx.lineWidth = 2;
    for (const sign of [1, -1]) {
      ctx.beginPath();
      for (let i = 0; i < R.length; i++) {
        const xArcsec = sign * R[i];
        const vKmS = vsys + sign * VROT[i] * Math.sin(incRad);
        const px = toPx(xArcsec), py = toPy(vKmS);
        i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
        ctx.save(); ctx.beginPath(); ctx.arc(px, py, 2.5, 0, 2 * Math.PI); ctx.fill(); ctx.restore();
      }
      ctx.stroke();
    }
  }

  ctx.strokeStyle = tickColor || '#ccc';
  ctx.strokeRect(x0, y0, w, h);

  // Axis tick marks + numeric labels -- only on the grid's outer edges
  // (left column gets V ticks, bottom row gets X ticks), matching the same
  // left/bottom-labeling convention already used in cornerPlot.js.
  ctx.fillStyle = tickColor || '#666';
  ctx.strokeStyle = tickColor || '#666';
  ctx.font = '9px sans-serif';
  ctx.lineWidth = 1;
  if (drawYTicks) {
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (const v of niceTicks(yMin, yMax, 4)) {
      const py = toPy(v);
      ctx.beginPath();
      ctx.moveTo(x0 - 3, py); ctx.lineTo(x0, py);
      ctx.stroke();
      ctx.fillText(v.toFixed(0), x0 - 5, py);
    }
  }
  if (drawXTicks) {
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (const v of niceTicks(xMin, xMax, 5)) {
      const px = toPx(v);
      ctx.beginPath();
      ctx.moveTo(px, y0 + h); ctx.lineTo(px, y0 + h + 3);
      ctx.stroke();
      ctx.fillText(v.toFixed(0), px, y0 + h + 4);
    }
  }
}

// ---------------------------------------------------------------------------
// drawPVGrid(canvas, report) -- the full 2x2 (major/minor x data/residual)
// layout, matching AvgPVPlots/DiffModelPVPlot's own panel arrangement.
// ---------------------------------------------------------------------------
function drawPVGrid(canvas, report) {
  const pv = report.pvDiagrams;
  const panelW = 320, panelH = 170, gap = 24, labelH = 18;
  const leftMargin = 34, bottomMargin = 16;
  canvas.width = leftMargin + panelW * 2 + gap;
  canvas.height = (panelH + labelH) * 2 + gap + bottomMargin;
  const ctx = canvas.getContext('2d');
  const styles = getComputedStyle(document.body);
  const panelBg = styles.getPropertyValue('--panel-2').trim() || '#fff';
  const ink = styles.getPropertyValue('--ink').trim() || '#000';
  const border = styles.getPropertyValue('--border').trim() || '#ccc';
  ctx.fillStyle = panelBg;
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = ink;
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'left';

  const vsys = report.VSYS[0];
  const incRad = (report.INCLINATION[0] * Math.PI) / 180;
  const rotCurve = { R: report.R, VROT: report.VROT, incRad };

  const commonOpts = {
    nSpatialPix: pv.nSpatialPix, nChan: pv.nChan, pixelSizeX: pv.pixelSizeX,
    channelVelsKmS: pv.channelVelsKmS, vsys, tickColor: ink,
  };

  const nRows = 2;
  function panel(col, row, axisData, label, isResidual, includeRotCurve) {
    const x0 = leftMargin + col * (panelW + gap), y0 = row * (panelH + labelH);
    ctx.fillStyle = ink;
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText(label, x0, y0 + 12);
    drawPVPanel(ctx, x0, y0 + labelH, panelW, panelH,
      isResidual ? axisData.diff : axisData.data,
      { ...commonOpts,
        modelFlat: axisData.model,
        noise: isResidual ? axisData.noiseDiff : axisData.noiseData,
        nLevels: isResidual ? 1 : 3,
        contourColor: isResidual ? '#9c27b0' : '#dc143c',
        diverging: isResidual,
        rotCurve: includeRotCurve ? rotCurve : null,
        drawYTicks: col === 0,
        drawXTicks: row === nRows - 1,
      });
  }

  panel(0, 0, pv.major, 'Major axis -- data', false, true);
  panel(1, 0, pv.minor, 'Minor axis -- data', false, false);
  panel(0, 1, pv.major, 'Major axis -- residual (data - model)', true, false);
  panel(1, 1, pv.minor, 'Minor axis -- residual (data - model)', true, false);

  ctx.fillStyle = ink;
  ctx.font = '10px sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'alphabetic';
  ctx.fillText('X (arcsec)', leftMargin + panelW, canvas.height - 2);
  ctx.save();
  ctx.translate(10, leftMargin + panelH / 2 + labelH);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText('V (km/s)', 0, 0);
  ctx.restore();
}

if (typeof module !== 'undefined') {
  module.exports = { drawPVGrid, drawPVPanel, coolwarmColor };
}
if (typeof globalThis !== 'undefined') {
  globalThis.PVDiagram = { drawPVGrid, drawPVPanel, coolwarmColor };
}
