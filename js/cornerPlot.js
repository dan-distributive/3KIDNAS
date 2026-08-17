'use strict';

// =============================================================================
// cornerPlot.js
// Port of BootstrapBoxPlot.py -- despite the Python filename, this is a
// lower-triangular corner/pair-plot (diagonal = step histograms,
// off-diagonal = scatter), NOT a box-and-whisker plot (no ax.boxplot() call
// anywhere in the source). See ~/.claude/plans/breezy-launching-nova.md
// Phase 2.
//
// Ships with all 8 parameters the reference GeoBoxPlot.png shows (XCENTER,
// YCENTER, INCLINATION, POSITIONANGLE, VSYS, VDISP, RHI, VHI) -- RHI/VHI
// computed per-realization via scalingParams.js's extractVHIForModel
// (Phase 4), the same per-model extraction computeScalingParams' own
// bootstrap-uncertainty loop uses internally.
//
// The Python source's own panel-placement arithmetic is a figure-fraction,
// Y-up-coordinate, per-column running-j-index scheme (BootstrapBoxPlot.py:
// 49-105) -- verified by hand that it's equivalent to the much simpler
// standard "lower-triangular grid, row >= col, diagonal at row===col"
// formulation used directly here for canvas's Y-down coordinates; the
// resulting arrangement (which parameter pairs land in which cell) is
// identical, just derived more directly.
// =============================================================================

// extractVHIForModel does a small linear scan per call (5-ish rings) -- run
// once per model and cached on the model object itself, so the RHI and VHI
// picks below (and any repeated draw) don't redo the extraction.
function scalingParamsFor(r) {
  if (!r.__scalingParams) r.__scalingParams = ScalingParams.extractVHIForModel(r);
  return r.__scalingParams;
}

const PARAMS = [
  { key: 'INCLINATION', label: 'Inc (deg)', pick: (r) => r.INCLINATION[0] },
  { key: 'POSITIONANGLE', label: 'PA (deg)', pick: (r) => r.POSITIONANGLE[0] },
  { key: 'VSYS', label: 'V_sys (km/s)', pick: (r) => r.VSYS[0] },
  { key: 'VDISP', label: 'V_disp (km/s)', pick: (r) => r.VDISP[0] },
  { key: 'XCENTER', label: 'X (pix)', pick: (r) => r.XCENTER[0] },
  { key: 'YCENTER', label: 'Y (pix)', pick: (r) => r.YCENTER[0] },
  { key: 'RHI', label: 'RHI (arcsec)', pick: (r) => scalingParamsFor(r).rhi },
  { key: 'VHI', label: 'VHI (km/s)', pick: (r) => scalingParamsFor(r).vhi },
];

// ---------------------------------------------------------------------------
// AdjustPA -- circular-wrap a bootstrap PA value into the same "winding" as
// the best-fit value, so a realization that landed just across the 0/360
// boundary doesn't show as a huge, physically-meaningless spread. Applied
// per-point, using the best-fit value as the reference every time (not a
// single global unwrap). Port of BootstrapBoxPlot.py:294-302, verbatim.
// ---------------------------------------------------------------------------
function adjustPA(pa, paBest) {
  const diff = paBest - pa;
  if (diff <= -180) return pa - 360;
  if (diff >= 180) return pa + 360;
  return pa;
}

function valuesFor(param, realizations, bestVal) {
  return realizations.map((r) => {
    const v = param.pick(r);
    return param.key === 'POSITIONANGLE' ? adjustPA(v, bestVal) : v;
  });
}

// ---------------------------------------------------------------------------
// drawCornerPlot(canvas, bootstrapResults, bestFit)
// bootstrapResults: array of runBootstrapRealization reports (galaxy-fit.html's
//   own `bootstrapResults` variable -- already exactly this shape).
// bestFit: the runInitialFit report (galaxy-fit.html's `initialFitResult`).
// ---------------------------------------------------------------------------
function drawCornerPlot(canvas, bootstrapResults, bestFit) {
  const ok = bootstrapResults.filter((r) => !r.sofiaFailed);
  const n = PARAMS.length;
  const cell = 140;
  const margin = 46; // room for axis labels on the outer edges
  canvas.width = margin + n * cell;
  canvas.height = margin + n * cell;
  const ctx = canvas.getContext('2d');
  const styles = getComputedStyle(document.body);
  const ink = styles.getPropertyValue('--ink').trim() || '#000';
  const panelBg = styles.getPropertyValue('--panel-2').trim() || '#fff';
  const border = styles.getPropertyValue('--border').trim() || '#ccc';
  ctx.fillStyle = panelBg;
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  // Best-fit values + per-parameter value arrays (PA already circularly
  // unwrapped relative to its own best-fit value) + shared axis ranges.
  const bestVals = PARAMS.map((p) => p.pick(bestFit));
  const allValues = PARAMS.map((p, k) => valuesFor(p, ok, bestVals[k]));
  const ranges = PARAMS.map((p, k) => {
    const vals = allValues[k].filter(Number.isFinite);
    let lo = Math.min(bestVals[k], ...vals);
    let hi = Math.max(bestVals[k], ...vals);
    if (!(hi > lo)) { lo -= 1; hi += 1; } // degenerate (all-equal) guard
    const pad = (hi - lo) * 0.1 || 1;
    return [lo - pad, hi + pad];
  });

  function cellRect(row, col) {
    return { x: margin + col * cell, y: margin + row * cell, w: cell, h: cell };
  }

  for (let row = 0; row < n; row++) {
    for (let col = 0; col <= row; col++) {
      const { x, y, w, h } = cellRect(row, col);
      ctx.strokeStyle = border;
      ctx.strokeRect(x, y, w, h);
      ctx.save();
      ctx.beginPath();
      ctx.rect(x, y, w, h);
      ctx.clip();

      const [yLo, yHi] = ranges[row];
      const toY = (v) => y + h - ((v - yLo) / (yHi - yLo || 1)) * h;

      if (row === col) {
        // Diagonal: horizontal step histogram, value on the (shared) Y
        // axis, count on X -- matches HistPlt's own horizontal orientation,
        // chosen there so the diagonal keeps the same Y-scale as the
        // scatter panels in its row. Silently renders empty on any error
        // (all-NaN, single-valued data, etc.) -- matches HistPlt's own
        // bare try/except: pass.
        try {
          const vals = allValues[row].filter(Number.isFinite);
          const bins = 20;
          const [lo, hi] = ranges[row];
          const counts = new Array(bins).fill(0);
          for (const v of vals) {
            let b = Math.floor(((v - lo) / (hi - lo || 1)) * bins);
            b = Math.max(0, Math.min(bins - 1, b));
            counts[b]++;
          }
          const maxCount = Math.max(1, ...counts);
          ctx.strokeStyle = '#f44336';
          ctx.lineWidth = 2;
          ctx.beginPath();
          for (let b = 0; b < bins; b++) {
            const v0 = lo + (b / bins) * (hi - lo);
            const v1 = lo + ((b + 1) / bins) * (hi - lo);
            const barW = (counts[b] / maxCount) * (w - 6);
            const py0 = toY(v0), py1 = toY(v1);
            ctx.moveTo(x, py0);
            ctx.lineTo(x + barW, py0);
            ctx.lineTo(x + barW, py1);
            ctx.lineTo(x, py1);
          }
          ctx.stroke();
        } catch (e) { /* degrade to an empty panel, matching HistPlt's own try/except: pass */ }

        // Best-fit value -- thick horizontal teal line.
        ctx.strokeStyle = '#008da9';
        ctx.lineWidth = 4;
        ctx.beginPath();
        ctx.moveTo(x, toY(bestVals[row]));
        ctx.lineTo(x + w, toY(bestVals[row]));
        ctx.stroke();
      } else {
        const [xLo, xHi] = ranges[col];
        const toX = (v) => x + ((v - xLo) / (xHi - xLo || 1)) * w;
        // Bootstrap realizations -- translucent red dots (density-by-overlap
        // is how spread reads visually, not a KDE/contour).
        ctx.fillStyle = 'rgba(244,67,54,0.5)';
        for (let k = 0; k < ok.length; k++) {
          const xv = allValues[col][k], yv = allValues[row][k];
          if (!Number.isFinite(xv) || !Number.isFinite(yv)) continue;
          ctx.beginPath();
          ctx.arc(toX(xv), toY(yv), 4, 0, 2 * Math.PI);
          ctx.fill();
        }
        // Best-fit point -- large teal star, drawn on top.
        drawStar(ctx, toX(bestVals[col]), toY(bestVals[row]), 9, '#008da9');
      }
      ctx.restore();
    }
    // Left-column Y label + bottom-row X label (PanelFmt/HistFmt's own
    // row===0/col===0 tick-label-suppression logic, simplified to just
    // showing the axis label once per row/column rather than replicating
    // full major-tick suppression).
    ctx.fillStyle = ink;
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    const leftCell = cellRect(row, 0);
    ctx.save();
    ctx.translate(margin - 6, leftCell.y + leftCell.h / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center';
    ctx.fillText(PARAMS[row].label, 0, 0);
    ctx.restore();
  }
  ctx.fillStyle = ink;
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  for (let col = 0; col < n; col++) {
    const bottomCell = cellRect(n - 1, col);
    ctx.fillText(PARAMS[col].label, bottomCell.x + bottomCell.w / 2, bottomCell.y + bottomCell.h + 4);
  }
}

function drawStar(ctx, cx, cy, r, color) {
  ctx.fillStyle = color;
  ctx.beginPath();
  for (let i = 0; i < 10; i++) {
    const ang = (Math.PI / 5) * i - Math.PI / 2;
    const rad = i % 2 === 0 ? r : r * 0.4;
    const px = cx + rad * Math.cos(ang), py = cy + rad * Math.sin(ang);
    i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
  }
  ctx.closePath();
  ctx.fill();
}

if (typeof module !== 'undefined') {
  module.exports = { drawCornerPlot, adjustPA, PARAMS };
}
if (typeof globalThis !== 'undefined') {
  globalThis.CornerPlot = { drawCornerPlot, adjustPA, PARAMS };
}
