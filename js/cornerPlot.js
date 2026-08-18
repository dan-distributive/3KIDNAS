'use strict';

// =============================================================================
// cornerPlot.js
// Port of BootstrapBoxPlot.py -- despite the Python filename, this is a
// lower-triangular corner/pair-plot (diagonal = histograms, off-diagonal =
// scatter), NOT a box-and-whisker plot (no ax.boxplot() call anywhere in
// the source). See ~/.claude/plans/breezy-launching-nova.md Phase 2.
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
//
// Visual style matches tools/run_both_report.html's SVG corner plot (real
// axis tick marks/values on the outer edges, filled-bar histograms, an
// outlined circle marking the mean) -- ported onto this file's existing
// Canvas architecture rather than switching to SVG, since Canvas is already
// wired into index.html's sizing/theming and there's only ever one dataset
// here (bootstrap realizations vs. the initial fit), not run_both_report's
// multi-leg overlay. Two markers, not one: the star is the initial/best
// fit (unchanged from before), the new outlined circle is the mean across
// bootstrap realizations -- distinct concepts, both worth seeing at once.
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

function cpMean(values) {
  const finite = values.filter(Number.isFinite);
  if (!finite.length) return NaN;
  return finite.reduce((a, b) => a + b, 0) / finite.length;
}

// ---------------------------------------------------------------------------
// niceTicks/fmtTick -- port of run_both_report.html's gbpNiceTicks/
// gbpFmtTick: pick ~targetCount round-number tick values covering [lo, hi],
// then format each to just enough decimals for that step size (e.g. step=5
// -> integers, step=0.5 -> 1 decimal), so ticks read as "240, 245, 250" not
// "240.0000, 245.0000, 250.0000".
// ---------------------------------------------------------------------------
function cpNiceTicks(lo, hi, targetCount) {
  if (!(hi > lo)) return { step: 1, values: [lo] };
  const rawStep = (hi - lo) / targetCount;
  const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const norm = rawStep / mag;
  const step = norm < 1.5 ? mag : norm < 3 ? 2 * mag : norm < 7 ? 5 * mag : 10 * mag;
  const start = Math.ceil(lo / step) * step;
  const values = [];
  for (let v = start; v <= hi + step * 1e-6; v += step) values.push(Math.round(v / step) * step);
  return { step, values };
}
function cpFmtTick(v, step) {
  const decimals = Math.max(0, Math.min(3, -Math.floor(Math.log10(step) + 1e-9)));
  return v.toFixed(decimals);
}

// ---------------------------------------------------------------------------
// drawCornerPlot(canvas, bootstrapResults, bestFit)
// bootstrapResults: array of runBootstrapRealization reports (index.html's
//   own `bootstrapResults` variable -- already exactly this shape).
// bestFit: the runInitialFit report (index.html's `initialFitResult`).
// ---------------------------------------------------------------------------
function drawCornerPlot(canvas, bootstrapResults, bestFit) {
  const ok = bootstrapResults.filter((r) => !r.sofiaFailed);
  const n = PARAMS.length;
  const cell = 140;
  const margin = 46; // room for axis labels + tick marks on the outer edges
  const bottomMargin = 40; // room for X-axis tick marks/values + the axis-name label below the last row
  canvas.width = margin + n * cell;
  canvas.height = margin + n * cell + bottomMargin;
  const ctx = canvas.getContext('2d');
  const styles = getComputedStyle(document.body);
  const ink = styles.getPropertyValue('--ink').trim() || '#000';
  const muted = styles.getPropertyValue('--muted').trim() || '#888';
  const panelBg = styles.getPropertyValue('--panel-2').trim() || '#fff';
  const border = styles.getPropertyValue('--border').trim() || '#ccc';
  // On-brand palette, pulled from the page's own theme tokens (so light/
  // dark mode both work automatically) instead of one-off hardcoded hex --
  // accent is the page's primary color (buttons, step numbers, eyebrow
  // text), warn is the existing secondary/amber token, reused here rather
  // than inventing a third color just for this chart.
  const bestColor = styles.getPropertyValue('--accent').trim() || '#008da9';  // initial fit -- star + its reference line
  const bootColor = styles.getPropertyValue('--warn').trim() || '#f44336';    // bootstrap realizations -- dots, histogram, mean circle
  ctx.fillStyle = panelBg;
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  // Best-fit + mean-of-bootstraps values, per-parameter value arrays (PA
  // already circularly unwrapped relative to the best-fit value) + shared
  // axis ranges.
  const bestVals = PARAMS.map((p) => p.pick(bestFit));
  const allValues = PARAMS.map((p, k) => valuesFor(p, ok, bestVals[k]));
  const meanVals = allValues.map((vals) => cpMean(vals));
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

  // Value-axis tick marks + labels along a cell's Y edge (left side, shared
  // by every cell in the row) or X edge (bottom side, shared by every cell
  // in the column) -- only drawn once per row/column, on the outermost
  // cell, matching run_both_report.html's row===n-1/col===0 gating so the
  // grid doesn't get cluttered with repeated numbers.
  function drawYAxisTicks(x, y, h, lo, hi) {
    const ticks = cpNiceTicks(lo, hi, 3);
    ctx.strokeStyle = border;
    ctx.fillStyle = muted;
    ctx.font = '8px sans-serif';
    ctx.textAlign = 'end';
    ctx.textBaseline = 'middle';
    for (const v of ticks.values) {
      const py = y + h - ((v - lo) / (hi - lo || 1)) * h;
      if (py < y - 1 || py > y + h + 1) continue;
      ctx.beginPath();
      ctx.moveTo(x - 4, py);
      ctx.lineTo(x, py);
      ctx.stroke();
      ctx.fillText(cpFmtTick(v, ticks.step), x - 6, py);
    }
  }
  function drawXAxisTicks(x, y, w, lo, hi) {
    const ticks = cpNiceTicks(lo, hi, 3);
    ctx.strokeStyle = border;
    ctx.fillStyle = muted;
    ctx.font = '8px sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (const v of ticks.values) {
      const px = x + ((v - lo) / (hi - lo || 1)) * w;
      if (px < x - 1 || px > x + w + 1) continue;
      ctx.beginPath();
      ctx.moveTo(px, y);
      ctx.lineTo(px, y + 4);
      ctx.stroke();
      ctx.fillText(cpFmtTick(v, ticks.step), px, y + 6);
    }
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
        // Diagonal: horizontal filled-bar histogram, value on the (shared)
        // Y axis, count on X -- horizontal orientation deliberately kept
        // (rather than run_both_report.html's vertical bars) so the
        // diagonal's value axis matches the scatter panels sharing its row.
        // Silently renders empty on any error (all-NaN, single-valued
        // data, etc.), matching the Python source's own bare try/except.
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
          ctx.fillStyle = bootColor;
          ctx.globalAlpha = 0.35;
          ctx.strokeStyle = bootColor;
          ctx.lineWidth = 0.75;
          for (let b = 0; b < bins; b++) {
            const v0 = lo + (b / bins) * (hi - lo);
            const v1 = lo + ((b + 1) / bins) * (hi - lo);
            const barW = (counts[b] / maxCount) * (w - 6);
            if (barW <= 0) continue;
            const py0 = toY(v0), py1 = toY(v1);
            const barH = Math.max(Math.abs(py1 - py0) - 0.5, 0.5);
            ctx.fillRect(x, Math.min(py0, py1), barW, barH);
            ctx.strokeRect(x, Math.min(py0, py1), barW, barH);
          }
          ctx.globalAlpha = 1;
        } catch (e) { /* degrade to an empty panel, matching the Python source's own try/except: pass */ }

        // Mean-of-bootstraps -- thin dashed line, drawn under the best-fit line.
        if (Number.isFinite(meanVals[row])) {
          ctx.strokeStyle = bootColor;
          ctx.lineWidth = 1.6;
          ctx.setLineDash([4, 3]);
          ctx.beginPath();
          ctx.moveTo(x, toY(meanVals[row]));
          ctx.lineTo(x + w, toY(meanVals[row]));
          ctx.stroke();
          ctx.setLineDash([]);
        }
        // Best-fit value -- thick solid teal line, on top.
        ctx.strokeStyle = bestColor;
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
        ctx.fillStyle = bootColor;
        ctx.globalAlpha = 0.55;
        for (let k = 0; k < ok.length; k++) {
          const xv = allValues[col][k], yv = allValues[row][k];
          if (!Number.isFinite(xv) || !Number.isFinite(yv)) continue;
          ctx.beginPath();
          ctx.arc(toX(xv), toY(yv), 4, 0, 2 * Math.PI);
          ctx.fill();
        }
        ctx.globalAlpha = 1;

        // Mean-of-bootstraps -- outlined circle, drawn under the star.
        if (Number.isFinite(meanVals[col]) && Number.isFinite(meanVals[row])) {
          cpDrawMeanCircle(ctx, toX(meanVals[col]), toY(meanVals[row]), 6, bootColor, panelBg);
        }
        // Best-fit point -- large teal star, drawn on top of everything.
        drawStar(ctx, toX(bestVals[col]), toY(bestVals[row]), 9, bestColor);
      }
      ctx.restore();

      // Tick marks/labels sit just outside the cell's own rect (below for
      // X, left for Y) -- must be drawn after ctx.restore() ends the clip
      // above, or they're clipped away along with everything outside the
      // cell bounds. Y-axis (value) ticks apply to any cell in the first
      // column, diagonal included -- the histogram's Y-axis is the same
      // value scale as its row's scatter panels. X-axis ticks only apply
      // to off-diagonal (scatter) cells in the last row -- the diagonal's
      // X-axis is bin count, not a value scale, so it isn't tick-labelled.
      if (row === n - 1 && row !== col) drawXAxisTicks(x, y + h, w, ranges[col][0], ranges[col][1]);
      if (col === 0) drawYAxisTicks(x, y, h, yLo, yHi);
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
    ctx.translate(margin - 22, leftCell.y + leftCell.h / 2);
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
    ctx.fillText(PARAMS[col].label, bottomCell.x + bottomCell.w / 2, bottomCell.y + bottomCell.h + 18);
  }

  // Legend: star = initial fit, outlined circle = mean of bootstrap realizations.
  const legendY = margin - 30;
  ctx.font = '11px sans-serif';
  ctx.textBaseline = 'middle';
  ctx.textAlign = 'left';
  drawStar(ctx, margin + 8, legendY, 7, bestColor);
  ctx.fillStyle = ink;
  ctx.fillText('Initial fit', margin + 20, legendY);
  cpDrawMeanCircle(ctx, margin + 130, legendY, 5, bootColor, panelBg);
  ctx.fillStyle = ink;
  ctx.fillText('Mean of bootstraps', margin + 142, legendY);
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

// Outlined circle -- fill in `color`, stroke in `outline` (the surrounding
// panel background) for contrast against overlapping bootstrap dots,
// matching run_both_report.html's .gbp-mean { stroke: var(--surface) }.
function cpDrawMeanCircle(ctx, cx, cy, r, color, outline) {
  ctx.beginPath();
  ctx.arc(cx, cy, r, 0, 2 * Math.PI);
  ctx.fillStyle = color;
  ctx.fill();
  ctx.strokeStyle = outline;
  ctx.lineWidth = 1.5;
  ctx.stroke();
}

if (typeof module !== 'undefined') {
  module.exports = { drawCornerPlot, adjustPA, PARAMS };
}
if (typeof globalThis !== 'undefined') {
  globalThis.CornerPlot = { drawCornerPlot, adjustPA, PARAMS };
}
