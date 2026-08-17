'use strict';

// =============================================================================
// radialProfilePlot.js
// The "R plots" from the bootstrap comparison tooling (run_all_three_results.
// html's gbpBuildProfile / "Rotation curve" & "Surface density" cards),
// adapted for galaxy-fit.html's single-galaxy bootstrap results: one V_rot-
// vs-R and one Sigma-vs-R panel, each showing every bootstrap realization's
// own radial profile as a thin translucent line (the spread IS the point --
// unlike run_all_three's multi-run comparison, there's no separate
// median/MAD summary needed on top of it) plus the initial fit's profile as
// a thick accent-colored line, matching the reference GeoBoxPlot.png's own
// upper-corner RC/Sigma spread panels. Styled to match cornerPlot.js (same
// translucent-red-for-realizations / teal-for-best-fit convention) rather
// than run_all_three's SVG styling, to stay consistent with this page's own
// look.
// =============================================================================

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

function drawRadialProfilePanel(ctx, x0, y0, w, h, realizations, bestFit, xKey, yKey, yLabel, colors) {
  const bestX = bestFit[xKey], bestY = bestFit[yKey];
  const allX = [], allY = [0];
  for (const r of realizations) {
    if (r[xKey]) allX.push(...r[xKey]);
    if (r[yKey]) allY.push(...r[yKey].filter(Number.isFinite));
  }
  if (bestX) allX.push(...bestX);
  if (bestY) allY.push(...bestY.filter(Number.isFinite));

  const xMin = 0, xMax = Math.max(1, ...allX) * 1.05;
  const yMin = Math.min(0, ...allY), yMax = Math.max(1, ...allY) * 1.08;
  const toPx = (x) => x0 + ((x - xMin) / (xMax - xMin || 1)) * w;
  const toPy = (y) => y0 + h - ((y - yMin) / (yMax - yMin || 1)) * h;

  // Per-realization profiles -- thin translucent lines + small dots at each
  // ring, same red used for the corner plot's own scatter points.
  ctx.strokeStyle = colors.realization;
  ctx.fillStyle = colors.realization;
  ctx.lineWidth = 1;
  for (const r of realizations) {
    const X = r[xKey], Y = r[yKey];
    if (!X || !Y || X.length < 2) continue;
    ctx.beginPath();
    for (let i = 0; i < X.length; i++) {
      const px = toPx(X[i]), py = toPy(Y[i]);
      i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
    }
    ctx.stroke();
  }

  // Best (initial) fit -- thick accent line + dots, drawn last/on top.
  if (bestX && bestY) {
    ctx.strokeStyle = colors.best;
    ctx.fillStyle = colors.best;
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    for (let i = 0; i < bestX.length; i++) {
      const px = toPx(bestX[i]), py = toPy(bestY[i]);
      i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
    }
    ctx.stroke();
    for (let i = 0; i < bestX.length; i++) {
      if (!Number.isFinite(bestY[i])) continue;
      ctx.beginPath();
      ctx.arc(toPx(bestX[i]), toPy(bestY[i]), 3.5, 0, 2 * Math.PI);
      ctx.fill();
    }
  }

  // Axis ticks + labels.
  ctx.fillStyle = colors.ink;
  ctx.strokeStyle = colors.ink;
  ctx.font = '9px sans-serif';
  ctx.lineWidth = 1;
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  for (const v of niceTicks(yMin, yMax, 4)) {
    const py = toPy(v);
    ctx.beginPath(); ctx.moveTo(x0 - 3, py); ctx.lineTo(x0, py); ctx.stroke();
    ctx.fillText(v.toFixed(v % 1 === 0 ? 0 : 1), x0 - 5, py);
  }
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  for (const v of niceTicks(xMin, xMax, 5)) {
    const px = toPx(v);
    ctx.beginPath(); ctx.moveTo(px, y0 + h); ctx.lineTo(px, y0 + h + 3); ctx.stroke();
    ctx.fillText(v.toFixed(0), px, y0 + h + 4);
  }

  ctx.strokeStyle = colors.border;
  ctx.strokeRect(x0, y0, w, h);

  ctx.fillStyle = colors.ink;
  ctx.font = '10px sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'alphabetic';
  ctx.fillText('R (arcsec)', x0 + w / 2, y0 + h + 16);
  ctx.save();
  ctx.translate(x0 - 30, y0 + h / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(yLabel, 0, 0);
  ctx.restore();
}

// ---------------------------------------------------------------------------
// drawRadialProfileGrid(canvas, bootstrapResults, bestFit) -- V_rot and
// Sigma panels side by side.
// ---------------------------------------------------------------------------
function drawRadialProfileGrid(canvas, bootstrapResults, bestFit) {
  const ok = bootstrapResults.filter((r) => !r.sofiaFailed);
  const panelW = 380, panelH = 220, gap = 30, leftMargin = 40, topLabelH = 18, bottomMargin = 24;
  canvas.width = leftMargin + panelW * 2 + gap;
  canvas.height = topLabelH + panelH + bottomMargin;
  const ctx = canvas.getContext('2d');
  const s = getComputedStyle(document.body);
  const colors = {
    ink: s.getPropertyValue('--ink').trim() || '#000',
    border: s.getPropertyValue('--border').trim() || '#ccc',
    best: s.getPropertyValue('--accent').trim() || '#008da9',
    realization: 'rgba(244,67,54,0.35)',
  };
  ctx.fillStyle = s.getPropertyValue('--panel-2').trim() || '#fff';
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  function panel(col, xKey, yKey, label, yLabel) {
    const x0 = leftMargin + col * (panelW + gap), y0 = topLabelH;
    ctx.fillStyle = colors.ink;
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText(label, x0, 12);
    drawRadialProfilePanel(ctx, x0, y0, panelW, panelH, ok, bestFit, xKey, yKey, yLabel, colors);
  }

  panel(0, 'R', 'VROT', 'Rotation curve', 'V (km/s)');
  panel(1, 'R_SD', 'SURFDENS', 'Surface density', 'Σ (M☉/pc²)');
}

if (typeof module !== 'undefined') {
  module.exports = { drawRadialProfileGrid, drawRadialProfilePanel };
}
if (typeof globalThis !== 'undefined') {
  globalThis.RadialProfilePlot = { drawRadialProfileGrid, drawRadialProfilePanel };
}
