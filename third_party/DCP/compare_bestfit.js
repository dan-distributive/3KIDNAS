'use strict';
// =============================================================================
// compare_bestfit.js
//
// Field-by-field %-difference table between two runs' *_BSModel.txt --
// the unified best-fit-plus-bootstrap-error summary BOTH UseDCP=0
// (fortran-local) and UseDCP=1 (js-dcp) write, via the SAME shared Python
// post-processing (Bootstrap_Error_Analysis/ExtractScalingParams/
// GeometryCorrection/BootstrapModelPlot -- only the anchor fit itself
// differs between the two paths, everything downstream is identical code).
// This is the "parameters that come out" comparison -- not just chi2 --
// covering geometry, the full rotation-curve/surface-density radial
// profiles, and the derived RHI/VHI scaling parameters, each with its own
// bootstrap-derived error bar from both sides.
//
// Usage:
//   node compare_bestfit.js <fortran_dir> <js_dcp_dir> [--md out.md]
//   (each <dir> is a leg's ObjName folder, e.g.
//    ../../3KIDNASTests/SingleGalaxyTest/TestFits_RunAllThree_FortranLocal/WALLABY_J103538-484832)
// =============================================================================

const fs = require('fs');
const path = require('path');

function arg(name, def) {
  const i = process.argv.indexOf('--' + name);
  return i >= 0 ? process.argv[i + 1] : def;
}

const [, , fortranDir, jsDir] = process.argv;
if (!fortranDir || !jsDir) {
  console.error('Usage: node compare_bestfit.js <fortran_leg_dir> <js_dcp_leg_dir> [--md out.md]');
  process.exit(2);
}

function findBSModelTxt(dir) {
  const files = fs.readdirSync(dir).filter((f) => f.endsWith('_BSModel.txt'));
  if (files.length !== 1) throw new Error(`expected exactly one *_BSModel.txt in ${dir}, found ${files.length}`);
  return path.join(dir, files[0]);
}

// ---------------------------------------------------------------------------
// parseBSModel -- tab/space-tolerant line-oriented parser for the format
// shown above (Geometry Parameters table, Rotation Curve block, Surface
// Density Profile block, RHI/VHI blocks). Deliberately simple/line-anchored
// rather than a generic tabular parser -- this file's own layout is fixed
// and hand-written by BootstrapModelPlot.py's WriteBSModelFile, not a
// format that needs to tolerate arbitrary reordering.
// ---------------------------------------------------------------------------
function parseBSModel(text) {
  const lines = text.split('\n');
  const out = { rings: [] };

  function numsOnLine(line) {
    return (line.match(/-?\d+\.?\d*(?:[eE][+-]?\d+)?/g) || []).map(Number);
  }

  const geomFieldOrder = ['X', 'Y', 'RA', 'DEC', 'Inc', 'PA', 'PA_g', 'VSys', 'VDisp'];
  let geomIdx = 0;
  let section = null;
  let ringSection = null; // 'rc' | 'sd'

  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    const trimmed = line.trim();

    if (/^RMS \(mJy\/beam\)/.test(trimmed)) { out.RMS = numsOnLine(trimmed)[0]; continue; }
    if (/^SN_Int\b/.test(trimmed)) { out.SN_Int = numsOnLine(trimmed)[0]; continue; }
    if (/^SN_Peak\b/.test(trimmed)) { out.SN_Peak = numsOnLine(trimmed)[0]; continue; }
    if (/^SN_Avg\b/.test(trimmed)) { out.SN_Avg = numsOnLine(trimmed)[0]; continue; }
    if (/^SN_Median\b/.test(trimmed)) { out.SN_Median = numsOnLine(trimmed)[0]; continue; }

    if (/^Geometry Parameters/.test(trimmed)) { section = 'geom'; geomIdx = 0; continue; }
    if (section === 'geom' && /_model/.test(trimmed)) {
      const nums = numsOnLine(trimmed);
      const key = geomFieldOrder[geomIdx++];
      if (key && nums.length >= 1) { out[key] = nums[0]; out[key + '_err'] = nums[1]; }
      if (geomIdx >= geomFieldOrder.length) section = null;
      continue;
    }

    if (/^Rotation Curve/.test(trimmed)) { ringSection = 'rc'; continue; }
    if (/^Surface Density Profile/.test(trimmed)) { ringSection = 'sd'; continue; }
    if (ringSection && /^nR=/.test(trimmed)) continue;
    if (ringSection && /^Rad/.test(trimmed)) continue;
    if (ringSection && /^\(''\)/.test(trimmed)) continue;
    if (ringSection && trimmed === '') { ringSection = null; continue; }
    if (ringSection === 'rc') {
      const nums = numsOnLine(trimmed);
      if (nums.length >= 3) out.rings.push({ R: nums[0], VROT: nums[1], VROT_err: nums[2] });
      continue;
    }
    if (ringSection === 'sd') {
      const nums = numsOnLine(trimmed);
      if (nums.length >= 3) {
        const idx = out.rings.findIndex((r) => Math.abs(r.R - nums[0]) < 1e-6);
        if (idx >= 0) { out.rings[idx].SD = nums[1]; out.rings[idx].SD_err = nums[2]; }
      }
      continue;
    }

    if (/^RHI and limits \(arcsec\)/.test(trimmed)) {
      const nums = numsOnLine(lines[i + 1]);
      [out.RHI_as, out.RHI_as_lo, out.RHI_as_hi] = nums;
      continue;
    }
    if (/^RHI and limits \(kpc\)/.test(trimmed)) {
      const nums = numsOnLine(lines[i + 1]);
      [out.RHI_kpc, out.RHI_kpc_lo, out.RHI_kpc_hi] = nums;
      continue;
    }
    if (/^VHI and error \(km\/s\)/.test(trimmed)) {
      const nums = numsOnLine(lines[i + 1]);
      [out.VHI, out.VHI_err] = nums;
      continue;
    }
  }
  return out;
}

function pctDiff(a, b) {
  if (!Number.isFinite(a) || !Number.isFinite(b)) return null;
  const denom = Math.max(Math.abs(a), Math.abs(b), 1e-12);
  return (Math.abs(a - b) / denom) * 100;
}

function fmt(v, digits = 4) {
  return Number.isFinite(v) ? v.toFixed(digits) : '--';
}

const fortranModel = parseBSModel(fs.readFileSync(findBSModelTxt(fortranDir), 'utf8'));
const jsModel = parseBSModel(fs.readFileSync(findBSModelTxt(jsDir), 'utf8'));

const SCALAR_FIELDS = [
  ['RMS (mJy/beam)', 'RMS'], ['SN_Integrated', 'SN_Int'], ['SN_Peak', 'SN_Peak'],
  ['SN_Avg', 'SN_Avg'], ['SN_Median', 'SN_Median'],
  ['X (pix)', 'X'], ['Y (pix)', 'Y'], ['RA (deg)', 'RA'], ['DEC (deg)', 'DEC'],
  ['Inclination (deg)', 'Inc'], ['PA (deg)', 'PA'], ['PA global (deg)', 'PA_g'],
  ['V_sys (km/s)', 'VSys'], ['V_disp (km/s)', 'VDisp'],
  ['R_HI (arcsec)', 'RHI_as'], ['R_HI (kpc)', 'RHI_kpc'], ['V_HI (km/s)', 'VHI'],
];

// Bootstrap-derived error bars -- unlike the central values above (fully
// deterministic given the same cube/mask/PA/Inc estimate, so identical
// across separate runs regardless of bootstrap seed/count), these DO
// depend on which realizations this specific run happened to draw --
// the part actually exercised fresh by a new dispatch.
const ERROR_FIELDS = [
  ['X error (pix)', 'X_err'], ['Y error (pix)', 'Y_err'],
  ['Inclination error (deg)', 'Inc_err'], ['PA error (deg)', 'PA_err'],
  ['V_sys error (km/s)', 'VSys_err'], ['V_disp error (km/s)', 'VDisp_err'],
];

const rows = [];
for (const [label, key] of SCALAR_FIELDS) {
  const a = fortranModel[key], b = jsModel[key];
  rows.push({ label, fortran: a, js: b, pct: pctDiff(a, b) });
}
for (const [label, key] of ERROR_FIELDS) {
  const a = fortranModel[key], b = jsModel[key];
  rows.push({ label, fortran: a, js: b, pct: pctDiff(a, b) });
}
for (let i = 0; i < Math.max(fortranModel.rings.length, jsModel.rings.length); i++) {
  const fr = fortranModel.rings[i], jr = jsModel.rings[i];
  const R = fr ? fr.R : jr ? jr.R : NaN;
  rows.push({ label: `V_rot @ R=${fmt(R, 2)}" (km/s)`, fortran: fr && fr.VROT, js: jr && jr.VROT, pct: pctDiff(fr && fr.VROT, jr && jr.VROT) });
}
for (let i = 0; i < Math.max(fortranModel.rings.length, jsModel.rings.length); i++) {
  const fr = fortranModel.rings[i], jr = jsModel.rings[i];
  const R = fr ? fr.R : jr ? jr.R : NaN;
  rows.push({ label: `Sigma @ R=${fmt(R, 2)}" (Msol/pc^2)`, fortran: fr && fr.SD, js: jr && jr.SD, pct: pctDiff(fr && fr.SD, jr && jr.SD) });
}

const colW = { label: Math.max(...rows.map((r) => r.label.length), 20), num: 14 };
function padL(s, w) { s = String(s); return s.length >= w ? s : ' '.repeat(w - s.length) + s; }
function padR(s, w) { s = String(s); return s.length >= w ? s : s + ' '.repeat(w - s.length); }

console.log(`\nFortran (${fortranDir})  vs  JS/DCP (${jsDir})\n`);
console.log(padR('Parameter', colW.label) + padL('Fortran', colW.num) + padL('JS/DCP', colW.num) + padL('% diff', colW.num));
console.log('-'.repeat(colW.label + colW.num * 3));
let maxPct = 0;
for (const r of rows) {
  const pctStr = r.pct == null ? '--' : (r.pct < 0.01 ? '<0.01%' : r.pct.toFixed(3) + '%');
  if (r.pct != null) maxPct = Math.max(maxPct, r.pct);
  console.log(padR(r.label, colW.label) + padL(fmt(r.fortran), colW.num) + padL(fmt(r.js), colW.num) + padL(pctStr, colW.num));
}
console.log('-'.repeat(colW.label + colW.num * 3));
console.log(`Largest % difference across all ${rows.length} compared parameters: ${maxPct.toFixed(3)}%`);

const mdPath = arg('md');
if (mdPath) {
  const md = [
    `| Parameter | Fortran | JS/DCP | % diff |`,
    `|---|---:|---:|---:|`,
    ...rows.map((r) => `| ${r.label} | ${fmt(r.fortran)} | ${fmt(r.js)} | ${r.pct == null ? '--' : r.pct.toFixed(3) + '%'} |`),
  ].join('\n');
  fs.writeFileSync(mdPath, md);
  console.log(`\nMarkdown table written to ${mdPath}`);
}
