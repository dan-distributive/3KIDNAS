/**
 * build_partial_report.js -- generates a run_both-shaped report JSON
 * purely by READING whatever's already on disk for each leg, without
 * running or wiping anything.
 *
 * Exists specifically for "one leg is done, the other is still running,
 * I want to look at what's finished so far" -- run_both.js itself can't
 * safely be used for that: its main() unconditionally rewrites both
 * run_both_*_config.py files (including the one a still-running leg's
 * Python process may reference) before it even looks at any --skip flag,
 * every single invocation. This script never touches those config files,
 * never wipes a folder, never spawns a subprocess -- it only reads each
 * leg's *_BootstrapFits.csv/*_BootstrapTimings.json if present. A leg
 * that's still running simply won't have those files yet (they're only
 * written at the very end, after all realizations + aggregation finish),
 * so it's correctly omitted rather than read partially/corrupted.
 *
 * Usage:
 *   node build_partial_report.js [--json PATH]
 */
'use strict';
const path = require('node:path');
const fs = require('node:fs');

const TEST_DIR = path.join(__dirname, '..', '..', '3KIDNASTests', 'SingleGalaxyTest');
const GALAXY = { ObjName: 'WALLABY_J103538-484832' };
const folders = {
  fortranLocal: 'TestFits_RunAllThree_FortranLocal',
  jsDcp: 'TestFits_RunAllThree_JSDcp',
};


// ---- Copied verbatim from run_both.js (kept in sync by hand -- small,
//      stable, read-only helpers) ----
function parseCsvLine(line) {
  const fields = [];
  let cur = '';
  let inQuotes = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (c === '"') { inQuotes = !inQuotes; continue; }
    if (c === ',' && !inQuotes) { fields.push(cur); cur = ''; continue; }
    cur += c;
  }
  fields.push(cur);
  return fields;
}
function readBootstrapCsv(objFolder) {
  const csvPath = path.join(objFolder, `${GALAXY.ObjName}_BootstrapFits.csv`);
  if (!fs.existsSync(csvPath)) return null;
  const lines = fs.readFileSync(csvPath, 'utf8').trim().split('\n');
  const header = parseCsvLine(lines[0]);
  return lines.slice(1).map((line) => {
    const fields = parseCsvLine(line);
    const row = {};
    header.forEach((h, i) => { row[h] = fields[i]; });
    return row;
  });
}
function readTimingsJson(objFolder) {
  const jsonPath = path.join(objFolder, `${GALAXY.ObjName}_BootstrapTimings.json`);
  if (!fs.existsSync(jsonPath)) return null;
  return JSON.parse(fs.readFileSync(jsonPath, 'utf8'));
}
// Reproducibility metadata (seed, nBootstraps, CPU count, cloud density)
// written by run_both.js into the leg's own folder -- null (not an error)
// for data predating this file, or for a leg run before this feature
// existed.
function readRunMeta(objFolder) {
  const p = path.join(objFolder, `${GALAXY.ObjName}_RunMeta.json`);
  if (!fs.existsSync(p)) return null;
  return JSON.parse(fs.readFileSync(p, 'utf8'));
}
// This script has no live process to ask for a real wall-clock duration --
// a leg's results may come from an entirely separate, already-finished
// invocation. BootstrapTimings.json carries no timestamps, only durations,
// so approximate wall time as the spread between the earliest and latest
// file mtimes in the leg's output folder (first fixture-fit output through
// the final BootstrapFits.csv/PNG write). Verified against a real run: e.g.
// TestFits_RunAllThree_JSDcp's earliest file (07:59:04) to its final CSV
// (08:58:05) brackets that leg's actual ~59min run closely.
function estimateWallSeconds(objFolder) {
  const files = fs.readdirSync(objFolder)
    .map((f) => path.join(objFolder, f))
    .filter((f) => fs.statSync(f).isFile());
  if (files.length < 2) return null;
  const mtimes = files.map((f) => fs.statSync(f).mtimeMs);
  return (Math.max(...mtimes) - Math.min(...mtimes)) / 1000;
}
const SCALAR_FIELDS = ['X_model', 'Y_model', 'Inc_model', 'PA_model', 'Vsys_model',
  'RA_model', 'DEC_model', 'Vdisp_model', 'RHI_AS', 'VHI'];
function compareBootstraps(rowsA, rowsB) {
  if (!rowsA || !rowsB) return { skipped: true, fields: [] };
  const n = Math.min(rowsA.length, rowsB.length);
  const rowCountMismatch = rowsA.length !== rowsB.length;
  const fields = [];
  for (const field of SCALAR_FIELDS) {
    let maxDiff = 0, sumDiff = 0, count = 0;
    for (let i = 0; i < n; i++) {
      const a = parseFloat(rowsA[i][field]);
      const b = parseFloat(rowsB[i][field]);
      if (!Number.isFinite(a) || !Number.isFinite(b)) continue;
      const d = Math.abs(a - b);
      maxDiff = Math.max(maxDiff, d);
      sumDiff += d;
      count += 1;
    }
    if (count === 0) continue;
    fields.push({ field, maxDiff, meanDiff: sumDiff / count, count });
  }
  return { skipped: false, rowCountMismatch, rowCountA: rowsA.length, rowCountB: rowsB.length, fields };
}
function summarizePerf(wallSeconds, timings) {
  const perf = { totalWallSeconds: wallSeconds, nRealizations: timings ? timings.length : 0 };
  if (!timings || timings.length === 0) return perf;
  const totals = timings.map((t) => t.totalMs).filter((v) => typeof v === 'number');
  if (totals.length) {
    perf.avgRealizationMs = totals.reduce((a, b) => a + b, 0) / totals.length;
    perf.minRealizationMs = Math.min(...totals);
    perf.maxRealizationMs = Math.max(...totals);
    perf.sumRealizationSeconds = totals.reduce((a, b) => a + b, 0) / 1000;
    if (wallSeconds != null) perf.overheadSeconds = wallSeconds - perf.sumRealizationSeconds;
  }
  const breakdown = {};
  for (const field of ['resampleMs', 'sofiaMs', 'fixtureFitMs', 'fitMs', 'convolveMs']) {
    const vals = timings.map((t) => t[field]).filter((v) => typeof v === 'number');
    if (vals.length) breakdown[field] = vals.reduce((a, b) => a + b, 0) / vals.length;
  }
  if (Object.keys(breakdown).length) perf.avgBreakdownMs = breakdown;
  const counts = {};
  for (const field of ['convolveCalls', 'evalCount']) {
    const vals = timings.map((t) => t[field]).filter((v) => typeof v === 'number');
    if (vals.length) counts[field] = vals.reduce((a, b) => a + b, 0) / vals.length;
  }
  if (Object.keys(counts).length) perf.avgCounts = counts;
  return perf;
}

function argVal(name, def) {
  const args = process.argv.slice(2);
  const i = args.indexOf('--' + name);
  return i >= 0 ? args[i + 1] : def;
}

function main() {
  const jsonPath = argVal('json', path.join(__dirname, 'run_both_report_partial.json'));

  const runs = {};
  for (const key of ['fortranLocal', 'jsDcp']) {
    const objFolder = path.join(TEST_DIR, folders[key], GALAXY.ObjName);
    const csvRows = fs.existsSync(objFolder) ? readBootstrapCsv(objFolder) : null;
    const timings = fs.existsSync(objFolder) ? readTimingsJson(objFolder) : null;
    // Every leg is "not run this invocation" here by definition -- this
    // script never runs anything, only reads. stale:true on whichever legs
    // have complete on-disk results; legs with nothing on disk (e.g. still
    // running, hasn't reached final aggregation/write yet) come back null
    // and are correctly omitted below, exactly like run_both.js's own
    // handling of a skipped leg with no prior results.
    let perf = null;
    if (csvRows) {
      const wallSeconds = estimateWallSeconds(objFolder);
      perf = summarizePerf(wallSeconds, timings);
      if (wallSeconds != null) perf.wallTimeSource = 'estimated-from-file-timestamps';
    }
    const runMeta = readRunMeta(objFolder);
    runs[key] = csvRows ? { exitCode: null, results: csvRows, timings, perf, stale: true, runMeta } : null;
    console.log(`${key}: ${csvRows ? csvRows.length + ' realization row(s) found' : 'no results on disk yet -- omitted'}`);
  }

  const comparisonPairs = [['fortranLocal', 'jsDcp']];
  const comparisons = {};
  for (const [a, b] of comparisonPairs) {
    if (!runs[a] || !runs[b]) continue;
    comparisons[`${a}_vs_${b}`] = compareBootstraps(runs[a].results, runs[b].results);
  }

  // Best-effort top-level summary: the first leg with a recorded RunMeta.
  // Legs can legitimately disagree (different invocations) -- runs[key]
  // .runMeta is the authoritative per-leg source; this is only a convenience
  // for the viewer's subtitle line.
  const firstMeta = ['fortranLocal', 'jsDcp'].map((k) => runs[k] && runs[k].runMeta).find(Boolean);
  const report = {
    seed: (firstMeta && firstMeta.seed) || null,
    nBootstraps: (firstMeta && firstMeta.nBootstraps) || null,
    runs,
    comparisons,
  };
  fs.writeFileSync(jsonPath, JSON.stringify(report, null, 2));
  console.log(`\nWrote ${jsonPath}`);
}

main();
