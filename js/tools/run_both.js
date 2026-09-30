/**
 * run_both.js -- runs TWO variants of the same bootstrap-fit config on the
 * same galaxy/seed and compares numerical results + performance:
 *
 *   - fortran-local: fully-native Fortran pipeline (UseDCP=0).
 *   - js-dcp: the JS/DCP pipeline (UseDCP=1), EITHER dispatched for real
 *     over the DCP network (needs --apiKey; --computeGroups/--slicePrice
 *     optional) OR run in-process with --local (no network, no
 *     credentials, no compute spend -- for when DCP is down/unreachable).
 *     Real dispatch costs real compute credits; wall time then includes
 *     DCP scheduling overhead, not just compute.
 *
 * (This script used to also run a third "js-local" leg as an ALWAYS-ON
 * third comparison point -- removed 2026-08 as a permanent leg (it was
 * consistently ~4x slower than js-dcp on this hardware, with no remaining
 * diagnostic value once that gap was understood), but --local below
 * revives the same underlying mode as an opt-in swap-in for the js-dcp
 * leg specifically, not a third leg -- for exactly the situation that
 * removal didn't anticipate: DCP being unreachable at all. Both routes
 * through the same bootstrap-realization-launcher.js `--local N` mode,
 * via RunBootstrapsDCP.py/RunInitialFitDCP.py's DCP_FORCE_LOCAL fallback.)
 *
 * Both legs run SEQUENTIALLY off the SAME BootstrapSeed (a matched-seed
 * diff, not independent-random noise), each into its own TargFolder so
 * results don't collide.
 *
 * Usage:
 *   node run_both.js --seed <idum> [--nBootstraps N]
 *     [--objName NAME [--cube PATH] [--mask PATH] --pa DEG --inc DEG]
 *     [--apiKey 0x... | --local]
 *     [--computeGroups joinKey[,joinSecret][:joinKey[,joinSecret]...]]
 *     [--slicePrice N] [--skip-fortran] [--skip-js-dcp]
 *     [--skip-wipe] [--json PATH]
 *
 * js-dcp needs --apiKey (or DCP_API_KEY in the environment) -- without it
 * AND without --local, that leg is skipped automatically, same spirit as
 * --skip-js-dcp. --local always wins over a real dispatch when both an
 * apiKey and --local are present (never sends DCP_API_KEY through in that
 * case, even if set in the environment for some unrelated reason).
 *
 *
node run_both.js \
  --seed 42 \
  --nBootstraps 1000 \
  --cloudDensity 500 \
  --apiKey 0xf1512793d2dcb94a0102d53e6ab55ac8b145982342eae999be826aed54533ec7 \
  --computeGroups bell,18be80 \
  --slicePrice 1.012 \
  --json 1000_bootstraps_bell.json \
  --skip-fortran
 *
 *
 */
'use strict';
const { spawn, execSync, exec } = require('node:child_process');
const path = require('node:path');
const fs = require('node:fs');
const os = require('node:os');

// ---------------------------------------------------------------------------
// getPerformanceCoreCount
// See the identical helper in bootstrap-realization-launcher.js for the full
// rationale/measurements -- os.cpus().length counts Apple Silicon's
// performance and efficiency cores as equivalent, but E-cores are much
// slower for this kind of sustained CPU-bound work, so sizing a worker pool
// to logical-cores-1 oversubscribes the fast cores. Same fix applied here
// for fortran-local's nProcessors (Python's multiprocessing.Pool, one
// subprocess per bootstrap realization -- FullSingleGalaxyFit.py:118),
// which had the identical logical-cores-1 sizing and therefore the same
// P/E oversubscription exposure as js-local did before that fix. Returns
// null (caller falls back to logical-cores-1) on non-Darwin platforms or
// Intel Macs, where cores are already homogeneous.
// ---------------------------------------------------------------------------
function getPerformanceCoreCount() {
  if (process.platform !== 'darwin') return null;
  try {
    const out = execSync('sysctl -n hw.perflevel0.physicalcpu', { stdio: ['ignore', 'pipe', 'ignore'] })
      .toString().trim();
    const n = parseInt(out, 10);
    return Number.isFinite(n) && n > 0 ? n : null;
  } catch (e) {
    return null;
  }
}

// ---------------------------------------------------------------------------
// writeAndOpenHtmlReport
//
// run_both_report.html (this same directory) is a static viewer template
// -- normally you'd drag a *_report.json onto it by hand. This writes a
// fresh copy of that template with the just-produced `report` baked in as
// a <script id="embedded-report-data" type="application/json"> (inserted
// right after the EMBEDDED_REPORT_DATA_MARKER comment near the end of the
// file, replacing whatever a prior run left there), then opens it directly
// as a file:// URL in the default browser -- no server needed.
//
// (A local HTTP server was tried first, to work around what looked like a
// file://-specific Chrome error -- "Unsafe attempt to load URL file://X
// from frame with URL file://X. 'file:' URLs are treated as unique
// security origins." Turned out to be a red herring: the REAL bug was
// run_both_report.html's own auto-load script running in a <script> block
// that appears BEFORE the embedded-data <script> tag in document order, so
// getElementById('embedded-report-data') returned null every time
// (confirmed directly) -- fixed there by deferring that check to
// DOMContentLoaded. Once that was actually fixed, file:// worked fine, so
// the server was removed again rather than kept as unneeded complexity.)
//
// Overwrites run_both_report.html in place (not a timestamped copy) --
// deliberate: this is a dev tool's "latest run" view, same spirit as the
// JSON report defaulting to one fixed filename unless --json says
// otherwise. Never throws: a report-viewing convenience failing (missing
// template, no GUI to open a browser on a headless box, etc.) shouldn't
// fail the run or its exit code.
// ---------------------------------------------------------------------------
function writeAndOpenHtmlReport(report) {
  try {
    const templatePath = path.join(__dirname, 'run_both_report.html');
    let html = fs.readFileSync(templatePath, 'utf8');

    const marker = '<!-- EMBEDDED_REPORT_DATA_MARKER';
    const markerIdx = html.indexOf(marker);
    if (markerIdx === -1) {
      console.log('[run_both] WARNING: run_both_report.html has no EMBEDDED_REPORT_DATA_MARKER -- skipping HTML report');
      return;
    }
    const markerEndIdx = html.indexOf('-->', markerIdx) + '-->'.length;

    // Strip any embedded-data script tag a PRIOR run of this function left
    // right after the marker, so re-running against the same file is
    // idempotent instead of stacking duplicate script tags.
    const afterMarker = html.slice(markerEndIdx);
    const staleScriptMatch = afterMarker.match(/^\s*<script id="embedded-report-data"[^>]*>[\s\S]*?<\/script>/);
    const afterMarkerClean = staleScriptMatch ? afterMarker.slice(staleScriptMatch[0].length) : afterMarker;

    const dataScript = `\n<script id="embedded-report-data" type="application/json">${
      JSON.stringify(report).replace(/</g, '\\u003c')
    }</script>`;
    html = html.slice(0, markerEndIdx) + dataScript + afterMarkerClean;

    fs.writeFileSync(templatePath, html);
    console.log(`[run_both] wrote HTML report to ${templatePath}`);

    const fileUrl = require('node:url').pathToFileURL(templatePath).href;
    const openCmd = process.platform === 'darwin' ? 'open'
      : process.platform === 'win32' ? 'start ""'
      : 'xdg-open';
    exec(`${openCmd} ${JSON.stringify(fileUrl)}`, (err) => {
      if (err) console.log(`[run_both] could not auto-open the HTML report (${err.message}) -- open it manually: ${templatePath}`);
    });
  } catch (e) {
    console.log(`[run_both] WARNING: failed to write/open the HTML report: ${e.message}`);
  }
}

const TEST_DIR = path.join(__dirname, '..', '..', '3KIDNASTests', 'SingleGalaxyTest');
const DRIVER = path.join(__dirname, '..', '..', 'WRKP_GalaxyFitDriver.py');
const BASE_FITTING_OPTIONS = path.join(__dirname, '..', '..', 'Inputs', 'SingleGalaxyTestFittingOptions_Base.txt');

// Default test galaxy -- overridable per-invocation via --objName (+
// --cube/--mask if the mask filename doesn't follow this galaxy's own
// WALLABY_<objName>_mask.fits convention, e.g. WALLABY_J100336-262923's
// mask is SoFiA_J100336-262923_mask.fits, a SoFiA-catalogue name, not a
// WALLABY_ one) and --pa/--inc. See main()'s arg parsing below.
const GALAXY = {
  CubeName: '../TestData/WALLABY_Test_sources/WALLABY_J103538-484832/WALLABY_J103538-484832_VelCube.fits',
  MaskName: '../TestData/WALLABY_Test_sources/WALLABY_J103538-484832/WALLABY_J103538-484832_mask.fits',
  ObjName: 'WALLABY_J103538-484832',
  PA_Estimate: 241.319,
  Inc_Estimate: 89.00,
};

// Fortran and the JS/DCP bootstrap payload both ultimately read the SAME
// fitting-options file (RunWRKP.LoadDefaultWRKPFiles for Fortran's own run,
// RunBootstrapsDCP.ParseFittingOptionsExtras replaying Fortran's own read
// order over those same lines for the payload) -- its path is a single
// hardcoded default in SetFileLocations.py (WRKP_GeneralOptionsIn), but
// GalaxyFitParameters.OverwriteDefaults already lets ANY GeneralDict default
// be overridden just by defining a same-named variable in the run's config
// .py file (that's how BootstrapSeed already works). So overriding the
// cloud density doesn't need touching SetFileLocations.py or Fortran at
// all: write a per-leg copy of the base options file with the density line
// swapped, then have writeConfig() point WRKP_GeneralOptionsIn at it --
// both legs pick up the same new value through the existing mechanism.
function writeCloudDensityOptionsFile(cloudDensity, outPath) {
  const base = fs.readFileSync(BASE_FITTING_OPTIONS, 'utf8');
  const linesArr = base.split('\n');
  const labelIdx = linesArr.findIndex((l) => l.replace(/^#\t*/, '').trim() === 'The base cloud surface density');
  if (labelIdx === -1 || labelIdx + 1 >= linesArr.length) {
    throw new Error(`writeCloudDensityOptionsFile: couldn't find the cloud-surface-density line in ${BASE_FITTING_OPTIONS}`);
  }
  const valueStr = Number.isInteger(cloudDensity) ? `${cloudDensity}.` : String(cloudDensity);
  linesArr[labelIdx + 1] = valueStr;
  fs.writeFileSync(outPath, linesArr.join('\n'));
}

// The value actually in effect when --cloudDensity isn't passed -- read
// fresh off the base file rather than hardcoding "400" here, so this stays
// correct if the base file's own default ever changes.
function readBaseCloudDensity() {
  const base = fs.readFileSync(BASE_FITTING_OPTIONS, 'utf8');
  const linesArr = base.split('\n');
  const labelIdx = linesArr.findIndex((l) => l.replace(/^#\t*/, '').trim() === 'The base cloud surface density');
  return parseFloat(linesArr[labelIdx + 1]);
}

// Reproducibility metadata, written into a leg's OWN output folder right
// alongside its BootstrapFits.csv/BootstrapTimings.json -- travels with the
// results regardless of which invocation produced them, so a later, separate
// invocation (or build_partial_report.js, reading a stale leg from disk)
// can still recover what this leg was actually run with. nProcessors/
// totalLogicalCores are null for js-dcp: that leg's work is dispatched to
// real DCP workers over the network, whose CPU counts aren't knowable here
// at all (unlike fortran-local, which runs as a process pool sized off THIS
// machine's own os.cpus()).
function writeRunMeta(objFolder, { seed, nBootstraps, nProcessors, cloudDensity }) {
  const meta = {
    seed: seed || null,
    nBootstraps,
    nProcessors: nProcessors != null ? nProcessors : null,
    totalLogicalCores: nProcessors != null ? os.cpus().length : null,
    cloudDensity,
    timestamp: new Date().toISOString(),
  };
  fs.mkdirSync(objFolder, { recursive: true });
  fs.writeFileSync(path.join(objFolder, `${GALAXY.ObjName}_RunMeta.json`), JSON.stringify(meta, null, 2));
}
function readRunMeta(objFolder) {
  const p = path.join(objFolder, `${GALAXY.ObjName}_RunMeta.json`);
  if (!fs.existsSync(p)) return null;
  return JSON.parse(fs.readFileSync(p, 'utf8'));
}

function writeConfig(filePath, { targFolder, nBootstraps, nProcessors, useDCP, seed, cloudDensity }) {
  const lines = [
    `CubeName="${GALAXY.CubeName}"`,
    `MaskName="${GALAXY.MaskName}"`,
    `ObjName="${GALAXY.ObjName}"`,
    `TargFolder="${targFolder}/"`,
    // .toFixed, not template-literal interpolation: GalaxyFitParameters.
    // CheckParamTypes requires these as Python float, and JS numbers don't
    // preserve a trailing ".00" (89.00 stringifies to "89", which Python's
    // ast/exec reads back as an int, failing the type check).
    `PA_Estimate= ${GALAXY.PA_Estimate.toFixed(3)}`,
    `Inc_Estimate=${GALAXY.Inc_Estimate.toFixed(2)}`,
    `nBootstraps= ${nBootstraps}`,
    `nProcessors_Bootstraps=${nProcessors}`,
    `UseDCP=${useDCP ? 1 : 0}`,
    `BootstrapSeed=${seed}`,
  ];
  if (cloudDensity != null) {
    // Same folder/name pattern as the leg's own config file, not the shared
    // Inputs/ location -- this is a per-leg, per-run copy, never the base
    // file other tests may still be reading.
    const optionsPath = filePath.replace(/\.py$/, '_fitting_options.txt');
    writeCloudDensityOptionsFile(cloudDensity, optionsPath);
    lines.push(`WRKP_GeneralOptionsIn="${optionsPath}"`);
  }
  lines.push('');
  fs.writeFileSync(filePath, lines.join('\n'));
}

function run(cmd, args, opts = {}) {
  return new Promise((resolve, reject) => {
    const t0 = Date.now();
    const child = spawn(cmd, args, { stdio: ['ignore', 'pipe', 'pipe'], ...opts });
    let out = '';
    let err = '';
    child.stdout.on('data', (d) => { out += d; process.stdout.write(d); });
    child.stderr.on('data', (d) => { err += d; process.stderr.write(d); });
    child.on('close', (code) => resolve({ code, seconds: (Date.now() - t0) / 1000, stdout: out, stderr: err }));
    child.on('error', reject);
  });
}

// Minimal CSV line parser respecting double-quoted fields that themselves
// contain commas (this project's array-valued columns, e.g. "7.49, 22.48, ...").
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

const SCALAR_FIELDS = ['X_model', 'Y_model', 'Inc_model', 'PA_model', 'Vsys_model',
  'RA_model', 'DEC_model', 'Vdisp_model', 'RHI_AS', 'VHI'];

// Pairwise diff between two runs' BootstrapFits.csv rows.
function compareBootstraps(rowsA, rowsB) {
  if (!rowsA || !rowsB) return { skipped: true, fields: [] };
  const n = Math.min(rowsA.length, rowsB.length);
  const rowCountMismatch = rowsA.length !== rowsB.length;
  const fields = [];
  for (const field of SCALAR_FIELDS) {
    let maxDiff = 0, sumDiff = 0, count = 0;
    // %diff per row: |a-b| / mean(|a|,|b|) -- symmetric, so it doesn't
    // matter which side is "reference." Rows where both sides are ~0 (e.g.
    // a field that's genuinely zero, like a fixed Vdisp_model) are excluded
    // from the % average/max (0/0 is meaningless, not "0% different") but
    // still count toward the absolute maxDiff/meanDiff above.
    let maxPctDiff = 0, sumPctDiff = 0, pctCount = 0;
    const ZERO_EPS = 1e-9;
    for (let i = 0; i < n; i++) {
      const a = parseFloat(rowsA[i][field]);
      const b = parseFloat(rowsB[i][field]);
      if (!Number.isFinite(a) || !Number.isFinite(b)) continue;
      const d = Math.abs(a - b);
      maxDiff = Math.max(maxDiff, d);
      sumDiff += d;
      count += 1;
      const denom = (Math.abs(a) + Math.abs(b)) / 2;
      if (denom > ZERO_EPS) {
        const pct = (d / denom) * 100;
        maxPctDiff = Math.max(maxPctDiff, pct);
        sumPctDiff += pct;
        pctCount += 1;
      }
    }
    if (count === 0) continue;
    fields.push({
      field, maxDiff, meanDiff: sumDiff / count, count,
      maxPctDiff: pctCount ? maxPctDiff : null,
      meanPctDiff: pctCount ? sumPctDiff / pctCount : null,
    });
  }
  return { skipped: false, rowCountMismatch, rowCountA: rowsA.length, rowCountB: rowsB.length, fields };
}

// Performance summary for one run: total wall time (the whole subprocess,
// including Python/Fortran/DCP-dispatch overhead) plus the average
// per-realization time and a breakdown, from the same
// {realizationIndex, resampleMs, sofiaMs, fixtureFitMs, fitMs, convolveMs,
// convolveCalls, evalCount, totalMs} timings schema both sides already
// write (Bootstrap_Outputs.StoreBootstrapTimings_JSON /
// RunBootstrapsDCP.SummarizeWorkerTimings).
function summarizePerf(wallSeconds, timings) {
  const perf = { totalWallSeconds: wallSeconds, nRealizations: timings ? timings.length : 0 };
  if (!timings || timings.length === 0) return perf;
  const totals = timings.map((t) => t.totalMs).filter((v) => typeof v === 'number');
  if (totals.length) {
    perf.avgRealizationMs = totals.reduce((a, b) => a + b, 0) / totals.length;
    perf.minRealizationMs = Math.min(...totals);
    perf.maxRealizationMs = Math.max(...totals);
    perf.sumRealizationSeconds = totals.reduce((a, b) => a + b, 0) / 1000;
    // Wall time minus sum of workers' own execution time: process/Python
    // overhead for fortran-local, real DCP scheduling/network overhead for
    // js-dcp -- the number that separates dispatch cost from compute cost.
    // Left unset (not NaN) when wallSeconds is null (a leg reported from a
    // prior invocation's on-disk results, not run now).
    if (wallSeconds != null) perf.overheadSeconds = wallSeconds - perf.sumRealizationSeconds;
  }
  const breakdown = {};
  for (const field of ['resampleMs', 'sofiaMs', 'fixtureFitMs', 'fitMs', 'convolveMs']) {
    const vals = timings.map((t) => t[field]).filter((v) => typeof v === 'number');
    if (vals.length) breakdown[field] = vals.reduce((a, b) => a + b, 0) / vals.length;
  }
  if (Object.keys(breakdown).length) perf.avgBreakdownMs = breakdown;

  // Raw counts, kept separate from avgBreakdownMs (which stays purely
  // ms-denominated). evalCount distinguishes "each eval is slower" from
  // "the optimizer needed more evals" -- either produces the same fitMs gap.
  const counts = {};
  for (const field of ['convolveCalls', 'evalCount']) {
    const vals = timings.map((t) => t[field]).filter((v) => typeof v === 'number');
    if (vals.length) counts[field] = vals.reduce((a, b) => a + b, 0) / vals.length;
  }
  if (Object.keys(counts).length) perf.avgCounts = counts;
  return perf;
}

async function main() {
  const args = process.argv.slice(2);
  const argVal = (name, def) => {
    const i = args.indexOf('--' + name);
    return i >= 0 ? args[i + 1] : def;
  };
  const seed = parseInt(argVal('seed', '0'), 10);
  const nBootstraps = parseInt(argVal('nBootstraps', '5'), 10);
  // Cap parallel workers independent of nBootstraps -- tying
  // nProcessors_Bootstraps 1:1 to nBootstraps let a large --nBootstraps spawn
  // that many concurrent SoFiA/fit worker processes with no relation to the
  // machine's actual core count, which can exhaust memory and thermally
  // throttle a fanless machine badly enough to crash WindowServer. Sized to
  // the performance-core count when known (see getPerformanceCoreCount
  // above) rather than logical-cores-1, which oversubscribes the fast cores
  // on Apple Silicon's heterogeneous P+E designs -- falls back to the
  // original logical-cores-1 sizing where that split doesn't apply.
  const pCores = getPerformanceCoreCount();
  const nProcessors = Math.max(1, Math.min(nBootstraps, pCores || (os.cpus().length - 1)));
  // Overrides Inputs/SingleGalaxyTestFittingOptions_Base.txt's cloud density
  // (400 by default) for both legs -- see writeCloudDensityOptionsFile for
  // why this doesn't need to touch Fortran or the shared base file.
  // Undefined (flag omitted) means "don't override," i.e. today's behavior.
  const cloudDensityFlag = argVal('cloudDensity');
  const cloudDensity = cloudDensityFlag !== undefined ? parseFloat(cloudDensityFlag) : undefined;
  // What's actually in effect this run, whether overridden or not -- recorded
  // into each leg's RunMeta.json below so a report always says what density
  // produced it, instead of leaving it implicit whenever --cloudDensity was
  // omitted.
  const effectiveCloudDensity = cloudDensity != null ? cloudDensity : readBaseCloudDensity();
  const skipWipe = args.includes('--skip-wipe');
  const skipFortran = args.includes('--skip-fortran');
  let skipJsDcp = args.includes('--skip-js-dcp');

  const apiKey = argVal('apiKey', process.env.DCP_API_KEY);
  const computeGroups = argVal('computeGroups', process.env.DCP_COMPUTE_GROUPS);
  const slicePrice = argVal('slicePrice', process.env.DCP_SLICE_PRICE);

  const jsonPath = argVal('json', path.join(__dirname, 'run_both_report.json'));

  // Galaxy override -- default is GALAXY (WALLABY_J103538-484832) above.
  // --objName switches to <objName>'s own folder under TestData/
  // WALLABY_Test_sources/, defaulting cube/mask to that galaxy's own
  // WALLABY_<objName>_VelCube.fits/_mask.fits naming -- override either
  // with --cube/--mask (relative to TEST_DIR, or absolute) when a mask
  // doesn't follow that convention. --pa/--inc override the estimates.
  const objNameFlag = argVal('objName');
  if (objNameFlag) {
    GALAXY.ObjName = objNameFlag;
    GALAXY.CubeName = `../TestData/WALLABY_Test_sources/${objNameFlag}/${objNameFlag}_VelCube.fits`;
    GALAXY.MaskName = `../TestData/WALLABY_Test_sources/${objNameFlag}/${objNameFlag}_mask.fits`;
  }
  const cubeFlag = argVal('cube');
  if (cubeFlag) GALAXY.CubeName = cubeFlag;
  const maskFlag = argVal('mask');
  if (maskFlag) GALAXY.MaskName = maskFlag;
  const paFlag = argVal('pa');
  if (paFlag !== undefined) GALAXY.PA_Estimate = parseFloat(paFlag);
  const incFlag = argVal('inc');
  if (incFlag !== undefined) GALAXY.Inc_Estimate = parseFloat(incFlag);

  if (!seed) {
    console.log('[run_both] WARNING: no --seed given -- both runs will use unseeded, '
      + 'time-based randomness. The comparison below will only be a loose statistical sanity '
      + 'check, not a matched-seed diff.');
  }

  // --local: run the js leg in-process (bootstrap-realization-launcher.js's
  // own --local N mode via RunBootstrapsDCP.py/RunInitialFitDCP.py's
  // DCP_FORCE_LOCAL fallback -- see their own comments) instead of real DCP
  // dispatch. No network, no credentials, no compute spend -- for when DCP
  // is down/unreachable or a --apiKey isn't available. Mutually exclusive
  // with real dispatch: forces DCP_FORCE_LOCAL=1 in the leg's own
  // subprocess env regardless of whether --apiKey/DCP_API_KEY is also set,
  // and (unlike a real dispatch) never needs an apiKey to avoid the
  // no-apiKey auto-skip below.
  const runLocal = args.includes('--local');

  if (!skipJsDcp && !runLocal && !apiKey) {
    console.log('[run_both] WARNING: no --apiKey given and DCP_API_KEY not set -- '
      + 'skipping the js-dcp leg (nothing to authenticate a real dispatch with). '
      + 'Pass --local to run it in-process instead (no dispatch, no credentials needed).');
    skipJsDcp = true;
  }

  // The js leg's key/folder/label all follow --local, so a local run and a
  // real-dispatch run never collide in the same folder or get mislabeled as
  // each other in the console/report (hit directly: an earlier version
  // always called this leg "jsDcp" even when --local was used, making a
  // local-only run look like it had been dispatched for real).
  const jsLegKey = runLocal ? 'jsLocal' : 'jsDcp';
  const jsLegFolderLabel = runLocal ? 'JSLocal' : 'JSDcp';
  const jsLegLogLabel = runLocal ? 'JS-LOCAL (in-process, no dispatch)' : 'JS-DCP (real network dispatch)';

  const folders = {
    fortranLocal: 'TestFits_RunAllThree_FortranLocal',
    [jsLegKey]: `TestFits_RunAllThree_${jsLegFolderLabel}`,
  };
  const configPaths = {
    fortranLocal: path.join(TEST_DIR, 'run_both_fortran_local_config.py'),
    [jsLegKey]: path.join(TEST_DIR, `run_both_${jsLegKey}_config.py`),
  };

  writeConfig(configPaths.fortranLocal, { targFolder: folders.fortranLocal, nBootstraps, nProcessors, useDCP: false, seed, cloudDensity });
  writeConfig(configPaths[jsLegKey], { targFolder: folders[jsLegKey], nBootstraps, nProcessors, useDCP: true, seed, cloudDensity });

  if (!skipWipe) {
    // Only wipe folders for legs actually running this invocation -- wiping
    // a skipped leg's folder unconditionally destroyed its still-valid
    // prior results for no reason (hit directly: re-running just js-dcp
    // after --skip-fortran deleted the fortran-local CSVs from the
    // immediately-preceding full run).
    const activeFolders = [
      !skipFortran && folders.fortranLocal,
      !skipJsDcp && folders[jsLegKey],
    ].filter(Boolean);
    console.log(`[run_both] wiping prior output for: ${activeFolders.join(', ') || '(none -- everything skipped)'}`);
    for (const f of activeFolders) {
      fs.rmSync(path.join(TEST_DIR, f), { recursive: true, force: true });
    }
  }

  console.log(`[run_both] seed=${seed || '(none)'} nBootstraps=${nBootstraps}`);

  const results = { fortranLocal: null, [jsLegKey]: null };

  if (!skipFortran) {
    console.log('\n[run_both] running FORTRAN-LOCAL (fully-native) pipeline...');
    results.fortranLocal = await run('python3', [DRIVER, path.basename(configPaths.fortranLocal)], { cwd: TEST_DIR });
    console.log(`[run_both] fortran-local: exit=${results.fortranLocal.code} ${results.fortranLocal.seconds.toFixed(1)}s`);
    writeRunMeta(path.join(TEST_DIR, folders.fortranLocal, GALAXY.ObjName), { seed, nBootstraps, nProcessors, cloudDensity: effectiveCloudDensity });
  }

  if (!skipJsDcp) {
    console.log(`\n[run_both] running ${jsLegLogLabel}...`);
    const dcpEnv = { ...process.env };
    if (runLocal) {
      dcpEnv.DCP_FORCE_LOCAL = '1';
      delete dcpEnv.DCP_API_KEY; // never needed for --local; don't leak one through if set for other reasons
    } else {
      dcpEnv.DCP_API_KEY = apiKey;
    }
    if (computeGroups) dcpEnv.DCP_COMPUTE_GROUPS = computeGroups;
    if (slicePrice) dcpEnv.DCP_SLICE_PRICE = slicePrice;
    results[jsLegKey] = await run('python3', [DRIVER, path.basename(configPaths[jsLegKey])], { cwd: TEST_DIR, env: dcpEnv });
    console.log(`[run_both] ${jsLegKey}: exit=${results[jsLegKey].code} ${results[jsLegKey].seconds.toFixed(1)}s`);
    // nProcessors omitted (null) for real dispatch: that leg's work runs on
    // real DCP workers over the network, not as a local process pool --
    // their CPU counts aren't ours to report. Reported normally for
    // --local, since that DOES run as a local worker_threads pool (see
    // bootstrap-realization-launcher.js's own --local sizing).
    writeRunMeta(path.join(TEST_DIR, folders[jsLegKey], GALAXY.ObjName), { seed, nBootstraps, nProcessors: runLocal ? nProcessors : null, cloudDensity: effectiveCloudDensity });
  }

  console.log('\n=== run_both report ===');

  const runs = {};
  for (const key of ['fortranLocal', jsLegKey]) {
    const r = results[key];
    const objFolder = path.join(TEST_DIR, folders[key], GALAXY.ObjName);
    const csvRows = fs.existsSync(objFolder) ? readBootstrapCsv(objFolder) : null;
    const timings = fs.existsSync(objFolder) ? readTimingsJson(objFolder) : null;
    // Read back rather than reusing this invocation's own seed/nBootstraps/
    // etc: for a leg NOT run this invocation, those don't apply to it at all
    // (it may be from an entirely earlier invocation with different
    // settings) -- RunMeta.json lives in the leg's own folder and was
    // written by whichever invocation actually produced these results, so
    // it's correct either way. null (older data, predating this file) is a
    // legitimate "unknown," not an error.
    const runMeta = readRunMeta(objFolder);

    if (!r) {
      // Not run this invocation -- if its folder still has results from an
      // earlier invocation (now preserved, see the wipe-only-active-folders
      // fix above), report those instead of dropping the leg from the
      // combined report entirely. No fresh wall-time to report for it.
      runs[key] = csvRows ? { exitCode: null, results: csvRows, timings, perf: summarizePerf(null, timings), stale: true, runMeta } : null;
      if (runs[key]) console.log(`\n-- ${key} -- (not run this invocation, showing prior results)`);
      continue;
    }
    const perf = summarizePerf(r.seconds, timings);
    runs[key] = {
      exitCode: r.code,
      results: csvRows,
      timings,
      perf,
      runMeta,
    };
    console.log(`\n-- ${key} --`);
    console.log(`  exit=${r.code}  totalWallSeconds=${r.seconds.toFixed(1)}`
      + (perf.avgRealizationMs != null ? `  avgRealizationMs=${perf.avgRealizationMs.toFixed(0)}` : ''));
    if (csvRows) console.log(`  ${csvRows.length} realization row(s) parsed from BootstrapFits.csv`);
  }

  console.log('\n=== Pairwise numerical comparison ===');
  const comparisonPairs = [
    ['fortranLocal', jsLegKey],
  ];
  const comparisons = {};
  for (const [a, b] of comparisonPairs) {
    if (!runs[a] || !runs[b]) continue;
    const cmp = compareBootstraps(runs[a].results, runs[b].results);
    comparisons[`${a}_vs_${b}`] = cmp;
    if (!cmp.skipped) {
      console.log(`\n${a} vs ${b}:`);
      console.log(`  ${'field'.padEnd(12)} ${'max|diff|'.padStart(12)} ${'mean|diff|'.padStart(12)} ${'max%diff'.padStart(10)} ${'mean%diff'.padStart(10)}`);
      for (const f of cmp.fields) {
        const maxPct = f.maxPctDiff != null ? f.maxPctDiff.toFixed(2) + '%' : 'n/a';
        const meanPct = f.meanPctDiff != null ? f.meanPctDiff.toFixed(2) + '%' : 'n/a';
        console.log(`  ${f.field.padEnd(12)} ${f.maxDiff.toFixed(6).padStart(12)} ${f.meanDiff.toFixed(6).padStart(12)} ${maxPct.padStart(10)} ${meanPct.padStart(10)}`);
      }
    }
  }

  const report = {
    seed: seed || null,
    nBootstraps,
    runs,
    comparisons,
  };
  fs.writeFileSync(jsonPath, JSON.stringify(report, null, 2));
  console.log(`\n[run_both] wrote report to ${jsonPath}`);

  writeAndOpenHtmlReport(report);

  const ranOk = (skipFortran || (results.fortranLocal && results.fortranLocal.code === 0))
    && (skipJsDcp || (results[jsLegKey] && results[jsLegKey].code === 0));
  process.exit(ranOk ? 0 : 1);
}

main().catch((e) => { console.error(e); process.exit(2); });
