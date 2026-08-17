/**
 * run_both.js -- runs TWO variants of the same bootstrap-fit config on the
 * same galaxy/seed and compares numerical results + performance:
 *
 *   - fortran-local: fully-native Fortran pipeline (UseDCP=0).
 *   - js-dcp: the JS/DCP pipeline (UseDCP=1), dispatched for real over the
 *     DCP network (needs --apiKey; --computeGroups/--slicePrice optional).
 *     Costs real compute credits; wall time includes DCP scheduling
 *     overhead, not just compute.
 *
 * (This script used to also run a third "js-local" leg -- the same JS/DCP
 * pipeline executed in-process on this machine, no network, no dcp-client.
 * Removed: it was consistently ~4x slower than js-dcp on this hardware,
 * with no remaining diagnostic value once that gap was understood.
 * bootstrap-realization-launcher.js's own `--local N` mode that leg used
 * still exists and still works for other local/manual testing -- only its
 * use as a leg in this comparison harness is gone.)
 *
 * Both legs run SEQUENTIALLY off the SAME BootstrapSeed (a matched-seed
 * diff, not independent-random noise), each into its own TargFolder so
 * results don't collide.
 *
 * Usage:
 *   node run_both.js --seed <idum> [--nBootstraps N]
 *     [--apiKey 0x...] [--computeGroups joinKey[,joinSecret][:joinKey[,joinSecret]...]]
 *     [--slicePrice N] [--skip-fortran] [--skip-js-dcp]
 *     [--skip-wipe] [--json PATH]
 *
 * js-dcp needs --apiKey (or DCP_API_KEY in the environment) -- without it,
 * that leg is skipped automatically, same spirit as --skip-js-dcp.
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
const { spawn, execSync } = require('node:child_process');
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

const TEST_DIR = path.join(__dirname, '..', '..', '3KIDNASTests', 'SingleGalaxyTest');
const DRIVER = path.join(__dirname, '..', '..', 'WRKP_GalaxyFitDriver.py');
const BASE_FITTING_OPTIONS = path.join(__dirname, '..', '..', 'Inputs', 'SingleGalaxyTestFittingOptions_Base.txt');

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

  if (!seed) {
    console.log('[run_both] WARNING: no --seed given -- both runs will use unseeded, '
      + 'time-based randomness. The comparison below will only be a loose statistical sanity '
      + 'check, not a matched-seed diff.');
  }
  if (!skipJsDcp && !apiKey) {
    console.log('[run_both] WARNING: no --apiKey given and DCP_API_KEY not set -- '
      + 'skipping the js-dcp leg (nothing to authenticate a real dispatch with).');
    skipJsDcp = true;
  }

  const folders = {
    fortranLocal: 'TestFits_RunAllThree_FortranLocal',
    jsDcp: 'TestFits_RunAllThree_JSDcp',
  };
  const configPaths = {
    fortranLocal: path.join(TEST_DIR, 'run_both_fortran_local_config.py'),
    jsDcp: path.join(TEST_DIR, 'run_both_js_dcp_config.py'),
  };

  writeConfig(configPaths.fortranLocal, { targFolder: folders.fortranLocal, nBootstraps, nProcessors, useDCP: false, seed, cloudDensity });
  writeConfig(configPaths.jsDcp, { targFolder: folders.jsDcp, nBootstraps, nProcessors, useDCP: true, seed, cloudDensity });

  if (!skipWipe) {
    // Only wipe folders for legs actually running this invocation -- wiping
    // a skipped leg's folder unconditionally destroyed its still-valid
    // prior results for no reason (hit directly: re-running just js-dcp
    // after --skip-fortran deleted the fortran-local CSVs from the
    // immediately-preceding full run).
    const activeFolders = [
      !skipFortran && folders.fortranLocal,
      !skipJsDcp && folders.jsDcp,
    ].filter(Boolean);
    console.log(`[run_both] wiping prior output for: ${activeFolders.join(', ') || '(none -- everything skipped)'}`);
    for (const f of activeFolders) {
      fs.rmSync(path.join(TEST_DIR, f), { recursive: true, force: true });
    }
  }

  console.log(`[run_both] seed=${seed || '(none)'} nBootstraps=${nBootstraps}`);

  const results = { fortranLocal: null, jsDcp: null };

  if (!skipFortran) {
    console.log('\n[run_both] running FORTRAN-LOCAL (fully-native) pipeline...');
    results.fortranLocal = await run('python3', [DRIVER, path.basename(configPaths.fortranLocal)], { cwd: TEST_DIR });
    console.log(`[run_both] fortran-local: exit=${results.fortranLocal.code} ${results.fortranLocal.seconds.toFixed(1)}s`);
    writeRunMeta(path.join(TEST_DIR, folders.fortranLocal, GALAXY.ObjName), { seed, nBootstraps, nProcessors, cloudDensity: effectiveCloudDensity });
  }

  if (!skipJsDcp) {
    console.log('\n[run_both] running JS-DCP (real network dispatch)...');
    const dcpEnv = { ...process.env, DCP_API_KEY: apiKey };
    if (computeGroups) dcpEnv.DCP_COMPUTE_GROUPS = computeGroups;
    if (slicePrice) dcpEnv.DCP_SLICE_PRICE = slicePrice;
    results.jsDcp = await run('python3', [DRIVER, path.basename(configPaths.jsDcp)], { cwd: TEST_DIR, env: dcpEnv });
    console.log(`[run_both] js-dcp: exit=${results.jsDcp.code} ${results.jsDcp.seconds.toFixed(1)}s`);
    // nProcessors omitted (null): this leg's work runs on real DCP workers
    // over the network, not as a local process pool -- their CPU counts
    // aren't ours to report.
    writeRunMeta(path.join(TEST_DIR, folders.jsDcp, GALAXY.ObjName), { seed, nBootstraps, nProcessors: null, cloudDensity: effectiveCloudDensity });
  }

  console.log('\n=== run_both report ===');

  const runs = {};
  for (const key of ['fortranLocal', 'jsDcp']) {
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
    ['fortranLocal', 'jsDcp'],
  ];
  const comparisons = {};
  for (const [a, b] of comparisonPairs) {
    if (!runs[a] || !runs[b]) continue;
    const cmp = compareBootstraps(runs[a].results, runs[b].results);
    comparisons[`${a}_vs_${b}`] = cmp;
    if (!cmp.skipped) {
      console.log(`\n${a} vs ${b}:`);
      console.log(`  ${'field'.padEnd(12)} ${'max|diff|'.padStart(12)} ${'mean|diff|'.padStart(12)}`);
      for (const f of cmp.fields) {
        console.log(`  ${f.field.padEnd(12)} ${f.maxDiff.toFixed(6).padStart(12)} ${f.meanDiff.toFixed(6).padStart(12)}`);
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

  const ranOk = (skipFortran || (results.fortranLocal && results.fortranLocal.code === 0))
    && (skipJsDcp || (results.jsDcp && results.jsDcp.code === 0));
  process.exit(ranOk ? 0 : 1);
}

main().catch((e) => { console.error(e); process.exit(2); });
