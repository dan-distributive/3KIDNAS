#!/usr/bin/env node
'use strict';
// =============================================================================
// run-galaxy-fit-cli.js
//
// A from-scratch, Python-free CLI entry point for fitting one galaxy through
// the JS/DCP pipeline: reads a cube + mask straight off disk, builds the
// payload itself (via src/PayloadBuilder/buildFitPayloads.js, the JS
// equivalent of RunInitialFitDCP.py/RunBootstrapsDCP.py's payload builders),
// and dispatches through the SAME, already-hardened bootstrap-realization-
// launcher.js CLI (as a subprocess) that FitDriverScripts/*.py already
// shells out to -- this reuses its dispatch mechanics (--local worker pool
// sizing, --apiKey/--computeGroups/--slicePrice real-DCP dispatch, the
// ENOPROGRESS/Buffer/payload-binding fixes) unchanged rather than
// duplicating any of it.
//
// estimate-geometry, plus two explicit dispatchable subcommands matching
// two separate DCP jobs (and the deliberate pause between them -- the
// bootstrap job's cost scales with nBootstraps, so this never auto-chains):
//
//   node run-galaxy-fit-cli.js estimate-geometry --cube <path> [--out <path>]
//     (always local, in-process -- see that branch's own comment for why)
//
//   node run-galaxy-fit-cli.js initial-fit \
//     --cube <path> --mask <path> --pa <deg> --inc <deg> \
//     [--objName <name>] [--out <path>] \
//     [--apiKey <key> | --local] [--computeGroups a,b] [--slicePrice 1.0]
//
//   node run-galaxy-fit-cli.js bootstrap \
//     --cube <path> --initialFitResult <path from the step above> \
//     --nBootstraps <n> [--seed <n>] \
//     [--objName <name>] [--out <path>] \
//     [--apiKey <key> | --local] [--computeGroups a,b] [--slicePrice 1.0]
//
// galaxy-fit.html is the browser equivalent of this same flow -- both call
// the same buildFitPayloads.js functions; only the dispatch entry point
// differs (this file uses dcp-client's Node identity.set(apiKey) path via
// the launcher's own CLI tail; the browser uses dcp.wallet.get()).
// =============================================================================

const fs = require('fs');
const path = require('path');
const { execFileSync } = require('child_process');
const { buildGeometryEstimatePayload, buildInitialFitPayload, buildBootstrapPayload } = require('./src/PayloadBuilder/buildFitPayloads');
const { defaultFittingOptions, SOFIA_PAR_TEMPLATE_RELATIVE_PATH } = require('./src/PipelineConfig/defaultFittingOptions');

const LAUNCHER = path.join(__dirname, 'bootstrap-realization-launcher.js');

function arg(name, def) {
  const eq = process.argv.find((a) => a.startsWith('--' + name + '='));
  if (eq) return eq.slice(name.length + 3);
  const i = process.argv.indexOf('--' + name);
  return i >= 0 ? process.argv[i + 1] : def;
}
function flagPresent(name) {
  return process.argv.includes('--' + name);
}

// Shells out to the existing launcher's own CLI (same one FitDriverScripts/
// *.py already invokes via subprocess) -- see this file's header for why.
function dispatch(extraArgs, outPath) {
  const args = [LAUNCHER, ...extraArgs, '--out', outPath];
  const apiKey = arg('apiKey') || process.env.DCP_API_KEY;
  if (apiKey && !flagPresent('local')) {
    args.push('--apiKey', apiKey);
  } else {
    if (!apiKey) console.log("No --apiKey/DCP_API_KEY set -- falling back to --local (direct in-process execution, no dispatch)");
    // --local expects a following value even for the single-job initial-fit
    // path (the launcher's shared arg() parser takes argv[i+1]
    // unconditionally) -- see RunInitialFitDCP.py's LaunchInitialFitJob for
    // the same gotcha, matched here.
    args.push('--local', arg('local', '1'));
  }
  const computeGroups = arg('computeGroups');
  if (computeGroups) args.push('--computeGroups', computeGroups);
  const slicePrice = arg('slicePrice');
  if (slicePrice) args.push('--slicePrice', slicePrice);
  execFileSync('node', args, { stdio: 'inherit', cwd: __dirname });
  return JSON.parse(fs.readFileSync(outPath, 'utf8'));
}

function usageAndExit(msg) {
  console.error(msg);
  console.error('\nUsage:');
  console.error('  node run-galaxy-fit-cli.js estimate-geometry --cube <path> [--out <path>]');
  console.error('  node run-galaxy-fit-cli.js initial-fit --cube <path> --mask <path> --pa <deg> --inc <deg> [--objName <name>] [--out <path>] [--apiKey <key> | --local] [--computeGroups a,b] [--slicePrice 1.0]');
  console.error('  node run-galaxy-fit-cli.js bootstrap --cube <path> --initialFitResult <path> --nBootstraps <n> [--seed <n>] [--objName <name>] [--out <path>] [--apiKey <key> | --local] [--computeGroups a,b] [--slicePrice 1.0]');
  process.exit(1);
}

const mode = process.argv[2];

if (mode === 'estimate-geometry') {
  // Always local, in-process -- no dispatch machinery at all. Measured well
  // under a second on a real test cube (SoFiA source-finding on one
  // un-resampled cube, no fit), so there's nothing to gain from a real DCP
  // dispatch here, only its scheduling/sandbox overhead.
  const { estimateGeometry } = require('./bootstrap-realization-launcher');
  const cubePath = arg('cube');
  const outPath = arg('out', 'geometry_estimate.json');
  if (!cubePath) {
    usageAndExit('estimate-geometry: --cube is required.');
  }

  const sofiaParTemplateText = fs.readFileSync(path.join(__dirname, SOFIA_PAR_TEMPLATE_RELATIVE_PATH), 'utf8');
  const payload = buildGeometryEstimatePayload({
    cubeBytes: new Uint8Array(fs.readFileSync(cubePath)),
    sofiaParTemplateText,
  });
  console.log('Estimating geometry from cube...');
  estimateGeometry(0, payload).then((result) => {
    fs.writeFileSync(outPath, JSON.stringify(result));
    if (result.sofiaFailed) {
      console.error('Estimate failed:', result.error || '(no error message)');
      process.exit(1);
    }
    console.log(`\nPA=${result.paEstDeg.toFixed(3)} deg  Inc=${result.incEstDeg.toFixed(3)} deg${result.objectName ? '  object=' + result.objectName : ''}`);
    if (result.velocityAxisWarning) {
      console.warn(`WARNING: ${result.velocityAxisWarning}`);
    }
    console.log(`Result written to ${outPath}`);
  }).catch((e) => { console.error(e); process.exit(1); });
} else if (mode === 'initial-fit') {
  const cubePath = arg('cube');
  const maskPath = arg('mask');
  const paEstDeg = parseFloat(arg('pa'));
  const incEstDeg = parseFloat(arg('inc'));
  const objName = arg('objName', 'galaxy');
  const outPath = arg('out', objName + '_initial_fit_result.json');
  if (!cubePath || !maskPath || !Number.isFinite(paEstDeg) || !Number.isFinite(incEstDeg)) {
    usageAndExit('initial-fit: --cube, --mask, --pa, and --inc are required.');
  }

  const options = defaultFittingOptions();
  const payload = buildInitialFitPayload({
    cubeBytes: new Uint8Array(fs.readFileSync(cubePath)),
    maskBytes: new Uint8Array(fs.readFileSync(maskPath)),
    paEstDeg, incEstDeg, options,
  });
  const payloadPath = outPath + '.payload.json';
  fs.writeFileSync(payloadPath, JSON.stringify(payload));
  console.log(`Dispatching initial fit for ${objName}...`);
  try {
    const result = dispatch(['--initialFit', '--payload', payloadPath, '--jobName', objName], outPath);
    if (!result.FITAchieved) {
      console.error('Initial fit did not converge / failed:', result.error || '(no error message)');
      process.exit(1);
    }
    console.log(`\nInitial fit converged: chi2=${result.chi2.toFixed(3)}  rings=${result.R.length}`);
    console.log(`  Inc=${result.INCLINATION[0].toFixed(2)} deg  PA=${result.POSITIONANGLE[0].toFixed(2)} deg  VSys=${result.VSYS[0].toFixed(1)} km/s`);
    console.log(`\nResult written to ${outPath}`);
    console.log(`Review it, then pass it as --initialFitResult to the 'bootstrap' subcommand.`);
  } finally {
    fs.rmSync(payloadPath, { force: true });
  }
} else if (mode === 'bootstrap') {
  const cubePath = arg('cube');
  const initialFitResultPath = arg('initialFitResult');
  const nBootstraps = parseInt(arg('nBootstraps'), 10);
  const seed = arg('seed') ? parseInt(arg('seed'), 10) : 0;
  const objName = arg('objName', 'galaxy');
  const outPath = arg('out', objName + '_bootstrap_results.json');
  if (!cubePath || !initialFitResultPath || !Number.isFinite(nBootstraps)) {
    usageAndExit('bootstrap: --cube, --initialFitResult, and --nBootstraps are required.');
  }

  const initialFitResult = JSON.parse(fs.readFileSync(initialFitResultPath, 'utf8'));
  const sofiaParTemplateText = fs.readFileSync(path.join(__dirname, SOFIA_PAR_TEMPLATE_RELATIVE_PATH), 'utf8');
  const options = defaultFittingOptions();
  const payload = buildBootstrapPayload({
    initialFitResult,
    cubeBytes: new Uint8Array(fs.readFileSync(cubePath)),
    sofiaParTemplateText, bootstrapSeed: seed, options,
  });
  const payloadPath = outPath + '.payload.json';
  fs.writeFileSync(payloadPath, JSON.stringify(payload));
  console.log(`Dispatching ${nBootstraps} bootstrap realization(s) for ${objName}...`);
  try {
    const results = dispatch(['--nBootstraps', String(nBootstraps), '--payload', payloadPath, '--jobName', objName], outPath);
    const nOk = results.filter((r) => !r.sofiaFailed).length;
    console.log(`\n${nOk}/${nBootstraps} realizations succeeded. Results written to ${outPath}.`);
  } finally {
    fs.rmSync(payloadPath, { force: true });
  }
} else {
  usageAndExit(`Unknown or missing subcommand: ${mode || '(none)'}`);
}
