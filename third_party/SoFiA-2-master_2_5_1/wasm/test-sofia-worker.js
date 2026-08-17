/**
 * @file        test-sofia-worker.js
 * @description Confirms the SoFiA-2 + wcslib wasm build (sofia-wasm.js /
 *              sofia-module.js) actually runs inside a real, remote DCP
 *              worker sandbox -- not just local Node. That's the one gap
 *              local testing (see ../wasm/README.md) couldn't close: the
 *              build has no Buffer, process.*, or fetch() usage, and it
 *              doesn't force a specific Emscripten ENVIRONMENT, the same
 *              profile that worked for the published duckdbwasm package,
 *              but inference from the code isn't the same as a real run.
 *
 *              Ships two small real WALLABY HI cubes (~260-290 KB each,
 *              from 3KIDNASTests/TestData) as the input set, one per
 *              worker, base64-encoded since compute.for's input set goes
 *              through DCP's normal (de)serialization. The parameter
 *              file text is a static argument: identical pipeline
 *              settings shared by every slice, the same distinction
 *              this project already draws in bootstrap-fit-launcher.js
 *              (INPUT SET = per-slice data, STATIC ARGS = shared config).
 *
 *              Each worker runs SoFiA-2's real S+C finder + linker + WCS
 *              parameterisation on its cube and returns the detected
 *              source count and catalogue text -- enough to eyeball that
 *              real science happened, not just "exit code 0". Local
 *              Node runs against these same two cubes already produced
 *              exactly one source each with sane sky coordinates; a
 *              passing worker run should match that.
 *
 * @usage       node test-sofia-worker.js --apiKey=0x<identity> [--computeGroup=key,secret]
 */
'use strict';

const fs = require('fs');
const path = require('path');

// -----------------------------------------------------------------------
// WORK FUNCTION -- module scope, no closures over outer variables: DCP's
// compute.for ships this function's source to remote worker sandboxes,
// which don't have access to this file's module scope (same rule
// bootstrap-fit-launcher.js documents for its own work function).
// -----------------------------------------------------------------------
async function testSofia({ name, cubeBase64 }, parText) {
  const setProgress = (p) => { if (typeof progress === 'function') progress(p); };
  setProgress(0.0);

  const report = { name };
  const t0 = Date.now();

  try {
    const sofia = require('./sofia-wasm');
    setProgress(0.1);

    const cubeBytes = Uint8Array.from(atob(cubeBase64), (c) => c.charCodeAt(0));
    setProgress(0.3);

    const { exitCode, log, files } = await sofia.run({ cube: cubeBytes, par: parText });
    setProgress(0.9);

    report.exitCode = exitCode;
    report.outputFiles = [...files.keys()];
    report.logTail = log.split('\n').slice(-15).join('\n');
    if (files.has('result_cat.txt')) {
      report.catalog = new TextDecoder('utf-8').decode(files.get('result_cat.txt'));
      // Real catalogue rows start with a quoted source name; comment/blank
      // lines don't. Cheap way to report a source count without parsing
      // the whole fixed-width table.
      report.sourceCount = report.catalog.split('\n').filter((l) => l.trim().startsWith('"')).length;
    }
    report.totalMs = Date.now() - t0;
  } catch (e) {
    report.error = (e && e.message) || String(e);
    report.totalMs = Date.now() - t0;
  }

  setProgress(1.0);
  return report;
}

// -----------------------------------------------------------------------
// DRIVER
// -----------------------------------------------------------------------
async function main() {
  const identity = require('dcp/identity');
  const compute = require('dcp/compute');

  // --apiKey=0x...            identity / API key (required)
  // --computeGroup=key,secret join key + secret (optional; default: public group)
  function getFlag(name) {
    const pfx = `--${name}=`;
    const hit = process.argv.find((a) => a.startsWith(pfx));
    return hit ? hit.slice(pfx.length) : undefined;
  }

  const apiKey = getFlag('apiKey');
  const cg = getFlag('computeGroup');
  if (!apiKey) {
    console.error('ERROR: --apiKey=0x... is required');
    process.exit(1);
  }
  let computeGroup;
  if (cg) {
    const [joinKey, joinSecret] = cg.split(',');
    computeGroup = joinSecret ? { joinKey, joinSecret } : { joinKey };
  }

  await identity.set(apiKey);

  // INPUT SET -- one real WALLABY cube per worker, base64-encoded.
  const TEST_DATA = path.join(__dirname, '../../../3KIDNASTests/TestData/WALLABY_Test_sources');
  const CUBES = [
    { name: 'WALLABY_J103554-475245', file: path.join(TEST_DATA, 'WALLABY_J103554-475245/WALLABY_J103554-475245_cube.fits') },
    { name: 'WALLABY_J103458-495128', file: path.join(TEST_DATA, 'WALLABY_J103458-495128/WALLABY_J103458-495128_cube.fits') },
  ];
  const inputSet = CUBES.map(({ name, file }) => ({
    name,
    cubeBase64: fs.readFileSync(file).toString('base64'),
  }));

  // STATIC ARG -- pipeline settings shared by every slice. See
  // ../template_par_file.par for the full parameter grammar.
  const PAR_TEXT = `
pipeline.verbose     = false
input.data           = /work/cube.fits

scfind.enable        = true
linker.enable        = true
parameter.enable     = true
parameter.wcs        = true

output.directory     = /work/out
output.filename      = result
output.writeCatASCII = true
output.overwrite     = true
`;

  const job = compute.for(inputSet, testSofia, [PAR_TEXT]);

  // REQUIRED MODULES -- job.requires(['./sofia-wasm']) alone is enough:
  // sofia-wasm.js's own require('./sofia-module.js') is a transitive
  // dependency the walker picks up from that one listed path (same
  // mechanism sofia-wasm.js's own header comment describes, itself
  // following duckdb-wasm.js's precedent in the Edequity project).
  job.requires(['./sofia-wasm']);

  job.computeGroups = computeGroup ? [computeGroup] : [{ joinKey: 'public' }];

  job.public = {
    name: '🔭 SoFiA-2 wasm worker test',
    description: 'Confirms the SoFiA-2 + wcslib wasm build runs inside a real DCP worker (source finding on WALLABY HI cubes)',
    link: 'https://github.com/SoFiA-Admin/SoFiA-2',
  };

  job.on('readystatechange', (ev) => console.log(`Ready state: ${ev}`));
  job.on('accepted', () => console.log(`  Job id: ${job.id}\n  Awaiting results...`));
  job.on('error', (error) => console.error('  Job error:', error));
  job.on('nofunds', (ev) => console.log(ev));
  job.on('result', (ev) => {
    const r = ev.result;
    if (r.error) {
      console.log(`  ${r.name}: FAILED -- ${r.error}`);
    } else {
      console.log(`  ${r.name}: exit=${r.exitCode}  sources=${r.sourceCount}  files=[${r.outputFiles.join(', ')}]  (${r.totalMs} ms)`);
    }
  });

  const results = await job.exec();

  console.log('\n=== SUMMARY ===');
  let allOk = true;
  for (const r of results) {
    const ok = !r.error && r.exitCode === 0 && r.sourceCount > 0;
    allOk = allOk && ok;
    console.log(`  ${r.name}: ${ok ? 'PASS' : 'FAIL'}`);
    if (!ok) console.log(`    ${r.error || `exitCode=${r.exitCode}, sourceCount=${r.sourceCount}`}`);
  }
  console.log(allOk ? '\nAll slices passed: SoFiA-2 wasm runs correctly in a real DCP worker.' : '\nAt least one slice failed -- see details above.');
  if (!allOk) process.exitCode = 1;
}

require('dcp-client').init().then(main).catch((e) => { console.error(e); process.exit(1); });
