/**
 * @file        test-cfitsio-worker.js
 * @description Confirms the cfitsio 4.6.3 + zlib wasm build
 *              (cfitsio-wasm.js / cfitsio-module.js) runs correctly
 *              inside a real, remote DCP worker sandbox -- not just
 *              local Node (see README.md's "Verifying it").
 *
 *              INPUT SET: two real WALLABY HI cubes (base64-encoded,
 *              same shard-by-file shape as SoFiA-2's own test job), one
 *              per worker. STATIC ARGUMENT: the shared header keyword
 *              every worker reads (a real use of the pattern: which
 *              keyword to check is shared policy, not per-file data).
 *              Each worker reads image dimensions, the shared keyword,
 *              and pixel min/max/mean, then does an exact round-trip
 *              write+read of a small synthetic image to confirm both
 *              read AND write paths work remotely, not just read.
 *
 * @usage       node test-cfitsio-worker.js --apiKey=0x<identity> [--computeGroup=key,secret]
 */
'use strict';

const fs = require('fs');
const path = require('path');

// -----------------------------------------------------------------------
// WORK FUNCTION -- module scope, no closures over outer variables: DCP's
// compute.for ships this function's source to remote worker sandboxes,
// which don't have access to this file's module scope.
// -----------------------------------------------------------------------
async function testCfitsio({ name, cubeBase64 }, keyName) {
  const setProgress = (p) => { if (typeof progress === 'function') progress(p); };
  setProgress(0.0);

  const report = { name };
  const t0 = Date.now();

  try {
    const cfitsio = require('./cfitsio-wasm');
    setProgress(0.1);

    const cubeBytes = Uint8Array.from(atob(cubeBase64), (c) => c.charCodeAt(0));
    setProgress(0.2);

    const info = await cfitsio.readImageInfo(cubeBytes);
    report.naxes = info.naxes;
    setProgress(0.4);

    try {
      report[keyName] = await cfitsio.readKeyDouble(cubeBytes, keyName);
    } catch (e) {
      report[keyName] = null; // keyword may legitimately not exist on every file
    }
    setProgress(0.5);

    const { data } = await cfitsio.readImageDouble(cubeBytes);
    let min = Infinity, max = -Infinity, sum = 0, finiteCount = 0;
    for (const v of data) {
      if (Number.isFinite(v)) { min = Math.min(min, v); max = Math.max(max, v); sum += v; finiteCount++; }
    }
    report.pixelMin = min;
    report.pixelMax = max;
    report.pixelMean = sum / finiteCount;
    report.finiteFraction = finiteCount / data.length;
    setProgress(0.8);

    // Exact round-trip: confirms the write path works remotely too, not just read.
    const wNaxes = [4, 3];
    const wData = [1.5, -2.25, 3, 0.125, 100, -0.5, 7.75, 42, -13.375, 0, 1, -1];
    const written = await cfitsio.writeImageDouble(wNaxes, wData);
    const readBack = await cfitsio.readImageDouble(written);
    report.roundTripExact = JSON.stringify(readBack.data) === JSON.stringify(wData);
    setProgress(0.95);

    report.pass = report.finiteFraction === 1 && report.roundTripExact;
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

  // INPUT SET -- two real WALLABY cubes, base64-encoded.
  const TEST_DATA = path.join(__dirname, '../../../3KIDNASTests/TestData/WALLABY_Test_sources');
  const CUBES = [
    { name: 'WALLABY_J103554-475245', file: path.join(TEST_DATA, 'WALLABY_J103554-475245/WALLABY_J103554-475245_cube.fits') },
    { name: 'WALLABY_J103458-495128', file: path.join(TEST_DATA, 'WALLABY_J103458-495128/WALLABY_J103458-495128_cube.fits') },
  ];
  const inputSet = CUBES.map(({ name, file }) => ({
    name,
    cubeBase64: fs.readFileSync(file).toString('base64'),
  }));

  // STATIC ARG -- header keyword shared by every worker.
  const KEY_NAME = 'CRPIX1';

  const job = compute.for(inputSet, testCfitsio, [KEY_NAME]);
  job.requires(['./cfitsio-wasm']);
  job.computeGroups = computeGroup ? [computeGroup] : [{ joinKey: 'public' }];
  job.public = {
    name: '\u{1F5C2}\u{FE0F} cfitsio wasm worker test',
    description: 'Confirms the cfitsio wasm build runs inside a real DCP worker (real FITS read/write, WALLABY HI cubes)',
    link: 'https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html',
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
      console.log(`  ${r.name}: naxes=${JSON.stringify(r.naxes)}  ${KEY_NAME}=${r[KEY_NAME]}  min=${r.pixelMin.toFixed(6)}  max=${r.pixelMax.toFixed(6)}  roundTrip=${r.roundTripExact}  (${r.totalMs} ms)`);
    }
  });

  const results = await job.exec();

  console.log('\n=== SUMMARY ===');
  let allOk = true;
  for (const r of results) {
    const ok = !r.error && r.pass;
    allOk = allOk && ok;
    console.log(`  ${r.name}: ${ok ? 'PASS' : 'FAIL'}`);
    if (!ok) console.log(`    ${r.error || JSON.stringify(r)}`);
  }
  console.log(allOk ? '\nAll slices passed: cfitsio wasm runs correctly in a real DCP worker.' : '\nAt least one slice failed -- see details above.');
  if (!allOk) process.exitCode = 1;
}

require('dcp-client').init().then(main).catch((e) => { console.error(e); process.exit(1); });
