/**
 * @file        test-published-package.js
 * @description Run this AFTER publishing package.dcp, to confirm
 *              cfitsio4wasm/cfitsio-wasm.js works as an actual published
 *              package -- job.requires(['cfitsio4wasm/cfitsio-wasm.js']) +
 *              require('cfitsio-wasm.js') (bare filename, no path), per
 *              the package manager's usage convention (see
 *              docs/patching-wasm-libraries-for-dcp.pdf in the Edequity
 *              repo). If this passes, "require it and just use it" is
 *              real for anyone else who adds the same job.requires()
 *              line -- they never need to know about Emscripten, MEMFS,
 *              or the byte-order fix in ../README.md.
 *
 *              This is deliberately a separate, narrower check than
 *              ../test-cfitsio-worker.js: that one already confirmed the
 *              unpublished build's logic locally. This one's only job is
 *              to confirm the *published* artifact -- the bravojs wrap,
 *              specifically -- resolves and runs the same way, for real,
 *              on a real worker.
 *
 * @usage       node test-published-package.js --apiKey=0x<identity> [--computeGroup=key,secret]
 */
'use strict';

const fs = require('fs');
const path = require('path');

async function testCfitsio({ cubeBase64 }) {
  const setProgress = (p) => { if (typeof progress === 'function') progress(p); };
  setProgress(0.0);

  const report = {};
  const t0 = Date.now();

  try {
    // The whole point of publishing: this is ALL a consumer needs to write.
    // No wasm/, no Emscripten awareness, no MEMFS/byte-order gotchas.
    const cfitsio = require('cfitsio-wasm.js');
    setProgress(0.1);

    const cubeBytes = Uint8Array.from(atob(cubeBase64), (c) => c.charCodeAt(0));
    setProgress(0.2);

    const info = await cfitsio.readImageInfo(cubeBytes);
    report.naxes = info.naxes;
    setProgress(0.5);

    const written = await cfitsio.writeImageDouble([4, 3], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
    const readBack = await cfitsio.readImageDouble(written);
    report.roundTripExact = JSON.stringify(readBack.data) === JSON.stringify([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
    setProgress(0.9);

    report.pass = report.roundTripExact === true;
    report.totalMs = Date.now() - t0;
  } catch (e) {
    report.error = (e && e.message) || String(e);
    report.totalMs = Date.now() - t0;
  }

  setProgress(1.0);
  return report;
}

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

  const CUBE = path.join(
    __dirname,
    '../../../../3KIDNASTests/TestData/WALLABY_Test_sources/WALLABY_J103554-475245/WALLABY_J103554-475245_cube.fits'
  );
  const cubeBase64 = fs.readFileSync(CUBE).toString('base64');

  const job = compute.for([{ cubeBase64 }], testCfitsio, []);
  job.requires(['cfitsio4wasm/cfitsio-wasm.js']);
  job.computeGroups = computeGroup ? [computeGroup] : [{ joinKey: 'public' }];
  job.public = {
    name: 'cfitsio4wasm published-package test',
    description: 'Confirms the published cfitsio4wasm package works via job.requires()',
  };

  job.on('readystatechange', (ev) => console.log(`Ready state: ${ev}`));
  job.on('accepted', () => console.log(`  Job id: ${job.id}\n  Awaiting result...`));
  job.on('error', (error) => console.error('  Job error:', error));
  job.on('nofunds', (ev) => console.log(ev));
  job.on('result', (ev) => {
    console.log('\n=== RESULT ===\n');
    console.log(JSON.stringify(ev.result, null, 2));
  });

  const [result] = await job.exec();
  const ok = !result.error && result.pass && JSON.stringify(result.naxes) === JSON.stringify([33, 35, 53]);
  console.log(ok ? '\nPASS: published cfitsio4wasm package runs correctly.' : '\nFAIL: see result above.');
  if (!ok) process.exitCode = 1;
}

require('dcp-client').init().then(main).catch((e) => { console.error(e); process.exit(1); });
