/**
 * @file        test-published-package.js
 * @description Run this AFTER publishing package.dcp, to confirm
 *              sofia2wasm/sofia-wasm.js works as an actual published
 *              package -- job.requires(['sofia2wasm/sofia-wasm.js']) +
 *              require('sofia-wasm.js') (bare filename, no path), per the
 *              package manager's usage convention (see
 *              docs/patching-wasm-libraries-for-dcp.pdf in the Edequity
 *              repo). If this passes, "require it and just use it" is real
 *              for anyone else who adds the same job.requires() line --
 *              they never need to know about Emscripten, MEMFS, callMain,
 *              or any of the rest of the porting story in ../README.md.
 *
 *              This is deliberately a separate, narrower check than
 *              ../test-sofia-worker.js: that one already confirmed the
 *              unpublished build runs correctly on a real worker (job id
 *              8cG9uHpYtEGFDcyeo6hzvM, both slices PASS). This one's only
 *              job is to confirm the *published* artifact -- the bravojs
 *              wrap, specifically -- resolves and runs the same way.
 *
 * @usage       node test-published-package.js --apiKey=0x<identity> [--computeGroup=key,secret]
 */
'use strict';

const fs = require('fs');
const path = require('path');

async function testSofia(_input, cubeBase64, parText) {
  const setProgress = (p) => { if (typeof progress === 'function') progress(p); };
  setProgress(0.0);

  const report = {};
  const t0 = Date.now();

  try {
    // The whole point of publishing: this is ALL a consumer needs to write.
    // No wasm/, no Emscripten awareness, no MEMFS/callMain plumbing.
    const sofia = require('sofia-wasm.js');
    setProgress(0.1);

    const cubeBytes = Uint8Array.from(atob(cubeBase64), (c) => c.charCodeAt(0));
    setProgress(0.3);

    const { exitCode, files } = await sofia.run({ cube: cubeBytes, par: parText });
    setProgress(0.9);

    report.exitCode = exitCode;
    report.outputFiles = [...files.keys()];
    if (files.has('result_cat.txt')) {
      report.catalog = new TextDecoder('utf-8').decode(files.get('result_cat.txt'));
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
  const parText = `
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

  const job = compute.for([0], testSofia, [cubeBase64, parText]);
  job.requires(['sofia2wasm/sofia-wasm.js']);
  job.computeGroups = computeGroup ? [computeGroup] : [{ joinKey: 'public' }];
  job.public = {
    name: 'sofia published-package test',
    description: 'Confirms the published sofia package works via job.requires()',
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
  const ok = !result.error && result.exitCode === 0 && result.sourceCount > 0;
  console.log(ok ? '\nPASS: published sofia package runs correctly.' : '\nFAIL: see result above.');
  if (!ok) process.exitCode = 1;
}

require('dcp-client').init().then(main).catch((e) => { console.error(e); process.exit(1); });
