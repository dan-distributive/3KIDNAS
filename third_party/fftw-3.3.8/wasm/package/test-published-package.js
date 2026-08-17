/**
 * @file        test-published-package.js
 * @description Run this AFTER publishing package.dcp, to confirm
 *              fftw3wasm-v3/fftw-wasm.js works as an actual published
 *              package -- job.requires(['fftw3wasm-v3/fftw-wasm.js']) +
 *              require('fftw-wasm.js') (bare filename, no path), per the
 *              package manager's usage convention (see
 *              docs/patching-wasm-libraries-for-dcp.pdf in the Edequity
 *              repo). If this passes, "require it and just use it" is real
 *              for anyone else who adds the same job.requires() line --
 *              they never need to know about Emscripten, ccall/cwrap
 *              buffer plumbing, or fdlibm, or any of the rest of the
 *              porting story in ../README.md.
 *
 *              This is deliberately a separate, narrower check than
 *              ../test-fftw-worker.js: that one already confirmed the
 *              unpublished build runs correctly on a real worker. This
 *              one's only job is to confirm the *published* artifact --
 *              the bravojs wrap, specifically -- resolves and runs the
 *              same way.
 *
 * @usage       node test-published-package.js --apiKey=0x<identity> [--computeGroup=key,secret]
 */
'use strict';

async function testFftw({ re, im }, tolerance) {
  const setProgress = (p) => { if (typeof progress === 'function') progress(p); };
  setProgress(0.0);

  const report = {};
  const t0 = Date.now();

  try {
    // The whole point of publishing: this is ALL a consumer needs to write.
    // No wasm/, no Emscripten awareness, no ccall/cwrap/heap plumbing.
    const fftw = require('fftw-wasm.js');
    setProgress(0.1);

    function naiveDft(reIn, imIn) {
      const n = reIn.length;
      const reOut = new Array(n).fill(0);
      const imOut = new Array(n).fill(0);
      for (let k = 0; k < n; k++) {
        let sr = 0, si = 0;
        for (let t = 0; t < n; t++) {
          const a = (-2 * Math.PI * k * t) / n;
          sr += reIn[t] * Math.cos(a) - imIn[t] * Math.sin(a);
          si += reIn[t] * Math.sin(a) + imIn[t] * Math.cos(a);
        }
        reOut[k] = sr;
        imOut[k] = si;
      }
      return { reOut, imOut };
    }
    function maxAbsDiff(a, b) {
      let m = 0;
      for (let i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
      return m;
    }

    const ref = naiveDft(re, im);
    setProgress(0.3);

    const c2c = await fftw.dft1d(re, im);
    setProgress(0.5);
    report.c2cMaxDiff = Math.max(maxAbsDiff(c2c.re, ref.reOut), maxAbsDiff(c2c.im, ref.imOut));

    const r2c = await fftw.r2c1d(re);
    setProgress(0.6);
    const nc = r2c.re.length;
    report.r2cMaxDiff = Math.max(
      maxAbsDiff(r2c.re, ref.reOut.slice(0, nc)),
      maxAbsDiff(r2c.im, ref.imOut.slice(0, nc))
    );

    report.pass = report.c2cMaxDiff < tolerance && report.r2cMaxDiff < tolerance;

    // ---- Repeated-call timing: this is the part that actually matters --
    // a single call (above) can't reveal per-call cwrap overhead. Uses the
    // sync API (warmUp() once, then N calls with no re-instantiation) --
    // the exact same shape CubeKernelConvolution.js's hot path uses, and
    // the thing the cwrap-caching fix (fftw-wasm.js's cachedDftFn/
    // cachedR2cFn) targets. Real DCP worker measured 677ms/call for this
    // BEFORE the fix (found via bootstrap-realization-launcher.js's actual
    // convolveMs/convolveCalls instrumentation) -- this isolates whether
    // that's fixed, without the cost of a full bootstrap realization run.
    await fftw.warmUp();
    const N_REPEAT = 500;
    const tLoop0 = Date.now();
    for (let i = 0; i < N_REPEAT; i++) {
      fftw.dft1dSync(re, im);
      fftw.r2c1dSync(re);
    }
    report.repeatedCallMs = Date.now() - tLoop0;
    report.nRepeat = N_REPEAT;
    report.msPerCallPair = report.repeatedCallMs / N_REPEAT;

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

  const n = 20;
  const re = Array.from({ length: n }, (_, t) => Math.sin((2 * Math.PI * 7 * t) / n) + 0.2);
  const im = new Array(n).fill(0);
  const TOLERANCE = 1e-8;

  const job = compute.for([{ re, im }], testFftw, [TOLERANCE]);
  job.requires(['fftw3wasm-v3/fftw-wasm.js']);
  job.computeGroups = computeGroup ? [computeGroup] : [{ joinKey: 'public' }];
  job.public = {
    name: 'fftw3wasm-v3 published-package test',
    description: 'Confirms the published fftw3wasm-v3 package works via job.requires()',
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
  const ok = !result.error && result.pass;
  console.log(ok ? '\nPASS: published fftw3wasm-v3 package runs correctly.' : '\nFAIL: see result above.');
  if (!ok) process.exitCode = 1;
}

require('dcp-client').init().then(main).catch((e) => { console.error(e); process.exit(1); });
