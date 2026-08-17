/**
 * @file        test-fftw-worker.js
 * @description Confirms the FFTW 3.3.8 + fdlibm wasm build (fftw-wasm.js /
 *              fftw-module.js) runs correctly inside a real, remote DCP
 *              worker sandbox -- not just local Node (see README.md's
 *              "Verifying it"). Unlike SoFiA-2's CLI/file-I/O shape, this
 *              is a pure numeric library call: each slice generates a
 *              signal, runs it through both the complex-to-complex and
 *              real-to-complex DFT entry points, and checks the result
 *              against an independent naive O(n^2) DFT computed in the
 *              same work function -- no dependency on FFTW's own
 *              correctness to check FFTW's own correctness.
 *
 *              INPUT SET: four signals of different, deliberately awkward
 *              sizes (16, 17, 64, 100 -- power-of-2 and not), one per
 *              worker. STATIC ARGUMENT: the shared numerical tolerance
 *              every worker checks its result against, the same
 *              INPUT-SET-vs-STATIC-ARGS distinction this project already
 *              draws in bootstrap-fit-launcher.js.
 *
 * @usage       node test-fftw-worker.js --apiKey=0x<identity> [--computeGroup=key,secret]
 */
'use strict';

// -----------------------------------------------------------------------
// WORK FUNCTION -- module scope, no closures over outer variables: DCP's
// compute.for ships this function's source to remote worker sandboxes,
// which don't have access to this file's module scope.
// -----------------------------------------------------------------------
async function testFftw({ name, re, im }, tolerance) {
  const setProgress = (p) => { if (typeof progress === 'function') progress(p); };
  setProgress(0.0);

  const report = { name, n: re.length };
  const t0 = Date.now();

  try {
    const fftw = require('./fftw-wasm');
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
    setProgress(0.6);
    report.c2cMaxDiff = Math.max(maxAbsDiff(c2c.re, ref.reOut), maxAbsDiff(c2c.im, ref.imOut));

    const r2c = await fftw.r2c1d(re);
    setProgress(0.9);
    const nc = r2c.re.length;
    report.r2cMaxDiff = Math.max(
      maxAbsDiff(r2c.re, ref.reOut.slice(0, nc)),
      maxAbsDiff(r2c.im, ref.imOut.slice(0, nc))
    );

    report.pass = report.c2cMaxDiff < tolerance && report.r2cMaxDiff < tolerance;
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
function makeSignal(n, freqs) {
  const re = Array.from({ length: n }, (_, t) =>
    freqs.reduce((sum, [k, amp]) => sum + amp * Math.sin((2 * Math.PI * k * t) / n), 0.3)
  );
  const im = new Array(n).fill(0);
  return { re, im };
}

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

  // INPUT SET -- four deliberately awkward sizes (power-of-2 and not).
  const inputSet = [
    { name: 'n16-single-tone', ...makeSignal(16, [[3, 1]]) },
    { name: 'n17-odd', ...makeSignal(17, [[5, 2]]) },
    { name: 'n64-two-tones', ...makeSignal(64, [[6, 1], [21, 0.5]]) },
    { name: 'n100-composite', ...makeSignal(100, [[9, 1], [30, 0.7], [44, 0.3]]) },
  ];

  // STATIC ARG -- numerical tolerance every worker checks against.
  const TOLERANCE = 1e-8;

  const job = compute.for(inputSet, testFftw, [TOLERANCE]);

  // REQUIRED MODULES -- job.requires(['./fftw-wasm']) alone is enough:
  // fftw-wasm.js's own require('./fftw-module') is a transitive
  // dependency the walker picks up from that one listed path (same
  // mechanism documented in sofia-wasm.js / duckdb-wasm.js).
  job.requires(['./fftw-wasm']);

  job.computeGroups = computeGroup ? [computeGroup] : [{ joinKey: 'public' }];

  job.public = {
    name: '\u{1F30A} FFTW3 wasm worker test',
    description: 'Confirms the FFTW 3.3.8 wasm build runs inside a real DCP worker (DFT correctness checked against an independent naive DFT, per slice)',
    link: 'http://fftw.org/',
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
      console.log(`  ${r.name} (n=${r.n}): pass=${r.pass}  c2cMaxDiff=${r.c2cMaxDiff.toExponential(3)}  r2cMaxDiff=${r.r2cMaxDiff.toExponential(3)}  (${r.totalMs} ms)`);
    }
  });

  const results = await job.exec();

  console.log('\n=== SUMMARY ===');
  let allOk = true;
  for (const r of results) {
    const ok = !r.error && r.pass;
    allOk = allOk && ok;
    console.log(`  ${r.name}: ${ok ? 'PASS' : 'FAIL'}`);
    if (!ok) console.log(`    ${r.error || `c2cMaxDiff=${r.c2cMaxDiff}, r2cMaxDiff=${r.r2cMaxDiff}`}`);
  }
  console.log(allOk ? '\nAll slices passed: FFTW3 wasm runs correctly in a real DCP worker.' : '\nAt least one slice failed -- see details above.');
  if (!allOk) process.exitCode = 1;
}

require('dcp-client').init().then(main).catch((e) => { console.error(e); process.exit(1); });
