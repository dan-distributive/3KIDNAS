'use strict';

// =============================================================================
// plan_sweep.js
// Empirical validation harness for Planner/chooseDecomposition.js: for a wide
// range of N, asks the REAL compiled FFTW3 (via ground_truth_harness's `plan`
// mode) what it actually chose, and checks whether chooseComplex(N)/
// chooseReal(N)'s predicted top-level solver appears in that ground-truth
// text. This is a token-membership check, not a full structural diff --
// ground truth for a given N0 sometimes wraps the exact same underlying
// decomposition in a `dft-buffered` layer depending on the outer 2D
// transform's vector count (confirmed by direct comparison: N=216's column
// pass showed identical dft-ct-dit/8(t2_8 ... t1_3 ... n1_9) structure both
// bare and dft-buffered-wrapped, depending on context) -- buffering is a
// vector-loop-level wrapping decision this port deliberately doesn't predict
// (see chooseDecomposition.js header), so token-membership is the right
// granularity: did we get the underlying algorithm/codelet choice right.
//
// Run: node plan_sweep.js [maxN]
// =============================================================================

const { execFileSync } = require('child_process');
const path = require('path');
const { chooseComplex, chooseReal } = require('../Planner/chooseDecomposition');

const HARNESS = path.join(__dirname, 'ground_truth_harness');

function getPlanText(n0, n1) {
  const out = execFileSync(HARNESS, ['plan', String(n0), String(n1)], { encoding: 'utf8' });
  const r2cMatch = out.match(/=== r2c plan for \d+x\d+ ===\n([\s\S]*?)\n\n=== c2r/);
  return r2cMatch ? r2cMatch[1] : '';
}

// Extract every codelet-name / solver-name token that appears in a ground
// truth plan block, e.g. "n1_37", "t1_6", "r2cf_13", "dft-generic-43",
// "dft-rader-37", "dft-bluestein-103".
function extractTokens(planText) {
  const tokens = new Set();
  let m;
  const quoted = /"([a-zA-Z0-9_]+)"/g;
  while ((m = quoted.exec(planText))) tokens.add(m[1]);
  const named = /\((dft-generic|dft-rader|dft-bluestein|rdft-generic-r2hc|rdft-generic-hc2r)-(\d+)/g;
  while ((m = named.exec(planText))) tokens.add(`${m[1]}-${m[2]}`);
  return tokens;
}

// Reduce our own prediction's describe() string to the same kind of token
// set, so membership-checking is symmetric.
function predictionTokens(describeStr) {
  const tokens = new Set();
  let m;
  const codelet = /\b(n1_\d+|t1_\d+|t2_\d+|r2cf_\d+|r2cb_\d+|hf_\d+|hb_\d+)\b/g;
  while ((m = codelet.exec(describeStr))) tokens.add(m[1]);
  const named = /\b(dft-generic-\d+|dft-rader-\d+|dft-bluestein-\d+\/nb=\d+|rdft-generic-r2hc-\d+|rdft-generic-hc2r-\d+)\b/g;
  while ((m = named.exec(describeStr))) tokens.add(m[1].replace(/\/nb=\d+/, ''));
  return tokens;
}

function checkComplex(n) {
  const planText = getPlanText(n, n);
  const groundTokens = extractTokens(planText);
  const pred = chooseComplex(n);
  const predTokens = predictionTokens(pred.describe());
  const hit = [...predTokens].some((t) => groundTokens.has(t));
  return { n, ok: hit, predicted: pred.describe(), groundTokens: [...groundTokens] };
}

function checkReal(n) {
  const planText = getPlanText(n, n);
  const groundTokens = extractTokens(planText);
  const pred = chooseReal('R2HC', n);
  const predTokens = predictionTokens(pred.describe());
  const hit = [...predTokens].some((t) => groundTokens.has(t));
  return { n, ok: hit, predicted: pred.describe(), groundTokens: [...groundTokens] };
}

function main() {
  const maxN = parseInt(process.argv[2] || '250', 10);
  const complexResults = [];
  const realResults = [];

  for (let n = 2; n <= maxN; n++) {
    try {
      complexResults.push(checkComplex(n));
    } catch (e) {
      complexResults.push({ n, ok: null, error: e.message });
    }
    try {
      realResults.push(checkReal(n));
    } catch (e) {
      realResults.push({ n, ok: null, error: e.message });
    }
  }

  for (const [label, results] of [['COMPLEX (column, dft-*)', complexResults], ['REAL R2HC (row, rdft2-*)', realResults]]) {
    const matched = results.filter((r) => r.ok === true).length;
    const mismatched = results.filter((r) => r.ok === false);
    const errored = results.filter((r) => r.ok === null);
    console.log(`\n=== ${label}: ${matched}/${results.length} matched, ${mismatched.length} mismatched, ${errored.length} errored ===`);
    for (const r of mismatched) {
      console.log(`  MISMATCH n=${r.n}: predicted "${r.predicted}" -- not found in ground truth tokens [${r.groundTokens.join(', ')}]`);
    }
    for (const r of errored) {
      console.log(`  ERROR n=${r.n}: ${r.error}`);
    }
  }
}

main();
