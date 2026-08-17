'use strict';
// =============================================================================
// compare_fixtures.js
//
// Field-by-field diff between a Fortran-written diskfit_fixture.json
// (SingleGalaxyTests.f's DumpFittingFixture) and a JS-written one
// (bootstrap-realization-launcher.js's runInitialFit, result.fixtureJson) --
// same schema, written by two independent implementations. Verifies that
// closing the "run without Fortran" gap didn't just happen to produce a
// working fit by coincidence: this checks the actual INTERMEDIATE state
// (pvIni, param limits, TR_FittingOptions, the observed cube's own header/
// flux array, the beam) field-by-field, not just the final fit result.
//
// Usage:
//   node compare_fixtures.js <fortran_fixture.json> <js_fixture.json> [--tol 1e-5]
//
// KNOWN, EXPLAINED EXCEPTION: observedBeam.kernel will show as entirely
// mismatched (fortran=0, js=nonzero) for the cells near the kernel's centre.
// This is NOT a bug: Fortran's DumpFittingFixture runs right after
// PreGalaxyAnalysis, BEFORE GalaxyFit_Simple builds the actual kernel --
// Beam%Kernel is still zero-initialized at dump time. runInitialFit's own
// kernel is already computed by then (built explicitly, before
// initialAnalysis()). Confirmed inconsequential: NEITHER
// runBootstrapRealization NOR runInitialFit ever reads observedBeam.kernel
// back off the payload -- both always independently recompute it via
// calculate2DBeamKernel. Every other field, including nRadialCells,
// beamSigmaVector, and the kernel array's own LENGTH, matches exactly.
// =============================================================================

const fs = require('fs');

function arg(name, def) {
  const i = process.argv.indexOf('--' + name);
  return i >= 0 ? process.argv[i + 1] : def;
}

const [, , fortranPath, jsPath] = process.argv;
if (!fortranPath || !jsPath) {
  console.error('Usage: node compare_fixtures.js <fortran_fixture.json> <js_fixture.json> [--tol 1e-5]');
  process.exit(2);
}
const tol = parseFloat(arg('tol', '1e-5'));

const fortranFixture = JSON.parse(fs.readFileSync(fortranPath, 'utf8'));
const jsFixture = JSON.parse(fs.readFileSync(jsPath, 'utf8'));

const mismatches = [];
let nCompared = 0;

function relDiff(a, b) {
  if (a === b) return 0;
  const denom = Math.max(Math.abs(a), Math.abs(b), 1e-30);
  return Math.abs(a - b) / denom;
}

// Recursively walk both structures in parallel. Arrays of numbers are
// compared element-wise with a tolerance; arrays of objects (radialProfiles)
// recurse per-element; scalars compare directly (numeric with tolerance,
// everything else with strict equality); booleans compare strictly.
function compare(path, a, b) {
  if (a === null || b === null || a === undefined || b === undefined) {
    if (a !== b) mismatches.push({ path, a, b, reason: 'null/undefined mismatch' });
    return;
  }
  if (Array.isArray(a) || Array.isArray(b)) {
    if (!Array.isArray(a) || !Array.isArray(b)) {
      mismatches.push({ path, a: typeof a, b: typeof b, reason: 'one side is not an array' });
      return;
    }
    if (a.length !== b.length) {
      mismatches.push({ path, a: `length=${a.length}`, b: `length=${b.length}`, reason: 'array length mismatch' });
      return;
    }
    for (let i = 0; i < a.length; i++) compare(`${path}[${i}]`, a[i], b[i]);
    return;
  }
  if (typeof a === 'object' && typeof b === 'object') {
    const keys = new Set([...Object.keys(a), ...Object.keys(b)]);
    for (const k of keys) compare(path ? `${path}.${k}` : k, a[k], b[k]);
    return;
  }
  if (typeof a === 'number' && typeof b === 'number') {
    nCompared++;
    const d = relDiff(a, b);
    if (d > tol) mismatches.push({ path, a, b, reason: `relDiff=${d.toExponential(3)} > tol` });
    return;
  }
  if (typeof a === 'boolean' && typeof b === 'boolean') {
    nCompared++;
    if (a !== b) mismatches.push({ path, a, b, reason: 'boolean mismatch' });
    return;
  }
  // strings / anything else
  nCompared++;
  if (a !== b) mismatches.push({ path, a, b, reason: 'value mismatch' });
}

compare('', fortranFixture, jsFixture);

console.log(`Compared ${nCompared} leaf values between:`);
console.log(`  Fortran: ${fortranPath}`);
console.log(`  JS:      ${jsPath}`);
console.log(`  Tolerance (relative): ${tol}\n`);

if (mismatches.length === 0) {
  console.log(`ALL MATCH -- 0 mismatches out of ${nCompared} compared values.`);
} else {
  console.log(`${mismatches.length} MISMATCH(ES) out of ${nCompared} compared values:\n`);
  for (const m of mismatches.slice(0, 100)) {
    console.log(`  ${m.path}: fortran=${m.a} js=${m.b} (${m.reason})`);
  }
  if (mismatches.length > 100) console.log(`  ... and ${mismatches.length - 100} more`);
}
process.exit(mismatches.length === 0 ? 0 : 1);
