'use strict';

// =============================================================================
// compare.js
// Byte-exact comparator for the FFTW3-to-JS port ground-truth fixtures.
// Loads a raw double-precision fixture written by ground_truth_harness.c
// (see that file's header for the exact binary layout) and diffs it against
// the JS port's own output, element by element, using DataView reads so no
// decimal-text rounding can hide a genuine 1-ULP mismatch.
//
// Usage (as a library): compareBuffers(expectedBuffer, actualFloat64Array)
// Usage (as a CLI): node compare.js <fixture.bin> <label>
//   -- reads the fixture and just reports its own length/first few values,
//   for sanity-checking the harness output. Real comparisons happen from
//   each phase's own test script, which calls compareBuffers directly
//   after computing the JS side's output in memory.
// =============================================================================

const fs = require('fs');

// Reads a raw-double fixture file into a Float64Array.
function readFixture(path) {
  const buf = fs.readFileSync(path);
  if (buf.length % 8 !== 0) {
    throw new Error(`${path}: byte length ${buf.length} is not a multiple of 8`);
  }
  const out = new Float64Array(buf.length / 8);
  for (let i = 0; i < out.length; i++) {
    out[i] = buf.readDoubleLE(i * 8);
  }
  return out;
}

// Byte-exact (bit-exact) comparison of two Float64Arrays (or plain arrays
// of numbers). Returns { matches, mismatches: [{index, expected, actual}] }.
// mismatches is capped at maxReport entries so a totally-wrong array
// doesn't flood the output -- the count is still accurate.
function compareBuffers(expected, actual, maxReport = 20) {
  if (expected.length !== actual.length) {
    throw new Error(`length mismatch: expected ${expected.length}, actual ${actual.length}`);
  }
  const mismatches = [];
  let mismatchCount = 0;
  for (let i = 0; i < expected.length; i++) {
    const e = expected[i];
    const a = actual[i];
    // Bit-exact: NaN-safe, sign-of-zero-safe comparison via DataView round
    // trip rather than `e === a` (which treats -0 === 0, masking a real
    // sign-bit difference that would show up as a genuine mismatch upstream).
    const eb = Buffer.alloc(8); eb.writeDoubleLE(e);
    const ab = Buffer.alloc(8); ab.writeDoubleLE(a);
    if (!eb.equals(ab)) {
      mismatchCount++;
      if (mismatches.length < maxReport) {
        mismatches.push({ index: i, expected: e, actual: a, ulpDiff: ulpDistance(e, a) });
      }
    }
  }
  return { matches: mismatchCount === 0, mismatchCount, mismatches, total: expected.length };
}

// Approximate ULP distance between two doubles (same sign, both finite) --
// useful for characterizing HOW far off a mismatch is (1-ULP rounding noise
// vs. a genuine algorithmic bug look very different here).
function ulpDistance(a, b) {
  if (a === b) return 0;
  if (!Number.isFinite(a) || !Number.isFinite(b)) return Infinity;
  const ab = Buffer.alloc(8); ab.writeDoubleLE(a);
  const bb = Buffer.alloc(8); bb.writeDoubleLE(b);
  const ai = ab.readBigInt64LE();
  const bi = bb.readBigInt64LE();
  const toOrdered = (x) => (x < 0n ? 0x8000000000000000n - x : x);
  const diff = toOrdered(ai) - toOrdered(bi);
  return diff < 0n ? Number(-diff) : Number(diff);
}

function reportComparison(label, result) {
  if (result.matches) {
    console.log(`[PASS] ${label}: ${result.total} values, bit-exact`);
  } else {
    console.log(`[FAIL] ${label}: ${result.mismatchCount}/${result.total} values differ`);
    for (const m of result.mismatches) {
      console.log(`  [${m.index}] expected=${m.expected} actual=${m.actual} ulpDiff=${m.ulpDiff}`);
    }
    if (result.mismatchCount > result.mismatches.length) {
      console.log(`  ... and ${result.mismatchCount - result.mismatches.length} more`);
    }
  }
  return result.matches;
}

module.exports = { readFixture, compareBuffers, ulpDistance, reportComparison };

// ---------------------------------------------------------------------------
// CLI: inspect a fixture file directly (sanity check, not a real comparison)
// ---------------------------------------------------------------------------
if (require.main === module) {
  const [, , path] = process.argv;
  if (!path) {
    console.error('usage: node compare.js <fixture.bin>');
    process.exit(1);
  }
  const data = readFixture(path);
  console.log(`${path}: ${data.length} doubles`);
  console.log('first 8:', Array.from(data.slice(0, 8)));
}
