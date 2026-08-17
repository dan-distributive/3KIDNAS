'use strict';

// =============================================================================
// index.js -- registry of ported complex-side codelets, keyed by radix.
// Composite1D.js consumes this to dispatch base-case (noTwiddle) and
// Cooley-Tukey combine (twiddle) candidates for a given radix, mirroring
// chooseDecomposition.js's own view of which codelets exist.
//
// Radix 6 lives in this directory (n1_6.js/t1_6.js) like every other radix
// -- it used to be hand-inlined in RaderSolver.js, moved here so
// RaderSolver.js can depend on Composite1D.js (for its own generalized N-1
// sub-transform) without a circular require back through this registry.
// =============================================================================

const { n1_2 } = require('./n1_2');
const { n1_3 } = require('./n1_3');
const { n1_4 } = require('./n1_4');
const { n1_5 } = require('./n1_5');
const { n1_6 } = require('./n1_6');
const { n1_7 } = require('./n1_7');
const { n1_8 } = require('./n1_8');
const { n1_9 } = require('./n1_9');
const { n1_10 } = require('./n1_10');
const { n1_11 } = require('./n1_11');
const { n1_12 } = require('./n1_12');
const { n1_13 } = require('./n1_13');
const { n1_14 } = require('./n1_14');
const { n1_15 } = require('./n1_15');
const { n1_16 } = require('./n1_16');
const { n1_20 } = require('./n1_20');
const { n1_32 } = require('./n1_32');
const { t1_2 } = require('./t1_2');
const { radix3Twiddle: t1_3 } = require('../../CompositeSolver1D');
const { t1_4 } = require('./t1_4');
const { t1_5 } = require('./t1_5');
const { t1_6 } = require('./t1_6');
const { t1_7 } = require('./t1_7');
const { t1_8 } = require('./t1_8');
const { t1_9 } = require('./t1_9');
const { t1_12 } = require('./t1_12');
const { t2_5 } = require('./t2_5');
const { t2_8 } = require('./t2_8');
const { t2_10 } = require('./t2_10');
const { t2_16 } = require('./t2_16');
const { t2_20 } = require('./t2_20');
const { t2_25 } = require('./t2_25');
const { t2_4 } = require('./t2_4');
const { t2_32 } = require('./t2_32');
const { t2_64 } = require('./t2_64');

const noTwiddle = {
  2: n1_2,
  3: n1_3,
  4: n1_4,
  5: n1_5,
  6: n1_6,
  7: n1_7,
  8: n1_8,
  9: n1_9,
  10: n1_10,
  11: n1_11,
  12: n1_12,
  13: n1_13,
  14: n1_14,
  15: n1_15,
  16: n1_16,
  20: n1_20,
  32: n1_32,
};

// NOTE: no twiddle[10] yet -- t1_10 deferred (large, intricate codelet,
// ~60+ intermediate variables; higher transcription risk than the rest of
// this batch warranted rushing, and t2_10 already covers every currently-
// known plan needing radix 10). isFullyPortedComplex() correctly reports
// false for any plan needing it, same as any other unported codelet.
const twiddle = {
  2: t1_2,
  3: t1_3,
  4: t1_4,
  5: t1_5,
  6: t1_6,
  7: t1_7,
  8: t1_8,
  9: t1_9,
  12: t1_12,
};

// t2 family -- alternate-codegen ("twiddle-log3/precompute-twiddles")
// siblings of t1_r: same math, different rounding, keyed separately so
// Composite1D.js can pick the exact codelet real FFTW's plan.codeletGroup
// says it would use, not just any radix-compatible one. Sparse on purpose
// -- only add an entry once that specific t2_r has been ported.
const twiddleT2 = {
  5: t2_5,
  8: t2_8,
  10: t2_10,
  16: t2_16,
  20: t2_20,
  4: t2_4,
  25: t2_25,
  32: t2_32,
  64: t2_64,
};

module.exports = { noTwiddle, twiddle, twiddleT2 };
