'use strict';

// =============================================================================
// RealEngine1D.js
// General recursive 1D real DFT engine (both R2HC forward and HC2R
// backward directions), the
// real-side counterpart to Composite1D.js. Executes a radix-r Cooley-Tukey
// decimation-in-time decomposition using the ported r2cf_*/hf_* codelets,
// operating throughout in FFTW3's own internal packed "halfcomplex" layout
// (see Codelets/real/r2cf_2.js's header) rather than converting to a split
// re/im array at every recursion level -- that packing is what hf_* natively
// combines (rdft/hc2hc.c's apply_dit: sub-transform writes into O, then the
// twiddle stage runs in place on O), so working in it throughout avoids a
// transcription-risky reformulation of the source arithmetic.
//
// SCOPE, deliberately narrow: only n where EVERY factor in the recursion
// stays odd. Traced through rdft/hc2hc-direct.c's mkcldw: for odd m (the
// sub-transform size at any recursion level), the "middle" cldm branch is a
// structural no-op (CLDMP(m,mstart,mcount) is false), so only the ordinary
// R2HC-kind cld0 (handling the r phase-DC values, itself just another
// r-point R2HC -- reuses r2cf_r directly) is ever needed. For EVEN m, cldm
// requires an R2HCII (shifted-kind) sub-transform, which is out of this
// port's stated scope (see FFTW3JS/README.md) -- so this file simply
// doesn't support that path yet, and callers must check isFullyPortedReal()
// before trusting a result, exactly like Composite1D.js's
// isFullyPortedComplex(). Since odd n forces every recursive factor to stay
// odd automatically (an even factor could never multiply out to an odd n),
// this scope restriction costs nothing for the sizes this pipeline's row
// pass has ever needed (PlanTable.js: odd N1 only).
//
// OUTPUT INDEX MAPPING: physical position IS the final packed-halfcomplex
// output index, with NO final gather/remap needed -- the same "physical
// position = output index" invariant already established for the complex
// side's dft-ct-dit (RaderSolver.js's dftSize36 header). Verified against a
// from-scratch hand-computed DFT for n=9=3x3: this file's stages 1-3
// (before any gather step existed) matched the hand-computed answer to
// machine epsilon.
//
// FALLBACK CONVENTION MISMATCH (found while verifying the above): the first
// attempt to cross-check against GenericSolver1D.rdftGenericR2hc appeared to
// show that function returning wrong data for composite odd N -- but that
// was this file's own bug, not rdftGenericR2hc's. That function's output is
// documented (see its header in GenericSolver1D.js) as a single INTERLEAVED
// complex array (O[2k]=Re,O[2k+1]=Im), fused with the hc2c bridge step --
// NOT FFTW's packed-halfcomplex format this file works in internally. Once
// called with the right (offset, stride=2) convention and its output
// converted to packed-HC, it's correct (re-verified to 1e-15 for n=7..13).
// genericFallback() below does that conversion; do not call rdftGenericR2hc
// directly elsewhere in this file with stride 1.
// =============================================================================

const { chooseReal } = require('./Planner/chooseDecomposition');
const { rdftGenericR2hc, rdftGenericHc2r } = require('./GenericSolver1D');
const { realCexp } = require('./Trig');
const {
  applyR2HCViaDht, applyHC2RViaDht, isFullyPortedDhtRader,
} = require('./DhtRaderSolver');

const r2cf = {
  2: require('./Codelets/real/r2cf_2').r2cf_2,
  3: require('./Codelets/real/r2cf_3').r2cf_3,
  4: require('./Codelets/real/r2cf_4').r2cf_4,
  5: require('./Codelets/real/r2cf_5').r2cf_5,
  6: require('./Codelets/real/r2cf_6').r2cf_6,
  7: require('./Codelets/real/r2cf_7').r2cf_7,
  8: require('./Codelets/real/r2cf_8').r2cf_8,
  9: require('./Codelets/real/r2cf_9').r2cf_9,
  10: require('./Codelets/real/r2cf_10').r2cf_10,
  11: require('./Codelets/real/r2cf_11').r2cf_11,
  12: require('./Codelets/real/r2cf_12').r2cf_12,
  13: require('./Codelets/real/r2cf_13').r2cf_13,
  14: require('./Codelets/real/r2cf_14').r2cf_14,
  15: require('./Codelets/real/r2cf_15').r2cf_15,
  16: require('./Codelets/real/r2cf_16').r2cf_16,
  20: require('./Codelets/real/r2cf_20').r2cf_20,
  25: require('./Codelets/real/r2cf_25').r2cf_25,
  32: require('./Codelets/real/r2cf_32').r2cf_32,
};
const hf = {
  2: require('./Codelets/real/hf_2').hf_2,
  3: require('./Codelets/real/hf_3').hf_3,
  4: require('./Codelets/real/hf_4').hf_4,
  5: require('./Codelets/real/hf_5').hf_5,
  6: require('./Codelets/real/hf_6').hf_6,
  7: require('./Codelets/real/hf_7').hf_7,
  9: require('./Codelets/real/hf_9').hf_9,
  10: require('./Codelets/real/hf_10').hf_10,
  12: require('./Codelets/real/hf_12').hf_12,
  15: require('./Codelets/real/hf_15').hf_15,
};
// hf2 -- alternate-codegen ("twiddle-log3/precompute-twiddles") sibling of
// hf_r, keyed separately so executePlanR2HC can dispatch the exact codelet
// real FFTW's plan.codeletGroup says it would use (same pattern as
// Composite1D.js's twiddleT2 on the complex side). chooseReal recurses
// through hf2_5 for essentially any composite involving a factor of 5 --
// confirmed empirically (176 of the N=25..2000 range use it somewhere) --
// so this, not hf_15/hf_25, is what actually unlocks multiples of 5.
const hf2 = {
  5: require('./Codelets/real/hf2_5').hf2_5,
  8: require('./Codelets/real/hf2_8').hf2_8,
  16: require('./Codelets/real/hf2_16').hf2_16,
  20: require('./Codelets/real/hf2_20').hf2_20,
  25: require('./Codelets/real/hf2_25').hf2_25,
  32: require('./Codelets/real/hf2_32').hf2_32,
};
const r2cb = {
  2: require('./Codelets/real/r2cb_2').r2cb_2,
  3: require('./Codelets/real/r2cb_3').r2cb_3,
  4: require('./Codelets/real/r2cb_4').r2cb_4,
  5: require('./Codelets/real/r2cb_5').r2cb_5,
  6: require('./Codelets/real/r2cb_6').r2cb_6,
  7: require('./Codelets/real/r2cb_7').r2cb_7,
  8: require('./Codelets/real/r2cb_8').r2cb_8,
  9: require('./Codelets/real/r2cb_9').r2cb_9,
  10: require('./Codelets/real/r2cb_10').r2cb_10,
  11: require('./Codelets/real/r2cb_11').r2cb_11,
  12: require('./Codelets/real/r2cb_12').r2cb_12,
  13: require('./Codelets/real/r2cb_13').r2cb_13,
  14: require('./Codelets/real/r2cb_14').r2cb_14,
  15: require('./Codelets/real/r2cb_15').r2cb_15,
  16: require('./Codelets/real/r2cb_16').r2cb_16,
  20: require('./Codelets/real/r2cb_20').r2cb_20,
  25: require('./Codelets/real/r2cb_25').r2cb_25,
  32: require('./Codelets/real/r2cb_32').r2cb_32,
};
const hb = {
  2: require('./Codelets/real/hb_2').hb_2,
  3: require('./Codelets/real/hb_3').hb_3,
  4: require('./Codelets/real/hb_4').hb_4,
  5: require('./Codelets/real/hb_5').hb_5,
  6: require('./Codelets/real/hb_6').hb_6,
  7: require('./Codelets/real/hb_7').hb_7,
  9: require('./Codelets/real/hb_9').hb_9,
  10: require('./Codelets/real/hb_10').hb_10,
  12: require('./Codelets/real/hb_12').hb_12,
};
// hb2 -- alternate-codegen sibling of hb_r, same reasoning as hf2 above.
// hb2_25 (Codelets/real/hb2_25.js) is intentionally NOT registered here
// (unlike hf2_25, which IS registered but separately excluded via
// isFullyPortedR2HC's guard) -- deliberately kept UNREACHABLE rather than
// "registered but excluded", because the evidence here is stronger than
// hf2_25's: this isn't a tiny few-ULP discrepancy, it's a confirmed,
// large, structural bug.
//
// Verification story: an automated semantic-evaluation diff against
// hb2_25.c (every JS assignment's dependency-graph value, evaluated with
// identical random inputs on both sides, in declaration order) found ZERO
// mismatches across 519 checked variables -- the JS is a faithful
// line-for-line transcription of the C source's expression tree. A
// SEPARATE, independent test (composing hb2_25 with hf2_25 -- which IS
// proven correct standalone, 0 discrepancies across all 625 output
// positions of a real N=625 R2HC ground-truth run -- and checking that
// hb2_25(hf2_25(x)) == 25*x, the expected exact property of a twiddle-
// combine/uncombine pair for FFTW's unnormalized transform convention,
// already independently confirmed to hold EXACTLY for the hf2_5/hb2_5
// pair) instead shows a clear, reproducible, and LARGE failure pattern:
// exactly the phases p where p%5 is 0 or 3 (excluding p=0 itself) come
// back wrong (by 20-300%, not a rounding-scale error); all other phases
// round-trip to machine precision. Confirmed independently via the
// r2hc1d_harness ground truth too: real N=625 HC2R output differs from
// this engine's by ~1300 (not ~1e-14) when hb2_25 is used.
//
// This p%5-in-{0,3} pattern is a real, actionable clue for a future pass
// (points at something in how the m=25 sub-recursion's own 5-groups
// interact with hb2_25's phase indexing, not a stray sign/macro typo --
// which the semantic diff already rules out) but wasn't root-caused this
// session. Given the semantic diff proves the JS matches the C source
// bit-for-bit in STRUCTURE, if this is confirmed not to be a calling-
// convention mismatch in a future pass, the C source's own cr[]/ci[]
// aliasing (cr = cr+ms, ci = ci-ms each loop iteration -- see hb2_25.c's
// header) or the "-sign 1" generation flag (present for hb2_25.c/hb2_5.c,
// absent for hf2_25.c/hf2_5.c) would be the next things to dig into.
const hb2 = {
  5: require('./Codelets/real/hb2_5').hb2_5,
  8: require('./Codelets/real/hb2_8').hb2_8,
  16: require('./Codelets/real/hb2_16').hb2_16,
  20: require('./Codelets/real/hb2_20').hb2_20,
  32: require('./Codelets/real/hb2_32').hb2_32,
};

// r2cfII / r2cbIII -- R2HCII/HC2RIII "cldm" direct codelets (rdft/hc2hc.c's
// middle-column combine, needed whenever a CT sub-transform size m is
// EVEN -- see the "Stage cldm" comment in executePlanR2HC/executePlanHC2R
// below for the full math/derivation). Keyed by the OUTER radix r (same
// radix as the hf_r/hb_r twiddle codelet it accompanies), NOT by m.
// Narrow radix set so far (matches Planner/chooseDecomposition.js's
// CLDM_OPS, which only offers `needsCldm` candidates for radixes this
// registry actually covers).
const r2cfII = {
  3: require('./Codelets/real/r2cfII_3').r2cfII_3,
  5: require('./Codelets/real/r2cfII_5').r2cfII_5,
  6: require('./Codelets/real/r2cfII_6').r2cfII_6,
  7: require('./Codelets/real/r2cfII_7').r2cfII_7,
  8: require('./Codelets/real/r2cfII_8').r2cfII_8,
  9: require('./Codelets/real/r2cfII_9').r2cfII_9,
  10: require('./Codelets/real/r2cfII_10').r2cfII_10,
  12: require('./Codelets/real/r2cfII_12').r2cfII_12,
  16: require('./Codelets/real/r2cfII_16').r2cfII_16,
  20: require('./Codelets/real/r2cfII_20').r2cfII_20,
  32: require('./Codelets/real/r2cfII_32').r2cfII_32,
};
const r2cbIII = {
  3: require('./Codelets/real/r2cbIII_3').r2cbIII_3,
  5: require('./Codelets/real/r2cbIII_5').r2cbIII_5,
  6: require('./Codelets/real/r2cbIII_6').r2cbIII_6,
  7: require('./Codelets/real/r2cbIII_7').r2cbIII_7,
  8: require('./Codelets/real/r2cbIII_8').r2cbIII_8,
  9: require('./Codelets/real/r2cbIII_9').r2cbIII_9,
  10: require('./Codelets/real/r2cbIII_10').r2cbIII_10,
  12: require('./Codelets/real/r2cbIII_12').r2cbIII_12,
  16: require('./Codelets/real/r2cbIII_16').r2cbIII_16,
  20: require('./Codelets/real/r2cbIII_20').r2cbIII_20,
  32: require('./Codelets/real/r2cbIII_32').r2cbIII_32,
};

// genericCombineR2HC / genericCombineHC2R -- faithful port of FFTW3's THIRD
// real-side CT combine solver, rdft/hc2hc-generic.c ("hc2hc-generic"),
// used whenever the outer CT radix r has a direct r2cf_r/r2cb_r codelet
// but NO dedicated twiddle codelet (hf_r/hf2_r/hb_r/hb2_r) -- the only two
// such radixes in real FFTW's own registered codelet set are 11 and 13
// (see Planner/chooseDecomposition.js's GENERIC_COMBINE_RADIXES comment).
// Unlike every OTHER codelet in this port, this is NOT a mechanical
// per-line transcription of hc2hc-generic.c's C source (that source
// operates via raw in-place pointer arithmetic over the FULL n=r*m array
// -- reading it as a literal translation gives a plausible-looking but
// WRONG algorithm, confirmed by trial: a naive "twiddle then two separate
// real transforms of the Re-parts/Im-parts across phases, recombined via
// even/odd DFT symmetry" formulation is mathematically VALID (verified
// against a naive O(r^2) complex-DFT reference) but does NOT match real
// FFTW's bit pattern -- large (not rounding-scale) errors resulted).
//
// What IS below was derived by manually tracing hc2hc-generic.c's THREE
// helper functions against this file's own (col, phase) array convention:
//   - bytwiddle(): confirmed to be an ordinary complex twiddle multiply of
//     (cr[p],ci[p]) by (Wc[p],Ws[p]) for p=1..r-1, p=0 untouched -- same
//     twiddle table this file already builds for every other codelet.
//   - the "cld" sub-plan: NOT a single r-point COMPLEX DFT (despite that
//     being the mathematically natural reading) -- it's the "col" and
//     "(m-col)" physical positions (across all r phases) each getting
//     their OWN, INDEPENDENT r-point REAL R2HC/HC2R transform (reusing
//     r2cf_r/r2cb_r) -- i.e. exactly the same "pack two reals into one
//     complex FFT" trick used throughout real-FFT implementations, but
//     the C source hides this behind generic rank-1 RDFT sub-problem
//     machinery rather than a visible "ReSpec/ImSpec" split.
//   - reorder_dit()/reorder_dif() + swapri(): the recombination step
//     that turns those two independent transforms' packed-halfcomplex
//     outputs into this file's final per-phase (outCr,outCi) pairs.
//     swapri() reverses the phase index specifically on the "(m-col)"
//     side (phase p <-> phase r-1-p, middle phase untouched) -- easy to
//     miss since it operates on a DIFFERENT physical position than the
//     main reorder loop it's paired with.
// Verified bit-exact against real FFTW ground truth (r2hc1d_harness) for
// BOTH r=11 (n=121) and r=13 (n=169), both R2HC and HC2R directions, not
// just one -- confirming the derivation generalizes across radix, not a
// coincidence specific to one size.
function genericCombineR2HC(r, cr, ci, Wc, Ws) {
  const rePart = new Float64Array(r), imPart = new Float64Array(r);
  rePart[0] = cr[0]; imPart[0] = ci[0];
  for (let p = 1; p < r; p++) {
    rePart[p] = cr[p] * Wc[p] + ci[p] * Ws[p];
    imPart[p] = ci[p] * Wc[p] - cr[p] * Ws[p];
  }
  // r's own r-point R2HC sub-transform ("cld0"/"cld" in hc2hc-generic.c's
  // own mkcldw -- a GENERIC rdft_1d problem, handed to the full planner,
  // NOT hardcoded to a direct codelet): matches real FFTW exactly even
  // when r itself has no direct r2cf_r (e.g. r=17, confirmed via
  // fftw_fprint_plan for n=323=17*19 -- "hc2hc-generic-dit-17-19" whose
  // OWN cld0/cld are "rdft-generic-r2hc-17"/"-19", not direct codelets).
  const ReSpec = executePlanR2HC(chooseReal('R2HC', r).plan, r, rePart);
  const ImSpec = executePlanR2HC(chooseReal('R2HC', r).plan, r, imPart);
  const outCr = new Float64Array(r), outCi = new Float64Array(r);
  outCr[0] = ReSpec[0];
  for (let a = 1; a < r; a++) outCr[a] = ReSpec[a] - ImSpec[r - a];
  for (let a = 1; a < r - 1; a++) outCi[a] = ImSpec[r - 1 - a] + ReSpec[a + 1];
  outCi[r - 1] = ImSpec[0];
  outCi[0] = ImSpec[r - 1] + ReSpec[1];
  return [outCr, outCi];
}

function genericCombineHC2R(r, cr, ci, Wc, Ws) {
  const ReSpecP = new Float64Array(r), ImSpecP = new Float64Array(r);
  ReSpecP[0] = cr[0];
  ImSpecP[0] = ci[r - 1];
  for (let a = 1; a < r; a++) {
    ReSpecP[a] = 0.5 * (cr[a] + ci[a - 1]);
    ImSpecP[a] = 0.5 * (ci[r - 1 - a] - cr[r - a]);
  }
  const localRe = executePlanHC2R(chooseReal('HC2R', r).plan, r, ReSpecP);
  const localIm = executePlanHC2R(chooseReal('HC2R', r).plan, r, ImSpecP);
  const outCr = new Float64Array(r), outCi = new Float64Array(r);
  outCr[0] = localRe[0]; outCi[0] = localIm[0];
  for (let p = 1; p < r; p++) {
    outCr[p] = localRe[p] * Wc[p] - localIm[p] * Ws[p];
    outCi[p] = localIm[p] * Wc[p] + localRe[p] * Ws[p];
  }
  return [outCr, outCi];
}

// correctR2HCNode / correctHC2RNode -- rewrite a SINGLE plan node where
// chooseReal's naive recursive cost model disagrees with real FFTW's
// actual codelet choice (root-caused, not a codelet bug -- see the long
// comment previously attached to isFullyPortedR2HC's exclusion, now
// superseded by this correction). Applied at the top of every recursive
// walk/execute call, so every depth self-corrects without a full-tree
// pre-pass and without ever mutating chooseReal's memoized cache.
//
// R2HC, radix===3 && m in {7,9}: chooseReal predicts radix-3 outer /
// radix-m inner; real FFTW actually picks radix-m outer / radix-3 inner,
// at EVERY nesting depth. Confirmed via fftw_fprint_plan for n=21, 27,
// 105, 147, 189, 243, 1029, 1701, 2187 (both direct top-level occurrences
// and nested-under-rdft-vrank>=1 occurrences). hf_7/hf_9/r2cf_3 are all
// already-ported, verbatim-correct codelets -- this is a pure structural
// swap, not new arithmetic.
function isPrime(n) {
  if (n <= 1) return false;
  if (n % 2 === 0) return n === 2;
  for (let i = 3; i * i <= n; i += 2) if (n % i === 0) return false;
  return true;
}

// GENERIC_SAFE_MAX -- largest PRIME n for which chooseReal's 'generic' leaf
// (rdftGenericR2hc/rdftGenericHc2r, faithful ports of rdft/generic.c's
// hartley+cdot algorithm) is PROVEN to match real FFTW's actual chosen
// solver bit-exactly. Gated to prime n ONLY -- chooseReal can ALSO predict
// a 'generic' leaf for a COMPOSITE n whose only prime factor has no
// registered twiddle (hf/hf2) codelet (e.g. n=121=11^2, n=143=11*13: no
// hf_11 exists in real FFTW at all, so our chooseReal has no CT option and
// falls to whole-n 'generic'). Confirmed via fftw_fprint_plan that real
// FFTW does NOT do the same thing there -- it uses "hc2hc-generic-dit-11-*"
// (n=121: "(rdft-ct-dit/11 (hc2hc-generic-dit-11-11 (r2cf_11)(r2cf_11))
// (r2cf_11))"), a THIRD, unported real-side solver family: a genuine CT
// split whose TWIDDLE-COMBINE stage itself falls back to a generic
// (non-codelet) algorithm, structurally different from (and NOT bit-exact
// with) treating the whole composite n as one opaque generic leaf -- this
// produced small (~1e-13/1e-14) mismatches for every n reachable through
// 121 or 143 before this guard was added (n=121, 143, 363, 429, 605, 715,
// 847, 1001, 1089, 1287, 1815, 2145, 2541, 3003, 3267, 3861, 4235, found by
// a full sweep). hc2hc-generic is UNPORTED -- parked alongside dht-rader
// below, not a small fix (it's a distinct combine algorithm, not a codelet
// substitution).
//
// Above this threshold (still restricted to prime n), real FFTW switches
// to an entirely different,
// UNPORTED solver family: "r2hc-dht" / "dht-rader" -- a Discrete Hartley
// Transform formulation of Rader's algorithm for the real domain,
// structurally unrelated to both rdft/generic.c AND to the existing
// complex-domain RaderSolver.js (confirmed via fftw_fprint_plan: n=167 is
// still plain "(rdft-generic-r2hc-167)"/"(rdft-generic-hc2r-167)", while
// n=173 -- the very next odd prime, nothing in between -- is
// "(r2hc-dht-173 (dht-rader-173/360 ...))" / "(hc2r-dht-173 ...)"). This
// is a genuinely new algorithm family (its own codelet types, its own
// "rdft-buffered"/"rdft-rank0-iter-ci" wrapper solvers visible in the
// dumped plans), not a small guard fix -- parked, same "complex stuff,
// circle back later" bucket as q1 and the 2D-embedding characterization.
//
// Below/at 167: exhaustively verified bit-exact (maxDiff=0) against real
// fftw_fprint_plan-confirmed ground truth for EVERY odd prime from 17
// (GENERIC_MAX_SLOW's own lower gate, see Planner/chooseDecomposition.js)
// through 167, both R2HC and HC2R directions, both bare (e.g. n=157) and
// nested inside a larger composite (e.g. n=489=3*163, n=34=2*17 -- the
// per-node solver choice is a pure function of that node's own size,
// confirmed identical whether the prime appears standalone or embedded).
const GENERIC_SAFE_MAX = 167;

// BAD_DHT_RADER_N -- prime n where Planner/chooseDecomposition.js's
// dhtRaderCost formula wrongly predicts the UNPADDED dht-rader variant
// (npad=n-1) as cheaper than the PADDED one (npad=smallest even 5-smooth
// >= 2*(n-1)-1), when real FFTW's actual plan uses padding (confirmed via
// direct fftw_fprint_plan cross-check for every entry here, e.g. n=269:
// real plan is "dht-rader-269/540", not "-269/268"). DhtRaderSolver.js
// only implements the UNPADDED execution path, so these n produce a
// mathematically-valid-but-not-bit-identical result if trusted -- same
// "cost model doesn't perfectly replicate FFTW's real crossover, fix by
// evidence not by re-deriving the exact formula" pattern used throughout
// this file's other BAD_*_M sets, just keyed by n instead of m since
// dht-rader is a leaf, not a CT (radix,m) node. The root cause (why the
// literal ops-formula transcription is imprecise here) traces to
// dht-rader.c's cld1/cld2 being planned with NO_SLOW forcibly set via
// mkplan_f_d(plnr,cldp,NO_SLOW,0,0) -- not yet fully isolated, see
// dhtRaderCost's header. Collected via a direct sweep of every prime
// 168..8000 whose OWN unpadded cost we predict as winning, cross-checked
// against fftw_fprint_plan one by one -- 87 of 312 candidates were wrong;
// extend only from further direct evidence, not by pattern-guessing.
const BAD_DHT_RADER_N = new Set([
  269, 647, 677, 773, 1087, 1283, 1447, 1453, 1543, 1811, 1931, 2003, 2111,
  2179, 2309, 2411, 2423, 2663, 2693, 2887, 2897, 2903, 3089, 3203, 3209,
  3253, 3259, 3373, 3389, 3433, 3463, 3547, 3593, 3607, 3847, 3863, 3877,
  3881, 4021, 4049, 4217, 4271, 4289, 4337, 4339, 4357, 4373, 4813, 5197,
  5237, 5279, 5297, 5303, 5393, 5431, 5779, 5783, 5791, 5821, 6007, 6011,
  6287, 6449, 6563, 6571, 6637, 6659, 6689, 6701, 6709, 6733, 6737, 6761,
  6823, 6833, 6949, 7151, 7207, 7211, 7213, 7219, 7237, 7243, 7283, 7481,
  7691, 7723,
  // Extended to N=15000 via the same direct-sweep-against-r2hc1d_harness
  // methodology (safety_sweep.js) -- these showed up as genuine
  // isFullyPortedR2HC/HC2R=true-but-wrong mismatches, each one confirmed
  // (by construction of the sweep) to use the unpadded "dht-rader-n/n-1"
  // shape when real FFTW's compiled output disagrees.
  8009, 8011, 8053, 8089, 8093, 8669, 8677, 8689, 8713, 8731, 8741, 8863,
  9283, 9649, 10037, 10099, 10337, 10613, 10639, 10861, 10891, 11617,
  11621, 11831, 12007, 12109, 12113, 12671, 12821, 12907, 12967, 12973,
  12979, 13003, 13009, 13033, 13339, 13829, 14437, 14461, 14731, 14771,
]);

// BAD_HF2_8_M -- m values where chooseReal('R2HC'|'HC2R', 8*m) predicts
// radix=8 outer (hf2_8) but real FFTW actually picks radix=16 or radix=32
// outer instead (unported). Collected from a direct 4000-N sweep against
// r2hc1d_harness -- see isFullyPortedR2HC's hf2/radix===8 guard for the
// full explanation. Not a derived formula; extend from evidence only.
const BAD_HF2_8_M = new Set([
  16, 22, 26, 28, 32, 44, 52, 56, 64, 68, 72, 76, 84, 88, 92, 96, 104, 112,
  116, 124, 128, 136, 152, 168, 176, 184, 208, 216, 224, 232, 248, 256, 264,
  296, 312, 328, 344, 352, 376, 392, 408, 416, 424, 448, 456, 472, 488,
]);

// BAD_HF2_16_M -- same idea as BAD_HF2_8_M, one radix tier up: m values
// where chooseReal predicts radix=16 outer (hf2_16) but real FFTW prefers
// radix=32 outer (hf2_32, unported) instead. Collected the same way (direct
// 4000-N sweep vs r2hc1d_harness).
const BAD_HF2_16_M = new Set([
  12, 22, 26, 28, 34, 38, 46, 58, 62, 64, 68, 72, 76, 84, 88, 92, 96, 104,
  112, 116, 124, 128, 136, 148, 152, 164, 168, 172, 184, 188, 192, 196, 212,
  216, 232, 236, 244, 248,
]);

// BAD_HF2_32_M -- same idea one radix tier up again: m values where
// chooseReal predicts radix=32 outer (hf2_32) but real FFTW prefers a
// DIFFERENT radix instead -- unlike the 8/16 tiers, this isn't always a
// bigger radix (n=256 loses to radix=16, not radix=64), confirming this is
// a genuine per-N cost tie-break gap, not a simple "always prefer larger
// power of 2" rule. Collected via the same direct 4000-N sweep. Radix=32
// itself (hf2_32/hb2_32/r2cfII_32/r2cbIII_32) is independently verified
// correct -- see n=1024 and n=1056, both bit-exact both directions.
// m=24 added after wiring hf_12/r2cfII_12 (this session): n=768's own
// m=24 sub-decomposition (radix=12,m=2+cldm via hf_12/r2cfII_12) is
// independently verified correct (n=24 alone is bit-exact and matches
// fftw_fprint_plan exactly) -- the bug is purely that chooseReal(R2HC,768)
// picks radix=32 outer, while real FFTW's actual n=768 plan uses radix=16
// outer (hf2_16) with m=48 wrapped in an "rdft-vrank>=1" solver this port
// doesn't model (see the parked 2D-embedding/vrank items) -- same "wrong
// outer radix" cost-model gap as every other entry in this set, not a
// codelet bug.
const BAD_HF2_32_M = new Set([8, 12, 22, 24, 26, 28, 64, 128]);

// BAD_HF_9_M / BAD_HF2_5_M -- newly-reachable false positives found
// immediately after wiring r2cfII_5/r2cfII_9 (this session): before that,
// EVERY radix=5-cldm and radix=9-cldm candidate was unconditionally
// excluded (no r2cfII[5]/r2cfII[9]), so these are new classes of
// false positive, not regressions in previously-shipped behavior.
// BAD_HF_9_M: n=18 (radix=9,m=2+cldm) -- fftw_fprint_plan confirms real
// FFTW's actual n=18 plan uses radix=3 (hf_3) with m=6 DIRECT via r2cf_6,
// not radix 9 at all.
const BAD_HF_9_M = new Set([2]);
// BAD_HF2_5_M: m=24,28,44,52,60,68,76 -- each confirmed via
// fftw_fprint_plan to have real FFTW picking a DIFFERENT outer radix
// (e.g. n=120 -> radix=20 outer with m=6 direct via r2cf_6, not radix=5
// with m=24). Collected via a direct 4000-N sweep against r2hc1d_harness;
// narrower/nested cases (e.g. n=90, n=162, n=450, n=810, n=1458 and
// others) resolve automatically once these two root-cause blacklists are
// applied, since chooseReal is memoized and those n's sub-decompositions
// recompute through the now-corrected chooseReal(R2HC,18)/(R2HC,m) calls
// -- confirmed by re-sweeping after adding this set, not assumed.
const BAD_HF2_5_M = new Set([24, 28, 44, 52, 60, 68, 76]);

// BAD_HF2_20_M -- same idea, radix 20. Collected via a direct 4000-N sweep
// against r2hc1d_harness after wiring hf2_20/hb2_20/r2cfII_20/r2cbIII_20:
// only m=64 and m=128 came back wrong (fftw_fprint_plan confirms real
// FFTW's n=1280/n=2560 plans use radix=20 outer -- matching our own
// top-level prediction exactly -- but with m's OWN sub-transform forced to
// a DIRECT r2cf_64/r2cf_128 codelet rather than the recursive hf2_8/hf2_16
// decomposition chooseReal(R2HC,64)/chooseReal(R2HC,128) picks in
// isolation; confirmed via fftw_fprint_plan for bare n=64/n=128 that the
// recursive decomposition IS what real FFTW picks standalone -- so this is
// context-dependent planner behavior specific to being a hc2hc "cld"
// sub-transform, not a wrong standalone prediction). Exactly the same
// m=64/m=128 values already independently found and blacklisted for
// radix 8 and radix 16 above (BAD_HF2_8_M, BAD_HF2_16_M both contain
// 64 and 128) -- this is evidently a general phenomenon at these two m
// values, not radix-specific, though kept as a separate named set per
// this file's established one-set-per-radix convention.
const BAD_HF2_20_M = new Set([64, 128]);

// BAD_HF_3_M -- m values where chooseReal('R2HC', 3*m) predicts radix=3
// outer (hf_3) WITH cldm engaged, but real FFTW actually picks a
// DIFFERENT, unported outer radix instead (fftw_fprint_plan confirms:
// n=24 -> real FFTW uses radix=12 (hf_12, the PLAIN non-alt-codegen family,
// still not ported) with m=2+cldm via r2cfII_12 (also unported); n=36 ->
// radix=6). Only surfaced once r2cfII_3/r2cbIII_3 were wired (this
// session) -- before that, EVERY radix=3+cldm candidate was
// unconditionally excluded (no r2cfII[3]), so this is a newly-reachable
// class of false positive, not a regression in previously-shipped
// behavior. Collected via a direct 4000-N sweep against r2hc1d_harness
// immediately after wiring r2cfII_3/r2cbIII_3 -- narrow, enumerated set
// (m=8, 12, 44), same evidence-only convention as every other BAD_*_M set
// in this file. NOTE: once hf_6/r2cfII_6 were ALSO wired (this session,
// immediately after), radix=6 started winning chooseReal's cost
// comparison over radix=3 for m=8 and m=44 (both divisible by 4, and
// hf_6+cldm's cost came out lower) -- so this set's m=8/44 entries are
// believed DEAD (radix=3 no longer reachable there) but kept for safety
// since "no longer reachable" isn't the same guarantee as "verified
// wrong"; m=12 may still be live (12 not divisible by... n=36=3*12 was
// NOT re-checked after the radix-6 wiring, so this is deliberately not
// pruned without fresh evidence).
const BAD_HF_3_M = new Set([8, 12, 44]);

// BAD_HF_6_M -- same idea, radix 6 (hf_6). Collected via a direct 4000-N
// sweep immediately after wiring hf_6/r2cfII_6: real FFTW's n=24/48 plans
// use radix=12 (hf_12, unported) outer instead, and n=264(=6*44)'s uses
// radix=... (same category, unported alternative) -- see BAD_HF_3_M's
// header, this is the SAME underlying n=24/48/264 gap now reached via a
// different (cheaper, once available) our-side candidate.
const BAD_HF_6_M = new Set([4, 8, 22]);

// BAD_HB_*_M / BAD_HB2_*_M -- HC2R (backward) analogues of the BAD_HF2_*_M
// sets above. Backward's needsCldm claims were originally EXCLUDED
// ENTIRELY (see git history / this file's prior revisions) because an
// early sweep -- run when only radix 7 was ported -- found every reachable
// case wrong. Re-investigated after the radix 8/16/32 even-radix work
// landed: a direct sweep bypassing the blanket exclusion found the
// MAJORITY (68 of 93) of needsCldm-somewhere HC2R cases are actually
// bit-exact -- the failures are the SAME "wrong outer radix chosen by our
// cost model" pattern as the forward direction, just with a DIFFERENT
// (not necessarily overlapping) set of bad (radix,m) pairs, confirming
// forward and backward genuinely do pick different structures for some N
// (e.g. n=28: forward correctly uses radix=7, backward's true structure
// uses radix=2 -- an even radix this port doesn't model at all -- so
// radix=7 is wrong for backward specifically). Collected via the same
// direct-sweep-against-r2hc1d_harness methodology, narrowed to exactly
// what's been verified; extend only from further direct evidence.
const BAD_HB_7_M = new Set([4, 8]);
// BAD_HB_3_M -- HC2R analogue of BAD_HF_3_M, radix=3 (hb_3). Collected the
// same way, immediately after wiring r2cbIII_3 -- only m=8 and m=44
// (NOT m=12 -- see isFullyPortedHC2R's radix=3/m=12 note: real FFTW's
// n=36 HC2R plan already picks a different, unrelated outer radix (9) on
// its own, so that (radix=3,m=12) shape is simply never reached here,
// not verified-and-safe).
const BAD_HB_3_M = new Set([8, 44]);
// BAD_HB_6_M -- HC2R analogue of BAD_HF_6_M, radix=6 (hb_6). Collected the
// same way, immediately after wiring r2cbIII_6 -- m=4 and m=22 confirmed
// bad; m=8 not confirmed bad on the HC2R side within the swept range (see
// BAD_HB_3_M's m=12 note for the same "not reached, not verified" caveat).
const BAD_HB_6_M = new Set([4, 22]);
// BAD_HB_10_M -- HC2R analogue of BAD_HF_6_M/BAD_HB_6_M, radix=10 (hb_10).
// n=20 (radix=10, m=2+cldm): real FFTW's actual HC2R plan is a single
// DIRECT "r2cb_20" codelet (n<=20 small-N direct-preferred, same category
// as the radix=7/m=2 exclusion elsewhere in this file), not any CT
// decomposition at all -- confirmed via fftw_fprint_plan. R2HC(20) does
// NOT have this problem (chooseReal('R2HC',20) already correctly picks
// direct r2cf_20) -- this is HC2R-specific cost-model asymmetry.
const BAD_HB_10_M = new Set([2]);
// BAD_HB_9_M / BAD_HB2_5_M -- HC2R analogues of BAD_HF_9_M/BAD_HF2_5_M,
// found the same way immediately after wiring r2cbIII_5/r2cbIII_9. m=12
// (n=108, root: real FFTW uses radix=12 (hb_12) with m=9 direct/no-cldm)
// and m=2 (n=576, nested inside a radix=32 nested sub-decomposition --
// independently confirmed bad, not just inherited from m=12) both found
// bad for radix=9. m=24,44,52,68,76 found bad for radix=5 -- same
// evidence-only convention, collected via direct 4000-N sweep.
const BAD_HB_9_M = new Set([2, 12]);
const BAD_HB2_5_M = new Set([24, 44, 52, 68, 76]);
const BAD_HB2_8_M = new Set([8, 22, 26, 28]);
const BAD_HB2_16_M = new Set([6, 8, 22, 26, 28, 34, 38, 46, 58, 62]);
const BAD_HB2_32_M = new Set([6, 8, 12, 22, 26, 28, 56]);

function correctR2HCNode(plan) {
  if (plan.type === 'ct' && plan.radix === 3 && (plan.m === 7 || plan.m === 9)) {
    const r = plan.m;
    return {
      type: 'ct', radix: r, m: 3,
      codeletName: `hf_${r}`, codeletGroup: 'hf',
      sub: { type: 'direct', radix: 3, codeletName: 'r2cf_3' },
    };
  }
  return plan;
}

// HC2R is ASYMMETRIC with R2HC here (same class of asymmetry already found
// for the separate radix=5/m=5 HC2R-only mispredict): only m===7 needs the
// swap (confirmed via fftw_fprint_plan for n=21, 147, 189, 1029, 1701 --
// real FFTW picks radix-7 outer / radix-3 inner, same shape as R2HC).
// m===9 does NOT need correction -- chooseReal('HC2R', 27)'s naive
// prediction (radix-3 outer/hb_3, m=9 sub/r2cb_9 direct) already matches
// real FFTW's actual plan exactly, confirmed via fftw_fprint_plan for
// n=27, 243, 1215, 2187. Both hb_3 and r2cb_9 are already ported, so this
// case needed no exclusion at all -- it was a false exclusion, not a
// mispredict.
function correctHC2RNode(plan) {
  if (plan.type === 'ct' && plan.radix === 3 && plan.m === 7) {
    return {
      type: 'ct', radix: 7, m: 3,
      codeletName: 'hb_7', codeletGroup: 'hb',
      sub: { type: 'direct', radix: 3, codeletName: 'r2cb_3' },
    };
  }
  return plan;
}

function isFullyPortedR2HC(n) {
  function walk(rawPlan) {
    const plan = correctR2HCNode(rawPlan);
    if (plan.type === 'direct') return plan.codeletGroup !== 'hf2' && !!r2cf[plan.radix];
    // Only the plain 'hf' family is ported -- 'hf2' is a different codelet
    // (alternate codegen) that executePlanR2HC would silently substitute
    // hf_r for, which is mathematically correct but NOT FFTW3's exact bit
    // path. plan.codeletGroup on a 'direct' node is undefined (r2cf has no
    // alternate-codegen sibling in scope here), hence the r2cf check above.
    //
    // radix===3, m in {7,9}: corrected (not excluded) via correctR2HCNode
    // above -- see that function's header for the fftw_fprint_plan
    // evidence. Was previously a hard exclusion here before hf_7/hf_9/
    // r2cf_3 were all ported and the swap was proven correct.
    //
    // hf2 group, radix===5, m===25 ALSO excluded -- found while sweeping
    // after hf2_5/hb2_5 landed: real FFTW's n=125 R2HC plan is
    // rdft-ct-dit/25(hf2_25, cld=rdft-ct-dit/5(hf2_5, r2cf_5)) -- radix 25
    // OUTER (needs the unported hf2_25) -- while chooseReal('R2HC', 125)
    // predicts radix 5 outer with m=25, recursing into ANOTHER hf2_5 layer
    // for the m=25 sub-problem instead of real FFTW's direct-25-flavored
    // nesting. Confirmed via fftw_fprint_plan. hf2_5.js itself is correct
    // (independently verified for n=25, where chooseReal's radix-5-outer
    // prediction DOES match real FFTW) -- this is scoped to the (5,25)
    // shape specifically, not "radix 5 is unsafe."
    // hf2 group, radix===25 EXCLUDED ENTIRELY (regardless of m) -- unlike
    // the guards above (which exclude a wrong PREDICTION), hf2_25.js's own
    // OUTPUT has a tiny (~2.8e-14) discrepancy from real FFTW at specific
    // (col, phase) positions, found while verifying n=625 (25^2, the first
    // N reachable through this codelet) -- reproduced identically at
    // n=675 too (same 2.842170943040401e-14 magnitude, a different N/cld
    // shape entirely), ruling out an accidental one-off and pointing at
    // hf2_25.js itself rather than orchestration or a chooseReal
    // mispredict. Re-verified hf2_25.js line-for-line against hf2_25.c
    // MULTIPLE times across multiple sessions (most recently: the entire
    // dependency chain feeding outCi[11] at col=6 -- every FMA/FNMS macro
    // expansion from raw cr[]/ci[] reads through T3d/T3e/T3g/T3h/T3v/T3w/
    // T5C/T5D/T5F/T5G/T5E/T5H/T5I/T5K/T5A/T5B/T5J down to the outCi[11]
    // write -- with zero discrepancies found; the JS is a byte-for-byte
    // faithful expression tree of the C source).
    //
    // Also investigated and RULED OUT: FFTW's twiddle generator
    // (kernel/trig.c) has a second, alternate mode -- AWAKE_SQRTN_TABLE --
    // that builds twiddle factors via a sqrt(n)-decomposition (two small
    // cached tables W0/W1, combined per-index via a complex MULTIPLY
    // rather than a fresh trig call), which would NOT be bit-identical to
    // this file's "always call realCexp fresh" assumption (matching
    // FFTW's OTHER mode, AWAKE_SINCOS). Confirmed via api/apiplan.c
    // (`sizeof(trigreal) > sizeof(R)` selects SQRTN_TABLE, else SINCOS)
    // and kernel/ifftw.h/config.h (FFTW_LDOUBLE and FFTW_QUAD both
    // undefined in this vendored build -> trigreal==double==R) that this
    // build ALWAYS uses AWAKE_SINCOS -- a real, plausible hypothesis,
    // but ruled out, not just unchecked.
    //
    // Root cause STILL not isolated after this. hf2_25.js is KEPT in the
    // tree (the code is otherwise extensively verified and this may well
    // be fixable), but excluded from every isFullyPortedR2HC claim until
    // root-caused -- same "park it, document it, don't guess" treatment
    // as q1. Candidate next steps for a future pass: dump FFTW's actual
    // internal twiddle table values (via a modified trig.c or an
    // instrumented plan) for n=625, col=6 and diff against realCexp's
    // values bit-for-bit, rather than re-deriving the angle-reduction
    // math by hand again; or write an independent (non-hf2_25-derived)
    // reference combine for radix 25 to compare every intermediate, not
    // just the final output.
    if (plan.type === 'ct') {
      if (plan.codeletGroup === 'hf2' && plan.radix === 5 && plan.m === 25) return false;
      if (plan.codeletGroup === 'hf2' && plan.radix === 25) return false;
      // radix===7, m===16: found NOT bit-exact (~2.7e-15) via the same
      // direct sweep once r2cf_16 was ported and this shape became
      // executable for the first time (n=112). fftw_fprint_plan confirms
      // real FFTW's actual n=112 R2HC plan uses radix=8 OUTER ("hf2_8", an
      // EVEN-radix alternate-codegen twiddle codelet this port doesn't
      // model at all) with m=14 (direct r2cf_14) -- a completely different,
      // unrelated structure, not a cldm/codelet bug. Same category as the
      // (radix=7,m=2) exclusion below: our odd-radix-only search finds A
      // valid decomposition, just not the one real FFTW's fuller (odd+even
      // radix) search actually picks.
      if (plan.codeletGroup === 'hf' && plan.radix === 7 && plan.m === 16) return false;
      // radix===7, m===32: same category, found once r2cf_32 made n=224
      // executable. fftw_fprint_plan confirms real FFTW's actual n=224
      // R2HC plan uses radix=16 OUTER ("hf2_16", another even-radix
      // alternate-codegen twiddle codelet) with m=14 (direct r2cf_14).
      if (plan.codeletGroup === 'hf' && plan.radix === 7 && plan.m === 32) return false;
      // hf2, radix===8: a 4000-N sweep found 47 m values where real FFTW
      // actually prefers radix=16 or radix=32 outer instead (hf2_16/hf2_32,
      // neither ported) -- confirmed via fftw_fprint_plan for several
      // (n=128,176,224 -> hf2_16; n=2048,3584,3904 -> hf2_32). The pattern
      // isn't a clean function of m alone (m=24 and m=12 both correctly use
      // radix=8, while m=16 and m=22 don't) so this is a plain enumerated
      // set from direct evidence, not a derived rule -- extend it if a
      // future sweep at larger N finds more. chooseReal(kind,m) is memoized
      // and this comparison is intrinsic to the 8*m subproblem regardless
      // of what encloses it, so blacklisting by m alone is safe (context-
      // independent) for any n where this exact (radix=8, m) pair recurs.
      if (plan.codeletGroup === 'hf2' && plan.radix === 8 && BAD_HF2_8_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hf2' && plan.radix === 16 && BAD_HF2_16_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hf2' && plan.radix === 32 && BAD_HF2_32_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hf2' && plan.radix === 20 && BAD_HF2_20_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hf' && plan.radix === 9 && plan.needsCldm && BAD_HF_9_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hf2' && plan.radix === 5 && plan.needsCldm && BAD_HF2_5_M.has(plan.m)) return false;
      // where Stage 3's ordinary twiddle loop never runs at all -- colLimit
      // = m/2-1 = 0 -- and ONLY Stage cldm executes): found NOT bit-exact
      // (~1.1e-16, a genuine mismatch, not FP noise) via a direct sweep
      // against real compiled FFTW3 (r2hc1d_harness) for n=14, propagating
      // to n=98=7*14 and n=686=7*98 through the SAME recursive sub-plan.
      // fftw_fprint_plan confirms real FFTW's actual n=14 R2HC plan is a
      // single DIRECT "r2cf_14" codelet -- NOT a CT decomposition at all --
      // so chooseReal's cost model is choosing the wrong STRUCTURE here
      // (a composite-direct-codelet gap this session also found affects
      // n=6,9,10,18,24,36,... -- see chooseDecomposition.js's TODO), not a
      // codelet-arithmetic bug. Every OTHER radix=7 cldm shape (m=4,8,12,
      // 16,20,... through n=1000) verified bit-exact in the same sweep, so
      // this is scoped narrowly to m===2 specifically, not radix 7 broadly.
      if (plan.codeletGroup === 'hf' && plan.radix === 7 && plan.m === 2) return false;
      // hf, radix===3, needsCldm: see BAD_HF_3_M's header -- only reachable
      // once r2cfII_3 was wired this session; real FFTW prefers a
      // different, unported outer radix (12, 6, ...) for these specific m.
      if (plan.codeletGroup === 'hf' && plan.radix === 3 && plan.needsCldm && BAD_HF_3_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hf' && plan.radix === 6 && plan.needsCldm && BAD_HF_6_M.has(plan.m)) return false;
      // See executePlanR2HC's identical guard -- r can be even now (hf2),
      // so n%2===0 is no longer a valid proxy for needsCldm; check the
      // flag directly.
      if (plan.needsCldm && !r2cfII[plan.radix]) return false;
      if (plan.codeletGroup === 'hf') return !!hf[plan.radix] && walk(plan.sub);
      if (plan.codeletGroup === 'hf2') return !!hf2[plan.radix] && walk(plan.sub);
      // generic-combine's radix r does NOT need its own direct r2cf_r --
      // genericCombineR2HC uses r's own FULL R2HC engine internally (see
      // that function's header), so trust it exactly whenever r's own
      // R2HC is itself already trusted (direct, CT, OR the prime<=167
      // generic fallback all qualify -- confirmed via fftw_fprint_plan
      // for n=323=17*19, "hc2hc-generic-dit-17-19" whose own cld0/cld are
      // plain "rdft-generic-r2hc-17/19", radix 17 having no direct
      // codelet at all).
      // radix > GENERIC_SAFE_MAX (i.e. a prime needing its OWN dht-rader/
      // generic fallback): only reachable once dht-rader made such primes
      // trusted at all (this session) -- both instances found up to
      // N=8000 (n=7141=193*37, n=7913=193*41) are the classic "wrong
      // outer radix chosen by cost model" pattern (fftw_fprint_plan
      // confirms real FFTW picks the SMALLER prime as outer, e.g. radix=37
      // with 193 nested inside via dht-rader -- the reverse of what we
      // predict), not a codelet/dht-rader bug -- real FFTW's hc2hc-generic
      // apparently never actually prefers a >167 prime as the OUTER
      // combine radix over using it as a plain nested leaf. General rule,
      // not a per-radix enumeration, since the underlying mechanism (a
      // radix that itself needs a slow fallback is a bad outer combine
      // choice) isn't specific to 193.
      if (plan.codeletGroup === 'generic-combine' && plan.radix > GENERIC_SAFE_MAX) return false;
      if (plan.codeletGroup === 'generic-combine') return isFullyPortedR2HC(plan.radix) && walk(plan.sub);
      return false;
    }
    // 'generic' leaf -- see GENERIC_SAFE_MAX's header for the evidence.
    if (plan.type === 'generic') return isPrime(plan.n) && plan.n <= GENERIC_SAFE_MAX;
    // 'dht-rader' leaf -- see DhtRaderSolver.js's header. Safe exactly when
    // the n-1 R2HC sub-transform it's built on (unit stride, no mismatch
    // risk -- unlike RaderSolver.js/BluesteinSolver.js's complex sub-plans)
    // is itself already trusted.
    // pad:true plans are only ever a STRUCTURAL cost-comparison winner
    // (see Planner/chooseDecomposition.js's dhtRaderCost) -- DhtRaderSolver.js
    // doesn't implement the zero-padded-convolution execution path yet, so
    // never trust one even if the underlying npad transform is itself fine.
    if (plan.type === 'dht-rader') {
      return !plan.pad && !BAD_DHT_RADER_N.has(plan.n) && isFullyPortedDhtRader(plan.n, isFullyPortedR2HC);
    }
    return false;
  }
  return walk(chooseReal('R2HC', n).plan);
}

function genericFallback(n, x) {
  if (n % 2 === 0) {
    throw new Error(`RealEngine1D: no ported codelet for n=${n} (even) -- rdftGenericR2hc is only valid for odd n.`);
  }
  // rdftGenericR2hc returns an INTERLEAVED complex array (O[2k]=Re,
  // O[2k+1]=Im), not packed-HC -- see this file's header. Convert.
  const half = Math.floor(n / 2);
  const interleaved = new Float64Array((half + 1) * 2);
  rdftGenericR2hc(n, x, 0, 1, interleaved, 0, 2);
  const O = new Float64Array(n);
  for (let k = 0; k <= half; k++) O[k] = interleaved[2 * k];
  for (let k = 1; k < Math.ceil(n / 2); k++) O[n - k] = interleaved[2 * k + 1];
  return O;
}

// Execute a chooseReal('R2HC', n)-shaped plan against real input x (length
// n). Returns O, length n, in FFTW's packed-halfcomplex format (matching
// r2cf_n's own convention exactly).
function executePlanR2HC(rawPlan, n, x) {
  const plan = correctR2HCNode(rawPlan);
  if (plan.type === 'direct') {
    const fn = r2cf[plan.radix];
    if (fn) return fn(x);
    return genericFallback(n, x);
  }

  if (plan.type === 'ct') {
    const { radix: r, m, sub, needsCldm } = plan;
    const fn = plan.codeletGroup === 'hf2' ? hf2[r]
      : plan.codeletGroup === 'generic-combine' ? ((cr, ci, Wc, Ws) => genericCombineR2HC(r, cr, ci, Wc, Ws))
        : hf[r];
    // r can now be EVEN too (hf2 family, e.g. hf2_8 as an outer radix --
    // see RealEngine1D.js's r2cfII/hf2 registries) -- the old "n%2===0 iff
    // needsCldm" shortcut assumed r was always odd (true for hf/generic-
    // combine, but not hf2) and is no longer a valid proxy. The real gate
    // is simply: fn must exist, and if this m needs a cldm combine, that
    // codelet must exist too. genericFallback's own even-n throw is the
    // correct backstop if somehow neither holds.
    const cldmFn = needsCldm ? r2cfII[r] : null;
    if (!fn || (needsCldm && !cldmFn)) return genericFallback(n, x);

    // Stage 1: decimate into r phases (element p + k*r for k=0..m-1), full
    // size-m R2HC sub-transform per phase, each phase's packed-HC output
    // placed at work[p*m .. p*m+m-1] (matches rdft/hc2hc.c mkplan's cld
    // tensor: input stride r*is per phase, output stride os within a
    // contiguous m-block per phase -- see this file's header).
    const work = new Float64Array(n);
    const xPhase = new Float64Array(m);
    for (let p = 0; p < r; p++) {
      for (let k = 0; k < m; k++) xPhase[k] = x[p + k * r];
      const localHC = executePlanR2HC(sub, m, xPhase);
      for (let k = 0; k < m; k++) work[p * m + k] = localHC[k];
    }

    // Stage 2 (cld0): the r phases' own local DC values (work[p*m], one per
    // phase) are themselves an r-point real sequence -- combine via r's own
    // R2HC (usually a direct r2cf_r, but not always -- generic-combine's r
    // can lack one entirely, e.g. r=17, so this goes through the full
    // recursive engine exactly like genericCombineR2HC's own cld0/cld;
    // reduces to a plain r2cf[r] call whenever one exists, unchanged from
    // before for every hf/hf2 radix), written back to the same r positions.
    const dcVals = new Float64Array(r);
    for (let p = 0; p < r; p++) dcVals[p] = work[p * m];
    const dcOut = executePlanR2HC(chooseReal('R2HC', r).plan, r, dcVals);
    for (let p = 0; p < r; p++) work[p * m] = dcOut[p];

    // Stage 3 (hf_r): for each non-DC, non-middle column col=1..colLimit,
    // gather the r phases' local (Re,Im) pairs at (col, m-col), twiddle-
    // combine, write back in place. When needsCldm (m even), column m/2 is
    // EXCLUDED here (floor(m/2)===m/2 would otherwise self-pair col with
    // itself, which fn's ordinary two-column formula can't handle) and
    // handled separately below by Stage cldm.
    const colLimit = needsCldm ? m / 2 - 1 : Math.floor(m / 2);
    const cr = new Float64Array(r), ci = new Float64Array(r);
    const Wc = new Float64Array(r), Ws = new Float64Array(r);
    for (let col = 1; col <= colLimit; col++) {
      for (let p = 0; p < r; p++) {
        cr[p] = work[p * m + col];
        ci[p] = work[p * m + (m - col)];
      }
      for (let p = 1; p < r; p++) {
        const [c, s] = realCexp((col * p) % n, n);
        Wc[p] = c; Ws[p] = s;
      }
      const [outCr, outCi] = fn(cr, ci, Wc, Ws);
      for (let p = 0; p < r; p++) {
        work[p * m + col] = outCr[p];
        work[p * m + (m - col)] = outCi[p];
      }
    }

    // Stage cldm (R2HCII, only when m is even): column m/2 is the m-sized
    // sub-transform's own Nyquist bin -- a lone real value per phase with no
    // Im-pair, incompatible with fn's two-column twiddle formula. Gather the
    // r phases' values at that column into an r-point real sequence and
    // combine via r2cfII[r] (see that codelet's header for the R2HCII math
    // and its out[k]/out[r-1-k] convention), writing straight back to the
    // same r physical positions -- direct phase-index<->array-index
    // correspondence, the same convention already established for Stage
    // 2's cld0 combine above.
    if (needsCldm) {
      const mid = new Float64Array(r);
      for (let p = 0; p < r; p++) mid[p] = work[p * m + m / 2];
      const midOut = cldmFn(mid);
      for (let p = 0; p < r; p++) work[p * m + m / 2] = midOut[p];
    }

    // No final gather/remap needed: physical position IS the final packed-
    // halfcomplex output index, verified empirically (see this file's
    // header) against a from-scratch hand-computed DFT for n=9 -- same
    // "physical position = output index" invariant already established for
    // the complex side's dft-ct-dit (RaderSolver.js's dftSize36 header).
    return work;
  }

  if (plan.type === 'dht-rader') return applyR2HCViaDht(n, x);

  return genericFallback(n, x);
}

// dftR2HC -- returns FFTW's native packed-halfcomplex format (matching
// r2cf_n's own convention, see Codelets/real/r2cf_2.js's header).
function dftR2HC(n, x, xOff, xs, O, oOff, os) {
  const xArr = new Float64Array(n);
  for (let i = 0; i < n; i++) xArr[i] = x[xOff + i * xs];
  const plan = chooseReal('R2HC', n).plan;
  const packed = executePlanR2HC(plan, n, xArr);
  for (let i = 0; i < n; i++) O[oOff + i * os] = packed[i];
}

function isFullyPortedHC2R(n) {
  function walk(rawPlan) {
    const plan = correctHC2RNode(rawPlan);
    if (plan.type === 'direct') return plan.codeletGroup !== 'hb2' && !!r2cb[plan.radix];
    // radix===3, m===7: corrected (not excluded) via correctHC2RNode above.
    // radix===3, m===9: NOT corrected -- chooseReal's naive prediction
    // already matches real FFTW exactly for HC2R (see correctHC2RNode's
    // header) -- was previously excluded here too, unnecessarily.
    //
    // hb2 group, radix===5, m===5 ALSO excluded -- a DIFFERENT (and
    // asymmetric) mispredict found while verifying hb2_5.js: real FFTW's
    // n=25 HC2R plan is a single DIRECT "r2cb_25" codelet (unported, and
    // NOT the CT decomposition at all), while chooseReal('HC2R', 25)
    // predicts rdft-ct-dif/5(hb2_5, cld=r2cb_5). Confirmed via
    // fftw_fprint_plan -- note this is asymmetric with the R2HC direction:
    // real FFTW's n=25 R2HC plan genuinely IS rdft-ct-dit/5(hf2_5,
    // cld=r2cf_5), matching chooseReal('R2HC', 25) exactly (independently
    // verified bit-exact) -- so hf2_5.js is correct and this exclusion is
    // HC2R-only. hb2_5.js itself is a verbatim-correct transcription
    // (re-verified line-for-line against hb2_5.c); the ~3.5e-15 HC2R(25)
    // discrepancy was chooseReal picking a real-but-not-FFTW's-actual
    // decomposition, not a codelet bug.
    if (plan.type === 'ct') {
      // dht-rader EMBEDDED as a CT sub-transform (plan.sub, i.e. Stage 1's
      // per-phase recursion) is HC2R-specific broken, even though the
      // SAME prime alone (top-level, standalone HC2R) is bit-exact --
      // e.g. n=433 alone is fine, but n=866=2*433's "hb_2 x m=433 via
      // hc2r-dht(...)" is not (~1.1e-12 off). Same category as hb2_20's
      // "standalone fine, embedded-as-cld genuinely different real
      // structure" pattern (see that guard's header) -- confirmed via
      // direct sweep (3 root primes found up to N=8000: 433, 1129, 1297,
      // all cascading into many composites), not yet root-caused to a
      // specific FFTW-source mechanism. R2HC does NOT have this problem
      // (a parallel sweep found zero embedded-dht-rader R2HC mismatches),
      // so this exclusion is HC2R-only, matching this file's established
      // forward/backward-asymmetry pattern elsewhere.
      if (plan.sub && plan.sub.type === 'dht-rader') return false;
      if (plan.codeletGroup === 'hb2' && plan.radix === 5 && plan.m === 5) return false;
      // hb2, radix===20: BLANKET exclusion, not a narrow per-m blacklist --
      // a direct sweep (found while wiring radix 10 exposed a new pruning
      // candidate that changed which structure "wins" for many n divisible
      // by 20) showed ~100% of reachable hb2_20 HC2R cases are wrong (131
      // mismatches across a 4000-N sweep, sampled as 66/66 in a focused
      // check). Root-caused via fftw_fprint_plan, NOT a codelet bug: real
      // FFTW's actual HC2R plan for EVERY n divisible by 20 checked (20,
      // 40, 60, 80, 100, 200, 300, 400, 500, 700, 900, 1100, 1300) wraps in
      // "rdft-buffered", and for radix 20 specifically that wrapper's own
      // "cld" sub-transform uses a FORWARD-type codelet (e.g. "r2cf_2", not
      // "r2cb_2") -- a genuinely different, more complex internal structure
      // this port doesn't replicate (contrast radix 12's n=24/48/156, where
      // the same buffered wrapper's cld sub-transform matches our own
      // prediction exactly -- "r2cb_2" -- and IS bit-exact; buffering
      // itself is not the problem, radix 20's specific internal shape is).
      // hb2_20.js remains registered and is trusted for HC2R execution
      // paths reached via OTHER routes (e.g. as a sub-plan under a
      // different, verified-correct outer radix) -- only claims where
      // hb2_20 is itself a plan node are excluded here.
      if (plan.codeletGroup === 'hb2' && plan.radix === 20) return false;
      // See isFullyPortedR2HC's identical (radix,m) exclusions -- this is
      // the HC2R analogue, using the BAD_HB*_M sets (see their header for
      // why this is no longer a blanket needsCldm exclusion -- a direct
      // sweep after the radix 8/16/32 even-radix work landed found the
      // MAJORITY of needsCldm HC2R cases are actually bit-exact).
      if (plan.codeletGroup === 'hb' && plan.radix === 7 && BAD_HB_7_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hb' && plan.radix === 3 && plan.needsCldm && BAD_HB_3_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hb' && plan.radix === 6 && plan.needsCldm && BAD_HB_6_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hb' && plan.radix === 10 && plan.needsCldm && BAD_HB_10_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hb' && plan.radix === 9 && plan.needsCldm && BAD_HB_9_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hb2' && plan.radix === 5 && plan.needsCldm && BAD_HB2_5_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hb2' && plan.radix === 8 && BAD_HB2_8_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hb2' && plan.radix === 16 && BAD_HB2_16_M.has(plan.m)) return false;
      if (plan.codeletGroup === 'hb2' && plan.radix === 32 && BAD_HB2_32_M.has(plan.m)) return false;
      // r can be even now (hb2) -- check needsCldm directly against
      // r2cbIII's actual registered radixes rather than an n%2 shortcut.
      if (plan.needsCldm && !r2cbIII[plan.radix]) return false;
      if (plan.codeletGroup === 'hb') return !!hb[plan.radix] && walk(plan.sub);
      if (plan.codeletGroup === 'hb2') return !!hb2[plan.radix] && walk(plan.sub);
      // See isFullyPortedR2HC's identical guard for why r's own direct
      // codelet isn't required, and its radix>GENERIC_SAFE_MAX guard
      // right above it for why that's excluded here too.
      if (plan.codeletGroup === 'generic-combine' && plan.radix > GENERIC_SAFE_MAX) return false;
      if (plan.codeletGroup === 'generic-combine') return isFullyPortedHC2R(plan.radix) && walk(plan.sub);
      return false;
    }
    // 'generic' leaf -- see GENERIC_SAFE_MAX's header for the evidence.
    if (plan.type === 'generic') return isPrime(plan.n) && plan.n <= GENERIC_SAFE_MAX;
    // 'dht-rader' leaf -- ALWAYS checked against isFullyPortedR2HC (not
    // isFullyPortedHC2R), even here in the HC2R walk: dht-rader.c's own
    // cld1/cld2 sub-transforms are hardcoded R2HC kind regardless of the
    // OUTER rdft-dht wrapper's kind (R2HC_ONLY_CONV=1 in the C source --
    // see DhtRaderSolver.js's header), so that's genuinely what this path
    // depends on at runtime.
    // pad:true plans are only ever a STRUCTURAL cost-comparison winner
    // (see Planner/chooseDecomposition.js's dhtRaderCost) -- DhtRaderSolver.js
    // doesn't implement the zero-padded-convolution execution path yet, so
    // never trust one even if the underlying npad transform is itself fine.
    if (plan.type === 'dht-rader') {
      return !plan.pad && !BAD_DHT_RADER_N.has(plan.n) && isFullyPortedDhtRader(plan.n, isFullyPortedR2HC);
    }
    return false;
  }
  return walk(chooseReal('HC2R', n).plan);
}

function genericFallbackHC2R(n, packedHC) {
  if (n % 2 === 0) {
    throw new Error(`RealEngine1D: no ported codelet for n=${n} (even) -- rdftGenericHc2r is only valid for odd n.`);
  }
  const half = Math.floor(n / 2);
  const interleaved = new Float64Array((half + 1) * 2);
  for (let k = 0; k <= half; k++) {
    interleaved[2 * k] = packedHC[k];
    interleaved[2 * k + 1] = k === 0 ? 0 : packedHC[n - k];
  }
  const x = new Float64Array(n);
  rdftGenericHc2r(n, interleaved, 0, 2, x, 0, 1);
  return x;
}

// Execute a chooseReal('HC2R', n)-shaped plan against packed-halfcomplex
// input (length n, same format executePlanR2HC produces). Returns x, length
// n, real, UNNORMALIZED (matches r2cb_n's own convention -- see
// Codelets/real/r2cb_2.js's header).
function executePlanHC2R(rawPlan, n, packedHC) {
  const plan = correctHC2RNode(rawPlan);
  if (plan.type === 'direct') {
    const fn = r2cb[plan.radix];
    if (fn) return fn(packedHC);
    return genericFallbackHC2R(n, packedHC);
  }

  if (plan.type === 'ct') {
    const { radix: r, m, sub, needsCldm } = plan;
    const fn = plan.codeletGroup === 'hb2' ? hb2[r]
      : plan.codeletGroup === 'generic-combine' ? ((cr, ci, Wc, Ws) => genericCombineHC2R(r, cr, ci, Wc, Ws))
        : hb[r];
    // See executePlanR2HC's identical guard -- r can be even now (hb2), so
    // n%2===0 is no longer a valid proxy for needsCldm.
    const cldmFn = needsCldm ? r2cbIII[r] : null;
    if (!fn || (needsCldm && !cldmFn)) return genericFallbackHC2R(n, packedHC);

    // DIF order (rdft/hc2hc.c's apply_dif): twiddle stage FIRST, in place
    // on the input; sub-transform SECOND. Mirror image of R2HC's DIT order
    // (sub-transform first, twiddle second) -- see this file's header.
    const work = Float64Array.from(packedHC.subarray ? packedHC.subarray(0, n) : packedHC.slice(0, n));

    // Stage 1 (hb_r): for each non-DC, non-middle column col=1..colLimit,
    // gather the r phases' OVERALL (Re,Im) pairs at physical positions
    // (col, m-col) (same "physical position = frequency index" property as
    // R2HC), twiddle-uncombine, write back -- now representing each phase's
    // own LOCAL halfcomplex block instead of the overall spectrum. Column
    // m/2 is excluded (see executePlanR2HC's identical Stage 3 comment) and
    // handled by Stage cldm below instead.
    const colLimit = needsCldm ? m / 2 - 1 : Math.floor(m / 2);
    const cr = new Float64Array(r), ci = new Float64Array(r);
    const Wc = new Float64Array(r), Ws = new Float64Array(r);
    for (let col = 1; col <= colLimit; col++) {
      for (let p = 0; p < r; p++) {
        cr[p] = work[p * m + col];
        ci[p] = work[p * m + (m - col)];
      }
      for (let p = 1; p < r; p++) {
        const [c, s] = realCexp((col * p) % n, n);
        Wc[p] = c; Ws[p] = s;
      }
      const [outCr, outCi] = fn(cr, ci, Wc, Ws);
      for (let p = 0; p < r; p++) {
        work[p * m + col] = outCr[p];
        work[p * m + (m - col)] = outCi[p];
      }
    }

    // Stage cldm (HC2RIII, only when m is even): the r phases' packed
    // R2HCII-format values at column m/2 (same physical positions
    // executePlanR2HC's Stage cldm wrote, in the same r2cfII out[k]/
    // out[r-1-k] convention -- r2cbIII[r]'s input convention is that exact
    // format, see that codelet's header) invert via r2cbIII[r] straight
    // back to each phase's own local real Nyquist value, written back to
    // the same r physical positions (direct phase-index<->array-index
    // correspondence, matching Stage 2's cld0 convention below). Must run
    // BEFORE Stage 3 recurses per-phase, since each phase's local block
    // needs its own middle value already resolved.
    if (needsCldm) {
      const mid = new Float64Array(r);
      for (let p = 0; p < r; p++) mid[p] = work[p * m + m / 2];
      const midOut = cldmFn(mid);
      for (let p = 0; p < r; p++) work[p * m + m / 2] = midOut[p];
    }

    // Stage 2 (cld0 inverse): the r overall DC-related bins (work[p*m],
    // p=0..r-1) are themselves a packed-HC r-point spectrum -- invert via
    // r's own HC2R (see executePlanR2HC's identical Stage 2 comment for
    // why this goes through the full recursive engine rather than a
    // hardcoded r2cb[r]), written back to the same r positions, now each
    // phase's own local DC value.
    const dcVals = new Float64Array(r);
    for (let p = 0; p < r; p++) dcVals[p] = work[p * m];
    const dcOut = executePlanHC2R(chooseReal('HC2R', r).plan, r, dcVals);
    for (let p = 0; p < r; p++) work[p * m] = dcOut[p];

    // Stage 3 (cld): each phase's now-complete local m-sized packed-HC
    // block inverts to m real values via the SAME recursive plan tree
    // chooseReal('HC2R', m) predicted (sub), written INTERLEAVED by r
    // (matches rdft/hc2hc.c mkplan's cld tensor for HC2R: output stride
    // r*os for the size dim, os for the vector/phase dim).
    const x = new Float64Array(n);
    const blockHC = new Float64Array(m);
    for (let p = 0; p < r; p++) {
      for (let k = 0; k < m; k++) blockHC[k] = work[p * m + k];
      const localReal = executePlanHC2R(sub, m, blockHC);
      for (let k = 0; k < m; k++) x[p + k * r] = localReal[k];
    }
    return x;
  }

  if (plan.type === 'dht-rader') return applyHC2RViaDht(n, packedHC);

  return genericFallbackHC2R(n, packedHC);
}

function dftHC2R(n, packed, pOff, ps, x, xOff, xs) {
  const packedArr = new Float64Array(n);
  for (let i = 0; i < n; i++) packedArr[i] = packed[pOff + i * ps];
  const plan = chooseReal('HC2R', n).plan;
  const real = executePlanHC2R(plan, n, packedArr);
  for (let i = 0; i < n; i++) x[xOff + i * xs] = real[i];
}

// dftHC2RFromInterleaved -- same computation, but takes a single
// INTERLEAVED complex array as input (X[2k]=Re,X[2k+1]=Im), matching
// GenericSolver1D.rdftGenericHc2r's calling convention exactly so it's a
// drop-in replacement at Rank2Orchestration.js's call site.
function dftHC2RFromInterleaved(n, X, xOff, xs, O, oOff, os) {
  const half = Math.floor(n / 2);
  const packedArr = new Float64Array(n);
  for (let k = 0; k <= half; k++) {
    packedArr[k] = X[xOff + k * xs];
    if (k > 0) packedArr[n - k] = X[xOff + k * xs + 1];
  }
  const plan = chooseReal('HC2R', n).plan;
  const real = executePlanHC2R(plan, n, packedArr);
  for (let i = 0; i < n; i++) O[oOff + i * os] = real[i];
}

// dftR2HCInterleaved -- same computation, but returns a single INTERLEAVED
// complex array (O[2k]=Re,O[2k+1]=Im for bin k=0..floor(n/2)), matching
// GenericSolver1D.rdftGenericR2hc's calling convention exactly so it's a
// drop-in replacement at Rank2Orchestration.js's call site.
function dftR2HCInterleaved(n, x, xOff, xs, O, oOff, os) {
  const xArr = new Float64Array(n);
  for (let i = 0; i < n; i++) xArr[i] = x[xOff + i * xs];
  const plan = chooseReal('R2HC', n).plan;
  const packed = executePlanR2HC(plan, n, xArr);
  const half = Math.floor(n / 2);
  for (let k = 0; k <= half; k++) {
    O[oOff + k * os] = packed[k];
    O[oOff + k * os + 1] = (k === 0 || (n % 2 === 0 && k === half)) ? 0 : packed[n - k];
  }
}

module.exports = {
  dftR2HC, dftR2HCInterleaved, executePlanR2HC, isFullyPortedR2HC,
  dftHC2R, dftHC2RFromInterleaved, executePlanHC2R, isFullyPortedHC2R,
};

// From-scratch hand-computed R2HC oracle (deliberately independent of
// rdftGenericR2hc -- see this file's header for why an earlier attempt to
// use that function as the oracle was itself miscalled, not buggy). Used
// only here, for the self-test below.
function handR2HC(x) {
  const n = x.length;
  const re = new Array(n).fill(0), im = new Array(n).fill(0);
  for (let k = 0; k < n; k++) {
    for (let j = 0; j < n; j++) {
      const ang = -2 * Math.PI * k * j / n;
      re[k] += x[j] * Math.cos(ang);
      im[k] += x[j] * Math.sin(ang);
    }
  }
  const O = new Array(n).fill(0);
  for (let k = 0; k <= Math.floor(n / 2); k++) O[k] = re[k];
  for (let k = 1; k < Math.ceil(n / 2); k++) O[n - k] = im[k];
  return O;
}

// ---------------------------------------------------------------------------
// Self-test (node RealEngine1D.js) -- compares against a from-scratch
// hand-computed DFT (not rdftGenericR2hc, which has its own bug for
// composite odd n -- see this file's header) across composite odd N.
// ---------------------------------------------------------------------------
if (require.main === module) {
  for (const n of [9, 15, 21, 25, 27, 45]) {
    const x = new Float64Array(n);
    for (let i = 0; i < n; i++) x[i] = Math.sin(i * 1.3 + 0.4);
    const plan = chooseReal('R2HC', n).plan;
    const fullyPorted = isFullyPortedR2HC(n);
    try {
      const got = executePlanR2HC(plan, n, x);
      const want = handR2HC(Array.from(x));
      let maxDiff = 0;
      for (let i = 0; i < n; i++) maxDiff = Math.max(maxDiff, Math.abs(got[i] - want[i]));
      console.log(`n=${n}: maxDiff=${maxDiff} ${maxDiff < 1e-9 ? 'OK' : 'FAIL'} fullyPorted=${fullyPorted} plan=${chooseReal('R2HC', n).describe()}`);
    } catch (e) {
      console.log(`n=${n}: THROWN (${e.message}) fullyPorted=${fullyPorted} plan=${chooseReal('R2HC', n).describe()}`);
    }
  }
}
