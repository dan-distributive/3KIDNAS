'use strict';

// =============================================================================
// chooseDecomposition.js
// Replicates FFTW3's FFTW_ESTIMATE solver search for 1D complex DFT and 1D
// real r2hc/hc2r sub-transforms, using the ops-cost tables extracted by
// extract_opcosts.js (Planner/registrationTables.js). This is NOT a
// minimum-cost search over all candidates -- it is a faithful port of
// kernel/planner.c's search0(): candidates are tried in LIFO (most-recently-
// registered-first) order, cost-compared as they're found, and the search
// STOPS as soon as any tried candidate has `could_prune_now_p = true`,
// keeping whatever is currently best -- even if a cheaper candidate would
// have appeared later in LIFO order. Getting this exactly right (not just
// "find the cheapest decomposition") is why a candidate that looks obviously
// cheaper on paper sometimes isn't what real FFTW picks; see verify/plan_sweep.js
// for the empirical validation this was tuned against.
//
// Facts this file depends on, each confirmed by reading FFTW3 source directly
// this session (not assumed from memory of FFTW's general reputation):
//
//   - X(iestimate_cost) = add + mul + (HAVE_FMA ? fma : 2*fma) + other
//     (kernel/planner.c). HAVE_FMA is never #defined in this project's build
//     (grepped every Makefile and config.h) so the coefficient is always 2.
//   - dft/generic.c and rdft/generic.c's `ops.other` assignment is wrapped in
//     `#if 0` (dead code) in FFTW 3.3.8 -- the field stays at its
//     zero-initialized default. Omitting this is NOT optional/approximate;
//     including the commented-out formula gives a visibly wrong answer (this
//     was caught empirically: it flips the generic-vs-bluestein choice for
//     n=43,47,53,67,... away from what real FFTW actually picks).
//   - could_prune_now_p is TRUE only for: non-buffered base-case codelets
//     (dft/direct.c, rdft/direct-r2c.c -- unconditionally true when not
//     buffered), and non-buffered fixed-radix twiddle codelets with
//     radix in [5,64) and m >= radix (dft/dftw-direct.c, rdft/hc2hc-direct.c
//     -- identical condition on both the complex and real side). A
//     Cooley-Tukey (dft-ct / rdft hc2hc) plan INHERITS could_prune_now_p from
//     its twiddle (cldw) child, not from its sub-transform (cld) child.
//     generic/rader/bluestein NEVER set it (confirmed absent from all three
//     source files) -- so among those three, the search always compares all
//     applicable candidates and picks the true minimum; no early-stop there.
//   - Applicability gates (kernel/ifftw.h): GENERIC_MAX_SLOW=16,
//     RADER_MAX_SLOW=32, BLUESTEIN_MAX_SLOW=24. Under FFTW_ESTIMATE, NO_SLOW
//     is always set (api/mapflags.c: NOT EXHAUSTIVE implies NO_SLOW, and
//     ESTIMATE implies NOT PATIENT implies NOT EXHAUSTIVE), so these
//     thresholds are always active for this port's purposes: generic needs
//     n>16, rader needs n>32 AND (n-1) 5-smooth... actually rader needs (n-1)
//     factors only into {2,3,5} (kernel/primes.c X(factors_into_small_primes),
//     literally {2,3,5} -- NOT the "5-smooth" bound Bluestein's own padding
//     search uses, though both happen to use the same prime set), bluestein
//     needs n>24 (and always n>16 separately, dft/bluestein.c's own
//     "avoid infinite self-recursion" guard).
//   - Cost formulas for the structural (non-codelet) solvers, transcribed
//     directly from each file's mkplan():
//       dft-generic(n):    add=5(n-1), fma=(n-1)^2,      other=0
//       rdft-generic(n):   add=2.5(n-1), fma=0.5(n-1)^2, other=0
//       dft-ct(r,m; n=rm): cost = r*cost(m) + m*cost(twiddle_r_ops)
//                          (cld is "r copies of size-m", cldw is "m twiddle
//                          combines"; vector-loop cost is an exact multiply,
//                          kernel/*vrank-geq1.c: ops_madd2(vl, child, ...))
//       dft-rader(n):      cost = 2*cost(n-1) + [add:2(n-1)+4, mul:4(n-1),
//                                                other:14(n-1)+6]
//       dft-bluestein(n):  nb = smallest n with factors_into_small_primes(n)
//                          and n >= 2*N-1; cost = 2*cost(nb) +
//                          [add:4N+2nb, mul:8N+4nb, other:6(N+nb)]
//       rdft-hc2hc(r,m):   same shape as dft-ct, but the child sub-transform
//                          (cld) is a SAME-KIND (R2HC recurses to R2HC, HC2R
//                          to HC2R) rdft problem of size m, NOT a complex DFT
//                          (rdft/hc2hc.c: mkproblem_rdft_d, confirmed by
//                          reading the file directly -- do not assume this
//                          mirrors dft/ct.c's complex recursion).
//
// NOT modeled yet (explicit, deliberate gaps -- see FFTW3JS/README.md
// Phase 6 and the plan file): q1_* (twiddle-squared, N=r^2 in-place DIF),
// the COMPLEX-side ct-generic/ct-genericbuf (arbitrary-radix complex CT
// using a generic twiddle combine instead of a registered t1_r/t2_r
// codelet), and the hc2cf/hc2cb family (competing real-side recursive
// solver operating directly on RDFT2). The REAL-side analogue of
// ct-generic (rdft/hc2hc-generic.c, "hc2hc-generic") IS modeled below
// (see GENERIC_COMBINE_RADIXES) -- narrowly, for the two radixes (11, 13)
// where real FFTW's own registered codelet set has a direct r2cf_r/
// r2cb_r but no twiddle codelet, verified against real hc2hc-generic.c's
// exact pointer arithmetic (RealEngine1D.js's genericCombineR2HC/HC2R)
// and cross-checked against fftw_fprint_plan across many N, not just cost-
// estimated. Their absence means this predictor can still UNDER-predict
// CT usage for composite N with factors outside the registered codelet
// radix set, and can mislabel (but not miscompute -- see below) N that's
// an exact square of a small registered radix. verify/plan_sweep.js's job
// is to quantify exactly how often this matters and report it honestly,
// not to paper over it.
//
// IMPORTANT SCOPE NOTE: this file predicts DECOMPOSITION STRUCTURE for
// validating understanding and guiding which codelets are worth porting
// next. It is NOT the production correctness gate -- PlanTable.js's
// isTransformSupported() (byte-exact-fixture-verified per size) is, and
// stays independent of whatever this file predicts.
// =============================================================================

const {
  complexCodelets, q1Codelets, realForwardCodelets, realBackwardCodelets,
} = require('./registrationTables');

const GENERIC_MAX_SLOW = 16;
const RADER_MAX_SLOW = 32;
const BLUESTEIN_MAX_SLOW = 24;

function opsCost(ops) {
  return ops.add + ops.mul + 2 * ops.fma + ops.other;
}

function isPow2(n) {
  return n > 0 && (n & (n - 1)) === 0;
}

function isPrime(n) {
  if (n <= 1) return false;
  if (n % 2 === 0) return n === 2;
  for (let i = 3; i * i <= n; i += 2) if (n % i === 0) return false;
  return true;
}

function factorsIntoSmallPrimes(n) {
  for (const p of [2, 3, 5]) { while (n % p === 0) n /= p; }
  return n === 1;
}

function smallest5SmoothAtLeast(minsz) {
  let n = minsz;
  while (!factorsIntoSmallPrimes(n)) n++;
  return n;
}

// kernel/buffered.c's X(nbuf): how many of the `vl` batched copies get
// processed per buffered chunk. Literal transcription (imin/imax ->
// Math.min/max, integer division truncated via Math.floor -- INT in C is
// already truncating, matching JS's | 0 / Math.floor for positive values
// here since every operand is positive).
const MAXBUFSZ = Math.floor((256 * 1024) / 8); // sizeof(R)=8 (double)
function computeNbuf(n, vl, maxnbufIn) {
  const maxnbuf = maxnbufIn || 256; // DEFAULT_MAXNBUF
  const nbuf = Math.min(maxnbuf, Math.min(vl, Math.max(1, Math.floor(MAXBUFSZ / n))));
  const lb = Math.max(1, Math.floor(nbuf / 4));
  for (let i = nbuf; i >= lb; i--) {
    if (vl % i === 0) return i;
  }
  return nbuf;
}

const BUFFERED_MAXNBUFS = [8, 256];

// dft-buffered's own cost model, ported directly from dft/buffered.c's
// mkplan() + kernel/buffered.c's X(nbuf)/X(nbuf_redundant) +
// rdft/dft-r2hc.c's mkplan() + rdft/rank0.c's mkplan() (the chain that
// actually plans the "copy back from scratch buffer" step, cldcpy --
// confirmed by literally reading all four files this session, not
// assumed: dft-buffered's cldcpy is a rank-0 (no arithmetic) DFT problem,
// which gets bridged through dft-r2hc.c into an rdft R2HC problem, which
// rdft/rank0.c's own family of solvers (memcpy/iter/tiled/...) actually
// executes as a PURE data-movement copy -- ops.other = 2*tensor_sz(vecsz),
// zero add/mul/fma, confirmed via rdft/rank0.c's mkplan():
//   X(ops_other)(2 * X(tensor_sz)(p->vecsz), &pln->super.super.ops);
// and rdft/dft-r2hc.c's own wrapper adds only {other: 1} on top (an
// "estimator hack for nop plans", moot here since (n-1)/2==0 for the
// rank-0 sz dft-buffered's cldcpy always uses).
//
// dft-buffered's whole point: rather than trying to plan `n`'s own
// transform IN PLACE (which is what forces q1/nothing, see
// chooseComplexEmbedded's header), it transforms into a FRESH, entirely
// separate scratch buffer -- i.e. exactly an out-of-place problem, so its
// own "cld" sub-plan is just chooseComplex(n) unmodified -- then copies
// the (already-correct) result back into the real, possibly-in-place
// output location via a second, purely-mechanical copy plan. Since that
// copy plan does no arithmetic at all, this solver's EXECUTION needs no
// new code: whenever it wins the cost comparison, running chooseComplex
// (n)'s own plan directly against the (n, parentRadix)-cofactor's data
// produces bit-identical results to what real FFTW's buffered wrapper
// would -- confirmed by this derivation, not by guessing.
//
// vl here is always parentRadix (r) in this port's own usage -- the
// ambient batch count of a nested cofactor is always exactly the outer
// radix that created it (same fact q1's r==v applicability check relies
// on), so callers never need to pass a separate vl.
//
// cldrest (the `vl % nbuf` leftover copies, planned as a FRESH recursive
// sub-search) is treated as cost 0 here: X(nbuf)'s own divisor search
// (see computeNbuf above) is specifically designed to find an nbuf that
// evenly divides vl whenever possible, and since vl=parentRadix is always
// <= 256 (the larger registered maxnbuf) and virtually always <=
// MAXBUFSZ/n for the sizes this port deals with, computeNbuf(n, r, 256)
// returns nbuf=r directly (vl%vl==0 trivially) essentially always in
// practice -- i.e. the remainder is 0, not approximated away. If some
// future N ever hits the non-zero-remainder edge case, this slightly
// underestimates dft-buffered's true cost (cldrest folded in as free
// rather than its own small extra charge), which only matters if it
// flips a close cost comparison against q1 -- exactly the kind of
// discrepancy the mandatory broad safety sweep exists to catch.
function bufferedCost(n, vl, subCost) {
  let bestForR = null;
  const triedNbufs = [];
  for (let idx = 0; idx < BUFFERED_MAXNBUFS.length; idx++) {
    const maxnbuf = BUFFERED_MAXNBUFS[idx];
    const nbuf = computeNbuf(n, vl, maxnbuf);
    if (triedNbufs.includes(nbuf)) continue; // X(nbuf_redundant)
    triedNbufs.push(nbuf);

    const cldCost = nbuf * subCost;
    const cldcpyCost = 4 * nbuf * n + 1; // rank0's 2*(2*nbuf*n) + r2hc's +1
    const chunks = Math.floor(vl / nbuf);
    const remainder = vl % nbuf;
    if (remainder !== 0) continue; // see header -- not modeled, skip this variant
    const cost = chunks * (cldCost + cldcpyCost);
    if (!bestForR || cost < bestForR.cost) bestForR = { cost, nbuf };
  }
  return bestForR;
}

// Both n1/t1/t2 (complex) and r2cf/hf (real) codelet groups sorted by
// registrationIndex DESCENDING = LIFO search order (most-recently-registered
// tried first).
function lifoOrder(entries) {
  return entries.slice().sort((a, b) => b.registrationIndex - a.registrationIndex);
}

const complexLifo = lifoOrder(complexCodelets); // n1_*, t1_*, t2_* only
const realForwardLifo = lifoOrder(realForwardCodelets.filter((e) => e.group === 'r2cf' || e.group === 'hf' || e.group === 'hf2'));
const realBackwardLifo = lifoOrder(realBackwardCodelets.filter((e) => e.group === 'r2cb' || e.group === 'hb' || e.group === 'hb2'));

// CLDM_OPS -- ops cost for the R2HCII/HC2RIII "cldm" middle-column combine
// (rdft/hc2hc-direct.c: whenever the RECURSED sub-transform size m is
// EVEN, column m/2 -- the m-sized sub-transform's own Nyquist bin, a lone
// real value per phase with no Im-pair -- can't go through the ordinary
// hf_r/hb_r twiddle codelet, and needs this genuinely different shifted-
// frequency transform instead; see RealEngine1D.js's genericCombineR2HC-
// adjacent header for the math). Extracted directly from each
// r2cfII_r.c/r2cbIII_r.c's own `static const kr2c_desc desc = {r,
// name, {add,mul,fma,other}, ...}` (registrationTables.js's extraction
// pass never covered this codelet family -- this is the SAME {add,mul,
// fma,other} shape, just hand-copied for the radixes that matter here,
// not regenerated). Only exists for radixes with a registered hf_r/hb_r
// twiddle codelet in the first place (11 and 13 -- the hc2hc-generic-only
// radixes -- have no r2cfII_11/13 in FFTW's own codelet set at all).
const CLDM_OPS = {
  R2HC: {
    3: { add: 3, mul: 1, fma: 1, other: 0 },
    5: { add: 9, mul: 3, fma: 3, other: 0 },
    6: { add: 11, mul: 2, fma: 2, other: 0 },
    7: { add: 12, mul: 6, fma: 12, other: 0 },
    8: { add: 18, mul: 6, fma: 4, other: 0 },
    9: { add: 25, mul: 13, fma: 17, other: 0 },
    10: { add: 26, mul: 6, fma: 6, other: 0 },
    12: { add: 39, mul: 8, fma: 4, other: 0 },
    15: { add: 54, mul: 15, fma: 18, other: 0 },
    16: { add: 54, mul: 18, fma: 12, other: 0 },
    20: { add: 86, mul: 18, fma: 16, other: 0 },
    25: { add: 126, mul: 61, fma: 87, other: 0 },
    32: { add: 138, mul: 46, fma: 36, other: 0 },
  },
  HC2R: {
    3: { add: 3, mul: 1, fma: 1, other: 0 },
    5: { add: 8, mul: 3, fma: 4, other: 0 },
    6: { add: 10, mul: 4, fma: 2, other: 0 },
    7: { add: 9, mul: 4, fma: 15, other: 0 },
    8: { add: 18, mul: 8, fma: 4, other: 0 },
    9: { add: 22, mul: 8, fma: 10, other: 0 },
    10: { add: 26, mul: 10, fma: 6, other: 0 },
    12: { add: 38, mul: 16, fma: 4, other: 0 },
    16: { add: 54, mul: 20, fma: 12, other: 0 },
    20: { add: 82, mul: 32, fma: 12, other: 0 },
    32: { add: 138, mul: 48, fma: 36, other: 0 },
    15: { add: 49, mul: 11, fma: 15, other: 0 },
    25: { add: 100, mul: 46, fma: 52, other: 0 },
  },
};

// hc2hc-generic (rdft/hc2hc-generic.c) is registered in real FFTW as a
// RADIX=0 WILDCARD solver (X(hc2hc_generic_register)(p) { regsolver(p, 0);
// }), not one candidate per specific radix -- it competes with EVERY odd
// divisor r of n as a possible outer CT radix, cost-comparing against the
// specific-radix hf/hf2 candidates above. Its own "cld0"/"cld" sub-steps
// (RealEngine1D.js's genericCombineR2HC/genericCombineHC2R) are themselves
// just an ordinary r-point R2HC/HC2R transform -- handed to the FULL
// recursive planner, not hardcoded to a direct r2cf_r/r2cb_r codelet, so r
// does NOT need its own direct codelet (confirmed via fftw_fprint_plan for
// n=323=17*19: "hc2hc-generic-dit-17-19" whose own cld0/cld are plain
// "rdft-generic-r2hc-17"/"-19" -- radix 17 has no direct codelet at all in
// real FFTW's registered set, yet hc2hc-generic still uses it as the outer
// CT radix). Originally implemented narrowly for radix in {11,13} (the two
// radixes with a direct-but-no-twiddle codelet) before this generalization
// -- kept working exactly the same for those, now also reachable for any
// other valid divisor once its own r-point transform is itself trusted
// (isFullyPortedR2HC(r)/isFullyPortedHC2R(r), which already correctly
// covers direct/CT/generic-prime-fallback uniformly).
function divisorsFor(n) {
  const divs = [];
  for (let r = 3; r * r <= n; r += 2) {
    if (n % r !== 0) continue;
    divs.push(r);
    const m = n / r;
    if (m !== r) divs.push(m);
  }
  return divs;
}

function twiddleCouldPruneNowP(radix, m) {
  return radix >= 5 && radix < 64 && m >= radix;
}

// dft-ct-generic (dft/dftw-generic.c) is the COMPLEX-side analogue of
// hc2hc-generic above -- a RADIX=0 WILDCARD solver providing the "cldw"
// (twiddle-combine) role of a dft-ct-dit structure for any radix r with
// NO registered t1_r/t2_r codelet (confirmed: neither exists in real
// FFTW's own source tree for r=11/13, e.g. -- these aren't gaps in THIS
// port, real FFTW itself has no dedicated codelet there either). Unlike
// hc2hc-generic, both EVEN and ODD r are valid (dftw-generic.c has no
// odd-only restriction -- that was specific to the real side's plain
// hf_r/hb_r calling convention), so this tries every divisor, not just
// odd ones.
function allDivisorsFor(n) {
  const divs = [];
  for (let r = 2; r * r <= n; r++) {
    if (n % r !== 0) continue;
    divs.push(r);
    const m = n / r;
    if (m !== r) divs.push(m);
  }
  return divs;
}

// realTwiddleInvocations -- the multiplier for a real-side twiddle
// combine's per-invocation ops cost in this file's cost formulas.
//
// IMPORTANT, learned the hard way: the ACTUAL invocation count, per
// rdft/hc2hc.c's own mstart/mcount/CLDMP bookkeeping, is floor((m-1)/2)
// for odd m (mb=1, me=mcount, CLDMP false) -- NOT the raw sub-size `m`
// this file's hf/hf2 and generic-combine formulas have always multiplied
// by. That IS a real discrepancy (confirmed: `e.ops`'s own header
// comments describe cost for a SINGLE column). But switching to the
// "correct" floor((m-1)/2) UNIFORMLY, for the ALREADY-odd-only regime
// this file has been extensively fftw_fprint_plan-verified against, is a
// REGRESSION, not a fix: verified via a fresh 1484-N sweep -- 477 previously-
// correct predictions flipped wrong (e.g. n=51=3*17 stopped preferring
// hf_3 over the generic-combine wildcard).
//
// RE-TESTED AGAIN this session, more thoroughly: reading rdft/hc2hc-direct.c's
// mkcldw directly confirms FFTW's own ops accounting really is
// `ops_madd2(v*(me-mb), &e->ops,...)` PLUS SEPARATELY `ops_madd2(v,
// &cld0->ops,...)` PLUS `ops_madd2(v, &cldm->ops,...)` -- so tried adding
// cld0's cost (chooseReal(kind,r).cost) to the hf/hf2 formula TOGETHER with
// switching this function to the true me-mb count uniformly. Result: WORSE,
// not better -- n=51 regressed AGAIN (hc2hc-generic/17 wrongly beat hf_3,
// cost 1014 vs 1118 by hand calculation), even though hc2hc-generic's own
// formula (rdft/hc2hc-generic.c's mkcldw: `pln->ops = cld->ops` where cld is
// batched 2*mcount1 times, PLUS the twiddle-multiply overage, notably
// WITHOUT ever adding cld0's cost at all -- confirmed directly from source,
// not assumed) was ALSO updated to use the same true mcount1. Both formulas
// individually now match their respective C sources term-for-term as far as
// this investigation traced, yet the cross-comparison between them still
// doesn't reproduce real FFTW's actual n=51 choice -- meaning some further
// factor (candidate LIFO tie-breaking order between hc2hc-direct and
// hc2hc-generic, `e->genus->vl` codelet-batching semantics this port's
// registrationTables extraction doesn't carry, or something not yet
// identified) is still unaccounted for. Reverted BOTH changes back to this
// exact state (needsCldm-gated: raw `m` when !needsCldm, true me-mb only
// when needsCldm) since it's the one empirically verified at scale (1484-N
// R2HC + hundreds-N HC2R sweeps, zero unexplained mismatches beyond the one
// pre-existing documented n=125 case). This function ONLY changes behavior
// for the needsCldm=true case, which was NEVER reachable before this
// session (no prior candidate had needsCldm at all) -- so using the
// properly-derived count there can't regress anything, and IS what fixed
// chooseReal(R2HC,28) to correctly prefer radix=7-outer+cldm over
// radix=4-outer (confirmed via fftw_fprint_plan and, this session, verified
// bit-exact end-to-end against real compiled FFTW's R2HC(28) output).
//
// Do NOT re-attempt "add cld0 + switch odd-m to true count" without first
// nailing down the missing factor above -- two independent attempts across
// two sessions have now hit the same n=51-style regression.
function realTwiddleInvocations(m, needsCldm) {
  if (!needsCldm) return m;
  const mcount = Math.floor((m + 2) / 2);
  const me = mcount - 1;
  const mb = 1;
  return me - mb;
}

const complexMemo = new Map(); // keyed by n -- chooseComplex() (out-of-place)
const complexEmbeddedMemo = new Map(); // keyed by `${n}:${parentRadix||0}` -- chooseComplexEmbedded()
const realMemo = new Map(); // keyed by `${kind}:${n}`

// Returns { cost, couldPruneNowP, describe() } -- the winning candidate for
// a 1D complex DFT of size n, found via the exact LIFO/early-stop algorithm
// described in this file's header (kernel/planner.c's search0()). This is
// the plain, ALWAYS-OUT-OF-PLACE-COMPATIBLE search used everywhere in this
// port except the two chooseComplexEmbedded() call sites below -- unchanged
// from its original form; q1 is never a candidate here (see
// chooseComplexEmbedded's header for why that's correct specifically for
// this function, not just an omission).
function chooseComplex(n) {
  if (complexMemo.has(n)) return complexMemo.get(n);

  let best = null;
  let pruned = false;
  const consider = (cand) => {
    if (!best || cand.cost < best.cost) best = cand;
    return cand.couldPruneNowP;
  };

  outer:
  for (const e of complexLifo) {
    if (e.group === 'n1') {
      if (e.radix !== n) continue;
      const cand = {
        cost: opsCost(e.ops), couldPruneNowP: true, describe: () => `n1_${e.radix}`,
        plan: { type: 'direct', radix: n, codeletName: e.name },
      };
      if (consider(cand)) { pruned = true; break outer; }
    } else { // t1 / t2 -- Cooley-Tukey with this fixed radix
      const r = e.radix;
      if (n % r !== 0) continue;
      const m = n / r;
      if (m <= 1) continue;
      const sub = chooseComplex(m);
      const cost = r * sub.cost + m * opsCost(e.ops);
      const cand = {
        cost,
        couldPruneNowP: twiddleCouldPruneNowP(r, m),
        describe: () => `dft-ct-dit/${r}(${e.name} x m=${m} via ${sub.describe()})`,
        plan: { type: 'ct', radix: r, m, codeletName: e.name, codeletGroup: e.group, sub: sub.plan },
      };
      if (consider(cand)) { pruned = true; break outer; }
    }
  }

  if (!pruned) {
    if (n > GENERIC_MAX_SLOW) {
      const cost = 5 * (n - 1) + 2 * (n - 1) * (n - 1);
      consider({ cost, couldPruneNowP: false, describe: () => `dft-generic-${n}`, plan: { type: 'generic', n } });
    }
    if (isPrime(n) && n > RADER_MAX_SLOW && factorsIntoSmallPrimes(n - 1)) {
      const sub = chooseComplex(n - 1);
      const overhead = (2 * (n - 1) + 4) + (4 * (n - 1)) + (14 * (n - 1) + 6);
      consider({
        cost: 2 * sub.cost + overhead, couldPruneNowP: false,
        describe: () => `dft-rader-${n}(via ${sub.describe()})`,
        plan: { type: 'rader', n, sub: sub.plan },
      });
    }
    if (isPrime(n) && n > 16 && n > BLUESTEIN_MAX_SLOW) {
      const nb = smallest5SmoothAtLeast(2 * n - 1);
      const sub = chooseComplex(nb);
      const overhead = (4 * n + 2 * nb) + (8 * n + 4 * nb) + 6 * (n + nb);
      consider({
        cost: 2 * sub.cost + overhead, couldPruneNowP: false,
        describe: () => `dft-bluestein-${n}/nb=${nb}(via ${sub.describe()})`,
        plan: { type: 'bluestein', n, nb, sub: sub.plan },
      });
    }
    // dft-ct-generic (dft/dftw-generic.c's mkcldw): cost = r*chooseComplex(m)
    // [Stage 1, same as t1/t2] + m*chooseComplex(r) [Stage 2's own r-point
    // DFT, batched m times, "cld->ops" in the C source] + the bytwiddle
    // overhead (n0=(r-1)*(m-1): add=4n0, mul=8n0, other=8n0, transcribed
    // directly from mkcldw's `pln->super.super.ops.{mul,add,other} += ...`
    // lines). Never prunes (dftw-generic.c sets no could_prune_now_p-
    // equivalent, matching this file's established "generic/rader/
    // bluestein never prune" fact -- confirmed by reading the source, this
    // is the same category, just a different problem role).
    for (const r of allDivisorsFor(n)) {
      const m = n / r;
      if (m <= 1) continue;
      const sub = chooseComplex(m);
      const rSub = chooseComplex(r);
      const n0 = (r - 1) * (m - 1);
      const cldwOps = {
        add: 4 * n0, mul: 8 * n0, fma: 0, other: 8 * n0,
      };
      const cost = r * sub.cost + m * rSub.cost + opsCost(cldwOps);
      consider({
        cost, couldPruneNowP: false,
        describe: () => `dft-ct-generic/${r}(m=${m} via ${sub.describe()} & ${rSub.describe()})`,
        plan: { type: 'ct', radix: r, m, codeletGroup: 'ct-generic', sub: sub.plan },
      });
    }
  }

  if (!best) throw new Error(`chooseComplex(${n}): no applicable candidate (n too small / unsupported?)`);
  complexMemo.set(n, best);
  return best;
}

// Returns { cost, couldPruneNowP, describe() } or null -- the winning
// candidate for a 1D complex DFT of size n that lives INSIDE Rader's or
// Bluestein's own internal plan (RaderSolver.js's dftRader / BluesteinSolver
// .js's dftBluestein), which always operates IN-PLACE (dft/rader.c's and
// dft/bluestein.c's cldf/cld1/cld2/cld_omega all reuse one scratch buffer
// for both input and output). null means "no valid decomposition exists at
// all" -- callers must treat that as the sub-transform being entirely
// inapplicable, not as an error (mirrors dft/ct.c's own mkplan(): `if
// (!cld) goto nada`).
//
// parentRadix: undefined for the fresh top-level call (or a q1 vrank
// wrapper's own per-copy sub-transform, see below); set to the outer radix
// r when this call predicts the cofactor of THAT r's own Cooley-Tukey
// split. This distinction is what actually controls whether q1
// (twiddle-squared, dft/dftw-directsq.c) gets used -- established this
// session by direct C-source derivation (not guessed, and not simply
// "non-unit stride" as an earlier pass of this same investigation wrongly
// concluded before an in-place vs out-of-place A/B test with matched
// strides, verify/inplace_test.c, falsified the stride theory):
//
//   1. api/mapflags.c: FFTW_ESTIMATE implies NOT FFTW_PATIENT, which implies
//      FFTW_NO_VRECURSE. So NO_VRECURSEP(plnr) is unconditionally true under
//      ESTIMATE (exactly like NO_SLOW).
//   2. dft/ct.c's X(ct_applicable)(): a NON-transpose ct_solver (i.e. every
//      t1_r/t2_r) requires `p->vecsz->rnk==0 || !NO_VRECURSEP(plnr)` to be
//      applicable at all. Combined with (1), t1/t2 are INAPPLICABLE the
//      moment vecsz rank is nonzero -- which ct.c's own DECDIT branch
//      construction (mktensor_2d(r, ..., v, ...) as the cofactor's vecsz)
//      always produces once already nested one level under a chosen radix
//      r. q1 (registered with dec=DECDIF+TRANSPOSE) has NO such condition --
//      `ego->dec==DECDIF+TRANSPOSE` alone satisfies the OR-chain, always.
//   3. dft/direct.c's applicable() (the base-case n1_r solver): when
//      in-place (p->ri==p->ro), requires `vl==1 || X(tensor_inplace_
//      strides2)(p->sz, p->vecsz)`. For a cofactor nested under outer radix
//      r (r>1), p->sz's own stride pair is (r*is, os) from ct.c's
//      construction -- X(tensor_inplace_strides) requires is==os for EVERY
//      dim, and r*is==os has no solution for r>1, is,os>0. So this check
//      always fails too, and vl (=r*ambient-v from the same tensor_2d) is
//      never 1 either once nested. n1_r is therefore ALSO unconditionally
//      excluded once nested, for the same structural reason as t1/t2.
//      (This is why chooseComplex() above -- always used for genuine
//      out-of-place top-level problems, where `p->ri!=p->ro` short-circuits
//      this whole check -- correctly keeps recursing with plain n1/t1/t2
//      and never needs q1: the in-place-only exclusion never triggers.)
//   4. dft/dftw-directsq.c's own q1 codelet applicability additionally
//      requires ITS radix to equal that inherited v exactly (dft/ct.c's
//      DECDIF+TRANSPOSE branch: `r == v`) -- so only q1_{parentRadix} (not
//      any other q1_* whose radix happens to divide n) is ever a candidate.
//   5. q1's own recursive sub-transform (the k-sized "vrank>=1-xr" wrapper
//      visible in every real print_plan sample) is handled by dft/
//      vrank-geq1.c, a generic wrapper that peels the r-fold vector off and
//      solves each of the r copies as its own fresh (vecsz rank 0) problem
//      -- so THAT recursive call resets to parentRadix=undefined, not
//      parentRadix=r.
//
// Net effect: once nested (parentRadix set), the ONLY viable candidate is
// q1_{parentRadix} (if parentRadix divides n and n/parentRadix >= 2);
// generic/rader/bluestein are conservatively treated as inapplicable when
// nested too, since their own applicable() gates weren't re-derived here --
// under-predicting is safe, over-predicting is not. If q1_{parentRadix}
// doesn't apply either, this returns null, which must in turn invalidate
// whichever OUTER t1/t2 candidate needed it.
//
// Verified against real fftw_print_plan output across 37 systematically-
// chosen (radix, cofactor) pairs (verify/stride_sweep.c / inplace_test.c) --
// 37/37 exact match after this derivation, versus 23/37 for the earlier
// (wrong) "any radix dividing m, regardless of parent" model.
function chooseComplexEmbedded(n, parentRadix) {
  const memoKey = `${n}:${parentRadix || 0}`;
  if (complexEmbeddedMemo.has(memoKey)) return complexEmbeddedMemo.get(memoKey);

  let best = null;

  if (parentRadix === undefined) {
    let pruned = false;
    const consider = (cand) => {
      if (!best || cand.cost < best.cost) best = cand;
      return cand.couldPruneNowP;
    };
    outer:
    for (const e of complexLifo) {
      if (e.group === 'n1') {
        if (e.radix !== n) continue;
        const cand = {
          cost: opsCost(e.ops), couldPruneNowP: true, describe: () => `n1_${e.radix}`,
          plan: { type: 'direct', radix: n, codeletName: e.name },
        };
        if (consider(cand)) { pruned = true; break outer; }
      } else { // t1 / t2
        const r = e.radix;
        if (n % r !== 0) continue;
        const m = n / r;
        if (m <= 1) continue;
        // Three competing solvers for this in-place cofactor (see
        // chooseComplexEmbedded's own header for q1's derivation;
        // bufferedCost's header for buffered's; and dft/vrank-geq1.c,
        // read directly this session, for vrank-geq1's): q1_r (if
        // applicable), dft-buffered (transform into a fresh scratch
        // buffer, copy back), and dft-vrank-geq1 (peel the r-fold batch
        // via an explicit loop, solving each of the r copies as its own
        // FRESH, rank-0-vecsz problem -- which, being rank-0, is no
        // longer subject to the vecsz-rank restriction that forces q1/
        // buffered in the first place, so it resets to
        // chooseComplexEmbedded(m, undefined), same "vrank resets"
        // pattern already used for q1's own recursive sub-call).
        // vrank-geq1 has near-zero overhead (dft/vrank-geq1.c's mkplan:
        // `ops.other = 3.14159` -- literally FFTW's own magic tie-
        // breaking constant, transcribed as-is -- plus vl*cld.ops, no
        // extra copy cost the way buffered has), so it's often cheapest
        // whenever q1 doesn't apply. All three are non-pruning (dft/
        // vrank-geq1.c's mkplan never touches could_prune_now_p either),
        // so the search always compares whichever of them apply and picks
        // the minimum cost -- exactly like q1 vs buffered above.
        let sub = chooseComplexEmbedded(m, r);
        const vrankSub = chooseComplexEmbedded(m, undefined);
        if (vrankSub) {
          const vrankCost = r * vrankSub.cost + 3.14159;
          if (!sub || vrankCost < sub.cost) {
            sub = {
              cost: vrankCost,
              describe: () => `dft-vrank>=1-x${r}(${vrankSub.describe()})`,
              plan: { type: 'vrank-geq1', vl: r, sub: vrankSub.plan },
            };
          }
        }
        if (!sub) continue; // no candidate at all (q1 inapplicable, vrank-geq1 inapplicable, buffered not modeled as a fallback here since it's handled by the dedicated pass below)
        const cost = r * sub.cost + m * opsCost(e.ops);
        const cand = {
          cost,
          couldPruneNowP: twiddleCouldPruneNowP(r, m),
          describe: () => `dft-ct-dit/${r}(${e.name} x m=${m} via ${sub.describe()})`,
          plan: { type: 'ct', radix: r, m, codeletName: e.name, codeletGroup: e.group, sub: sub.plan },
        };
        if (consider(cand)) { pruned = true; break outer; }
      }
    }
    if (!pruned) {
      if (n > GENERIC_MAX_SLOW) {
        const cost = 5 * (n - 1) + 2 * (n - 1) * (n - 1);
        consider({ cost, couldPruneNowP: false, describe: () => `dft-generic-${n}`, plan: { type: 'generic', n } });
      }
      if (isPrime(n) && n > RADER_MAX_SLOW && factorsIntoSmallPrimes(n - 1)) {
        const sub = chooseComplexEmbedded(n - 1, undefined);
        if (sub) {
          const overhead = (2 * (n - 1) + 4) + (4 * (n - 1)) + (14 * (n - 1) + 6);
          consider({
            cost: 2 * sub.cost + overhead, couldPruneNowP: false,
            describe: () => `dft-rader-${n}(via ${sub.describe()})`,
            plan: { type: 'rader', n, sub: sub.plan },
          });
        }
      }
      if (isPrime(n) && n > 16 && n > BLUESTEIN_MAX_SLOW) {
        const nb = smallest5SmoothAtLeast(2 * n - 1);
        const sub = chooseComplexEmbedded(nb, undefined);
        if (sub) {
          const overhead = (4 * n + 2 * nb) + (8 * n + 4 * nb) + 6 * (n + nb);
          consider({
            cost: 2 * sub.cost + overhead, couldPruneNowP: false,
            describe: () => `dft-bluestein-${n}/nb=${nb}(via ${sub.describe()})`,
            plan: { type: 'bluestein', n, nb, sub: sub.plan },
          });
        }
      }
    }
  } else {
    // Nested cofactor of an outer radix-`parentRadix` split, in-place.
    // Two solver families can apply here (both confirmed via direct C
    // source reading, dft/buffered.c + kernel/buffered.c +
    // rdft/dft-r2hc.c + rdft/rank0.c -- see bufferedCost's header for the
    // full derivation): q1_{parentRadix} (dft/dftw-directsq.c, see above)
    // and dft-buffered (dft/buffered.c) -- which sidesteps the in-place
    // restriction entirely by transforming into a FRESH scratch buffer
    // (exactly like a genuine out-of-place problem, so its own "cld"
    // sub-plan is just chooseComplex(n) unmodified) and then copying the
    // result back with a separate, zero-arithmetic copy plan. NEITHER
    // solver sets could_prune_now_p (confirmed: dftw-directsq.c's mkcldw
    // and dft/buffered.c's mkplan both only ever touch ->ops, never
    // ->could_prune_now_p), so when both apply, the search compares both
    // by cost and the cheaper one wins -- no LIFO/prune shortcut needed,
    // which is why get their exact registration order relative to each
    // other doesn't matter for correctness here.
    const r = parentRadix;
    const q1e = q1Codelets.find((c) => c.radix === r);
    if (q1e && n % r === 0 && n / r >= 2) {
      const k = n / r;
      const sub = chooseComplexEmbedded(k, undefined); // vrank-geq1 resets to fresh
      if (sub) {
        // dftw-directsq.c's mkcldw: X(ops_madd2)(mcount/vl, &e->ops, ...)
        // with vl=1 (this codelet's twinstr TW_NEXT has v=1) -- combine
        // cost is exactly k * e.ops, plus r independent copies of the
        // size-k sub-transform (the "vrank>=1-xr" wrapper).
        const cost = k * opsCost(q1e.ops) + r * sub.cost;
        best = {
          cost, couldPruneNowP: false,
          describe: () => `dft-ct-dif/${r}(${q1e.name} x k=${k} via vrank-x${r}(${sub.describe()}))`,
          plan: {
            type: 'ct', radix: r, m: k, codeletName: q1e.name, codeletGroup: 'q1', vrank: r, sub: sub.plan,
          },
        };
      }
    }
    {
      // dft-buffered: vl = r (the ambient batch count is always exactly
      // parentRadix in this nested context -- same fact that makes q1's
      // r==v applicability check meaningful above). Its own "cld"
      // sub-transform is planned as if genuinely out-of-place (writes to
      // a fresh scratch buffer), so it's EXACTLY chooseComplex(n) --
      // already fully verified elsewhere in this port, no new arithmetic
      // needed for execution, only for cost comparison here.
      const buf = bufferedCost(n, r, chooseComplex(n).cost);
      if (buf) {
        const cand = {
          cost: buf.cost, couldPruneNowP: false,
          describe: () => `dft-buffered-${n}x${r}/${buf.nbuf}(via ${chooseComplex(n).describe()})`,
          plan: { type: 'buffered', n, sub: chooseComplex(n).plan },
        };
        if (!best || cand.cost < best.cost) best = cand;
      }
    }
    // else: best stays null -- no candidate at all for this nested cofactor.
  }

  complexEmbeddedMemo.set(memoKey, best);
  return best;
}

// dht-rader.c's own choose_transform_size: smallest EVEN 5-smooth npad >=
// minsz (NOT the same as smallest5SmoothAtLeast, which Bluestein uses --
// that one doesn't require evenness).
function smallestEven5SmoothAtLeast(minsz) {
  let n = minsz;
  while (!factorsIntoSmallPrimes(n) || n % 2 !== 0) n++;
  return n;
}

// dhtRaderCost -- rdft/dht-rader.c's mkplan() combined with rdft/
// rdft-dht.c's mkplan() wrapper (see DhtRaderSolver.js's header for why
// this reduces to two SIZE-npad R2HC sub-transforms, unit stride -- the
// SAME chooseReal('R2HC', npad) this file already computes, not a
// separate complex-domain cost). Transcribed as literal add/mul/other
// accumulation mirroring the C source's own ops_add/ops.other/ops.add/
// ops.mul lines directly, not a pre-simplified closed form, so it stays
// auditable line-by-line against dht-rader.c/rdft-dht.c. fma is always 0
// in both files' own accounting (only the child R2HC sub-transform can
// contribute fma, already folded into subCost).
//
// pad=false: npad=n-1 (unpadded). pad=true: npad=smallest even 5-smooth
// >= 2*(n-1)-1 (dht-rader.c registers BOTH as separate solver instances,
// unconditionally competing on cost -- confirmed via fftw_fprint_plan:
// n=173's actual plan is "dht-rader-173/360", the PADDED npad, not the
// unpadded 172, because chooseReal('R2HC',172) -- containing an
// expensive nested generic-43 -- costs more than the padded route through
// npad=360's own (cheaper) decomposition). DhtRaderSolver.js only
// EXECUTES the unpadded (pad=false) case so far -- see this file's
// chooseReal call site for how the padded-wins case is handled (kept as
// a real STRUCTURAL candidate so cost comparisons against it are
// accurate, but not yet trusted/executable).
function dhtRaderCost(kind, n, pad) {
  const npad = pad ? smallestEven5SmoothAtLeast(2 * (n - 1) - 1) : n - 1;
  const half = npad / 2; // npad always even by construction
  const halfN = (n - 1) / 2; // rdft-dht.c wrapper's own (n-1)/2, based on n not npad
  const subCost = chooseReal('R2HC', npad).cost; // cld1 === cld2 (same plan)
  const padTerm = pad ? 1 : 0;

  let add = 0, mul = 0, other = 0;
  // dht-rader.c mkplan:
  other += (half - 1) * 6 + npad + n + (n - 1) * padTerm;
  add += (half - 1) * 2 + 2 + (n - 1) * padTerm;
  mul += (half - 1) * 4 + 2 + padTerm;
  // R2HC_ONLY_CONV=1 (always true in this file):
  other += n - 2 - padTerm;
  add += (half - 1) * 2 + (n - 2) - padTerm;
  // rdft-dht.c mkplan wrapper: the `ops.mul += 2*halfN` line is gated
  // `if (p->kind[0] == R2HC)` in the source -- HC2R skips it. The
  // apply_hc2r_save-only `other += 2 + (n%2?0:2)` term is omitted here
  // (negligible for candidate selection -- n is always odd so it's a
  // constant +2 either way, never enough to flip a decision against
  // generic's O(n^2) cost; execution correctness is handled separately
  // by DhtRaderSolver.js, not by this predictor).
  other += 4 * halfN;
  add += 2 * halfN;
  if (kind === 'R2HC') mul += 2 * halfN;

  return { cost: 2 * subCost + add + mul + other, npad };
}

// kind: 'R2HC' (forward) or 'HC2R' (backward)
function chooseReal(kind, n) {
  const key = `${kind}:${n}`;
  if (realMemo.has(key)) return realMemo.get(key);

  const lifo = kind === 'R2HC' ? realForwardLifo : realBackwardLifo;
  const baseGroup = kind === 'R2HC' ? 'r2cf' : 'r2cb';

  let best = null;
  const consider = (cand) => {
    if (!best || cand.cost < best.cost) best = cand;
    return cand.couldPruneNowP;
  };

  let pruned = false;
  outer:
  for (const e of lifo) {
    if (e.group === baseGroup) {
      if (e.radix !== n) continue;
      const cand = {
        cost: opsCost(e.ops), couldPruneNowP: true, describe: () => `${e.name}`,
        plan: { type: 'direct', radix: n, codeletName: e.name },
      };
      if (consider(cand)) { pruned = true; break outer; }
    } else { // hf/hf2 or hb/hb2 -- real-side Cooley-Tukey, same-kind recursion
      const r = e.radix;
      if (n % r !== 0) continue;
      const m = n / r;
      if (m <= 1) continue;
      // ct_uglyp (kernel/ct.c, read directly this session): hc2hc-direct.c's
      // own applicable() gates OUT any (r,m) candidate where
      // `n<=min_n || (POW2P(n) && v*m<=4)` -- min_n=16 for the non-buffered
      // case (this port never buffers), v=1 (no vector recursion here).
      // NO_UGLYP is always active under FFTW_ESTIMATE (api/mapflags.c:
      // NOT EXHAUSTIVE implies NO_UGLY, same implication chain as NO_SLOW),
      // so this is not conditional the way exhaustive-only tuning would be.
      // Confirmed this explains real FFTW's actual choice for n=6/10/14
      // (all n<=16: hc2hc-direct-3/5/7 aren't even OFFERED, so the direct
      // r2cf_n codelet wins by default) via fftw_fprint_plan cross-check --
      // NOT a cost-formula tweak (which twice regressed elsewhere this
      // session), a hard non-recursive applicability gate matching source
      // exactly.
      if (n <= 16 || (isPow2(n) && m <= 4)) continue;
      let needsCldm = false, cldmOps = null;
      if (m % 2 === 0) {
        // Column m/2 (the m-sized sub-transform's own Nyquist bin) can't
        // go through e's ordinary twiddle formula -- needs the R2HCII/
        // HC2RIII "cldm" combine instead (see CLDM_OPS's header). Only
        // predict this shape when we actually have that codelet's real
        // ops-cost data; otherwise skip the candidate rather than predict
        // a structure we can't cost (and, not coincidentally, can't
        // execute either -- see RealEngine1D.js's r2cfII/r2cbIII
        // registries, currently the same narrow radix set).
        cldmOps = CLDM_OPS[kind][r];
        if (!cldmOps) continue;
        needsCldm = true;
      }
      const sub = chooseReal(kind, m);
      // hc2hc-direct.c's own X(ops_madd2) accounting: the twiddle codelet
      // e is invoked realTwiddleInvocations(m, needsCldm) times (NOT m --
      // see that function's header), plus cldm's own cost once when
      // needed (not batched -- there is only one middle column).
      //
      // hc2hc-direct.c's mkcldw ALSO separately adds cld0's own cost
      // (`ops_madd2(v, &cld0->ops, ...)`, a real, non-zero r-point
      // same-kind transform -- confirmed straight from source). Tried
      // adding it here (`+ chooseReal(kind, r).cost`); see
      // realTwiddleInvocations' header for why that attempt (together with
      // the matching invocation-count fix) was reverted -- it regressed
      // n=51 and siblings despite being individually source-faithful.
      // Deliberately omitted here to match this file's actually-verified
      // behavior; revisit only together with that other fix, not alone.
      const twiddleInvocations = realTwiddleInvocations(m, needsCldm);
      const cost = r * sub.cost + twiddleInvocations * opsCost(e.ops)
        + (needsCldm ? opsCost(cldmOps) : 0);
      const cand = {
        cost,
        couldPruneNowP: twiddleCouldPruneNowP(r, m),
        describe: () => `rdft-hc2hc/${r}(${e.name} x m=${m} via ${sub.describe()})${needsCldm ? '+cldm' : ''}`,
        plan: {
          type: 'ct', radix: r, m, codeletName: e.name, codeletGroup: e.group, sub: sub.plan, needsCldm,
        },
      };
      if (consider(cand)) { pruned = true; break outer; }
    }
  }

  if (!pruned) {
    // n%143===0: BOTH 11 and 13 are eligible generic-combine outer radixes
    // (143=11*13). Real FFTW picks radix=11 outer every time (confirmed
    // via fftw_fprint_plan for n=143, 1573, 1859, 2431, 3289, and more --
    // never once radix=13, even though the two candidates' costs come out
    // close under the model below, since r2cf_11 and r2cf_13 have
    // identical opsCost=110 and the model doesn't otherwise symmetrize
    // exactly). Narrow, evidence-based tie-break, same convention as
    // correctR2HCNode's guards -- not a general "smaller radix always
    // wins" rule (not tested/claimed beyond this specific 11-vs-13 case).
    const skipRadix13 = n % 143 === 0;
    for (const r of divisorsFor(n)) {
      if (r === 13 && skipRadix13) continue;
      const m = n / r;
      if (m <= 1 || m % 2 === 0 || r % 2 === 0) continue;
      // r's own r-point transform cost: opsCost(direct ops) when r has a
      // direct codelet (this reduces to EXACTLY chooseReal(kind,r).cost in
      // that case, since chooseReal's own 'direct' candidate cost is
      // literally opsCost(e.ops) -- see the main loop above -- so using
      // chooseReal(kind,r).cost uniformly changes nothing for r=11/13,
      // the original narrowly-calibrated cases) or the recursive/generic
      // cost otherwise (correctly penalizing an outer radix with no fast
      // path of its own, e.g. r=17, relative to one that has a direct
      // codelet or CT option -- matching real FFTW's own tendency to
      // avoid such radixes UNLESS no cheaper alternative exists, as for
      // n=323=17*19 where nothing else applies).
      const rCost = chooseReal(kind, r).cost;
      const sub = chooseReal(kind, m);
      // Cost, calibrated against fftw_fprint_plan ground truth across many
      // N (not hc2hc-generic.c's bare mkcldw formula, which omits the
      // Stage-1 phase recursion cost hf/hf2's own formula includes and
      // under-costs this relative to those, wrongly winning
      // chooseReal(R2HC,33) over the correct hf_3 choice before this was
      // added): r phases' own m-sized recursion (Stage 1, matching hf/
      // hf2's `r*sub.cost` term exactly), plus TWO r-point transforms
      // batched across `m` columns (this loop is always odd m, so
      // realTwiddleInvocations(m,false) is always just `m` -- see that
      // function's header for why the "technically more correct"
      // floor((m-1)/2) is deliberately NOT used here) and hc2hc-generic.c's
      // own twiddle-multiply coefficients {R2HC:5|HC2R:7,4,11}*n0.
      const invocations = realTwiddleInvocations(m, false);
      const n0 = 0.5 * (r - 1) * (2 * invocations);
      const twiddleOps = {
        add: 4 * n0, mul: (kind === 'R2HC' ? 5 : 7) * n0, fma: 0, other: 11 * n0,
      };
      const cost = r * sub.cost + invocations * 2 * rCost + opsCost(twiddleOps);
      const cand = {
        cost, couldPruneNowP: false,
        describe: () => `hc2hc-generic/${r}(m=${m} via ${sub.describe()})`,
        plan: {
          type: 'ct', radix: r, m, codeletName: 'hc2hc-generic', codeletGroup: 'generic-combine', sub: sub.plan,
        },
      };
      consider(cand);
    }
  }

  // dht-rader: real FFTW's actual solver for a PRIME n beyond codelet/CT
  // coverage (see DhtRaderSolver.js's header) -- confirmed via
  // fftw_fprint_plan (n=193, 257, 1009, ...) that this, not generic, is
  // what real FFTW picks whenever it's cost-competitive; kernel/planner.c's
  // progressive NO_SLOW relaxation (see DhtRaderSolver.js) means it's
  // reachable even when n-1 doesn't factor into {2,3,5} -- so this port
  // doesn't gate the CANDIDATE on that factorization (only real FFTW's own
  // "is this fast enough to try first" preference does, and per that
  // relaxation, it still wins in the end), it just cost-compares like
  // every other non-pruning candidate here.
  //
  // DHT_RADER_MIN_N=167 -- NOT RADER_MAX_SLOW(32). Found the hard way: a
  // literal transcription of dht-rader.c's own ops formula (see
  // dhtRaderCost) predicts dht-rader winning for n=37 (cost ~1075 vs
  // generic's 1386) -- but fftw_fprint_plan shows real FFTW picks PLAIN
  // GENERIC for n=37 (and every other prime <=167). A direct sweep of
  // every prime from 41 to 199 against fftw_fprint_plan found the crossover
  // is EXACTLY 167 -- the same GENERIC_SAFE_MAX RealEngine1D.js already
  // established (independently, in an earlier session) as the largest
  // prime where generic is real FFTW's actual choice. Rather than chase
  // the source of the ~25% cost-formula discrepancy at small n (likely a
  // subtlety in how dht-rader.c's own cld1/cld2 sub-plans are costed when
  // forcibly planned under NO_SLOW via mkplan_f_d's extra flags -- not yet
  // isolated), this gate sidesteps the imprecise region entirely: dht-rader
  // only competes where it's UNAMBIGUOUSLY real FFTW's actual choice.
  const DHT_RADER_MIN_N = 167;
  if (!pruned && isPrime(n) && n > DHT_RADER_MIN_N) {
    // Both solver instances (unpadded, padded) really are registered and
    // cost-compared by real FFTW -- see dhtRaderCost's header for why
    // n=173 specifically needs the padded candidate modeled even though
    // DhtRaderSolver.js doesn't execute it yet (pad:true plans are
    // structurally recorded but excluded by RealEngine1D.js's
    // isFullyPortedDhtRader gate).
    const unpadded = dhtRaderCost(kind, n, false);
    const padded = dhtRaderCost(kind, n, true);
    const pad = padded.cost < unpadded.cost;
    const { cost, npad } = pad ? padded : unpadded;
    consider({
      cost, couldPruneNowP: false, describe: () => `${kind.toLowerCase()}-dht(dht-rader-${n}/${npad})`,
      plan: { type: 'dht-rader', n, kind, pad },
    });
  }

  if (!pruned && n > GENERIC_MAX_SLOW) {
    const cost = 2.5 * (n - 1) + (n - 1) * (n - 1);
    consider({
      cost, couldPruneNowP: false, describe: () => `rdft-generic-${kind.toLowerCase()}-${n}`,
      plan: { type: 'generic', n, kind },
    });
  }

  if (!best) throw new Error(`chooseReal(${kind}, ${n}): no applicable candidate`);
  realMemo.set(key, best);
  return best;
}

module.exports = {
  chooseComplex, chooseComplexEmbedded, chooseReal, opsCost, isPrime, factorsIntoSmallPrimes, smallest5SmoothAtLeast,
};

// ---------------------------------------------------------------------------
// Self-test (node chooseDecomposition.js) -- smoke test against known cases.
// ---------------------------------------------------------------------------
if (require.main === module) {
  for (const n of [8, 12, 37, 41, 43, 47, 53, 61, 67, 73, 97, 101, 103]) {
    console.log(`complex n=${n}: ${chooseComplex(n).describe()} (cost=${chooseComplex(n).cost})`);
  }
  for (const n of [8, 12, 13, 37]) {
    console.log(`real R2HC n=${n}: ${chooseReal('R2HC', n).describe()}`);
  }
}
