'use strict';

// =============================================================================
// Composite1D.js
// General recursive 1D complex DFT engine: executes whatever decomposition
// Planner/chooseDecomposition.js's chooseComplex(n) predicts FFTW3 would
// choose, using the ported codelets in Codelets/complex/. This generalizes
// CompositeSolver1D.js's hardcoded radix-3-only driver and RaderSolver.js's
// hardcoded radix-6x6 split, WITHOUT replacing either -- both stay as
// verified, standalone back-compat entry points; this file is the new,
// general path.
//
// Byte-exact fidelity against real FFTW3 depends on using the EXACT SAME
// decomposition FFTW's planner would pick, not just any mathematically
// valid one -- different radix choices are all "correct" to double
// precision but round differently, which is visible at the bit level. So
// this file does NOT implement its own factor-picking heuristic; it always
// asks chooseComplex(n) and executes its .plan verbatim.
//
// FALLBACK SAFETY: if chooseComplex(n) picks a plan node this file can't
// execute yet (an unported codelet radix, or a rader/bluestein sub-call
// whose applicability this file doesn't implement), execution falls back to
// GenericSolver1D.dftGeneric for that subtree -- ONLY WHEN n IS ODD.
// dftGeneric is a faithful port of dft/generic.c, whose `for (i=1; i+i<n;
// ++i)` loop never writes the self-paired middle bin for EVEN n (confirmed
// against the real FFTW3 source: same bound, same gap). Real FFTW3 never
// notices because it never actually selects dft-generic for even N -- a
// radix-2-based candidate always wins the cost comparison first -- but that
// means dftGeneric is NOT a universally-correct fallback the way it is for
// odd n; calling it for even n silently returns wrong (zero) data at the
// middle bin. So the fallback here throws for even n instead of guessing
// wrong -- "false is always safe" per PlanTable.js's philosophy extends to
// "loud failure beats a wrong answer that looks like a value." Callers
// that need to know whether the EXACT FFTW3 path was used (for production
// gating) should use isFullyPortedComplex() below, not a try/catch.
// =============================================================================

const { chooseComplex, chooseComplexEmbedded } = require('./Planner/chooseDecomposition');
const { dftGeneric } = require('./GenericSolver1D');
const { noTwiddle, twiddle, twiddleT2 } = require('./Codelets/complex/index');
const { realCexp } = require('./Trig');
const { Q1_RADIXES } = require('./Codelets/complex/q1');

const Q1_RADIX_SET = new Set(Q1_RADIXES);

function isPrime(n) {
  if (n <= 1) return false;
  if (n % 2 === 0) return n === 2;
  for (let i = 3; i * i <= n; i += 2) if (n % i === 0) return false;
  return true;
}

// See isFullyPortedComplex's 'generic' leaf comment for the evidence.
const GENERIC_COMPLEX_SAFE_MAX = 167;

// True iff chooseComplexEmbedded(n, parentRadix)'s predicted plan is fully
// executable -- the trust-check counterpart to isFullyPortedComplex, but for
// the IN-PLACE, possibly-q1-using decomposition RaderSolver.js/
// BluesteinSolver.js's own internal sub-transform needs (see
// chooseComplexEmbedded's header in Planner/chooseDecomposition.js and
// Codelets/complex/q1.js's header for why this is a genuinely different
// prediction than isFullyPortedComplex's out-of-place one).
//
// q1 codelet coverage: radixes 2, 3, 5, 8 were independently verified
// bit-exact end-to-end against real fftw_execute() for real embedded cases
// (n=24, n=27, n=125, n=128 respectively, all in-place). Radixes 4 and 6
// were NOT independently numerically verified (no simple, non-bluestein-
// buried real N was found to exercise them standalone) -- they're included
// here on the strength of the row-DFT-reuse pattern being verified
// byte-identical to noTwiddle[r] for r=3 AND r=4 via direct C source
// comparison, and the wrapper logic (transpose+twiddle) being completely
// generic (zero radix-specific code) across all six radixes already proven
// correct. This is exactly the kind of gap the mandatory broad safety sweep
// (see verify/, this session's established methodology) is supposed to
// catch -- if r=4 or r=6 turn out to be wrong, the sweep will surface it as
// a mismatch, at which point EXCLUDE those two radixes from Q1_RADIX_SET
// (or from Codelets/complex/q1.js's Q1_RADIXES directly) rather than
// trying to patch this comment's confidence claim.
function isFullyPortedComplexEmbedded(n, parentRadix) {
  const candidate = chooseComplexEmbedded(n, parentRadix);
  if (!candidate) return false;

  function walk(plan) {
    if (plan.type === 'direct') return !!noTwiddle[plan.radix];
    if (plan.type === 'ct') {
      // Same radix-25/64 outer/cofactor swap mispredict as
      // isFullyPortedComplex's walk() below (see that check's comment for
      // the full derivation) -- found here too via n=113's bluestein
      // (nb=225): chooseComplexEmbedded(225, undefined) picks
      // dft-ct-dit/25(t2_25, n1_9), but real FFTW's in-place n=225 plan is
      // dft-ct-dit/5(t2_5, dft-ct-dif/5(q1_5, vrank-x5(n1_9))) -- same
      // shape of cost-model divergence (radix 25 "eating" what should be
      // two nested radix-5 layers), just now surfacing in the embedded/
      // q1-aware predictor instead of the plain out-of-place one.
      if ((plan.radix === 25 || plan.radix === 64) && (twiddle[plan.m] || twiddleT2[plan.m])) return false;
      if (plan.codeletGroup === 'q1') {
        return Q1_RADIX_SET.has(plan.radix) && !!noTwiddle[plan.radix] && walk(plan.sub);
      }
      if (plan.codeletGroup === 't1') return !!twiddle[plan.radix] && walk(plan.sub);
      if (plan.codeletGroup === 't2') return !!twiddleT2[plan.radix] && walk(plan.sub);
      return false;
    }
    // dft-buffered (dft/buffered.c, see Planner/chooseDecomposition.js's
    // bufferedCost header): the sub-transform it wraps is planned as if
    // genuinely out-of-place, so plan.sub here is exactly chooseComplex
    // (plan.n)'s own plan -- trust it via isFullyPortedComplex, the SAME
    // check already used everywhere else in this port for that plan
    // shape, not a new predicate.
    if (plan.type === 'buffered') return isFullyPortedComplex(plan.n);
    if (plan.type === 'generic') return isPrime(plan.n) && plan.n <= GENERIC_COMPLEX_SAFE_MAX;
    if (plan.type === 'rader') return require('./RaderSolver.js').isFullyPortedRader(plan.n);
    if (plan.type === 'bluestein') return require('./BluesteinSolver.js').isFullyPortedBluestein(plan.n);
    return false;
  }

  return walk(candidate.plan);
}

// True iff every node of chooseComplex(n)'s chosen plan is executable by
// this file's codelet registry (no unported codelet, no rader/bluestein).
// Use this -- not a try/catch around dftComposite -- to decide whether a
// size is safe to promote in PlanTable.js, since a caught fallback still
// produces a CORRECT answer, just not necessarily FFTW3's exact bit
// pattern; PlanTable's production gate needs the latter.
function isFullyPortedComplex(n) {
  // n=4096=64^2 mispredict -- found empirically, not derived: real
  // fftw_print_plan for n=4096 is dft-ct-dit/64(t2_64, ..., n1_64), but
  // chooseComplex(4096) doesn't even consider radix 64 as the top-level
  // choice -- it picks dft-ct-dit/32(t2_32, ..., dft-ct-dit/8(t1_8,
  // n1_16)) instead. Unlike the radix-25/64 swap check below (which
  // catches chooseComplex picking the WRONG SIDE of a radix it did
  // consider), this is chooseComplex's cost search never reaching radix
  // 64 as a candidate at all for this specific N, so no node-level check
  // on its own predicted tree can catch it. Moot for execution anyway --
  // the real decomposition needs n1_64, which isn't ported -- but the
  // exclusion still matters so this doesn't get falsely marked verified.
  // Scoped to this single confirmed value; a broader 4..5000 sweep found
  // no other case of this shape.
  if (n === 4096) return false;
  function walk(plan) {
    if (plan.type === 'direct') return !!noTwiddle[plan.radix];
    // Only the plain 't1' family is ported here -- 't2' (radix 4,5,8,...;
    // an alternate-codegen twiddle-batching strategy, same math, different
    // instruction scheduling) would be silently substituted with t1_r by
    // executePlan, which is mathematically correct but NOT FFTW3's exact
    // bit path (same class of gap as RealEngine1D.js's hf-vs-hf2 check --
    // found by that same audit, applied here too).
    //
    // m===radix*radix (a perfect square of the SAME radix used one level
    // up) is excluded too -- found empirically, not derived: n=27's real
    // fftw_print_plan is dft-ct-dit/3(t1_3, dft-ct-dif/3(q1_3, n1_3)), i.e.
    // real FFTW prefers the twiddle-SQUARED q1_3 codelet (unported, see
    // README.md's Phase 6) for the m=9 middle layer specifically because
    // 9=3^2 shares the outer radix -- even though n=9 in TOP-LEVEL isolation
    // correctly uses the plain n1_9 direct codelet chooseComplex(9) already
    // predicts. chooseComplex is a size-only predictor and has no way to
    // know a sub-transform is embedded under a same-radix parent, so this
    // check is a deliberately blunt, conservative stand-in: exclude ANY
    // node whose sub-size m is a perfect square, rather than exactly
    // replicating FFTW's real applicability condition for q1.
    if (plan.type === 'ct') {
      if (plan.m === plan.radix * plan.radix) return false;
      // radix=25 outer/cofactor SWAP mispredict -- found empirically, not
      // derived: n=250's real fftw_print_plan is dft-ct-dit/10(t2_10,
      // ..., n1_25), but chooseComplex(250) predicts
      // dft-ct-dit/25(t2_25, ..., n1_10) -- same math, wrong bit pattern.
      // The naive recursive cost model here doesn't replicate whatever
      // makes real FFTW prefer swapping outer/cofactor roles for this
      // specific shape. Confirmed via a 16-point sweep of real
      // fftw_print_plan output across n=125..2000 (all multiples of
      // 125=5^3) that this ONLY happens when the outer radix is 25 AND
      // the cofactor m is ITSELF a registered twiddle radix (n=250: m=10,
      // n=500: m=20 both mispredicted; n=350/750/800/1000/... , where m
      // isn't a registered twiddle radix, verified bit-exact correct as
      // chooseComplex predicts). The SAME shape recurs at radix=64: real
      // fftw_print_plan for n=2048 is dft-ct-dit/32(t2_32, ..., n1_64),
      // but chooseComplex(2048) predicts dft-ct-dit/64(t2_64, ..., n1_32)
      // -- confirmed via fftw_fprint_plan. Scoped to radixes 25 and 64 --
      // the only two shapes with real ground-truth evidence -- rather
      // than assumed to generalize to every t2 radix without separate
      // verification.
      if ((plan.radix === 25 || plan.radix === 64) && (twiddle[plan.m] || twiddleT2[plan.m])) return false;
      // A NESTED (non-top) ct sub-layer whose own radix equals its own
      // cofactor (radix===m, a "self-squared" layer -- e.g. n=625's outer
      // radix=25 layer sitting over an inner radix=5,m=5 layer) is unsafe
      // for the same reason as the m===radix*radix check above: real FFTW
      // diverges from the naive recursive cost estimate specifically when
      // a same-radix-squared shape is EMBEDDED under another CT layer.
      // Confirmed empirically via ground truth: n=625 (outer radix=25,m=25
      // over inner radix=5,m=5), n=1280 (outer radix=20,m=64 over inner
      // radix=8,m=8), and n=1600 (outer radix=25,m=64 over the same inner
      // radix=8,m=8) all mispredict -- while the SAME radix=8,m=8 shape
      // used as N=64's OWN top-level (single-layer, not nested)
      // decomposition verifies bit-exact. So this only applies to an
      // embedded sub-layer, not a standalone top-level one (which the
      // m===radix*radix check and this check both leave alone).
      //
      // Scoped to codeletGroup t1/t2 ONLY -- does NOT apply when the
      // nested self-squared layer is 'ct-generic' (dft/dftw-generic.c has
      // no hand-tuned twiddle codelet to exhibit the SAME real-FFTW
      // cost-model quirk this restriction guards against). Directly
      // re-verified via dft1d_harness for every self-squared-ct-generic
      // nested case up to n=3000 (51 total, e.g. n=242=2*11^2, n=1331=
      // 11^3) -- all bit-exact.
      if (plan.sub.type === 'ct' && plan.sub.radix === plan.sub.m
        && (plan.sub.codeletGroup === 't1' || plan.sub.codeletGroup === 't2')) return false;
      if (plan.codeletGroup === 't1') return !!twiddle[plan.radix] && walk(plan.sub);
      if (plan.codeletGroup === 't2') return !!twiddleT2[plan.radix] && walk(plan.sub);
      // 'ct-generic' (dft/dftw-generic.c) -- radix r does NOT need its own
      // direct n1_r codelet or twiddle codelet; ctGeneric() uses r's own
      // FULL chooseComplex(r) engine internally (see that function's
      // header), so trust it exactly whenever r's own complex DFT is
      // itself already trusted -- same pattern as RealEngine1D.js's
      // generic-combine guard on the real side.
      if (plan.codeletGroup === 'ct-generic') return isFullyPortedComplex(plan.radix) && walk(plan.sub);
      return false;
    }
    // 'rader': delegate to RaderSolver.js's own applicability+sub-transform
    // check (lazy require -- RaderSolver.js requires this file for its
    // dftComposite(n-1,...) sub-transform, so a top-level require here would
    // be circular; by the time this FUNCTION actually runs, both modules
    // have finished loading, so a require() here just returns the cached,
    // fully-populated module).
    //
    // plan.n, NOT the outer n -- walk() is a closure over isFullyPortedComplex's
    // OWN top-level n parameter, fixed for the entire recursive walk, so for
    // a rader/bluestein plan NESTED under a ct-generic (or t1/t2) layer,
    // that outer n is a completely different, larger size than this plan
    // node's own. Bug found this session while investigating why some
    // embedded-rader structures were mispredicted -- this check was
    // silently validating the WRONG size the whole time. Every case found
    // so far happened to still resolve to the correct (false) answer by
    // coincidence (the wrong size also failed its own gate), which is how
    // this went undetected through extensive prior sweeps; do not assume
    // that coincidence holds for every possible n without re-sweeping.
    if (plan.type === 'rader') {
      return require('./RaderSolver.js').isFullyPortedRader(plan.n);
    }
    // 'bluestein': delegate to BluesteinSolver.js's own applicability+
    // sub-transform check (lazy require -- BluesteinSolver.js requires this
    // file for dftComposite(nb,...), same circular-require situation as
    // RaderSolver.js above). plan.n, not outer n -- see 'rader' comment above.
    if (plan.type === 'bluestein') {
      return require('./BluesteinSolver.js').isFullyPortedBluestein(plan.n);
    }
    // 'generic' leaf: safe IFF the generic leaf's OWN size (plan.n --
    // which may be a NESTED sub-transform's size, not the top-level n --
    // e.g. n=51=3*17 is "dft-ct-dit/3(t1_3 x m=17 via dft-generic-17)",
    // generic leaf size 17) is PRIME. A prime has no factors, so real
    // FFTW's own planner ALSO has no CT candidate for it and genuinely
    // falls to generic (or Rader/Bluestein beyond their own thresholds,
    // a separate case handled above) -- dftGeneric matches bit-for-bit
    // there. A COMPOSITE generic leaf (e.g. bare n=121=11^2, or n=143=
    // 11*13) means chooseComplex couldn't find ANY registered-codelet CT
    // radix for that composite and gave up to generic -- but real FFTW
    // may still have an UNPORTED twiddle codelet (e.g. t1_11, which this
    // port's Codelets/complex/ registry doesn't have) that lets it build
    // a genuine CT structure there instead, which does NOT numerically
    // match plain generic.
    //
    // Rigorously reverified this session after an EARLIER attempt at this
    // exact fix turned out to be built on a broken test harness that
    // silently vacuously passed every comparison (a parsing bug produced
    // NaN arrays; a max-tracking loop that only updates on `d > m` never
    // updates when d is NaN, so it always "found" maxDiff=0 no matter
    // what). Once the harness was fixed: a direct dft1d_harness sweep of
    // every odd n up to 3000 whose chooseComplex(n) plan contains a
    // generic leaf anywhere (573 total) found EXACTLY the prime/composite
    // split predicts -- all 185 prime-generic-leaf cases bit-exact, all
    // 388 composite-generic-leaf cases wrong (by ~1e-13 to ~1e-14, not
    // rounding noise). GENERIC_COMPLEX_SAFE_MAX mirrors RealEngine1D.js's
    // GENERIC_SAFE_MAX=167 (O(n^2) becomes a real DCP-worker-budget
    // concern well before correctness is ever in question at large N;
    // 167 is a conservative, evidence-backed starting point -- this
    // session's sweep covered generic leaves up to 3000, comfortably
    // beyond it). Any future attempt to touch this MUST first prove the
    // verification harness actually detects a known-bad case (e.g.
    // n=121, generic leaf 121) before trusting a broad "all clear" sweep
    // result -- that check would have caught the earlier broken harness
    // immediately.
    if (plan.type === 'generic') return isPrime(plan.n) && plan.n <= GENERIC_COMPLEX_SAFE_MAX;
    return false;
  }
  return walk(chooseComplex(n).plan);
}

// True iff chooseComplex(n)'s plan is 'direct', or a SINGLE 'ct' layer over
// a 'direct' sub-transform (radix r, base-case m) -- i.e. no nested CT
// recursion. This matters ONLY for sub-transforms embedded inside Rader's/
// Bluestein's OWN internal plan (RaderSolver.js/BluesteinSolver.js), which
// -- confirmed empirically via verify/stride_test.c -- always builds its
// cldf/cld_omega sub-problem with a NON-UNIT stride (dft/rader.c's and
// dft/bluestein.c's mktensor_1d(..., 2, 2)-shaped tensors). Real FFTW's
// planner can choose a DIFFERENT decomposition for a non-unit-stride
// problem than for the unit-stride one chooseComplex(n) predicts -- found
// by comparing fftw_print_plan for a stride-1 vs stride-2 size-225 complex
// DFT: stride-1 uses dft-ct-dit/5(t2_5, dft-ct-dit/5(t2_5, n1_9)) (matching
// chooseComplex(225) exactly), but stride-2 uses dft-ct-dit/5(t2_5,
// dft-ct-dif/5(q1_5, ...)) for the SAME logical decomposition -- the nested
// (depth>=2) inner CT layer gets replaced with the unported q1 family.
// n=36/40/60 (single-level: radix-r CT over a DIRECT base case) verified
// bit-exact under both strides; n=225 (two-level nesting) did not. Until
// FFTW's planner is replicated stride-aware (a real, separate research
// task -- see README.md's Phase 6 scope), isFullyPortedRader/
// isFullyPortedBluestein require this single-level check IN ADDITION TO
// isFullyPortedComplex, rather than trusting the stride-1 predictor alone.
function isSingleLevelOrDirect(n) {
  const plan = chooseComplex(n).plan;
  if (plan.type === 'direct') return true;
  return plan.type === 'ct' && plan.sub.type === 'direct';
}

// dftGeneric is only a valid stand-in for a missing codelet when n is odd
// (see this file's header) -- centralize that guard so every fallback site
// fails loudly on even n instead of silently returning zeroed bins.
function genericFallback(n, reIn, imIn) {
  if (n % 2 === 0) {
    throw new Error(
      `Composite1D: no ported codelet for n=${n} (even) and dftGeneric is `
      + `not valid for even n (misses the self-paired middle bin -- see `
      + `this file's header). Port the missing codelet, or route this size `
      + `through a different solver.`
    );
  }
  const ro = new Float64Array(n), io = new Float64Array(n);
  dftGeneric(n, reIn, 0, imIn, 0, 1, ro, 0, io, 0, 1);
  return [ro, io];
}

// dft-ct-generic's "cldw" role (dft/dftw-generic.c's bytwiddle + batched
// r-point DFT, see Planner/chooseDecomposition.js's cost-model comment for
// the derivation) -- a drop-in replacement for the twiddle[r]/twiddleT2[r]
// codelet call in Stage 2 below, same (br,bi,Wc,Ws)->[outR,outI] signature.
// br[0]/bi[0] and the Wc[0]/Ws[0] "twiddle" (always (1,0), untouched by the
// caller) need no multiply; p=1..r-1 get the ordinary twiddle-multiply this
// file already uses everywhere (FFT_SIGN=-1: multiply by (Wc[p]-i*Ws[p])),
// then the twiddled r values go through r's own FULL chooseComplex(r)
// engine -- not a hardcoded codelet, exactly like RealEngine1D.js's
// genericCombineR2HC/HC2R trusts r's own engine on the real side.
function ctGeneric(r, br, bi, Wc, Ws) {
  const tr = new Float64Array(r), ti = new Float64Array(r);
  tr[0] = br[0]; ti[0] = bi[0];
  for (let p = 1; p < r; p++) {
    const xr = br[p], xi = bi[p];
    tr[p] = xr * Wc[p] + xi * Ws[p];
    ti[p] = xi * Wc[p] - xr * Ws[p];
  }
  return executePlan(chooseComplex(r).plan, r, tr, ti);
}

// Execute a chooseComplex()-shaped plan against (reIn, imIn) (plain arrays
// or typed arrays, length n). Returns [reOut, imOut] as Float64Array(n).
function executePlan(plan, n, reIn, imIn) {
  if (plan.type === 'direct') {
    const fn = noTwiddle[plan.radix];
    if (fn) return fn(reIn, imIn);
    return genericFallback(n, reIn, imIn);
  }

  if (plan.type === 'ct') {
    const { radix: r, m, sub } = plan;
    const fn = plan.codeletGroup === 'ct-generic' ? ((br, bi, Wc, Ws) => ctGeneric(r, br, bi, Wc, Ws))
      : plan.codeletGroup === 't2' ? twiddleT2[r] : twiddle[r];
    if (!fn) return genericFallback(n, reIn, imIn);

    // Stage 1: decimate into r phases (element p + k*r for k=0..m-1), full
    // size-m sub-transform per phase, executed via the SAME plan tree
    // chooseComplex(m) predicted (sub), not re-chosen here.
    //
    // EXCEPT when sub is a 'q1' node (only ever produced by
    // Planner/chooseDecomposition.js's chooseComplexEmbedded(), used by
    // RaderSolver.js/BluesteinSolver.js for their own in-place internal
    // sub-transform, never by plain chooseComplex()) -- q1 genuinely reads
    // and writes ACROSS what would otherwise be independent phases (see
    // Codelets/complex/q1.js's header), so it can't be executed one phase
    // at a time. executeQ1Embedded takes the FULL (un-gathered) buffer
    // directly and produces workRe/workIm in the exact same p*m+col layout
    // the per-phase loop below would have, so no other change is needed.
    let workRe;
    let workIm;
    if (sub && sub.codeletGroup === 'q1') {
      const { executeQ1Embedded } = require('./Codelets/complex/q1');
      [workRe, workIm] = executeQ1Embedded(
        sub.radix, sub.m,
        (k, reK, imK) => executePlan(sub.sub, k, reK, imK),
        reIn, imIn,
      );
    } else {
      workRe = new Float64Array(n);
      workIm = new Float64Array(n);
      const riP = new Float64Array(m);
      const iiP = new Float64Array(m);
      for (let p = 0; p < r; p++) {
        for (let k = 0; k < m; k++) {
          riP[k] = reIn[p + k * r];
          iiP[k] = imIn[p + k * r];
        }
        const [roP, ioP] = executePlan(sub, m, riP, iiP);
        for (let k = 0; k < m; k++) {
          workRe[p * m + k] = roP[k];
          workIm[p * m + k] = ioP[k];
        }
      }
    }

    // Stage 2: twiddle-combine, one column (fixed k, all r phases) at a
    // time. Column position in the size-m sub-transform IS the twiddle
    // index (same "physical position == output index" property used
    // throughout this port -- see RaderSolver.js's dftSize36 header).
    const br = new Float64Array(r), bi = new Float64Array(r);
    const Wc = new Float64Array(r), Ws = new Float64Array(r);
    const outRe = new Float64Array(n), outIm = new Float64Array(n);
    for (let col = 0; col < m; col++) {
      for (let p = 0; p < r; p++) {
        br[p] = workRe[p * m + col];
        bi[p] = workIm[p * m + col];
      }
      for (let p = 1; p < r; p++) {
        const [c, s] = realCexp((col * p) % n, n);
        Wc[p] = c; Ws[p] = s;
      }
      const [outR, outI] = fn(br, bi, Wc, Ws);
      for (let p = 0; p < r; p++) {
        outRe[p * m + col] = outR[p];
        outIm[p * m + col] = outI[p];
      }
    }
    return [outRe, outIm];
  }

  if (plan.type === 'rader') {
    // Lazy require -- see isFullyPortedComplex's identical comment on why.
    const { dftRader } = require('./RaderSolver.js');
    const ro = new Float64Array(n), io = new Float64Array(n);
    dftRader(n, reIn, 0, imIn, 0, 1, ro, 0, io, 0, 1);
    return [ro, io];
  }

  if (plan.type === 'bluestein') {
    const { dftBluestein } = require('./BluesteinSolver.js');
    const ro = new Float64Array(n), io = new Float64Array(n);
    dftBluestein(n, reIn, 0, imIn, 0, 1, ro, 0, io, 0, 1);
    return [ro, io];
  }

  if (plan.type === 'buffered') {
    // dft-buffered's cldcpy step is pure data movement (see
    // isFullyPortedComplexEmbedded's 'buffered' comment) -- executing
    // plan.sub (chooseComplex(n)'s own, already-verified plan) directly
    // is bit-identical to what real FFTW's buffer-then-copy-back actually
    // produces.
    return executePlan(plan.sub, n, reIn, imIn);
  }

  // generic: FFTW3's own exact choice when n is odd (dftGeneric matches
  // exactly there); this is the final fallback for anything else.
  return genericFallback(n, reIn, imIn);
}

// dftComposite -- matches GenericSolver1D.dftGeneric's call signature
// (n, ri, riOff, ii, iiOff, is, ro, roOff, io, ioOff, os) so it drops into
// Rank2Orchestration.js's dft1D dispatcher the same way.
function dftComposite(n, ri, riOff, ii, iiOff, is, ro, roOff, io, ioOff, os) {
  const reIn = new Float64Array(n), imIn = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    reIn[i] = ri[riOff + i * is];
    imIn[i] = ii[iiOff + i * is];
  }
  const plan = chooseComplex(n).plan;
  const [reOut, imOut] = executePlan(plan, n, reIn, imIn);
  for (let i = 0; i < n; i++) {
    ro[roOff + i * os] = reOut[i];
    io[ioOff + i * os] = imOut[i];
  }
}

// dftComposeEmbedded -- same call signature/shape as dftComposite, but uses
// chooseComplexEmbedded(n, undefined) instead of chooseComplex(n) -- i.e.
// the IN-PLACE-aware (possibly q1-using) prediction, not the out-of-place
// one. This is what RaderSolver.js/BluesteinSolver.js's own internal
// sub-transform should call instead of dftComposite, since that sub-
// transform really is planned in-place by real FFTW (see
// isFullyPortedComplexEmbedded's header) even though this JS port always
// works through fresh (never literally aliased) typed arrays internally --
// the "in-place-ness" that matters here is about which SOLVER FFTW's own
// planner would pick, not about this file's own memory management.
function dftComposeEmbedded(n, ri, riOff, ii, iiOff, is, ro, roOff, io, ioOff, os) {
  const reIn = new Float64Array(n), imIn = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    reIn[i] = ri[riOff + i * is];
    imIn[i] = ii[iiOff + i * is];
  }
  const plan = chooseComplexEmbedded(n, undefined).plan;
  const [reOut, imOut] = executePlan(plan, n, reIn, imIn);
  for (let i = 0; i < n; i++) {
    ro[roOff + i * os] = reOut[i];
    io[ioOff + i * os] = imOut[i];
  }
}

module.exports = {
  dftComposite, isFullyPortedComplex, executePlan, isSingleLevelOrDirect, isFullyPortedComplexEmbedded,
  dftComposeEmbedded,
};

// ---------------------------------------------------------------------------
// Self-test (node Composite1D.js) -- smoke test only, NOT a substitute for
// verify/run_phase1_tests.js's real cross-language byte-exact comparison.
// Checks: DC term (X[0] must equal the plain sum of the input) across a
// spread of N built from ported radices {2,3,4,5,6}, plus a couple of
// generic-fallback sizes.
// ---------------------------------------------------------------------------
if (require.main === module) {
  for (const n of [2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 30, 36, 13, 17]) {
    const ri = new Float64Array(n), ii = new Float64Array(n);
    for (let i = 0; i < n; i++) { ri[i] = i + 1; }
    const sum = Array.from(ri).reduce((a, b) => a + b, 0);
    const fullyPorted = isFullyPortedComplex(n);
    try {
      const ro = new Float64Array(n), io = new Float64Array(n);
      dftComposite(n, ri, 0, ii, 0, 1, ro, 0, io, 0, 1);
      const ok = Math.abs(ro[0] - sum) < 1e-9;
      console.log(`n=${n}: DC=${ro[0]} (expect ${sum}) ${ok ? 'OK' : 'FAIL'} -- fullyPorted=${fullyPorted} plan=${chooseComplex(n).describe()}`);
    } catch (e) {
      console.log(`n=${n}: THROWN (${e.message.split('\n')[0]}) -- fullyPorted=${fullyPorted} plan=${chooseComplex(n).describe()}`);
    }
  }
}
