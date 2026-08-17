'use strict';

// =============================================================================
// PlanTable.js
// Coverage predicate for the FFTW3JS engine (Rank2Orchestration.js): which
// (N0, N1) 2D r2c/c2r sizes this hand-port can currently handle bit-exactly
// against real FFTW3, versus which must still fall back to ndarray-fft.
//
// Current coverage (see GenericSolver1D.js / RaderSolver.js /
// CompositeSolver1D.js headers) -- ROW and COLUMN axes have DIFFERENT
// coverage, because FFTW3 uses genuinely different solver families for the
// real (r2hc/hc2r) row pass versus the full-complex column pass:
//   - Row pass (r2hc/hc2r, size N1): any ODD PRIME N1 -- rdft-generic-r2hc/
//     hc2r, confirmed via fftw_print_plan for every prime probed so far
//     (37, 43, 47, 53). No composite-size real (RDFT) decomposition is
//     ported yet -- e.g. N=57 as a ROW size would need FFTW3's real-specific
//     radix-3 codelets (r2cf_3/hf_3), which are DIFFERENT from the complex
//     radix-3 codelet (t1_3) CompositeSolver1D.js ports; not attempted, out
//     of scope until a real galaxy size actually needs it as a row size.
//   - Column pass (full complex DFT, size N0): any ODD PRIME N0 -- dft-generic
//     for most primes; dft-rader-37 for N0=37 (36=6*6 factors favorably);
//     dft-ct-dit/3 (radix-3 Cooley-Tukey) for N0=3*p, p an odd prime (e.g.
//     N0=57=3*19, a real production padded size). All three dispatched
//     automatically inside Rank2Orchestration.js's dft1D, transparent here.
// Both N0 and N1 must be covered for the *whole* 2D transform to use this
// engine -- the row and column passes share one in-place working array
// (Rank2Orchestration.js), so there is no sound way to run one pass on this
// engine and the other on ndarray-fft.
//
// NOT YET covered: even sizes; composite row (N1) sizes; column (N0) sizes
// with no useful small-radix decomposition (e.g. products of two large
// primes). Extend incrementally as real galaxy padded sizes require it --
// anything not covered falls back to ndarray-fft, never silently produces
// wrong output.
//
// PRODUCTION-GATING RULE (binding, see FFTW3JS/README.md):
// isRowSizeSupported/isColumnSizeSupported/isTransformSupported returning
// `true` is a CORRECTNESS CLAIM, not an optimistic guess -- a size may only
// be added here after it has a corresponding byte-exact-verified case in
// verify/fixtures/ and passes verify/run_phase1_tests.js. Never widen these
// predicates as an incidental side effect of refactoring or generalizing the
// underlying solvers; widen them as a deliberate last step, once the new
// case is proven. `false` is always the safe default -- it only costs speed
// (fallback to ndarray-fft), whereas an unearned `true` would silently ship
// wrong output to production.
//
// UPDATE (wiring the general engines in -- column axis SEVERELY RESTRICTED,
// see below): the odd-prime / 3*oddprime column check and odd-prime-only
// row check above predate Composite1D.js's general recursive complex engine
// and RealEngine1D.js's general recursive real engine, which cover far more
// sizes than this file claimed -- but isFullyPortedComplex(n)/chooseComplex
// are STANDALONE 1D predictions, and the column pass here never runs
// standalone: it's embedded in a 2D r2c/c2r transform, invoked as a
// vector-looped sub-solver (FFTW's own plan text: "dft-vrank>=1-xK/1 (...)"
// wrapping the column transform). FFTW_ESTIMATE implies FFTW_NO_VRECURSE,
// which can make FFTW pick a COMPLETELY DIFFERENT decomposition when
// embedded than when the same N is planned standalone -- this is the exact
// same NO_VRECURSE mechanism documented (and left unresolved) for q1's
// in-place trigger in Composite1D.js/README.md, just surfaced here via
// vector-loop embedding instead of in-place-ness. Confirmed via
// fftw_fprint_plan on a real fftw_plan_dft_r2c_2d: standalone N=25 and N=64
// each predict a Cooley-Tukey recursion (dft-ct-dit/5 and dft-ct-dit/8
// respectively), but EMBEDDED in a 2D r2c plan, FFTW picks a single DIRECT
// codelet instead (n1_25 / n1_64) -- a different codelet sequence, silently
// wrong bit pattern if executed via the standalone-predicted CT structure.
// N=128 diverges even more severely, embedding as a q1-based decomposition
// (q1 is entirely unported). Composites like N=40/50, whose standalone plan
// is ALSO 'ct', stayed byte-identical when embedded (FFTW just wraps the
// SAME codelet sequence in a memory-layout-only "directbuf" adapter) --
// so this is NOT a blanket "CT is unsafe" rule, it's specific to whichever
// N NO_VRECURSE's cost adjustments cause an outright DIFFERENT solver
// choice for, and that set isn't characterized well enough yet to trust in
// general (same open problem as q1 -- see Composite1D.js's isFullyPortedComplex
// header). The one PROVABLY safe subset: N whose standalone chooseComplex(N)
// plan is 'direct' (a single base-case codelet, no CT recursion at all) --
// there's no recursive candidate for NO_VRECURSE to exclude in the first
// place, so embedding cannot change the choice. Confirmed via
// fftw_fprint_plan across N=12,16,20,32 (all 'direct' standalone) staying
// exactly "dft-direct-N-xK" when embedded. This is the ONLY new column
// widening applied here -- composite CT sizes (even ones already proven
// bit-exact standalone, like 40/44/50/64/100 from this session's sweep)
// stay UNsupported as column sizes until the embedded-vs-standalone
// divergence is properly characterized, not assumed safe by analogy.
// The row-axis widening (isFullyPortedR2HC/isFullyPortedHC2R) doesn't have
// this problem YET -- every size RealEngine1D.js currently fully-ports
// (2,3,4,5) is ALSO a 'direct' base case (chooseReal confirmed), and
// fftw_fprint_plan on an embedded 2D case (37x4) confirms the row codelet
// stays identical, just vector-looped. Revisit this comment if
// RealEngine1D.js ever grows a 'ct'-recursive path.
// =============================================================================

const { isFullyPortedComplex } = require('./Composite1D.js');
const { isFullyPortedR2HC, isFullyPortedHC2R } = require('./RealEngine1D.js');
const { chooseComplex } = require('./Planner/chooseDecomposition');

// True iff n's STANDALONE chooseComplex plan is a single direct codelet --
// the only column shape proven immune to the embedded/vector-loop
// NO_VRECURSE divergence documented above (no CT recursion exists to be
// excluded, so there is nothing for embedding to change). n<2 guarded
// separately -- chooseComplex/chooseReal throw for degenerate sizes rather
// than returning an "unsupported" plan, and this predicate must stay a
// safe `false`, never an exception, for any n a caller might pass in.
function isDirectOnlyComplex(n) {
  return n >= 2 && isFullyPortedComplex(n) && chooseComplex(n).plan.type === 'direct';
}

function isPrime(n) {
  if (n <= 1) return false;
  if (n % 2 === 0) return n === 2;
  for (let i = 3; i * i <= n; i += 2) {
    if (n % i === 0) return false;
  }
  return true;
}

// Odd prime specifically -- n=2 is prime but not reachable/useful here (row
// pass needs the "+i < n" hartley-fold structure, which assumes odd n; even
// n=2 would need its own direct codelet path, not attempted).
function isOddPrime(n) {
  return n % 2 === 1 && isPrime(n);
}

// Row pass (r2hc/hc2r): odd primes (rdft-generic-r2hc/hc2r, always correct
// there) OR whatever RealEngine1D.js's own general recursive engine has a
// ported codelet path for (both directions required, since a round-trip
// convolution needs forward AND inverse on the same size -- see header).
// n<2 guarded for the same reason as isDirectOnlyComplex above.
function isRowSizeSupported(n) {
  return isOddPrime(n) || (n >= 2 && isFullyPortedR2HC(n) && isFullyPortedHC2R(n));
}

// Column pass (full complex DFT): odd primes (dft-generic / dft-rader-37),
// n = 3*p for odd prime p (dft-ct-dit/3, CompositeSolver1D.js), OR n whose
// standalone decomposition is a single direct codelet (isDirectOnlyComplex
// -- see header for why ONLY this subset of Composite1D.js's broader
// coverage is safe to claim here).
function isColumnSizeSupported(n) {
  if (isOddPrime(n)) return true;
  if (n % 3 === 0) {
    const p = n / 3;
    if (p > 1 && isOddPrime(p)) return true;
  }
  return isDirectOnlyComplex(n);
}

// Backwards-compatible alias -- historically both axes had identical
// (odd-prime-only) coverage, so callers that don't care about the row/column
// distinction can still ask "is n supported at all, on either axis".
function isSizeSupported(n) {
  return isRowSizeSupported(n) || isColumnSizeSupported(n);
}

function isTransformSupported(N0, N1) {
  return isColumnSizeSupported(N0) && isRowSizeSupported(N1);
}

module.exports = {
  isPrime, isOddPrime,
  isSizeSupported, isRowSizeSupported, isColumnSizeSupported,
  isTransformSupported,
};

// ---------------------------------------------------------------------------
// Self-test (node PlanTable.js) -- smoke test only.
// ---------------------------------------------------------------------------
if (require.main === module) {
  const cases = [37, 43, 47, 53, 57, 36, 38, 40, 46, 2, 1, 0, 9];
  console.log('=== isRowSizeSupported / isColumnSizeSupported ===');
  for (const n of cases) {
    console.log(`n=${n}: row=${isRowSizeSupported(n)} column=${isColumnSizeSupported(n)}`);
  }
  console.log('=== isTransformSupported (N0=column, N1=row) ===');
  console.log('(37,37):', isTransformSupported(37, 37), '(expect true)');
  console.log('(37,43):', isTransformSupported(37, 43), '(expect true)');
  console.log('(57,53):', isTransformSupported(57, 53), '(expect true -- realization-200 case)');
  console.log('(53,57):', isTransformSupported(53, 57), '(expect true -- 57=3*19 became a covered ROW size once radix-3 (hf_3/hb_3) and prime<=167 generic-fallback trust landed)');
  console.log('(36,38):', isTransformSupported(36, 38), '(expect false)');
  console.log('(43,40):', isTransformSupported(43, 40), '(expect false)');
}
