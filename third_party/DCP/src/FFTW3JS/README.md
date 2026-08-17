# FFTW3JS

A hand-ported JavaScript re-implementation of FFTW3 3.3.8, built to be
bit-for-bit identical to real compiled FFTW3 output under `FFTW_ESTIMATE`
planning -- not just numerically close, and not just *a* correct FFT, but
*FFTW3's own chosen decomposition and arithmetic*, reproduced exactly. See
`NOTICE` for licensing (GPLv2+, inherited from FFTW3 -- distinct from the
MIT license covering the rest of this repository).

## Scope

General in transform size N. Fixed to:

- **Rank 2** (2D transforms only)
- **Unit stride** (no arbitrary strides, no vector loops)
- **Real <-> complex** (`r2c`/`c2r`) only -- no complex-to-complex top-level
  entry point, no DCT/DST (`REDFT`/`RODFT`)

This matches what the 3KIDNAS beam-convolution pipeline actually needs.
FFTW3 itself implements DCT/DST as a separate re-indexing layer (`reodft/`)
over R2HC rather than as dedicated codelets, so extending this port to cover
them later is additive, not a redesign. Arbitrary rank and non-unit strides
are a similar story -- structurally separate from what's built here, not
blocked by it.

## Ground truth

`verify/ground_truth_harness.c` links real compiled FFTW3 and can print its
actual chosen plan (`fftw_print_plan`) or execute a transform and dump raw
output, for any N. This is the source of truth this port is checked against
-- never hand-derived expectations.

## Production-gating rule

`PlanTable.js`'s `isTransformSupported(N0, N1)` controls whether production
code (`ConvolveCube/CubeKernelConvolution.js`) routes a given size through
this engine or falls back to the existing (slower, but independently
correct) `ndarray-fft`-based path. A size may only be flipped from
unsupported to supported *after* it has a corresponding byte-exact-verified
case in `verify/fixtures/` and passes `verify/run_phase1_tests.js` --
never as an incidental side effect of a coverage predicate becoming more
permissive for unrelated reasons. `false` here is always safe (slower
fallback); `true` is a correctness claim and must be earned by a passing
fixture.

## Layout

- `Trig.js` -- twiddle-factor / trig generation (`kernel/trig.c` port)
- `GenericSolver1D.js` -- universal O(n) fallback DFT/r2hc/hc2r (`rdft/generic.c`, `dft/generic.c` ports); always correct, used whenever no faster codelet path applies
- `RaderSolver.js` -- Rader's algorithm for prime N (complex side)
- `BluesteinSolver.js` -- Bluestein's algorithm for prime N the Rader path doesn't cover (complex side)
- `CompositeSolver1D.js` -- small Cooley-Tukey composite cases (pre-`Composite1D.js`, kept as a proven fallback)
- `Composite1D.js` -- general recursive complex 1D engine: executes whatever `Planner/chooseDecomposition.js` predicts, via the ported codelet set in `Codelets/complex/`. This is where most of the complex-side coverage now lives (radixes 2-16, 20, 25, 32, 64, both direct and twiddle families).
- `RealEngine1D.js` -- the real-side analogue of `Composite1D.js`, currently only at its initial "thin slice" (radixes 2-5) -- the codelet-family widening `Composite1D.js` went through has NOT yet happened here. This is the main remaining engine gap, not the complex side.
- `Rank2Orchestration.js` -- 2D r2c/c2r driver combining row + column 1D sub-transforms
- `PlanTable.js` -- production support-gating predicate (see above)
- `Planner/` -- FFTW's `FFTW_ESTIMATE` solver-selection logic, replicated (ops-cost tables + LIFO search + pruning)
- `Codelets/` -- one file per ported FFTW codelet, mirroring FFTW's own per-file layout
- `verify/` -- ground-truth harness, fixtures, and the regression test suite

## Known gap: embedded (2D) solver choice can differ from standalone

`Composite1D.js`'s `isFullyPortedComplex(n)` and `chooseComplex(n)` predict
what FFTW would choose for `n` as a **standalone** 1D transform. But
`Rank2Orchestration.js`'s column pass never runs standalone -- it's embedded
in a 2D r2c/c2r transform as a vector-looped sub-solver, and
`FFTW_ESTIMATE`'s implied `FFTW_NO_VRECURSE` can make real FFTW pick a
**different** decomposition once embedded than it would standalone (same
underlying mechanism as the in-place `q1` codelet family staying unported --
see `Composite1D.js`'s own comments on that). Confirmed via
`fftw_fprint_plan`: standalone N=25 and N=64 each predict a Cooley-Tukey
recursion, but embedded in a real `fftw_plan_dft_r2c_2d`, FFTW picks a
single direct codelet instead -- a different, unported codelet path, so
executing the standalone-predicted structure there would silently produce
the wrong bit pattern. `PlanTable.js`'s `isColumnSizeSupported` only widens
past the pre-existing odd-prime/3\*oddprime check for sizes whose standalone
plan is a single direct codelet (`isDirectOnlyComplex`) -- the one shape
provably immune to this divergence, since there's no recursion for
`NO_VRECURSE` to exclude. Composite (`ct`-recursive) column sizes are NOT
claimed here even when `Composite1D.js` proves them bit-exact standalone
(e.g. 40, 44, 50, 100) -- see `PlanTable.js`'s header for the full writeup
and the specific `fftw_fprint_plan` evidence. Properly characterizing which
composite sizes ARE safe when embedded is unstarted, open work, in the same
vein as the `q1` investigation.
