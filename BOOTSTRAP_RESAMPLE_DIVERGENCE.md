# Bootstrap-resampling Fortran/JS divergence — investigation log

**Status as of 2026-08-18 (later pass): a genuine, unresolved paradox.**
The `gasdev` cache-reset fix (below) is now verified, via fine-grained
bisection, to hold **perfectly across all 5 rings' entire particle loops**
(every 200th particle plus the last one — `idum` and velocity bit-exact,
every checkpoint). And yet the pre-convolution model cube, checked via a
separate, independently-verified diagnostic print, still shows real,
repeatable, substantial divergence at the same tracked pixels. These two
findings are individually solid but contradict each other, and that
contradiction is not yet resolved — see "The bisection paradox" near the
end. This doc exists so the investigation can be picked up cold, by anyone
(including a future session with no memory of how we got here).

## The problem, precisely

The **single initial fit** (no bootstrapping, `nBootstraps=0`) achieves
near-perfect Fortran/JS agreement at `cdens=400`: Inc/PA/VSys/chi2 all agree
to ~1e-7 relative (the accepted, unavoidable float32 cross-platform noise
floor for this whole project). This was validated repeatedly and is not in
question.

**Real bootstrap realizations do not show this same agreement.** A
`run_both.js` run of 20 real bootstrap realizations (same seed=42, matched
`cdens`) shows mean Inc divergence ~3-4°, PA ~3-4°, with individual
realizations up to 10-13° off — while the single initial fit underneath it
matches almost perfectly. This was confirmed at both `cdens=500` and
`cdens=400` (ruling out `cdens=500` being simply "untested territory" as
the sole explanation), and confirmed not to be caused by:
- A stale/cached published DCP package (republished under a bumped
  version, then under a brand-new package name `3kidnas-test2`; identical
  divergence both times — ruled out caching definitively).
- `cdens` failing to propagate to one leg (verified directly: the actual
  dispatched JSON payload contains the correct `cloudBaseSurfDens`, and the
  worker-side code requires that exact field, no silent fallback).

## Where it's isolated to

A full pixel-by-pixel diff of the **resampled bootstrap cube itself**
(before it even reaches the fit) between Fortran and JS, same seed, same
`cdens=400`, single realization (`nBootstraps=1`), shows:
- 81,501-81,585 of 179,520 cells (~45%) are **bit-identical**.
- 63 cells have absolute flux difference > 1e-4 (max ~2.14e-4), clustered
  in small (~3-5 pixel) patches within individual channels, spread across
  ~36 of 102 channels, **never at the cube's spatial edges**.
- This magnitude (up to 2500% relative on some cells) is ~1000x+ larger
  than the established float32 noise floor — not explainable by simple
  continuous ULP-level amplification through interpolation.

This pattern was traced (see "Diagnostic infrastructure" below) down to:
**the *model* cube (each platform's own independently-generated best-fit
model, built by Monte-Carlo particle binning) differs substantially and
systematically at these locations, while the *observed* cube (read
directly from the same FITS file, no randomness involved) is
bit-identical.** This is now proven, not inferred — see SRCTRACE results
below.

So the divergence is definitively NOT in: the resampling transform code,
the coordinate flip/rotation math, trilinear interpolation, boundary
checks, or corner-pixel selection. It is definitively IN: how each
platform's own model cube gets populated with particles.

## Hypotheses tested, with hard evidence (four disproven, one confirmed)

### 1. Discrete corner-selection flip in `GetFluxAtPoint`/`getFluxAtPoint`
**Theory:** `int(Pt)` truncation picks which 8 corners get trilinear-interpolated;
if `Pt` lands near a pixel boundary, the two platforms could pick different
corners.
**Fix applied:** `RoundForInterpStability`/`roundForInterpStability` — masks
the low 12 mantissa bits (same technique as the particle-count fix earlier
this session) before the truncation.
Files: `src/BootstrapSampler/GenerateBootstrap.f`,
`js/src/BootstrapSampler/GenerateBootstrap.js`.
**Result: zero effect.** Byte-identical pixel diff before/after (same 63
cells, same max diff to 17 significant figures).
**Why it didn't work:** direct tracing (`PXTRACE`) showed the coordinate
and corner selection already matched between platforms for the tracked
pixel — there was no boundary to fix.

### 2. Native trig instead of fdlibm-forced trig in `PhysCoordTransform`
**Theory:** this file was missed by the earlier project-wide fdlibm-forcing
sweep (`SingleRingGeneration.f`, `EstimateShape.f`, etc. were fixed;
`BootstrapSampler/PhysCoordTransform.f` was not) — it used plain
`cos`/`sin`/`atan2` instead of `fd_cos`/`fd_sin`/`FullCircATan`, which is
exactly the class of cross-platform trig ULP mismatch this project already
solved everywhere else.
**Fix applied:** routed `GetPhysCoords`/`getPhysCoords` and
`GetCubeCoords`/`getCubeCoords` through `fd_cos`/`fd_sin`/`FullCircATan`
(Fortran) and `fdCos`/`fdSin`/`fullCircATan` (JS).
Files: `src/BootstrapSampler/PhysCoordTransform.f`,
`js/src/BootstrapSampler/PhysCoordTransform.js`.
**Result: zero effect.** Same 63 cells, same max diff (0.00021399115212261677,
identical to 17 sig figs) as before the fix.
**Note:** this is still a real, legitimate correctness fix (worth keeping
regardless) — it just isn't the cause of this specific divergence.

### 3. Upstream geometry mismatch between the two independent initial fits
**Theory:** each bootstrap realization resamples using its OWN platform's
independently-computed best-fit geometry (X/Y/PA/Inc/VSys) — if those
differ non-trivially, everything downstream would too.
**Test:** compared full-precision geometry directly:
`Inc_kin=84.055496°` (Fortran) vs `84.05549621582031°` (JS),
`PA_kin=236.757477°` vs `236.75750732421875°`,
`VSys=5740.654297` vs `5740.654296875`.
**Result: ruled out.** Agrees to ~1e-7 relative — the same accepted noise
floor as everywhere else, not a meaningfully different input.

### 4. Discrete per-particle binning-cell flip in `FindParticleCellLocation`
**Theory:** every generated particle's flux gets added (not interpolated)
to whichever single cell `int(Pos+0.5)` rounds it into. With 50,000-70,000+
particles per ring, even a tiny per-particle chance of landing near a
`.5` boundary becomes a near-certainty that *some* particles round
differently between platforms — and since flux is added wholesale (not
blended), one straddling particle yanks real flux from one cell into its
neighbor.
**Fix applied:** `RoundForBinStability`/`roundForBinStability` — same
mantissa-masking technique, applied to the position immediately before
`int(pos+0.5)`.
Files: `src/TiltedRingToDataCube/FillDataCubeByTiltedRing.f`,
`js/src/TiltedRingToDataCube/FillDataCubeByTiltedRing.js`.
**Result: zero effect.** Same 63 cells, same max diff again.
**Also disproved a stronger version of this idea directly:** forcing BOTH
platforms to resynthesize the model from an *identical, externally-supplied
15-parameter vector* (`TRACE_OVERRIDE_PARAMS`, see below) — removing
parameter differences from the equation entirely — left the divergence at
the tracked pixels essentially unchanged (same magnitude, same sign, same
cells). This ruled out parameter differences as the cause and strongly
implied `idum` state was the remaining variable... which led to:

### 5. `idum` divergence from differing Nelder-Mead search paths — initially (wrongly) marked disproven, then confirmed via a corrected test
**Theory:** even when two independent fits converge to nearly the same
final answer, the *sequence* of evaluations each platform's optimizer
takes to get there is not guaranteed identical (this project's own
well-documented Nelder-Mead instability). Since `idum` is a running total
of every random draw consumed since the fit began, and the model cube is
resynthesized *after* convergence using whatever `idum` happens to be at
that point, the two platforms could easily be drawing from completely
different, uncorrelated positions in their own (individually bit-exact)
`ran2`/`gasdev` streams.

**First test (flawed) — appeared to disprove this.** Added
`TRACE_OVERRIDE_IDUM` on both platforms — an env var read right before the
final model-cube resynthesis call (`OutputCube` in `src/Outputs/FitOutput.f`
via a new `MaybeOverrideIdum` subroutine; `runInitialFit` in
`js/bootstrap-realization-launcher.js`, replacing `state.rng` with a fresh
`makeRng(overrideIdum)` right before `tiltedRingModelComparison`). Ran both
platforms with the same fixed idum (-12345) AND the same fixed 15-param
vector (Fortran's own converged values). The divergence persisted
essentially unchanged, which looked like a clean disproof.

**The test itself was confounded — found by tracing individual particles.**
`Ring_ParticleGeneration`'s first 5 particles per ring are already dumped by
pre-existing `PARTTRACE` instrumentation (gated on `WRKP_TRACE_DEBUG=1`/
`TRACE_DEBUG=1`). Running it under the override test showed: **particle 0's
*position* matched bit-exactly between platforms, but its *velocity*
(`Ring_CalcParticle_VSys`, which draws one `gasdev` call) did not, and
`idum`'s value right after particle 0 had already diverged.** That
pinpointed the mismatch to `gasdev` specifically, not `ran2` or particle
positions.

**Root cause: `gasdev`'s cached-spare-value state is not reset the same way
`idum` is.** The classic Numerical-Recipes `gasdev` (Box-Muller transform)
generates a *pair* of Gaussian values per expensive computation, returns one,
and caches the other (`gset`) for the very next call — governed by a toggle
flag (`iset`). Fortran's `iset`/`gset` were `SAVE`'d local variables,
persisting for the *entire program run*, only cleared by the line
`if (idum.lt.0) iset=0` — i.e. only when `gasdev` observes `idum` as
negative *at the moment it's called*. But `ran2`'s own reinitialization block
(`if(idum.le.0) then idum=max(-idum,1) ... end if`) converts a negative seed
to a positive working value on its *very first call* — and position
generation (`Ring_ParticlePosSelect`, 3 `ran2` calls) always runs before
`Ring_CalcParticle_VSys`'s `gasdev` call. So by the time `gasdev` is ever
reached for particle 0, `idum` is already positive, and the reset never
fires — Fortran's `gasdev` used whatever cache state was left over from
**the fit's entire prior evaluation history**, completely unrelated to the
override. JS's `makeRng(idum)`, by contrast, allocates a brand-new
`gasdevState` object every time it's called, starting at `iset=0`
unconditionally — always "fresh," regardless of the override.

**This means the "disproof" was itself an artifact of the test, not a real
result about the underlying theory.**

**Fix:** moved `gasdev`'s `iset`/`gset` from function-local `SAVE` variables
to module-level variables in `BasicRanNumGen` (`src/StandardMath/random.f`),
and added `ResetGasdevCache()`, called from `MaybeOverrideIdum` right after
setting the overridden `idum` — giving the diagnostic override a real,
complete RNG reset on both platforms for the first time.
Files: `src/StandardMath/random.f`, `src/Outputs/FitOutput.f`.

**Result, re-tested with the fix: real, measured, but partial improvement.**
Re-ran the identical-idum + identical-params test:
- Particle-level (`PARTTRACE`, first 5 particles of **all 5 rings**):
  bit-exact match on both platforms — position, velocity, and `idum` state,
  every single one. Confirms the fix works correctly at the point it's
  applied.
- Full model-cube diff: max abs diff dropped from 0.00115 (broken cache) to
  **0.00049** (fixed cache) — real, substantial, reproducible improvement —
  but not zero. 1743 cells still exceed 1e-4 (down from 2683).

**So this is a real, confirmed bug, now fixed** — but it does not fully
explain the divergence by itself. Since only the first 5 particles per ring
are traced (not all 7,837-70,529), the residual could be a *similar* desync
that resurfaces later in a ring's particle loop (not yet localized), or
could be unrelated noise (e.g. the already-characterized FFTW native-vs-WASM
convolution precision difference from earlier this session, since beam
convolution runs after particle binning and hasn't been isolated out of
either cube compared here).

## Where this leaves us

**Ruled out, all confirmed bit-identical between platforms in the
gasdev-fixed controlled test:**
- The model cube's own header fields (`start`, `channelSize`, `refVal`) —
  `HEADERTRACE`.
- Every field `generalizedParamVectorToTiltedRing` produces for ring 0,
  both actively-fitted (`Rmid`, `Inclination`, `PositionAngle`, `VSys`,
  `VRot`, `SigUse`, `CentPos`) and fixed/constant (`VRad`, `VDisp`, `Vvert`,
  `dvdz`, `z0`, `zGradiantStart`) — `RINGTRACE`/`RINGTRACE2`.
- `X**2.` vs `X*X` in `gasdev`'s rejection-sampling loop (`rsq=v1**2.+v2**2.`)
  — a real, plausible-sounding lead (Fortran real-exponent power vs JS's
  plain multiply aren't guaranteed identical by the language standard) —
  directly tested with a standalone gfortran program across 10
  representative values spanning the `[-1,1]` range `ran2` draws produce;
  bit-identical (verified via hex dump of the IEEE754 representation) in
  every case at `-O0`. Not the cause.
- Other Fortran `SAVE` state: audited every `SAVE` statement in `src/`
  (`grep -rniE "^\s*save\b"`) — only 4 exist. `ran2`'s own `iv/iy/idum2`
  resets reliably (its own `if(idum.le.0)` check is the first thing that
  runs, and `ran2` is always the first RNG call for a fresh `idum`, unlike
  `gasdev`'s problem). `ran0`/`ran1`/`ran3` have their own `SAVE`'d state
  but are never called anywhere outside `random.f` itself (confirmed via
  grep) — dead code, not a concern. `FullModelComparison.f`'s
  `TraceCallCounter` only gates `print` statements, confirmed no effect on
  any computed value.
- JS module-level persistent state: audited every top-level `let`/`var` in
  `js/src/` — all are either load-once caches for pure/deterministic
  functions (fdlibm-wasm module handles, FFT twiddle-factor memoization in
  `RaderSolver.js`'s `_omegaCache`) or purely-diagnostic counters
  (`traceCallCounter`, `evalCount`, `convolveMs`/`convolveCalls`, all
  reset once per realization via their own `resetXStats()` calls). None
  hold sequence-dependent state the way `gasdev`'s cache did.

**The bisection paradox — genuinely unresolved.** Added a coarse
checkpoint trace (`BISECTTRACE`, every 200th particle plus the last one)
across each ring's *entire* particle loop, not just `PARTTRACE`'s first 5.
Under the same controlled test (matched idum, matched params): **every
single checkpoint, across all 5 rings, matches bit-exactly between
platforms** — both `idum` and the per-particle velocity. Since `ran2`/
`gasdev` are pure functions of `idum` (given the now-fixed cache), and
`idum` matching at every 200-particle boundary should mathematically force
every particle in between to also match, this looks like proof that the
*entire* particle-generation stage is now correct, ring to ring, top to
bottom.

**And yet** the pre-convolution model cube — checked via `PRECONVTRACE`/
`PRECONVPIX`, a simple, direct, `console.error`/`print` style checksum and
8-sample-pixel dump, re-verified consistent and reproducible across
multiple independent runs of the identical controlled test — still shows
real, substantial divergence at the same tracked pixels (e.g. Fortran
`0.00025334` vs JS `0.00032994` at one cell — not noise-floor, a genuine
~30% difference), unchanged from before this session's bisection work.

These two results cannot both be true under the model of the code this
investigation has built up. One of them is wrong, or there's a mechanism
neither test is sensitive to. **Not yet resolved.**

**A dead end, noted so it isn't repeated:** attempted to resolve the
paradox with a full pixel-by-pixel dump of the pre-convolution cube (FITS
via `WriteDataCubeToFITS` on the Fortran side, a raw synchronous
`Float32Array` dump on the JS side). The dumped cubes' aggregate sum and
max did not match *either platform's own* `PRECONVTRACE` print from the
exact same run (Fortran: FITS-dump sum `0.718` vs its own print `0.620`)
— meaning `WriteDataCubeToFITS` (or something in that path) applies some
transform not present in a direct in-memory read, making it unsuitable for
this specific unscaled, pre-convolution diagnostic. This code was added,
found unreliable, and **removed** rather than left in the tree — don't
rebuild it without first understanding what `WriteDataCubeToFITS` does
differently from a raw read (check `WriteFITSHeader`/`WriteFITSBody` for
implicit scaling, e.g. BSCALE/BZERO or a beam-area-style conversion, before
trusting a FITS round-trip for cubes that aren't already in the final
output's unit convention).

Concrete next steps, in rough priority order:

1. **Resolve the paradox directly.** Since `BISECTTRACE` only checks
   velocity (a `gasdev` output) at 200-particle intervals, add a similar
   checkpoint for *position* (`Pos(0)`/`Pos(1)`, the `ran2`-only outputs)
   at the *same* checkpoints, and — more importantly — narrow the
   checkpoint interval (e.g. every particle, or binary-search within a
   200-particle window) specifically around whichever ring turns out to
   contribute to the tracked pixel region, to find the first individual
   particle where Fortran and JS actually disagree. If none is found even
   at full resolution, the bug is downstream of particle generation after
   all (binning or convolution) and the aggregate checksum divergence
   needs a *correctly-scaled* full-cube comparison (fix the `WriteDataCubeToFITS`
   scaling issue above, or find another synchronous, unscaled dump path).
2. Which ring(s) actually contribute flux to the tracked pixel region
   (`x∈{14,15}, y∈{23,24}, ch∈{89,90}`) has never been directly confirmed
   — only estimated from radius. Multiple rings' particles can land in the
   same output cell. Confirming this would narrow where to look.
3. Not yet ruled out: some other input to `Ring_ParticleGeneration` besides
   idum/Sigma/ring-geometry (e.g. `CloudBaseSurfDens`, `cmode` — simple
   scalars, unlikely but not directly re-verified in this specific
   controlled test).

## Diagnostic infrastructure now in place (kept in the codebase, TRACE-gated)

All of the following are cheap, permanent, gated behind existing
`WRKP_TRACE_DEBUG=1` (Fortran/Python side) / `TRACE_DEBUG=1` (JS side) or
their own dedicated env vars — safe to leave in place, zero cost when unset.

| Trace | What it prints | Where |
|---|---|---|
| `"TRACE resampled cube"` | shape/sum/min/max/corner-pixel checksum of the resampled cube | `Bootstrap_Error_Analysis.py`, `bootstrap-realization-launcher.js` (pre-existing) |
| `TRACE_DUMP_RESAMPLE_PATH` | writes the JS resampled cube's exact FITS bytes to a given path | `bootstrap-realization-launcher.js` (pre-existing) |
| `PXTRACE` | `FlipType, ChanID, j, i, CubeCoordFlip(1:3), BoundCheck` for a **hardcoded** pixel (currently `ChanID=14, j=15, i=29`) | `src/BootstrapSampler/FlippingBootstrap.f`, `js/src/BootstrapSampler/FlipBootstrap.js` |
| `CORNERTRACE` / `CORNERTRACE_OUT` | the 8 corner flux values and interpolated output, for a **hardcoded** `CurrIndx` (currently `14,23,89`) | `src/BootstrapSampler/GenerateBootstrap.f`, `js/src/BootstrapSampler/GenerateBootstrap.js` |
| `SRCTRACE` | observed/model/diff flux values for a **hardcoded** pixel range (currently `i∈{14,15}, j∈{23,24}, k∈{89,90}`) | `src/BootstrapSampler/CubeDifference.f`, `js/src/BootstrapSampler/CubeDifference.js` |
| `NPTRACE` | per-ring `Rmid, Sigma, DensMultiplications, PixelRing, nParticles` | `src/TiltedRingModelGeneration/SingleRingGeneration.f`, `.js` (pre-existing, from earlier session work) |
| `TRACE_OVERRIDE_PARAMS` | forces the 15-param vector used for the post-convergence model resynthesis (JSON array of 15 numbers) | `js/bootstrap-realization-launcher.js` only (pre-existing) — **no Fortran equivalent exists**; Fortran's own converged vector is read from its `CONVERGED_VECTOR` trace line instead |
| `TRACE_OVERRIDE_IDUM` | forces `idum` used for that same resynthesis call, on **both** platforms, AND now also resets `gasdev`'s cached-spare-value state (via `ResetGasdevCache`) — a true full RNG reset, not just the LCG portion | `src/Outputs/FitOutput.f` (`MaybeOverrideIdum`), `js/bootstrap-realization-launcher.js` (added during this investigation) |
| `HEADERTRACE` | model cube header: `nPixels, nChannels, Start(0:2), RefVal(0:2), ChannelSize` — fires once, gated on `TRACE_OVERRIDE_IDUM` | `src/TiltedRingToDataCube/FillDataCubeByTiltedRing.f`, `js/src/TiltedRingToDataCube/FillDataCubeByTiltedRing.js` |
| `RINGTRACE` / `RINGTRACE2` | ring 0's full field set (`Rmid, Rwidth, Inclination, PositionAngle, VSys, VRot, VRad, VDisp, Vvert, dvdz, SigUse, z0, zGradiantStart, CentPos`) right after param-vector deserialization — gated on `TRACE_OVERRIDE_IDUM`, fires every evaluation (grab the **last** occurrence for the override-resynthesis call) | `src/Outputs/FitOutput.f` (`MaybeTraceRingFields`), `js/src/CompareCubes/FullModelComparison.js` |
| `PARTTRACE` | per-particle `Rmid, i, idum, Pos(0:2), ProjectedPos(0:1), ProjectedVel(2)` for the first 5 particles of **every ring** | `src/TiltedRingModelGeneration/SingleRingGeneration.f`, `js/src/TiltedRingModelGeneration/SingleRingGeneration.js` (pre-existing, gated on `WRKP_TRACE_DEBUG=1`/`TRACE_DEBUG=1`, not `TRACE_OVERRIDE_IDUM`) — this is what actually caught the gasdev-cache bug |
| `BISECTTRACE` | `Rmid, i, idum, ProjectedVel(2)` for **every 200th particle plus the last one, across every ring's full loop** (not just the first 5) — gated on `TRACE_OVERRIDE_IDUM`. Verified: matches bit-exactly at every checkpoint, all 5 rings, in the gasdev-fixed controlled test | `src/TiltedRingModelGeneration/SingleRingGeneration.f`, `js/src/TiltedRingModelGeneration/SingleRingGeneration.js` (added this session) |
| `PRECONVTRACE` / `PRECONVPIX` | model-cube `sum/min/max` checksum and 8 tracked pixel values (`x∈{14,15}, y∈{23,24}, ch∈{89,90}`), taken right **before** beam convolution — gated on `TRACE_OVERRIDE_IDUM`. Verified reproducible across independent runs; this is what caught the still-unresolved residual divergence | `src/Outputs/FitOutput.f` (`MaybeTracePreConvChecksum`), `js/src/CompareCubes/FullModelComparison.js` (added this session) |

**Tried and removed — do not rebuild without reading this first:** a full
pixel-by-pixel pre-convolution cube dump (`MaybeDumpPreConvCube` in
Fortran via `WriteDataCubeToFITS`, plus a matching raw `Float32Array`
dump in JS, both gated on a `TRACE_DUMP_PRECONV_PATH` env var). Removed
after its own aggregate sum/max didn't match `PRECONVTRACE`'s print from
the identical run — see "A dead end, noted so it isn't repeated" above.
Neither the Fortran subroutine nor the JS block exist in the tree anymore.

**To retarget the hardcoded pixel/index values above to a different cell**,
just edit the literal `if` conditions at each trace site — they're
intentionally narrow (avoid flooding output across ~180k cells/tens of
thousands of particles) rather than parameterized via env var. If you need
a different pixel, grep for `PXTRACE`/`CORNERTRACE`/`SRCTRACE` and change
the hardcoded indices, matching a divergent cell from a fresh pixel diff.

## Fixes kept from this investigation

Three protective fixes that turned out not to explain this specific
divergence, but are real, legitimate improvements worth keeping regardless
(same class of bug as the particle-count fix from earlier this session):

1. **`RoundForInterpStability`/`roundForInterpStability`** — protects
   `GetFluxAtPoint`'s corner-selection truncation against the general class
   of float32-boundary-flip bug.
2. **`PhysCoordTransform` fdlibm fix** — closes a real gap in the
   project's fdlibm-forcing coverage; `BootstrapSampler/` was simply missed
   by the earlier sweep.
3. **`RoundForBinStability`/`roundForBinStability`** — same protective
   technique applied to per-particle nearest-pixel binning.

One fix that **is** confirmed to matter, with a measured, reproducible
effect on the actual divergence (see hypothesis 5 above):

4. **`gasdev` cache-reset fix** — moved `iset`/`gset` from function-local
   `SAVE` variables to module-level (`src/StandardMath/random.f`), added
   `ResetGasdevCache()`, wired into `MaybeOverrideIdum`
   (`src/Outputs/FitOutput.f`). This specifically fixes the diagnostic
   override mechanism (`TRACE_OVERRIDE_IDUM`) to do a *true* full RNG
   reset — it does not change behavior for any normal (non-diagnostic) run,
   since `iset`/`gset` are still functionally the same persistent state,
   just reachable from outside `gasdev` now. **Whether an analogous
   mismatch (not from this specific override mechanism, but from the real
   bootstrap pipeline's own gasdev-call-count history) explains any part of
   the original real-dispatch divergence is still an open question** — this
   fix only touches the diagnostic path; the real bootstrap resampling
   flow was never confirmed to have a matching cache-desync issue,
   only proven capable of one in a controlled, artificial test.

**Audit for other state-preservation bugs of the same class (requested
explicitly, completed this session, came up empty):** every Fortran
`SAVE` statement in `src/` and every module-level `let`/`var` in
`js/src/` was enumerated and individually checked (see "Ruled out" list
above for the specifics). No second instance of gasdev's
reset-timing-asymmetry pattern — a cache whose correctness depends on
call order, initialized via a "reset on negative seed" convention that a
downstream caller's earlier RNG call can silently skip — was found
anywhere else in the codebase. If a similar bug is suspected elsewhere,
the pattern to grep for is: (1) a `SAVE`/module-level variable, (2) whose
reset is gated on a condition that's only true on the "true first" call
of a sequence, (3) where another stateful call can run first and disturb
that condition before the guarded code ever sees it.

## How to reproduce / continue this investigation

Isolated single-bootstrap-realization test, scratchpad-based (doesn't
touch the tracked test-fixture directories):

```bash
# Fortran side — nBootstraps=1, same seed, cdens defaults to 400 (Inputs/SingleGalaxyTestFittingOptions_Base.txt)
cat > /tmp/fortran_config.py << 'EOF'
CubeName="/Users/dandesjardins/DCP/3KIDNAS/3KIDNASTests/TestData/WALLABY_Test_sources/WALLABY_J103538-484832/WALLABY_J103538-484832_VelCube.fits"
MaskName="/Users/dandesjardins/DCP/3KIDNAS/3KIDNASTests/TestData/WALLABY_Test_sources/WALLABY_J103538-484832/WALLABY_J103538-484832_mask.fits"
ObjName="WALLABY_J103538-484832"
TargFolder="/tmp/fortran_out/"
PA_Estimate= 241.319
Inc_Estimate=89.00
nBootstraps= 1
nProcessors_Bootstraps=1
UseDCP=0
BootstrapSeed=42
EOF
cd /tmp && WRKP_TRACE_DEBUG=1 python3 /Users/dandesjardins/DCP/3KIDNAS/WRKP_GalaxyFitDriver.py fortran_config.py > fortran_run.log 2>&1
# The resampled cube appears at TargFolder/WALLABY_J103538-484832/BootstrapCubes/WALLABY_J103538-484832_Bootstrap_0.fits
# and gets cleaned up shortly after — poll for it and copy it out immediately:
#   for i in $(seq 1 90); do f=$(find /tmp/fortran_out/.../BootstrapCubes -iname "*.fits"); [ -n "$f" ] && cp "$f" /tmp/caught.fits && break; sleep 0.3; done

# JS side — matched seed, --local (no real dispatch needed/wanted for this)
cd /Users/dandesjardins/DCP/3KIDNAS
node js/run-galaxy-fit-cli.js initial-fit --cube ... --mask ... --pa 241.319 --inc 89.00 \
  --objName WALLABY_J103538-484832 --local --out /tmp/js_initial_fit.json
TRACE_DUMP_RESAMPLE_PATH=/tmp/js_resampled.fits node js/run-galaxy-fit-cli.js bootstrap \
  --cube ... --initialFitResult /tmp/js_initial_fit.json --nBootstraps 1 --seed 42 \
  --objName WALLABY_J103538-484832 --local --out /tmp/js_bs.json
```

Then diff with `astropy`/`numpy` (see the pattern used throughout this
investigation — `fits.open(path)[0].data.astype(np.float64)`, compare
element-wise, look at `max`/`mean`/count-above-threshold, and specifically
inspect the tracked corner region `x∈{14,15}, y∈{23,24}, ch∈{89,90}` for a
quick sanity check against the numbers recorded in this doc).

For the `TRACE_OVERRIDE_PARAMS`/`TRACE_OVERRIDE_IDUM` controlled tests,
grab Fortran's `CONVERGED_VECTOR` line from its trace log (15 numbers, drop
the trailing chi2) and pass it as a JSON array via `TRACE_OVERRIDE_PARAMS`
to JS's `initial-fit --local` run; add `TRACE_OVERRIDE_IDUM=<any integer>`
to both the Fortran run and the JS run to force a matched RNG state.

## Test galaxy / config used throughout

- Galaxy: `WALLABY_J103538-484832`
- `PA_Estimate=241.319`, `Inc_Estimate=89.00`
- `BootstrapSeed=42`, realization 0
- `cdens=400` (the validated, matched-on-both-platforms default;
  `cmode=0`) — see the separate, resolved investigation in
  `UPSTREAM_SYNC.md` for why this differs from upstream's new
  particle-count formula, which remains deliberately dormant.
- Reference converged parameter vector (Fortran, this config, this seed):
  `[21.9886761, 19.1209297, 1.46704519, 5.70299435, 5740.65430,
  47.9252739, 107.418015, 118.114029, 149.004059, 124.014374,
  1.65046062E-04, 2.95489468E-03, 3.70019116E-03, 7.05593266E-04,
  1.64893863E-05]` — chi2 179994.469.
