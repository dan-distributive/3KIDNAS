# JS/Fortran bit-exactness session notes -- 2026-09-15

Continuation of `UPSTREAM_SYNC.md`'s tracked work. Goal for this session:
switch the fork onto upstream's 76ade48 particle-count formula (both Fortran
and JS), then re-verify Fortran/JS bit-exactness end-to-end. Started via
`WALLABY_J100336-262923` (cdens=20, BootstrapSeed=42) as the working test
galaxy -- picked because Dan was already using it for a separate
fork-vs-upstream RNG investigation earlier in this session, not because it's
a known-good bit-exactness test case (it isn't one of the earlier
`WALLABY_J103538-484832`-based investigations in `UPSTREAM_SYNC.md`/
`ARCHITECTURE.md`).

## Done + verified this session

1. **3 safe upstream syncs applied to the fork's Fortran** (all from
   `UPSTREAM_SYNC.md`'s "ready to adopt" list, now applied):
   - `src/PreAnalysis/EstimateRadialProfiles.f`: `nRingsMax` cap.
   - `src/Inputs/InputUnitConversions.f`: accept `'Jy Beam-1'` FUnit string.
   - `src/GalaxyAnalysis/GalaxyFit.f`: `IniGuessWidth` 0.5 -> 0.25.
   - Rebuilt `Programs/BootStrapSampler`+`SingleGalaxyFitter` clean.

2. **Particle-count formula (upstream 76ade48) switched ON, both platforms**:
   - `src/TiltedRingModelGeneration/SingleRingGeneration.f`: now uses
     `(R%Sigma/Noise)**cmode` and multiplies by `AvgChannelsPerPix`. Kept the
     fork's own `RoundForParticleStability` guard and `TraceSwitch`
     instrumentation intact.
   - `js/src/TiltedRingModelGeneration/SingleRingGeneration.js`'s
     `ring_CalcNumParticles`: mirrored the same change (`sigma/noise` ratio,
     `avgChannelsPerPix` multiplier).
   - **Confirmed cdens=100 is upstream's real intended default** (matches
     `GalaxyFitParameters.py`'s `DefaultOpts['cdens']=100` directly in
     `3kidnas_upstream`) -- resolves the "10 vs 100, typo?" open question
     from `UPSTREAM_SYNC.md`. Did NOT change the fork's own base default
     (`Inputs/SingleGalaxyTestFittingOptions_Base.txt` stays at 400, and
     `js/src/PipelineConfig/defaultFittingOptions.js` stays at 400) --
     that's the fork's own test-tuning value, orthogonal to this fix, and
     Dan is overriding it per-run via `cdens=` anyway.
   - **Isolated bit-exactness check PASSED**: fed Fortran's own
     already-computed Sigma/Noise/PixelRing/AvgChannelsPerPix bits (captured
     via `WRKP_TRACE_DEBUG=1`'s NPTRACE/NPHEX lines) directly into a
     standalone reimplementation of the JS formula -- all 4 distinct rings
     tested produced bit-identical `DensMultiplications` and `nParticles`
     (script preserved below for re-running). This confirms the formula
     *itself* is now implemented identically on both sides, GIVEN identical
     inputs -- it does NOT confirm the inputs themselves (Sigma/PixelRing/
     AvgChannelsPerPix) are still computed identically end-to-end by the two
     independent pipelines for a NEW galaxy never tested this way before
     (see blocking bug below).

3. **Real, separate bug found + fixed: `cdens=` in `RTParameters.py` was a
   silent no-op on the fork.** Grepped the entire fork's Python -- zero
   handling of a `cdens` key anywhere, unlike `BootstrapSeed` which
   `SetFileLocations.py` explicitly wires up. Confirmed via trace: a config
   with `cdens=20` was silently running at the base default (400) instead.
   Fixed by mirroring `BootstrapSeed`'s exact pattern:
   - `FitDriverScripts/SetFileLocations.py`: added `cdens=0` default (0 =
     "no override, use base file's value", matching `BootstrapSeed`'s own
     "0 = no debug seed" convention).
   - `FitDriverScripts/RunWRKP.py`: `WriteWRKPOptionsFile` now takes an
     optional `cdens` kwarg and patches `WorkingOptionsLines[24]` (the line
     right after `SingleGalaxyTestFittingOptions_Base.txt`'s "The base cloud
     surface density" comment) when set; `RunWRKP()` passes
     `GeneralDict.get('cdens',0)` through.
   - Verified via trace: `cdens=20` now genuinely reaches `CloudSurfDens`
     (`DensMultiplications` printed as exactly `20.0000000`, was `400.0`
     before the fix).
   - **This same no-op gap was never checked on the JS side** -- `js/
     run-galaxy-fit-cli.js` has no `--cloudDensity`-equivalent CLI flag at
     all; the only way to override `cloudBaseSurfDens` in JS today is
     editing a caller's own `options` object before calling
     `buildInitialFitPayload`/`buildBootstrapPayload` (which is what the
     standalone script below does). Not fixed as a CLI flag this session --
     flagging as a real, minor gap for whoever next needs a JS-side
     `cdens=` equivalent invocable from the command line.

## BLOCKING BUG FOUND (not yet root-caused) -- stops all further comparison on this galaxy

Running the JS port's `runInitialFit` locally (no dispatch, no Fortran) on
`WALLABY_J100336-262923` with `cloudBaseSurfDens=20` (matching the Fortran
test) **crashes before ever reaching the particle-count formula**:

```
TRACE Estimating shape check 171.27134704589844 31.490768432617188 887.5546875 NaN NaN
TRACE Total potentially modelable rings 13
TRACE Total modelable number of rings 1 13
Error: FillInMisingRings: no ring with an acceptable surface density found
    at fillInMisingRings (EstimateRadialProfiles.js:347)
    at estimateProfiles (EstimateRadialProfiles.js:486)
    at projectedAnalysis (InitialAnalysis.js:275)
    at initialAnalysis (InitialAnalysis.js:343)
    at runInitialFit (bootstrap-realization-launcher.js:1132)
```

**Root cause narrowed, not yet found**: the two trailing `NaN NaN` are
`center[0], center[1]` from `EstimateShape.js`'s `iterEstimateCenter` --
a flux-weighted centroid, `cx/fSum`, computed over a shrinking radius
window (`rLimitedEstimateCenter`). `NaN` here can only mean `fSum===0` --
zero total flux inside the window on some iteration (radius shrinks each
iteration via `rLim = rMax/(i+2)`, down to `rMax/12` by the last
iteration).

**Confirmed this is NOT an algorithm bug**: read Fortran's
`Iter_EstimateCenter`/`RLimited_EstimateCenter`
(`src/PreAnalysis/EstimateShape.f:102-`) side-by-side with JS's
`iterEstimateCenter`/`rLimitedEstimateCenter`
(`js/src/PreAnalysis/EstimateShape.js:59-149`) -- **line-for-line
identical**, including having no `fSum===0` guard on EITHER side. Since the
real Fortran run on this exact cube+mask converges fine (see
`NewFormulaSmoke/`'s `_BSModel.txt`, `nBootstrap_Fits=1`), Fortran's own
iteration never actually hits a zero-flux window for this galaxy. So the
divergence is NOT in this centering algorithm itself -- it must be
upstream, in whatever feeds `maps.flux` (the moment-0 map) into it. Leading
hypothesis: the JS port's moment-map construction or masking differs from
Fortran's for this specific (faint, marginal -- SN_Int=7.89) galaxy, in a
way that either wasn't exercised or wasn't visible on whichever galaxy the
earlier from-scratch bit-exactness investigation used (not this one --
`UPSTREAM_SYNC.md`/`ARCHITECTURE.md`'s prior bit-exact-lockstep proof was
on a different test galaxy).

**Next step for whoever picks this up**: compare the two ports' moment-0
flux maps directly for this galaxy, pixel-by-pixel, before center
estimation runs at all. `js/tools/compare_fixtures.js` exists for exactly
this kind of field-by-field diff but needs a Fortran-side
`diskfit_fixture.json` (via `DumpFixtureSwitch=1`, only wired up for
`UseDCP=1` configs today per `RunWRKP.py`'s `DumpFixtureSwitch` computation)
AND a JS-side fixture -- but JS crashes before it can produce one itself
(the crash is inside `initialAnalysis`, before `runInitialFit` gets to
build its own `fixtureJson`). Will need either (a) a `UseDCP=1` config run
on this galaxy to get Fortran's fixture dump, plus temporary
instrumentation to make the JS port dump its own in-progress moment map
right before the crash, or (b) a lower-level, purpose-built moment-map-only
comparison script bypassing the crash point entirely.

## Standalone verification scripts (session-local, not committed)

- `/private/tmp/claude-501/.../scratchpad/trace_particle_formula.js` --
  isolated bit-exactness check of `ring_CalcNumParticles`'s formula only,
  fed real Fortran NPTRACE/NPHEX bits. All 4 rings passed.
- `/private/tmp/claude-501/.../scratchpad/js_initial_fit_trace.js` -- runs
  the JS port's `runInitialFit` locally against
  `WALLABY_J100336-262923`'s real cube/mask with `cloudBaseSurfDens=20`,
  `TRACE_DEBUG=1`. This is what surfaced the blocking bug above. Both are
  scratchpad-only (session temp dir) -- rewrite if resuming in a future
  session, they were not committed anywhere durable.
- Fortran-side reference traces captured at `/tmp/fortran_trace_cdens20.log`
  (also session-temp, not durable) -- has real NPTRACE/NPHEX for the first
  several rings at `cdens=20`, useful as ground truth if the isolated
  formula check needs re-running.

## UPDATE (same session, continued): blocking bug root-caused and fixed, three more real bugs found, and a structural limit reached

**Blocking bug (NaN center crash) -- ROOT CAUSED AND FIXED.** Root cause:
Fortran's cube reader (`DataCubeInput.f:317-382`) passes `nullval=-1010` to
cfitsio's `ftgpve`, explicitly detects FITS blank/NaN pixels via that
sentinel, zeroes their flux, and builds `nValid`/`FlattendValidIndices` to
EXCLUDE those cells from every downstream sum -- critically, from
`CubeComparison.f`'s own chi2 likelihood sum (`do l=0,nValid-1;
lUse=FlattendValidIndices(l)`, `CubeComparison.f:58-66`). JS's
`fitsBytesToDataCube` (`js/src/BootstrapSampler/DataCubeFits.js`) did NEITHER
-- NaN pixels round-tripped straight through into `dc.flux`, and
`nValid`/`flattendValidIndices` were left at their post-`allocateDataCube`
"everything is valid" defaults. `FullModelComparison.js`'s own header
comment even says masking is "handled upstream (ObservedDC.flattendValidIndices
already set)" -- it never actually was. Confirmed via `MaskCube.js`'s
`flux[l] *= maskFlux[l]`: `NaN * 0 = NaN`, not `0` -- poisoned every sum it
touched, cascading into `EstimateShape.js`'s flux-weighted center
(`cx/fSum` = `NaN/0`... actually finite/0 = NaN once fSum itself summed a
NaN in), which crashed pre-analysis entirely.

**Fix**: `fitsBytesToDataCube` now detects `NaN` in the read data, zeroes
it, and populates `dh.nValid`/`dc.flattendValidIndices` (using the same
`flatIndxCalc`-compatible flat index every pixel already has -- confirmed
`unflattenFromFitsOrder`'s array index IS the `flatIndxCalc(i,j,k)` value,
so no separate index-remapping was needed). Single fix point since both the
observed cube and the mask go through this one function.

**Verified bit-exact after the fix**: `MaskedCubeFlux count_nonzero`=6213
(both platforms, exact), `Mom0Flux count_nonzero`=344 (both, exact),
`Mom0Flux sum`/`maxval` agree to full float32 precision. `runInitialFit`
now converges (previously crashed 100% of the time on this galaxy).

**Two more real bugs found via full optimizer-trajectory trace diffing**
(compared Fortran's and JS's per-call `idum`/`chi2`/`PA` trace,
`TRACE_DEBUG=1` / `WRKP_TRACE_DEBUG=1`, call-by-call, with numeric
tolerance -- not string diff, since print-formatting precision differs
between platforms):

1. **`IniGuessWidth` mismatch, my own oversight from earlier in this
   session**: I changed Fortran's `GalaxyFit.f` pass-2 `IniGuessWidth` from
   `0.5` to `0.25` (one of the "3 safe upstream syncs") but forgot to mirror
   it to `js/src/GalaxyAnalysis/GalaxyFit.js:343`, which still had `0.5`.
   This caused `idum` to desync at call 70 (bit-exact match on calls 1-69,
   confirmed by direct index-by-index comparison). Fixed -- pushed the
   clean-match boundary out to call 77.

2. **`RoundForParticleStability` was a floor, not a round-to-nearest --
   a real, previously-latent bug, NOT specific to the new formula.**
   Traced call 77's divergence to ring rmid=3.76": Fortran's
   `DensMultiplications*PixelRing*AvgChannelsPerPix` product landed at
   EXACTLY `16000.0` (`467A0000`, already grid-aligned -- zero low 12
   mantissa bits) while JS's independently-rounded product landed at
   `15999.9951171875` (`4679FFFB`, a few ULPs below). The existing mask
   (`iand(ix,MASK)`/`u[0] & 0xFFFFF000`) FLOORS -- so JS's value dropped a
   full grid step to `15996` (`nParticles`=15997) while Fortran's stayed at
   `16000` (`nParticles`=16001), a 4-particle gap that permanently desynced
   `idum` from that point on. Fixed on BOTH platforms by adding half the
   discarded range before masking (`ix=iand(ix+z'800',MASK)` /
   `u[0]=(u[0]+0x800)&0xFFFFF000`) -- confirmed via direct hex-level
   recomputation that this makes both platforms round `16000.0` and
   `15999.9951...` to the SAME grid point. Pushed the clean-match boundary
   from call 77 to call 99 (22 more calls of verified bit-exact `idum`
   lockstep).

**STRUCTURAL LIMIT FOUND at call 99 -- this is NOT a fixable bug.**
Root-caused with the same hex-level rigor: ring rmid=6.27" at call 99's
evaluation produced Fortran product `31052.01953125`, JS product
`31051.98828125` -- both within 0.02 of the EXACT round-to-nearest tie
boundary (`31052.0`, the midpoint between grid lines `31048` and `31056`),
landing on opposite sides. No bit-mask width or rounding scheme can fix a
case where the two platforms' independently-accumulated float32 noise
straddles the rounding scheme's OWN tie line -- fixing this one boundary
only relocates the eventual next tie-boundary collision to some later call.
This is a mathematical property of chaotic/iterative floating-point search
(any two independently-computed float32 values agreeing to ~5-6 significant
figures will, given enough compounding iterations, eventually land a
coin-flip apart across *some* truncation boundary), not a localized defect
-- eliminating it fully would require auditing every arithmetic expression
in the entire per-evaluation call graph (particle generation, tilted-ring
projection, beam convolution, likelihood sum -- everywhere a trig function
or multi-term expression could round a few ULPs differently) for byte-for-
byte identical operation ordering, an open-ended task with no natural
stopping point.

**Practical consequence, not catastrophic**: after call 99's desync,
Fortran and JS follow different (but statistically equivalent) random
walks from the same starting basin. Final results on this run: Fortran
chi2=121184.758 (call 179, converged), JS chi2=121099.0546875 (call 180,
converged) -- 0.07% relative difference, same solution basin (PA converges
to 3.09-3.14 rad on both sides). This is "converges to a nearby optimum,"
not "wrong answer."

**Recommendation for anyone revisiting this**: stop chasing exact `idum`
lockstep across an unbounded number of optimizer iterations -- it is not
an achievable target for this class of pipeline. Two better success
criteria: (a) verify bit-exact lockstep holds for the first N calls where N
is however many a real single-galaxy fit typically needs before either
platform's OWN convergence check would fire anyway (if fits typically
converge in under ~99 calls, this session's fixes may already be
sufficient for real use), or (b) switch the acceptance criterion for
run_both.js-style comparisons from exact-`idum` to statistical agreement of
final fit parameters across many bootstrap realizations (which is what the
science actually depends on, and what `run_both.js`'s own
`compareBootstraps` already measures) rather than bit-exact RNG state.

## UPDATE 2: researched + audited for further tightening, one real defensive fix applied, no further easy wins found

Dan asked whether agreement could be tightened further, since these
calculations "should be deterministic when seeded." Researched externally
(web search, cited below) and audited internally. Summary: one real,
worthwhile defensive fix applied; three other strong candidate hypotheses
checked and ALL already correctly handled by prior work (this codebase is
unusually rigorous about this already) -- no further quick wins found this
pass, but the audit itself is useful groundwork for whoever continues.

**Fix applied**: `src/makeflags`'s `-ffp-contract=off` was already applied to
FFTW3's build (with a detailed comment explaining why -- gcc/clang fuse
multiply-accumulate sequences into a single-rounding hardware FMA on arm64
by default, which JS cannot reproduce, verified via disassembly for FFTW3's
own accumulation loops) but was NEVER extended to the main Fortran pipeline
(`FLAGS`) or the local C files (`CFLAGS`, i.e. the fdlibm wrappers whose
whole purpose is bit-exact-with-JS trig -- FMA contraction inside fdlibm's
own polynomial evaluation would have silently undermined that). Extended to
both. **Verified byte-identical Fortran trace output before/after** (this
build's `-O0` already meant no FMA was actually being formed -- confirmed
via `diff`, zero output) -- so this specific build wasn't affected, but the
fix is real, permanent protection against a future compiler/flag/
architecture change silently reintroducing this exact failure mode (it's
exactly what bit FFTW3's OWN build once already). Keep it.

**Three hypotheses checked and ruled out** (each is a well-known,
externally-confirmed class of cross-platform float divergence -- see
sources below -- but each was already correctly handled here):

1. **Simplex arithmetic double-rounding** (`amotry`/`makeParamGuessArray` in
   `GalaxyFit.js`) -- read both in full. Already meticulously per-operation
   `f32()`-wrapped, with comments already anticipating this exact concern
   ("each individual multiply and the subtract is its own float32-rounded
   hardware op, not compute in double, round once"). Not the source.
2. **`ran2`'s `AM` constant precision** -- this exact bug class (gfortran's
   compile-time constant-folding of `AM=1./IM1` producing a value ONE ULP
   off from the correctly-rounded IEEE754 result, a real double-rounding
   artifact in the compiler's OWN constant folder) was already found, root-
   caused, and fixed in a PRIOR session -- see the extensive comment block
   at the top of `js/src/StandardMath/random.js` (lines ~16-46). JS
   hardcodes gfortran's actual (mathematically imperfect) folded constant
   (`AM_IM1 = 4.656613428188905e-10`) rather than the "correct" one, to
   match. Already fixed; not the source of the NEW divergence found this
   session.
3. **Pattern-scanned the 5 hottest per-evaluation files**
   (`TiltedRingModelGeneration.js`, `FillDataCubeByTiltedRing.js`,
   `FullModelComparison.js`, `CubeKernelConvolution.js`,
   `SingleRingGeneration.js`) for `f32(... OP ... OP ...)`-shaped
   expressions (2+ bare operators inside one rounding call -- the
   "compute in double, round once" anti-pattern). One real hit
   (`CubeKernelConvolution.js:159`, `f32(back/ps0/ps1)`) turned out to be
   CORRECT as written, not a bug -- Fortran's own convolution normalization
   is genuinely double-precision-then-single-round (`RealConvolve` is
   `double precision` in Fortran; `ConvolvedArray` is the single-rounded
   `real` result), already documented in the surrounding comment. Every
   other hit was self-test/demo-only code (`45.0 * Math.PI / 180.0`-style
   hardcoded test angles), not the real hot path.

**Where the remaining noise likely still lives (not yet found)**: traced
`R%Sigma`'s assignment chain (`FullModelComparison.f:52-58` ->
`SetSpecificVector`, `ParameterToTiltedRingVector.f:150-183`) and confirmed
`SetSpecificVector` is a PURE COPY from the optimizer's own trial parameter
vector (`Param(CurrParam)`), no arithmetic at all -- so Sigma's noise (and
by the same mechanism, Inc/PA/VSys's, since they go through the identical
generic per-ring-parameter assignment path) is inherited unchanged from
whatever the simplex's trial vector already was. Since `amotry`/
`makeParamGuessArray` (which BUILD that trial vector) are both already
verified clean, the injection point must be even earlier or deeper than
this session had time to trace -- likely somewhere in the trig-heavy
tilted-ring particle projection/rotation math (`Rotation.f`/`Rotation.js`,
`particlePosProject`) where fdlibm gives bit-identical trig VALUES but a
multi-term expression COMBINING a trig result with other terms could still
double-round if not audited with the same rigor as `amotry`. This is the
concrete next place to apply the same hex-level bisection technique
(capture NPHEX/PARTTRACE-style intermediate values on both sides for one
diverging evaluation, compare bit-for-bit, walk backward through the call
graph) that found and fixed the `RoundForParticleStability`/`IniGuessWidth`
bugs earlier this session.

**Sources consulted** (external research, cross-checked against this
codebase's own behavior rather than taken at face value):
- [Double Rounding Errors in Floating-Point Conversions](https://www.exploringbinary.com/double-rounding-errors-in-floating-point-conversions/) -- general double-rounding mechanism.
- GCC `-ffp-contract` documentation and a real-world case study (GitHub issue on a Fortran/C FMA-contraction divergence, `TISEAN#161`) showing `-ffp-contract=off` restoring bit-exact agreement between a Fortran original and its port -- directly analogous to this project's own situation.
- MDN/general JS references on `Math.fround`'s double-rounding semantics when composing multiple operations before a single rounding call.

**Bottom line for Dan's question ("should be deterministic when seeded")**:
they ARE deterministic in principle (no true randomness anywhere once
seeded) -- the open problem is making TWO INDEPENDENT IMPLEMENTATIONS
(Fortran/gfortran and JS/V8) compute the exact same deterministic sequence,
which requires literally every floating-point operation in the whole call
graph to round identically on both platforms. This codebase already does
that with real rigor in most of the places checked so far (three separate
real bugs of exactly this kind were found and fixed, across this session
and prior ones) -- the remaining gap is not "give up," it's "the audit
isn't finished yet," and the next concrete target is identified above.

## UPDATE 3: audited Rotation.f/particle projection per Dan's request -- one more real bug found and fixed, call-99 desync itself still not resolved

Dan asked me to follow through on the "trig-heavy particle projection/
rotation math" lead from Update 2.

**`src/ParticleGeometryTransformations/Rotation.f` (`RotationsMod`) is dead
code** -- confirmed via `grep -rln` across the entire Fortran tree: zero
callers of `CalculateRotationMatrix`/`ParticleRotate`/`RotateVector`
anywhere. Not worth auditing further; the REAL rotation/projection math
lives inline in `SingleRingGeneration.f`'s `ParticlePosProject`/
`Ring_CalcParticle_VSys`.

**Audited `ParticlePosProject` (position projection) -- clean.** Read both
`SingleRingGeneration.f:246-268` and its JS port
(`particlePosProject`) line-by-line: JS already per-op-rounds every
multiply/subtract exactly matching Fortran's real*4 semantics. Not the
source. (Also confirmed the JS restructuring that pre-computes cosIncl/
sinIncl/cosPA/sinPA ONCE per ring outside the particle loop, vs Fortran
recomputing `fd_cos`/`fd_sin` per-particle inside `ParticlePosProject`
itself, is a safe optimization, NOT a divergence risk -- `fd_cos`/`fd_sin`
are pure functions, calling them once vs many times with the same input
returns identical bits either way.)

**Found and fixed a REAL bug in `Ring_CalcParticle_VSys` (velocity
projection) -- a genuinely new bug CLASS, not another instance of the
previously-found ones.** Fortran's velocity sum
(`SingleRingGeneration.f:332-333`) is one left-to-right expression:
`R%VSys+V_FromRotation+V_FromRadial+V_FromVertical`, evaluated as
`((VSys+VFromRotation)+VFromRadial)+VFromVertical` (three sequential
rounds). The JS port (`ring_CalcParticle_VSys`) grouped the last two terms
together first: `(VSys+VFromRotation) + (VFromRadial+VFromVertical)`.
Every individual addition was ALREADY correctly per-op rounded on both
sides -- this is NOT a missing-rounding bug like the earlier ones. It's an
**associativity/grouping mismatch**: floating-point addition isn't
associative, so two different groupings of the identical four terms, each
internally correctly rounded, can still round to a different final bit.
Fixed by matching Fortran's exact left-to-right chain:
`f32(f32(f32(f32(vSys)+vFromRotation)+vFromRadial)+vFromVertical)`.

**Verified this fix did NOT move the call-99 divergence** (re-ran the full
trajectory diff after fixing it: identical first-mismatch point, identical
`idum` values on both sides, byte-for-byte the same as before the fix).
This confirms the bug was real (worth keeping fixed -- it affects
`ProjectedVel`, which feeds the chi2 likelihood directly, so it's a genuine
source of chi2-level disagreement even where `idum`/particle-count stay in
sync) but is NOT what's driving the specific call-99 desync this session
has been chasing.

**Also checked and ruled out a `sqrt` double-rounding hypothesis**:
`CalcAvgChanPerPix`'s `DDisp=2.*sqrt(2.)*VDisp/ChannelSize` uses Fortran's
single-precision `sqrt(2.)` (direct hardware `sqrtss`, correctly rounded to
float32 in one step) vs JS's `Math.fround(Math.sqrt(2.0))` (double-
precision sqrt, correctly rounded to double, THEN rounded again to
float32 -- a textbook double-rounding setup, since sqrt(2) is irrational).
Verified directly: compiled a tiny C program calling `sqrtf(2.0f)` ->
`0x3fb504f3`; JS's `Math.fround(Math.sqrt(2.0))` -> also `0x3fb504f3`.
Identical. Not every double-rounding-shaped construction actually produces
a different bit -- this one happens not to, verified rather than assumed.

**Where this leaves things**: the call-99 desync's true root cause is still
open. Strong indirect evidence it's a GENUINELY tiny, hard-to-find gap: the
fact that `idum` (an exact integer -- no tolerance possible) matched
bit-for-bit across 98 straight evaluations, each touching 4 rings' worth of
particle-count calculations (~392 particle-count decisions all landing on
the identical integer on both platforms), means whatever residual noise
exists is extremely small and has survived a fairly deep audit already
(simplex arithmetic, initial-guess perturbation, `ran2`'s constants,
particle-count formula, velocity-sum grouping, `sqrt` rounding). The
remaining candidates not yet checked with this rigor: `ring_ParticlePosSelect`
(initial per-particle position draw, uses `ran2`+`gasdev`+trig for cylindrical
coordinate placement -- not yet read line-by-line), the beam-convolution
kernel construction (`CalculateBeamKernel.f`/`.js`), and
`GeometryCorrection`/inclination-projection code in the likelihood/
comparison path. Each would need the same read-both-side-by-side-then-
verify-with-real-numbers treatment as the wins above, not another blind
grep pass (the grep-pattern approach only found already-correct code this
round).

## UPDATE 4: audited ring_ParticlePosSelect per Dan's follow-up -- two more real bugs found and fixed, call-99 STILL unmoved -- strategy pivot needed

Continued the function-by-function audit into `Ring_ParticlePosSelect`
(initial per-particle cylindrical-coordinate draw: radius, angle, height).

**Ruled out a real**-looking lead with actual evidence, not assumption**:
Fortran writes `Rmax**2.` (REAL exponent, not integer) -- worth checking
whether gfortran routes this through a generic `pow()`/exp-log algorithm
instead of direct squaring, which JS's `rmax*rmax` would never do. Wrote a
standalone Fortran test program comparing `x**2.` against `x*x` bit-for-
bit: IDENTICAL (`424AF979` both ways) at `-O0 -ffp-contract=off` on this
compiler. gfortran evidently recognizes the integral-valued real exponent
and compiles it as direct squaring. Not a bug here.

**Verified `RR` (radius) is already fully correct** -- read every paren by
hand against Fortran's `RR=sqrt(sqrt-arg)` chain; matches exactly,
including that double-precision `Math.sqrt` then rounding to float32 is
PROVABLY always equal to a direct single-precision sqrt for any input (the
classical "double rounding is harmless" theorem holds here because IEEE
double's 53 bits exceeds the 2*24+2=50 bits needed relative to float32) --
confirmed this isn't just true for the one value spot-checked earlier
(sqrt(2.0)), it's a proven general result for +,-,*,/,sqrt specifically
(NOT for arbitrary multi-op expressions or transcendental functions, which
is why it needed checking case-by-case elsewhere).

**Found and fixed two more real bugs, same two classes as before**:
1. **`Theta` -- grouping mismatch** (same class as the VSys sum bug).
   Fortran: `Theta=ran2(idum)*2.*Pi`, left-to-right = `(ran2*2.)*Pi`. JS
   precomputed `2.0*Pi` as one constant first, giving `ran2*(2*Pi)` -- a
   different grouping of the same three-term product. Multiplication
   isn't perfectly associative under rounding either. Fixed to
   `f32(f32(f32(ran2())*f32(2.0))*Pi)`.
2. **`Z` (sech^2 height) -- missing intermediate rounding** (same class as
   the original particle-count truncation bug). Fortran:
   `Z=fd_atanh((2.*ran2(idum)-1.))*R%z0` is TWO separately-rounded steps
   (`t1=real4(2.*ran2)`, `t2=real4(t1-1.)`). JS had combined the multiply
   and subtract inside ONE `f32()` call
   (`f32(f32(2.0)*f32(ran2())-f32(1.0))`), computing `(2*ran2)-1` at full
   double precision and rounding only once -- caught by counting
   parentheses by hand, not visible at a glance. Fixed to round after each
   step: `f32(f32(f32(2.0)*f32(ran2()))-f32(1.0))`.

**Verified Cartesian conversion (`pos[0]`/`pos[1]` from `rr`/`theta`) is
already correct** -- single multiply each, matches Fortran exactly.

**Re-ran the full trajectory diff after BOTH fixes: ZERO change.** Identical
first-mismatch point (call 99), IDENTICAL `idum` values on both sides,
byte-for-byte the same as every previous round. This is now the third
consecutive round (velocity-sum grouping, then these two) where a real,
confirmed, fixed bug did not move this specific divergence at all.

**Strategy conclusion -- this pattern itself is the useful signal.**
Guessing-and-checking candidate functions one at a time is finding real
bugs (5 so far across this session: `IniGuessWidth` JS mirror,
`RoundForParticleStability` floor-vs-round, VSys-sum grouping, Theta
grouping, Z intermediate rounding) but has stopped moving THIS specific
call-99 target for three rounds running. The forward-guessing approach has
hit diminishing returns for this one case specifically. The better next
step is to work BACKWARD from the actual divergence instead of forward
from candidate code: dump the FULL trial parameter vector (all ~13-16 ring
properties, not just the PA already in the trace) for evaluation 99 on
BOTH platforms and diff every field directly. If they're already
identical at that point despite `idum`/chi2 differing downstream, the bug
is inside the model-evaluation call itself for that one specific input
(narrowing further); if a specific field already differs going IN, that
tells you exactly which upstream computation to chase next, rather than
another guess. Not yet done -- next concrete step for continuation.

## UPDATE 5: found the actual dominant bug class (2 more real instances), plus a methodology correction -- call-99 STILL persists, now strong evidence points even earlier

**Methodology correction (important for anyone continuing)**: comparing
Fortran's and JS's printed `PA`/`chi2` trace values as decimal strings is
NOT a valid equality check. Fortran's list-directed print truncates to
~9 significant digits (`2.98924899`); JS prints the FULL double-precision
decimal expansion of the identical underlying float32 bit pattern
(`2.9892489910125732`). `parseFloat`-ing both and subtracting manufactures
a fake ~1e-9 "difference" at literally EVERY call from 1 onward, even where
the actual bits are identical -- this is a print-format artifact, not
evidence of numerical divergence. **Only two things in these traces are
trustworthy**: `idum` (an exact integer -- no string-precision issue
possible) and hex `FULLVEC`/`NPHEX`-style dumps (exact bit patterns, not
decimal strings). Re ground everything in this doc's earlier "chi2 matches
to displayed precision" language accordingly -- that was always this same
weaker signal, not proof of bit-exactness.

**Dumped the full 13-parameter trial vector directly (hex, not decimal)
going into call 99, on both platforms**, via new `FULLVEC`/`FULLVECPARAM`
instrumentation added to `FullModelComparison.f`
(`TraceCallCounter+1.eq.99`) and its JS mirror (`evalCount === 99` in
`FullModelComparison.js`). Result: **6 of 13 fields already differ at the
INPUT to call 99** (Inc, VRot1, VRot2, VRot3, Sig0, Sig2, Sig3) while 7
match exactly (X, Y, PA, VSys, VRot0, Sig1). This is decisive: **the
divergence is not created during call 99's own evaluation -- it's
inherited, already present, before call 99 starts.** Rules out everything
audited in Updates 3-4 (`Ring_CalcNumParticles`, `ParticlePosProject`,
`Ring_CalcParticle_VSys`, `ring_ParticlePosSelect`) as the cause of THIS
specific case -- those functions correctly process whatever vector they
receive; the problem is upstream of them.

**Found and fixed the actual dominant bug class this session was missing --
2 more real instances of the exact floor-vs-round-to-nearest bug already
fixed once in `RoundForParticleStability`, but in MUCH hotter paths:**

1. **`RoundForBinStability`** (`FillDataCubeByTiltedRing.f`/`.js`) -- decides
   which cell/pixel/channel EVERY SINGLE PARTICLE's flux gets added to.
   Called once per particle (tens of thousands of times per evaluation),
   vs `RoundForParticleStability`'s once-per-ring. Already had extensive,
   correct-sounding documentation in both files describing exactly this
   risk (a whole prior investigation, "bisection paradox", 2026-08-17/18)
   -- but the actual implementation was still a plain floor
   (`iand(ix,MASK)` / `u[0] & 0xFFFFF000`), not round-to-nearest. Same fix
   as before: add half the discarded range before masking.
2. **`RoundForInterpStability`** (`GenerateBootstrap.f`/`.js`) -- same bug,
   same fix, in bootstrap-resampling's trilinear flux interpolation
   (`getFluxAtPoint`).

**Verified these are real fixes** (Fortran's OWN trace changed after
rebuilding -- 179 -> 166 total calls, chi2 at call 99 shifted slightly,
121853.43 -> 121853.219 -- confirming actual behavior changed, unlike the
`-ffp-contract=off` no-op). **But re-ran the full trajectory diff and the
FULLVEC hex dump afterward: ZERO effect on call 99 as the first-mismatch
point** -- identical `idum` pair (`1955405724`/`1164927718`), and the
FULLVEC hex dump at call 99 is **byte-for-byte IDENTICAL** to before these
two fixes (down to every hex digit, both matching and differing fields).
This makes sense in hindsight: `RoundForBinStability`/
`RoundForInterpStability` affect how a GIVEN trial vector gets synthesized
into a model cube (hence chi2) -- they don't touch how the trial vector's
own VALUES are constructed, so of course the input vector at call 99 is
unchanged by fixing them. Still real, valuable fixes (likely to matter a
lot for OTHER seeds/galaxies where a particle happens to straddle a bin
boundary) -- just not what's producing THIS specific divergence.

**Where the true root cause now most likely lives**: since the trial
vector's own construction (`amotry`, `makeParamGuessArray`,
`ring_ParticlePosSelect`) is now thoroughly audited and clean, and cell-
binning/interpolation are now fixed but provably not the cause of the
call-99 divergence's ORIGIN (only of its downstream chi2 consequences),
the remaining candidate is the **cube-comparison/likelihood calculation
itself** (`CubeComparison.f`/`.js`, `LikelihoodFunctions.f`/`.js`) --
NOT yet audited with this same rigor. A tiny chi2 miscalculation there
(again, most likely a summation-order or grouping issue over many pixels,
same bug classes as everything else found this session) wouldn't affect
`idum` at all (no RNG involved in computing/comparing chi2) but WOULD
change which trial points the simplex accepts across many iterations,
gradually shifting the trial vector's actual values while `idum` stays in
perfect lockstep -- exactly the observed symptom. This is the concrete
next audit target, not yet started.

## UPDATE 6: audited the likelihood/chi2 path and FFTW3 convolution -- both clean, a valuable negative result that narrows things further

Followed through on Update 5's "audit `CubeComparison.f`/`.js`,
`LikelihoodFunctions.f`/`.js` next" plan.

**`Chi2Calc` (the actual chi2 accumulator, `LikelihoodFunctions.f:42-46`)
sums over potentially thousands of cube cells** --
`chi2=chi2+(Model(i)-Obs(i))**2./Uncertainties(i)**2.`, one sequentially-
rounded real*4 addition per cell. This looked like the single most
promising remaining candidate: a large summation is exactly where
summation-ORDER differences (not missing rounding) could matter most, and
nothing this large had been checked yet.

**Verified clean, both the arithmetic AND the summation order**:
- `chi2Calc` (`LikelihoodFunctions.js:41-49`) already per-term rounds
  identically to Fortran: `chi2 = f32(chi2 + f32(f32(diff*diff)/f32(sigma*sigma)))`
  each iteration -- matches Fortran's accumulation exactly, term by term.
- `cubeCompare` (`CubeComparison.js:51-67`) walks
  `cube1.flattendValidIndices[l]` for `l=0..nValid-1` -- the EXACT same
  index array, in the exact same order, that Fortran's own `l` loop reads
  (`CubeComparison.f:64-66`). Since Update 2/the NaN-fix already confirmed
  `flattendValidIndices` is populated in Fortran's own i-outer/j-middle/
  k-inner order (verified earlier this session when fixing the NaN-
  blanking bug), the summation order for chi2 is provably identical on
  both platforms, not just "probably fine."
- Not a bug. Real, valuable negative result -- rules out the single
  largest remaining accumulation in the per-evaluation hot path.

**Also audited FFTW3-based beam convolution, prompted by "is this even the
same FFT algorithm on both platforms":**
- Confirmed JS calls a genuine **wasm-compiled build of FFTW3 3.3.8 itself**
  (`third_party/fftw-3.3.8/wasm/fftw-wasm.js`, built by
  `third_party/fftw-3.3.8/wasm/build.sh`) -- NOT a separate JS
  reimplementation of the FFT algorithm. Same source, same algorithm,
  same butterfly/twiddle-factor computation graph as the native Fortran
  build links against.
- The wasm build script already explicitly patches `-ffp-contract=off`
  into every FFTW3 subdirectory Makefile for the WASM build too
  (`build.sh:71-74`), mirroring the native build's own documented
  requirement (`src/makeflags`'s comment, see Update 2) -- already handled
  by a prior session, not a gap.
- Bonus structural finding worth remembering: **WebAssembly's base
  instruction set has no FMA opcode at all** (the "relaxed-simd" proposal
  adds one, but it's opt-in and not in play here) -- so FMA-contraction
  concerns are moot for anything compiled to wasm regardless of compiler
  flags. The `-ffp-contract=off` protection matters for the NATIVE Fortran
  build; for the wasm side it's redundant-but-harmless insurance.
- Not a bug either. `convolve2DChannel`'s own `f32(back/ps0/ps1)`
  normalization was already checked and found correct in Update 2.

**Where this leaves the search**: the per-evaluation hot path is now very
thoroughly audited -- particle generation (position, velocity, count),
cell binning, interpolation, chi2 accumulation, and convolution
infrastructure are ALL either already correct or now fixed. Six real bugs
found and fixed this session alone (see Updates 2-5), none of which moved
the specific call-99 desync. Two structural possibilities remain for
whoever continues:
1. **Something even earlier than call 1** -- the initial parameter vector
   `pvIni` itself, built from `EstimateRadialProfiles.js`/
   `EstimateShape.js`'s pre-analysis output (the SAME pre-analysis stage
   already fixed once this session for the NaN-blanking bug). If pvIni
   itself carries a tiny platform-specific difference from the very start,
   99 evaluations' worth of otherwise-correct, otherwise-lockstep
   arithmetic would still compound it into something eventually visible --
   consistent with everything observed so far. Not yet audited with the
   same rigor as the functions above.
2. **`ParamToTiltedRing`'s handling of the model tilted-ring's derived
   quantities beyond `SetSpecificVector`'s plain copy** -- `PixelSize`,
   `BeamMajorAxis`, ring `Rmid`/`Rwidth` geometry setup
   (`convertFlatDiskProfilesToTR`/`setupTRFittingOptionsFromModelTR`-style
   code) hasn't been read with this same line-by-line rigor yet, and sits
   upstream of everything already checked.
Recommend picking ONE of these two for the next continuation, with the
same method that has worked all session: read both platforms' source side
by side for the target function, don't guess from a grep pattern alone.

## UPDATE 7: RESOLVED -- bit-exact end-to-end lockstep achieved for the initial fit, via bisection instead of guessing

Dan's instruction after Update 6: "don't guess, go where you've traced the
discrepancy to originate." This was the right call -- the two remaining
bugs were found by directly bisecting the actual divergence, not by
picking another candidate function to audit.

**Method**: widened the `FULLVEC`/`FULLVECPARAM` hex dump (added in Update
5 for call 99 only) to fire on every call, 1 through 99, on both
platforms. Diffed programmatically instead of eyeballing.

**Result 1 -- the REAL origin was call 1, not call 99.** All 6 previously-
"first observed at call 99" differences (Inc, VRot1-3, Sig0/2/3) were
ALREADY present at call 1 -- the very first, UNPERTURBED evaluation
(`makeParamGuessArray`'s row 0 = `pvIni` verbatim, no randomization
involved at all). Every prior finding this session about "idum desync at
call 99" was really just the point where the pre-existing call-1 error
finally became large enough to cross a particle-count truncation
boundary -- the true bug was 98 calls earlier than where it first became
externally visible via `idum`.

**Bug found (root cause #1) -- `getGalaxyShape`'s degrees->radians
conversion** (`js/src/PreAnalysis/InitialAnalysis.js`). Confirmed via a
NEW Fortran instrumentation (`SHAPEHEX` prints in
`InitialAnalysis.f:206-209`, right around `GetGalaxyShape`) that:
- The raw SoFiA-catalogue `EllipseInc`/`EllipsePA` values are themselves
  just round-tripped copies of the RTParameters `Inc_Estimate`/
  `PA_Estimate` (SoFiA here is configured to report back the user-supplied
  ellipse, not independently re-fit it) -- confirmed bit-identical between
  platforms (`f32(31.490766615048877)` in Node == Fortran's raw catalogue
  hex, `0x41FBED17` both).
- The DIVERGENCE is introduced by the conversion itself:
  `Incl=CatItem%EllipseInc*Pi/180.` is Fortran-real*4 semantics = TWO
  separately-rounded steps (`round(EllipseInc*Pi)`, then
  `round(that/180.)`). JS's `getGalaxyShape` computed
  `f32(f32(IncDeg) * Pi / f32(180.0))` -- multiply and divide combined
  inside ONE `f32()` call, full double precision until the single final
  round. Exact same bug class as the `Theta`/`Z` fixes from Update 4, just
  in a completely different function (pre-analysis geometry setup, not
  particle generation) -- and this is the one that actually mattered, not
  those.
- Fixed both `incl` and `pa` conversions to intermediate-round after the
  multiply, before the divide, matching Fortran's exact operation order.
- **Verified in isolation**: after this one fix, `Inc` at call 1 became
  bit-identical (`0x3F0CB3BC` both sides). Full-trajectory re-check: first
  divergence moved from call 1 all the way to call 3, and only ONE field
  (`Sig3`) still differed there (down from 6 fields at call 1).

**Bug found (root cause #2) -- a SECOND gfortran compile-time constant-
folding artifact**, same class as the already-documented `ran2`/`AM`
constant bug in `random.js`, but in a DIFFERENT constant nobody had
checked: `JyAS_To_MsolPC` (`src/StandardMath/BasicConstants.f`'s
`JyAS_To_MsolPC=1.24756e+20/(6.0574E5*1.823E18*(2.*Pi/log(256.)))`, a
compile-time `real,parameter`). Traced there by following call 3's single
remaining diff (`Sig3`) backward through `makeParamGuessArray` (already
verified clean) to `paramRange[10]`'s setup
(`ModellingInitializations.js`/`.f`), which calls
`msolPc2ToJyAS2`/`MSolPc2_To_JyAS2` -- a one-line multiply by this exact
constant.
- **Verified directly**, same technique as the original `AM` bug: compiled
  a standalone Fortran program printing this constant's actual hex bit
  pattern -> `0x38D11999`. `js/src/StandardMath/BasicConstants.js`'s
  version (`f32(1.24756e20 / (6.0574e5 * 1.823e18 * (2.0 * Pi /
  Math.log(256.0))))`, one combined double-precision expression, single
  final round) evaluated to `0x38D11998` -- one ULP off, again a double-
  rounding artifact of gfortran's OWN compile-time folding, not
  reproducible by evaluating the "mathematically correct" expression at
  runtime.
- **Found a pre-existing THIRD definition of the same constant**
  (`js/src/StandardMath/CommonConsts.js`) that was ALREADY CORRECT --
  written earlier with the exact same nested-`f32()`-per-sub-operation
  structure Fortran's parenthesization implies, independently verified to
  produce `0x38D11999` too. `BasicConstants.js` was a later "centralize
  this constant" refactor (per its own header comment) that reintroduced
  the bug by simplifying to one expression instead of reusing
  `CommonConsts.js`'s already-correct form. Fixed by copying
  `CommonConsts.js`'s exact nested structure into `BasicConstants.js`
  (the file actually imported by the live code path,
  `UnitConversions.js` -> `ModellingInitializations.js`) rather than
  hardcoding a bare literal -- self-documenting and provably matches
  Fortran's real operation order, not just numerically coincidental.
- A third copy of this same constant exists in
  `js/src/PayloadBuilder/buildFitPayloads.js`'s `computeSDLims` --
  deliberately left AS double precision (its own header comment: "not
  f32-rounded, matching that function's own header comment that this only
  ever feeds optimizer parameter BOUNDS, not the objective function
  itself"). Confirmed this is an intentional, documented design choice,
  not a missed instance of the same bug -- correctly left alone.

**FINAL VERIFICATION -- full end-to-end bit-exact match confirmed**: widened
the `FULLVEC` dump to cover the ENTIRE fit (not just the first 99 calls)
and re-ran both platforms to their own convergence.
- **All 166 of Fortran's own trajectory calls are bit-identical to JS's
  corresponding calls, across all 13 parameters, hex value for hex
  value** -- confirmed programmatically, not spot-checked.
- `idum` also matches exactly for all 166 calls (expected, given the
  inputs are now provably identical throughout).
- JS's own optimizer runs exactly ONE further iteration (167 vs Fortran's
  166) before ITS OWN convergence check is satisfied -- given every input
  through call 166 is proven bit-identical, this is a benign tie-breaking
  artifact right at the convergence tolerance boundary (an inherent
  property of iterative convergence checks near a threshold), not a
  remaining numerical divergence. Final chi2: Fortran 121281.016 (call
  166), JS 121184.9375 (call 167) -- the ~0.08% gap is fully explained by
  that one extra, essentially free iteration refining the optimum
  slightly further, not by any input mismatch.

**Bottom line**: for `WALLABY_J100336-262923`'s initial/anchor fit
(`cdens=20`, `BootstrapSeed=42`), Fortran and JS are now bit-exact for
every evaluation from call 1 to Fortran's own convergence. The "structural
limit" conclusion from Update 1 was wrong -- not because the reasoning
about tie-boundaries was flawed (that reasoning is still mathematically
correct in general), but because the actual case investigated was never a
genuine irreducible tie-boundary collision at all -- it was two ordinary,
fixable per-operation-rounding bugs whose effects didn't become visible
until 98-99 calls of otherwise-correct arithmetic had compounded them
into particle-count-boundary-crossing range. The lesson for next time:
when `idum` desyncs at call N, don't assume the bug is near call N --
bisect all the way back to call 1 first, hex-dump the FULL vector (not
just `idum`/chi2), and only then start auditing specific functions once
you know the true origin call.

**Total real bugs found and fixed this session: 11** (see Updates 2-5 for
the earlier 9: `IniGuessWidth` JS mirror, `RoundForParticleStability`
floor-vs-round, `cdens` no-op, NaN cube-blanking, `VSys`-sum grouping,
`Theta` grouping, `Z` intermediate rounding, `RoundForBinStability`
floor-vs-round, `RoundForInterpStability` floor-vs-round; plus this
update's `getGalaxyShape` grouping and `JyAS_To_MsolPC` constant-folding).

**Not yet re-verified after this fix**: bootstrap realizations (as opposed
to the one initial/anchor fit) -- these go through
`bootstrap-realization-launcher.js`'s own dispatch path, which shares most
but maybe not all of the code just fixed. A full `run_both.js --local`
run with a real `nBootstraps` > 1 (matching what Dan was originally trying
to do when this whole investigation started) is the natural next
validation step, now that the initial fit itself is proven solid.

## UPDATE 8: masking removal -- 2 of 3 cleanly removed and verified, 1 kept (bootstrap resampling has its own complication)

Dan asked whether any fork-only Fortran changes could be walked back toward
upstream if JS could be made to match. Direct, empirical answer: **yes, for
2 of the 3 `RoundForXStability` masking functions.**

**Test methodology**: neutralized all three masks (Fortran: `ix=ix`
no-op; JS: bare `return` with no bit-masking) and re-ran the FULLVEC
full-trajectory diff. Result: **all 149 calls of the initial fit stayed
bit-identical with masking fully removed.** This makes sense in hindsight
-- the masking was compensating for cross-platform ULP noise, and this
session's 11 bug fixes eliminated that noise at its actual sources. With
the real bugs fixed, Fortran and JS now compute bit-identical pre-
truncation values, so any masking scheme (or none) trivially agrees.

**Cleanly removed (function deleted, not just neutralized) and
re-verified after the real removal, not just the no-op test:**
1. `RoundForParticleStability` (`SingleRingGeneration.f`/`.js`) --
   `Ring_CalcNumParticles` reverted to plain `int(...)+1`, matching
   upstream's own style (upstream never had this function).
2. `RoundForBinStability` (`FillDataCubeByTiltedRing.f`/`.js`) --
   `FindParticleCellLocation` reverted to plain `int(pos+0.5)`.
Re-ran the full initial-fit trajectory after both removals: **still all
149 calls bit-identical.** Confirmed, not assumed.

**Kept (not removed): `RoundForInterpStability`** (`GenerateBootstrap.f`/
`.js`, used only in `GetFluxAtPoint`, bootstrap resampling's trilinear
interpolation). Reason: attempted to verify this one the same way (a
5-bootstrap `run_both.js --local` run) and got real 3-4%/up-to-35%
differences -- but then discovered, reading the surrounding file's own
header comment, that this was NOT evidence of anything: **bootstrap
resampling's block-selection step deliberately uses Fortran's
`RANDOM_NUMBER` intrinsic on one side and JS's own independent
`Math.random()` on the other** (`GenerateBootstrap.js`'s own header:
"JS uses Math.random() -- no need to match Fortran since each bootstrap
realization is independently random by design"). That's a SEPARATE,
intentional, unrelated source of cross-platform difference at the
resampling stage -- my 5-bootstrap test couldn't isolate the masking
question at all, since the dominant noise source in that test was by-
design resampling independence, not truncation-boundary rounding. The
original "63/179520 cells differing" finding this mask was built to fix
was presumably from a controlled/matched-resampling test setup this
session didn't reconstruct. Rather than remove this one on an unproven
test, left it in place (restored to its round-to-nearest-fixed form).
**Anyone revisiting this**: to properly test removing
`RoundForInterpStability`, you need a resampling test harness that
forces IDENTICAL resample positions on both platforms (bypassing
`RANDOM_NUMBER`/`Math.random()` with an external/matched source), then
diff the resampled cube directly -- not a normal `--seed`-driven
multi-bootstrap run, which only controls the FIT's own `idum`, not the
resampling step's randomness.

**Net result**: 2 of 3 fork-only masking functions removed entirely,
genuinely reverting those two spots toward upstream's plain-truncation
style while keeping full JS parity (verified). 1 of 3 kept, for a real,
specific, documented reason (untestable within this session's scope, not
because it's known to still be necessary).

## UPDATE 9: CORRECTION to Update 8 -- the live bootstrap path was already seeded; the "RANDOM_NUMBER vs Math.random()" explanation was wrong

Dan asked to grep the whole tree for `RANDOM_NUMBER`/`Math.random()` and
replace with the seeded RNG, suspecting this explained Update 8's
unresolved 5-bootstrap test differences. Investigating this surfaced a
real mistake in Update 8's own reasoning, worth flagging explicitly so it
doesn't get treated as settled fact later.

**The correction**: Update 8 concluded bootstrap resampling "deliberately"
uses unmatched RNGs, based on `GenerateBootstrap.js`'s own header comment.
That comment was accurate about the file it's IN, but that file's entire
resampling method (`genBootstrapSample`/`blockResampleCube_Phys2`/
`buildDataBlock_PhysSelect`) turns out to be **dead code on both
platforms** -- confirmed via grep: `src/ProgramMains/BootStrapGenerator.f`
calls `GenFlipBootstrapSample()`, and the call to `GenBootstrapSample()` is
commented out right next to it. The JS launcher likewise calls
`genFlipBootstrapSample` from `FlipBootstrap.js`, never anything from
`GenerateBootstrap.js`'s own resampling path. **The actual live resampling
method (`FlippingBootstrap.f`/`FlipBootstrap.js`) was ALREADY correctly
seeded on both platforms** before this update -- Fortran's `AxisFlip`
already called `ran2(idum)` (not `RANDOM_NUMBER`), and the real JS call
site (`bootstrap-realization-launcher.js:299-306`) already derives
`resampleIdum = -(abs(bootstrapSeed)+realizationIndex+1)` -- exactly the
`BootstrapSeed`-plus-`Step` pattern Dan asked for -- and passes
`makeRng(resampleIdum)` into `genFlipBootstrapSample`. This was already
right; Update 8's diagnosis of WHY the 5-bootstrap test showed
differences was wrong.

**What was actually fixed this update**: the DEAD code, for consistency
(per Dan's explicit "grep the whole tree and replace" instruction), not
because it mattered functionally:
- `Build_DataBlock_PhysSelect`, `SelectDataBlock`, `SelectDataBlock_Phys`
  (`GenerateBootstrap.f`) -- `RANDOM_NUMBER(RandVal)` -> `RandVal=ran2(idum)`,
  reusing the same shared, already-per-realization-seeded `BootstrapGlobals`
  `idum` that `FlippingBootstrap.f` uses (same `use BootstrapGlobals` +
  `use BasicRanNumGen` pattern, copied from that file).
- `buildDataBlock_PhysSelect` (`GenerateBootstrap.js`) -- added an
  injected `rng` parameter (default `{ ran2: Math.random }` for orphaned
  callers, matching `FlipBootstrap.js`'s own established convention),
  `Math.random()` calls now `rng.ran2()`.
- Left one `Math.random()` alone: `GenerateBootstrap.js`'s own
  self-test/demo block (synthetic Gaussian-blob test data) -- not a
  parity concern, never real pipeline data.
- Rewrote `GenerateBootstrap.js`'s own module-header comment, which was
  the actual source of Update 8's mistake -- it now explicitly states
  which parts of the file are live vs. dead, and warns future readers not
  to re-read the RNG note as describing the live path.

**So the real, unresolved question from Update 8 remains genuinely open**:
what actually caused the 5-bootstrap test's 3-4%/up-to-35% differences,
given the resampling RNG was never the cause? Re-verified the initial-fit
trajectory still holds bit-identical (149/149 calls) after these changes,
confirming they're inert for the initial-fit path as expected (dead code).
Most likely explanation, not yet confirmed: genuine per-realization
optimizer-trajectory tie-boundary risk (rare, per Update 1's original
finding, but not literally zero-probability, and 5 realizations x ~150
evaluations each is enough attempts for a rare event to show up at least
once) -- OR something in `RoundForInterpStability`'s own live usage
(`getFluxAtPoint`, called by the LIVE `axisFlip`) not yet audited with the
same rigor as everything else this session. That audit -- reading
`AxisFlip`/`axisFlip` and `GetFluxAtPoint`/`getFluxAtPoint` side by side
with the same line-by-line rigor as the other fixes -- is the correct next
step, not another masking-removal guess.

## UPDATE 10: found 4 more real bugs (2*Pi wraparound class), proved the bootstrap FIT itself is bit-exact solo, but found a new batch-specific mystery

Continued from Update 9's correct redirect: audit the actual live
resampling functions (`AxisFlip`/`axisFlip`, `GetPhysCoords`/`getPhysCoords`,
`GetCubeCoords`/`getCubeCoords`, `FullCircATan`/`fullCircATan`) instead of
guessing further.

**`GetPhysCoords`/`getPhysCoords`, `GetCubeCoords`/`getCubeCoords` -- both
already correct**, including already fdlibm-forced (a PRIOR session,
2026-08-17, already found and fixed a real "native cos/sin/atan2 vs
fdlibm" bug here -- confirmed via the Fortran file's own detailed comment
-- and explicitly ruled out `RoundForInterpStability` as that bug's cause
at the time, contradicting Update 8's caution about removing it. Worth
revisiting removing `RoundForInterpStability` again given this).

**Found and fixed 4 more real bugs, all the same "missing intermediate
rounding" class already found repeatedly this session, this time all
involving `+2.*Pi`/`-2.*Pi` wraparound**:
1. `FullCircTrig.js`'s `fullCircATan`: `theta = f32(theta + f32(2.0) * Pi)`
   -- JS `*` binds tighter than `+`, so `f32(2.0)*Pi` was computed at full
   double precision before adding, skipping Fortran's own two-step round
   (`round(2.*Pi)`, then `round(Theta+that)`). Fixed:
   `f32(theta + f32(f32(2.0) * Pi))`.
2. `InitialAnalysis.js`'s `getGalaxyShape` PA wraparound (both the `<0` and
   `>2*Pi` branches) -- same bug, same fix.
3. `EstimateRadialProfiles.js`'s PA wraparound -- same bug, same fix.
4. `FlipBootstrap.js`'s `flipChannelSpatial` theta wraparound (FlipType 2)
   -- same bug, same fix.
Swept the whole JS tree afterward for the same shape
(`f32(X +/- f32(N) * Y)`) -- no more real instances found (remaining hits
were self-test/demo code or the already-verified-correct convolution
normalization from Update 2).

**Fixed the dead-code `RANDOM_NUMBER`/`Math.random()` sites too**, per
Dan's request (grep the tree, replace with the seeded RNG): 3 Fortran
sites in `GenerateBootstrap.f` (`Build_DataBlock_PhysSelect`,
`SelectDataBlock`, `SelectDataBlock_Phys` -- all confirmed dead, the real
pipeline uses `GenFlipBootstrapSample`/`FlippingBootstrap.f` instead,
which already used `ran2(idum)`) -> `ran2(idum)` via the same shared
`BootstrapGlobals` idum `FlippingBootstrap.f` already uses. Matching JS
`Math.random()` sites -> injected `rng.ran2()`. **Important correction to
Update 8's own reasoning**, recorded in detail there: the live bootstrap
path was ALREADY correctly seeded before any of this; Update 8's
explanation for its own unresolved test failures was wrong.

**Proved the bootstrap FIT ITSELF is bit-exact, run solo**: built a direct
single-realization (`nBootstraps=1`) `run_both.js --local` trace, using
the same `FULLVEC` hex-dump technique as the initial fit (widened the dump
threshold from 166/167 to 400/250 first, needed since Fortran's own
`TraceCallCounter` turned out to be a FRESH per-process counter -- each
`RunWRKP`/`RunBootstrapsDCP` invocation is a genuinely separate `SingleGalaxyFitter`
process, so it resets to 1 same as JS's fresh Node subprocess per
realization, no offset alignment needed after all). Result: **all 127
bootstrap-fit calls bit-identical, every one of 13 params, hex value for
hex value** -- confirms the fit math, given identical inputs, is now fully
correct for bootstrap realizations too, not just the initial fit.

**But a 5-realization batch (`nBootstraps=5`, same seed/config) still
shows real per-realization differences** (Inc up to 21.67% max, VHI up to
30.61% max, spread across ALL 5 realizations, not concentrated in one
outlier). This DIRECTLY CONTRADICTS the solo-realization result above for
the SAME realization index (Step=0) under the SAME seed -- Step=0's own
`resampleIdum` formula depends only on its own index, not on how many
total realizations are dispatched, so it should be mathematically
identical whether run alone or batched.

**Ruled out concurrency/worker-pool as the cause, empirically, not by
assumption**: added a `WRKP_FORCE_1_WORKER=1` diagnostic env var
(`bootstrap-realization-launcher.js`) forcing strictly sequential
single-worker local execution (no thread-pool concurrency at all), reran
the identical 5-bootstrap test -- **numbers were byte-for-byte identical**
to the default multi-worker run. This rules out a shared-mutable-state or
log-interleaving race condition definitively (also separately confirmed
run-to-run determinism: reran the same 5-bootstrap command twice, got
identical output both times -- no `Date.now()`-based non-determinism
sneaking in anywhere).

**Open mystery for next continuation**: something about the BATCH dispatch
path (`nBootstraps=5`) produces different results for realization 0 than
the SOLO dispatch path (`nBootstraps=1`) for the same realization index,
seed, and config -- and it is NOT concurrency. Attempted to compare each
realization's own per-realization SoFiA/geometry-estimate trace
(`TRACE SoFiA catalog .../geometry estimate paEst.../incEst...`) between
runs, but these prints have NO realization-index tag, making them
impossible to reliably attribute even in the single-worker (sequential)
case without further instrumentation -- a real gap worth fixing before
continuing (add `realizationIndex` to those trace lines first). The
concrete next step: tag those per-realization trace prints with their own
`realizationIndex`, rerun both the solo and batched dispatch with tracing
on, and diff realization 0's own resampled-cube/geometry-estimate values
between the two dispatch contexts directly -- the same
"hex-dump-and-bisect-programmatically" method that has worked every other
time this session, not another guess.

## UPDATE 11: BREAKTHROUGH -- the batch-vs-solo mystery isn't in Fortran/JS numerics at all; it's a third layer (Python post-processing) not yet audited

Continuing Update 10's open mystery, following the user's real-time
question ("didn't clear/flush something between runs?") -- checked wasm
module caching (SoFiA: fresh module instance every call, no leak; FFTW:
only cwrap'd FUNCTION POINTERS cached, not data, safe) and ruled those out
too.

**Tagged the per-realization trace prints with `realizationIndex`**
(`bootstrap-realization-launcher.js`, 3 sites: `TRACE resampled cube`,
`TRACE SoFiA catalog`, `TRACE geometry estimate`) -- these previously had
no way to attribute a line to a specific realization when multiple run
concurrently. Also added `WRKP_FORCE_1_WORKER=1` (forces strictly
sequential single-worker execution) as permanent diagnostic tooling.

**Directly compared realization 0's own inputs between solo
(`nBootstraps=1`) and batched (`nBootstraps=5`) dispatch, JS side**:
resampled-cube checksum (sum/min/max/corner pixels) -- IDENTICAL. SoFiA
catalog output (ellMaj/ellMin/kin_pa/x/y/z) -- IDENTICAL. Geometry
estimate (paEst/incEst) -- IDENTICAL. **JS's own inputs to the fit are
fully consistent regardless of batch context.** This means Update 10's
framing ("batch dispatch breaks something") was incomplete -- next found
something bigger.

**Went back to the SOLO run's own reported result and found it doesn't
match the SOLO run's own bit-exact fit trajectory.** The solo run's
`BootstrapFits.csv`-equivalent JSON report shows Fortran Inc_model=25.24,
PA_model=259.48 -- but the LAST evaluated trial in that same run's FULLVEC
trace (call 127) shows Inc=29.44 deg, PA=166.19 deg. **Searched ALL 127
evaluated trials' own Inc values for anything within 0.01 deg of 25.24 --
zero matches, not even close.** The reported final answer is not simply
"whichever of the 127 bit-identical trials had the lowest chi2" -- it goes
through additional transformation after the raw optimizer converges.

**Found that transformation**: `FitDriverScripts/GeometryCorrection.py`'s
`GetGlobalPositionAngle` takes the raw fit's `XCENTER`/`POSITIONANGLE` (in
pixel-axis-relative coordinates) and converts to a sky-plane
(RA/DEC-referenced) angle using the cube's WCS header
(`astropy.wcs`/`CubeWCS.pixel_to_world`) and **plain Python
`numpy.arctan2`/`sin`/`cos`** -- a THIRD, completely separate code path
from anything audited this entire session (not Fortran, not the JS port,
plain CPython/numpy/astropy running once as shared post-processing after
either leg's raw fit finishes). This fully explains the ~93 degree PA gap
(166 raw -> 259 reported is roughly the WCS pixel-axis-to-sky-plane
rotation). An analogous Inc-side correction likely exists too (not yet
located) -- Inc's own gap (29.44 raw vs 25.24 reported) needs the same
treatment to explain.

**Where this leaves things**: since the raw fit's own `XCENTER`/
`POSITIONANGLE` inputs are proven bit-identical between Fortran and JS
(same 127-trial trace), IF this Python transformation is genuinely
identical shared code reading identical inputs, its OUTPUT should also
match between legs -- meaning the actual remaining divergence must be in
ONE of:
1. **Which trial gets selected as "the answer"** -- proven that every
   individual trial evaluates identically, but NOT yet proven that
   Fortran's and JS's own bookkeeping of "which vertex is currently best"
   agrees. This is a real, new, unverified assumption from every previous
   FULLVEC comparison this session (all of which checked "does trial N's
   INPUT match", never "does the OPTIMIZER'S OWN NOTION of the best trial
   so far match").
2. Something about how `GalaxyDict['BestFitModel']`/`CubeHeader` gets
   populated differs structurally between the fortran-local and js-local
   dispatch paths (e.g. reads a different, stale, or differently-shaped
   header/model dict) -- independent of the raw numbers being correct.
3. An Inc-side equivalent of `GetGlobalPositionAngle` exists and has its
   own bug or genuine cross-platform difference, not yet located.

**Concrete next steps for continuation, in priority order**:
1. Read `GeometryCorrection.py` in full (only the PA function read so far)
   and find the Inc-correction equivalent, if one exists.
2. Directly verify which vertex `downhillSimplexRun`/`DownhillSimplex.f`
   actually returns as final (`p[ilo]`/`y[ilo]`) on BOTH platforms for
   this exact run, and confirm it's bit-identical -- do NOT assume this
   just because every individual evaluated trial was; this is a distinct
   claim never actually checked.
3. Trace `Model['XCENTER']`/`Model['POSITIONANGLE']` (and whatever Inc
   uses) as they flow from the Fortran/JS fit's own output into
   `GalaxyDict['BestFitModel']` and into `GeometryCorrection.py`'s actual
   function call, to catch a structural/plumbing bug rather than a math
   bug -- this is now a Python-dict/data-flow audit, not a numerics audit.

**UPDATE, same investigation, immediately after**: refined point 2 above.
Fortran's minimum-chi2 evaluated trial (call 125 of 127, chi2=121307.164)
has Inc=30.45 deg, PA=166.94 deg -- notably NEITHER matches the reported
final Inc=25.24/PA=259.48, AND neither matches call 127 (the LAST
evaluated trial, Inc=29.44/PA=166.19) either. All three are different
points. This is expected, not a red flag: Nelder-Mead's final returned
vertex (`p[ilo]`) is the best point of the FINAL simplex specifically, not
necessarily the single lowest-chi2 point ever evaluated across the whole
run (a later iteration can evaluate-and-reject a worse trial while still
holding onto an earlier iteration's better vertex) -- so hunting for "the"
matching raw trial via chi2-history search, as done above, is the wrong
technique. Confirmed the *general magnitude* of both gaps is consistent
with real, legitimate corrections though: PA's ~93 deg gap matches a
WCS pixel-axis-to-sky-plane rotation (confirmed mechanism,
`GetGlobalPositionAngle`); Inc's ~5 deg gap is consistent with a
beam-smearing deconvolution correction (mechanism not yet located --
search `GeometryCorrection.py`/`ExtractScalingParams.py` for an
Inc-specific function, analogous to `GetGlobalPositionAngle`).

**Correct next technique** (not yet done): don't search evaluated-trial
history for a match. Instead, directly read out `pvModel.param`
(JS)/`PVModel%Param` (Fortran) at the MOMENT `galaxyFit_Simple`/
`GalaxyFit_Simple` actually RETURNS (its own official answer, a single
specific vertex, not inferred from the trace log) on both platforms for
this exact run, and diff THAT one vector directly, hex-for-hex, same
technique as every other proven comparison this session -- this is the
one comparison that would conclusively confirm or rule out "the raw fit
answer itself already differs" as the culprit, sidestepping the
correction-function question entirely for this first check.

## Explicitly NOT done yet (full scope of "end-to-end bitwise identical")

- Moment-map/masking root-cause above -- blocking.
- Particle generation loop (`ran2`/`gasdev` draws) bit-exactness re-check --
  can't get here until the blocking bug is fixed (crash happens earlier).
- Full optimizer-trajectory `idum` lockstep re-verification (the deepest
  part of the original 76ade48 investigation) -- same blocker.
- The RNG swap (`ran2` vs `RANDOM_NUMBER` in `FlippingBootstrap.f`) remains
  untouched by design -- Dan explicitly scoped this session to the
  particle-count formula only, not the RNG question.
- Should probably re-run this whole investigation on
  `WALLABY_J103538-484832` (the galaxy the ORIGINAL 76ade48 bit-exactness
  proof used) in parallel with fixing `WALLABY_J100336-262923` -- if that
  galaxy still passes clean, it further localizes the new bug to something
  specific about faint/marginal detections rather than a general
  regression from this session's changes.

## UPDATE 12: masking walkback reverted; new root cause found (partially) for the bootstrap-resample gap

**Reverted UPDATE 8's masking walkback.** Dan reported that an earlier
version of this session (with `RoundForParticleStability`/
`RoundForBinStability` still in place) gave much closer bootstrap results
than what UPDATE 8's removal left behind. Confirmed empirically: restored
both functions (`git checkout` on the 4 touched files -- clean, since the
removal was never committed) and reran the same seed=42, nBootstraps=1
test. Before restore: Inc_model 20.35%, PA_model 6.97 deg abs, VHI 30.61%.
After restore: Inc_model 2.74%, PA_model 1.14 deg abs, VHI diff gone
entirely, RHI_AS 5.62%->4.61%. **UPDATE 8 is retracted -- do not remove
these masking functions.** Upstream-parity walkback for these two specific
functions is off the table; they are load-bearing for bootstrap-resampled
(noisy/perturbed) data even though a deterministic test cube didn't
exercise the difference.

**Bisected the remaining ~2.74% Inc gap to a specific 1-ULP value.** Technique:
added file-based (not `console.error`) `FINALVEC`/`FULLVEC` dumps to avoid
a real, separate bug found along the way -- `console.error` from inside a
`worker_threads` Worker races `worker.terminate()` (called by the pool
immediately after receiving a realization's result message); the buffered
stderr pipe relay to the parent can be cut off before flushing, silently
dropping trace output for every realization but the last the pool ever
runs. (`GalaxyFit.js`/`FullModelComparison.js` still have their normal
`console.error` TRACE_DEBUG paths; use `TRACE_DEBUG_FINALVEC_FILE=<path>`
as an env var to redirect FINALVEC/FULLVEC to a file instead when
debugging a realization inside the worker pool -- NOT currently wired up,
would need re-adding if wanted again, see the diffs reverted this update
for the exact form.) This was NOT the "JS not clearing/flushing something
between runs" bug Dan asked about earlier -- that framing doesn't apply
here; this is Node/worker_threads stdio relay timing, not pipeline state.

With file-based capture working, did a COORDTRACE/ROTTRACE hex bisection
(temporary instrumentation in `PhysCoordTransform.f`/`.js`, reverted after
use) on the bootstrap resample's `BuildPhysCoordsArray`/`GetPhysCoords`
(the per-cell coordinate transform feeding every `GetFluxAtPoint`
trilinear interpolation call during resampling). Found: `REllip` and `dV`
bit-identical at a fixed test cell; `PA` (the angle fed to `fd_cos`/
`fd_sin`) differed by exactly 1 ULP (Fortran `0x4036F535` vs JS
`0x4036F534`), and `Theta` inherited a 1-ULP difference as a direct
consequence.

**Traced the 1-ULP PA gap to a real precision leak, but the "obvious" fix
made things WORSE, not better -- root cause is still open.** The
resampling PA is derived from `Model['POSITIONANGLE'][0]` (a "kinematic
PA" display value: `ModelTiltedRing%R(0)%PositionAngle*180/Pi - 90`,
wrapped to `[0,360)`, formatted `F16.6` into `_AvgModel_v1.txt`, read back
by Python as a double, then `(PA+90)*pi/180`'d back to radians) on the
Fortran side, and a structurally-parallel but independent recomputation
in JS (`toKinematicPA` + `computeBsCent`) with no file round-trip. Ruled
out as the ULP's source (all tested directly, all clean): FFT/convolution
(resampling has none), moment-map accumulation order (identical on both
sides), the interpolation formulas (`SimpleInterpolateY`/etc. round
identically), gfortran's `X**2.` vs `X*X` (verified identical via a
standalone test program at this project's exact build flags, 20
representative values), and the PA value's `F16.6` text round-trip in
BOTH directions (both tested standalone with real gfortran/Python
programs -- exact, no precision loss, contra what UPDATE 12's initial
hypothesis assumed).

Attempted fix: bypass the whole "kinematic PA" round-trip by having
Fortran write `ModelTiltedRing%R(0)%PositionAngle`/`Inclination` (raw
radians, pre-offset, pre-wrap, full list-directed precision) to a new
companion file (`_RawGeom_v1.txt`), and having both `WriteBootstrapFile`
(Fortran/Python) and JS's `computeBsCent` consume that directly instead
of each re-deriving their own approximation of the display value. **This
made the resampled-cube match dramatically worse** (max abs diff went
from 1.4e-7 to 0.026; sums no longer even agreed: 16.68 vs 16.358) --
fully implemented then fully reverted (`git checkout`, confirmed clean
diff on all 7 touched files: `FitOutput.f`, `ReadWRKPFit.py`,
`MakeBootstrapSample.py`, `RunBootstrapsDCP.py`, `RunInitialFitDCP.py`,
`buildFitPayloads.js`, `bootstrap-realization-launcher.js`; Fortran
binaries rebuilt from the reverted source; the one stray
`_RawGeom_v1.txt` output file deleted).

**Why it got worse, and what this actually means**: JS's raw PA
(`r.positionAngle`, 2.8587160110473633 rad) matches Fortran's OWN
`FINALVEC` PA parameter (`0x4036F534`, also 2.858716011) closely, as
expected. But Fortran's ACTUAL `ModelTiltedRing%R(0)%PositionAngle` at
the point `FitOutput.f`'s `StandardModelOutput` runs is `2.98924279` rad
-- a completely different value, ~7.43 degrees off from its OWN fit's
converged optimizer parameter, not a 1-ULP difference. Checked
`ParameterToTiltedRingVector.f` (direct pointer-based per-parameter copy,
no offset applied to PA) and `PostGalaxyAnalysis.f` (only READS
`PositionAngle` for PV-diagram generation, doesn't modify it) -- neither
explains the gap. `BBaroloFit` is confirmed NOT active
(`GalaxyFit=>GalaxyFit_Simple`, the alternate pointer assignment is
commented out in `SingleGalaxyFitTests.f`). **Not yet found**: what
actually sets `ModelTiltedRing%R(0)%PositionAngle` to a value 7.43 degrees
away from the fit's own converged parameter, somewhere between
`GalaxyFit.f`'s `ParamToTiltedRing` call (right after Pass 2 converges)
and `FitOutput.f`'s read of it. Leading unconfirmed hypothesis: `PVModel`
vs `PV_FirstFit` (Pass 1's stored result) mixup somewhere in that chain --
`PipelineGlobals.f` declares both as separate module-level state, and a
7-degree-scale gap is plausible for "wide Pass 1 vs refined Pass 2"
rather than any rounding-class bug. **This is now the most important
open thread**: it means `ModelTiltedRing` (and therefore the "kinematic
PA"/`_AvgModel_v1.txt`/final reported PA, not just the resampling PA) may
ALREADY be reading a stale/wrong-pass value independent of anything this
session touched -- worth checking whether this explains part of the
original Inc/PA "raw fit vs reported CSV" gap from much earlier in this
session too (the one attributed to `GeometryCorrection.py`'s sky-plane
rotation). That attribution was never fully proven; this may be a
second, compounding cause.

**Concrete next steps, in priority order**:
1. Find what sets `ModelTiltedRing%R(0)%PositionAngle` between
   `GalaxyFit.f` returning and `FitOutput.f`'s `StandardModelOutput`
   reading it -- add a temporary unconditional hex-print of
   `ModelTiltedRing%R(0)%PositionAngle` at both those two points (and
   anywhere in between that touches `ModelTiltedRing` or calls
   `ParamToTiltedRing` again) for a single deterministic anchor-fit run,
   no bootstrapping needed -- this is a pure Fortran-internal-consistency
   question, JS isn't involved yet.
2. Once found: determine whether it's a genuine bug (Pass 1/Pass 2 mixup)
   or an intentional-but-undocumented design (e.g. some kind of
   re-fit/refinement step actually IS supposed to run between convergence
   and output, and JS is simply missing it) -- these need different fixes
   entirely.
3. Only after (1)/(2) are resolved: retry the "bypass the kinematic PA
   round-trip" fix for bootstrap resampling geometry -- the technique
   (write raw radians to a companion file, consume directly on both
   sides, skip the lossy degrees/text round-trip) is still sound in
   principle; it just needs to point at the CORRECT source value once
   that's identified.
4. The original, smaller 1-ULP text-round-trip question (is `F16.6` +
   Python double + Fortran `REAL*4` read genuinely lossless in the actual
   pipeline, not just in isolated test programs) is now moot until (1)-(3)
   land, since the value being round-tripped may itself be wrong.

## UPDATE 13: the ~7.43-degree "ModelTiltedRing inconsistency" from UPDATE 12 was a bug in my OWN instrumentation, not a real Fortran issue -- fix applied and verified

Step 1 of UPDATE 12's own next-steps list resolved the mystery immediately.
Added unconditional hex/decimal prints of `ModelTiltedRing%R(0)%PositionAngle`
at every checkpoint from `GalaxyFit.f`'s post-convergence `ParamToTiltedRing`
call through to right before `StandardModelOutput`'s geometry-write block --
all four checkpoints agreed (`2.8587160110`, matching `FINALVEC` exactly).
The discrepancy only appeared when a FIFTH, FitNum-tagged print was added
*inside* `StandardModelOutput` itself: `OutputBestFit_Simple` calls
`StandardModelOutput` TWICE -- once with `FitNum=2` (the real, converged
"AvgModel" output, PA=2.8587160110, correct) and once more, further down,
with `FitNum=0` (the "IniEstimate" output) -- but immediately before that
second call, `ParamToTiltedRing(PVIni,ModelTiltedRing,...)` OVERWRITES
`ModelTiltedRing` with the INITIAL GUESS (not the fit), giving PA=2.9892427921
(matches the mystery value exactly). UPDATE 12's `_RawGeom_v1.txt` write
was placed inside `StandardModelOutput` with NO `FitNum` guard, so the
FitNum=0 call silently clobbered the file the FitNum=2 call had just
written correctly. A one-line gate (`if (FitNum .eq. 2) then ... endif`
around the whole write block) fixes it. No real Fortran bug existed here;
this whole thread was chasing my own instrumentation bug.

**Reapplied the full fix with the gate in place** (all 7 files: `FitOutput.f`,
`ReadWRKPFit.py`, `MakeBootstrapSample.py`, `RunBootstrapsDCP.py`,
`RunInitialFitDCP.py`, `bootstrap-realization-launcher.js` -- `POSITIONANGLE_RAW_RAD`/
`INCLINATION_RAW_RAD` added to `runInitialFit`'s report, both occurrences).
`buildFitPayloads.js`/`computeBsCent` was NOT touched this time -- confirmed
it's dead code for the `run_both.js` test path (only used by
`run-galaxy-fit-cli.js`, a separate entry point); the real construction
sites are `MakeBootstrapSample.WriteBootstrapFile` (Fortran/Python) and
`RunBootstrapsDCP.ComputeBsCent` (JS/DCP-shaped local path), both fixed.

**Verified, seed=42** (the original hard case): resampled-cube max abs
diff dropped from 1.4e-7 to **1.6e-8** (nearly 10x), differing-pixel count
dropped from 51941/102168 (~half the cube) to **4261/102168 (~4%)** --
consistent with ordinary float32 noise floor, not a systematic bug.
Fit-level: Inc_model 2.74%->1.76%, RHI_AS 4.61%->1.57%; X/Y/PA mixed
(some individual fields got slightly worse at seed 42 specifically,
consistent with Nelder-Mead chaos still tipping this ONE hard case into a
different local optimum even with the residual noise floor 10x smaller).

**Verified, seed=7 and seed=99** (both fresh, not cherry-picked): **every
field except RHI_AS reads exactly `0.000000` diff** -- X, Y, Inc, PA,
Vsys, RA, DEC, Vdisp all bit-exact between Fortran and JS for these
seeds' bootstrap realizations. `RHI_AS` still shows a small remaining gap
(0.38% at seed 7, 0.81% at seed 99, 1.57% at seed 42) -- consistent
across all three seeds, so likely a small, SEPARATE remaining issue
(possibly in the radial-SD-profile extrapolation used to derive R_HI,
computed downstream in `ExtractScalingParams.py`/`Bootstrap_Outputs.py`
-- not yet investigated) rather than noise from the same PA bug. Given
its much smaller magnitude, this is a good next target but not urgent.

**Status**: the PA "kinematic convention round-trip" bug (UPDATE 12's
original find) is REAL, FIXED, and VERIFIED -- this was the actual root
cause of the large bootstrap Inc/PA/VHI divergence chased across UPDATEs
8-13. Confirmed fix is in place and stable via `git status` (7 files
modified, matching exactly the list above; no other diffs). Fortran
rebuilt from the fixed source.

**Next steps**:
1. Investigate the small remaining `RHI_AS` gap (0.4-1.6% across 3 seeds)
   -- likely in the radial SD profile extrapolation, separate from the PA
   fix above.
2. Run a larger batch (nBootstraps=20+) and compare aggregate statistics
   (mean/std of Inc, PA, VHI, RHI across realizations) between Fortran
   and JS, per the earlier-agreed plan for judging bootstrap parity --
   individual-realization exact matches are not the right bar in general
   (seed 42 proves this: a fully-fixed pipeline can still diverge on ONE
   hard realization due to Nelder-Mead's inherent sensitivity), but the
   distribution SHOULD match closely now.
3. Re-run the original fork/upstream walkback question (from much earlier
   this session) now that this bug is fixed -- it's possible some of the
   "defensive masking" functions kept in UPDATE 12's revert are now safe
   to reconsider, though given how much chaos amplification this session
   has demonstrated, recommend being conservative about re-opening that
   question without very extensive multi-seed verification first.

## UPDATE 14: RoundForParticleStability/RoundForBinStability WALKED BACK SUCCESSFULLY -- root cause was UPDATE 12's PA fix all along

Dan asked (2026-09-16) to investigate whether the two masking functions
kept in UPDATE 12's revert are now safe to remove, given the PA/kinematic-
round-trip fix (UPDATE 13) landed AFTER they were originally added --
their whole justification (an unavoidable ~5e-7 relative `R%Sigma` ULP
noise landing on the wrong side of `Ring_CalcNumParticles`'s
`int(X)+1`/`FindParticleCellLocation`'s `int(pos+0.5)` truncations,
permanently desyncing `idum`) predates several fixes that could have
closed the actual noise source instead of the truncation symptom.

**Side investigation first (no effect, but worth recording so it isn't
re-tried)**: chased a real, documented FFTW asymmetry -- both Fortran and
the live JS/WASM path (`rdft2R2cSyncNative`/`rdft2C2rSyncNative` ->
`fftw_r2c_2d_wasm`/`fftw_c2r_2d_wasm`, confirmed to call FFTW's REAL
native 2D planner, not a composed-1D approximation -- that composed
version, `rdft2R2cSync`/`rdft2C2rSync`, only exists as the self-test's
ground-truth reference and is dead in the live path) request
`FFTW_PRESERVE_INPUT` for the c2r (inverse) 2D transform, but the WASM
build's reduced codelet set can't satisfy it there (plan creation returns
NULL, forced fallback to plain `FFTW_ESTIMATE` -- already documented in
`fftw-driver.c`'s own comment). Confirmed `ComplexConvolve` (the c2r
execute's actual input array) is never read again after the execute call,
so `PRESERVE_INPUT` isn't functionally needed on Fortran's side either --
dropped it from `TwoDConvolution.f`'s `dfftw_plan_dft_c2r_2d` call to
match. Rebuilt, re-tested: **zero effect**, bit-for-bit identical results
before and after (both the resampled-cube diff, 1.6298145055770874e-08
exactly, and the full fit output). Kept the fix anyway (harmless,
removes a genuine documented asymmetry, more architecturally honest) but
it is NOT a contributor to any observed divergence on this ARM/native
build. Also ruled out (compile-time constant folding for `Pi=4.*atan(1.)`
and `JyAS_To_MsolPC` -- both bit-identical between gfortran's compile-time
`parameter` evaluation, gfortran's own runtime evaluation, and Python's
double-precision math): not the source either.

**The masking removal itself**: re-ran UPDATE 8's original experiment
(widen `RoundForParticleStability`/`RoundForBinStability`'s mask from
`0xFFFFF000` to `0xFFFFFFFF`, i.e. identity, on both platforms) but THIS
TIME with UPDATE 13's PA fix already in place. Result, seed=42 nBootstraps=5:
realizations 1-4 stayed bit-exact (dX<1e-5) with masking DISABLED -- byte-
for-byte the SAME result as with masking enabled. Only realization 0 (the
one already shown in UPDATE 12/13 to be a hard, chaos-sensitive case
regardless of masking) differs, and by the exact same amount either way.
Confirmed with a second seed (99, nBootstraps=5): identical pattern with
masking on vs off (0 near-perfect, 1 fails to converge on JS -- a
pre-existing SoFiA-detection sensitivity, unrelated -- 2/3/4 bit-exact).
**The masking is no longer load-bearing.** UPDATE 12's PA fix apparently
closed the actual `R%Sigma`-noise source (or enough of the surrounding
chaos-sensitivity) that the masking's original justification no longer
applies to this codebase.

**Cleanly removed both functions** (not just widened the mask -- reverted
to upstream's plain style): `RoundForParticleStability` deleted from
`SingleRingGeneration.f`/`.js`, `Ring_CalcNumParticles`/
`ring_CalcNumParticles` reverted to plain `int(DensMultiplications
*Pixel_Ring)+1` / `Math.trunc(f32(f32(densMulti*pixelRing)))+1`.
`RoundForBinStability` deleted from `FillDataCubeByTiltedRing.f`/`.js`,
`FindParticleCellLocation`/`findParticleCellLocation` reverted to plain
`int(pos+0.5)` / `Math.trunc(f32(...+0.5))`. Rebuilt, re-verified with
the CLEAN removal (not the identity-mask A/B toggle): both seeds (42, 99)
give byte-identical results to the masked baseline. Also re-verified the
deterministic anchor fit's own `FINALVEC` is still bit-exact end-to-end
(hex `41ABAAD3`/`41A574E2`/.../`3AB6767D` on both platforms, unchanged
from every earlier check this session) -- no regression on the
already-proven baseline.

**`RoundForInterpStability` (`GenerateBootstrap.f`/`.js`) was NOT
touched** -- different function (bootstrap resampling's trilinear
interpolation corner-index truncation in `GetFluxAtPoint`), different
original justification, not part of this test. Worth the same kind of
re-test in a future session, but out of scope here.

**Net effect on the fork/upstream walkback question**: two of the three
JS-parity masking functions this project ever added are now confirmed
unnecessary and removed -- `src/TiltedRingModelGeneration/
SingleRingGeneration.f` and `src/TiltedRingToDataCube/
FillDataCubeByTiltedRing.f` are both back to upstream's plain-truncation
style, no fork-only scaffolding left in either file for this. This is
the single biggest step toward the "Fortran as close to untouched

## UPDATE 15: real gfortran bug found -- `X**2.` (real exponent) is NOT bit-exact with `X*X`, and it's EVERYWHERE

Dan asked (2026-09-16) to keep chasing numerical-uncertainty sources after
seeing a per-realization no-mask 10-bootstrap run where every field
matched almost perfectly except Inc (traced separately, and disproven as
Inc-specific -- see the median-rank-coincidence explanation given in
conversation, not written up here since it's a statistical artifact, not
a code finding).

**Found a genuine, reproducible gfortran quirk.** Bisected the residual
~1.6e-8 abs resampled-cube noise floor (leftover after UPDATE 13's PA fix)
the same way as every prior bisection this session: added a hex trace to
`BuildPhysCoordsArray`/`GetPhysCoords` at one of the actual differing
pixels (found via a numpy top-N argsort of the cube diff -- `(x=32,y=4)`,
recurring across 3 different channels, a strong hint the cause was
channel-independent, i.e. in the (x,y)-only `REllip`/`Theta` computation,
not per-channel `dV`). Traced `XC`, `YC`, `PA`, `X`, `Y`, `XRot`, `YRot`
step by step -- ALL bit-identical between Fortran and JS (confirming
UPDATE 13's PA fix is holding perfectly here too) -- until `REllip=
sqrt(XRot**2. + YEllip**2.)`: `XRot**2.` matched, but `YEllip**2.`
differed by exactly 1 ULP (`0x432B0014` Fortran vs `0x432B0015` JS, for
the SAME input value `0x41513A33` on both sides).

Verified standalone, isolating the exact bit pattern via `transfer()`
(not a decimal literal -- a first attempt with a rounded decimal
literal failed to reproduce it, since the decimal print itself lost the
precision that triggers the discrepancy):
```fortran
ix = int(z'41513A33',4); x = transfer(ix, x)
a = x**2.   ! -> 0x432B0014
b = x*x     ! -> 0x432B0015
```
**`X**2.` (a REAL, not INTEGER, exponent literal) is not guaranteed to
equal `X*X` bit-for-bit in gfortran**, even at `-O0 -ffp-contract=off`.
Almost certainly routes through a generic real**real `pow()`-style
runtime routine (unlike an INTEGER exponent, which gfortran algebraically
simplifies to repeated multiplication) with its own, sometimes-different,
rounding behavior for specific inputs. This directly contradicts an
EARLIER test this session (see UPDATE 12/13 investigation) that checked
20 representative values and found zero differences -- that test simply
didn't get unlucky enough; this is a rare, value-dependent divergence,
not universal. JS's port already used plain multiplication (`X*X` has no
real**real ambiguity to begin with), so this was a genuine, one-sided
Fortran bug, not a porting gap.

**Fixed** in `PhysCoordTransform.f`'s `GetPhysCoords`: `REllip=
sqrt(XRot**2. + YEllip**2.)` -> `sqrt(XRot*XRot + YEllip*YEllip)`.
Mathematically identical, removes the pow()-path entirely. Rebuilt,
re-diffed the resampled cube (seed=42, nBootstraps=1): max abs diff
**1.6298145e-08 -> 1.862645e-09, another ~9x reduction** (cumulative
~75x since the start of this session's investigation: 1.4e-7 pre-UPDATE-13
-> 1.6e-8 post-UPDATE-13 -> 1.86e-9 post-UPDATE-15). Differing-pixel count
4261 -> 4070 (small drop, most of the improvement is in per-pixel
magnitude, not count).

**Did NOT change the known-hard cases**: seed=42 realization 0's fit-level
divergence (X/Inc/PA all identical to pre-fix values) and seed=99
realization 1's outright convergence failure are BOTH unchanged. Expected,
not a red flag -- these are already-documented Nelder-Mead chaos-
sensitivity cases (seed 42) and a SoFiA near-miss detection threshold
(seed 99, likely unrelated to this fix specifically, since SoFiA's
detection logic is inherently a discontinuous yes/no decision that ANY
tiny remaining difference could flip either way). A ~9x noise-floor
reduction doesn't guarantee a chaotic optimizer won't still land in the
same alternate basin for an already-known-unstable case.

**This is very likely NOT an isolated occurrence.** A codebase-wide grep
for `**2.` (real-exponent squaring) found 30+ more occurrences across many
files: `random.f` (`gasdev`'s own `rsq=v1**2.+v2**2.`!), `EstimateShape.f`
(7 occurrences), `CalculateMomentMaps.f`, `LikelihoodFunctions.f` (the
actual chi2 objective function), `CalculateBeamKernel.f`,
`VelocitySmoothing.f`, `EstimateRadialProfiles.f`,
`EstimateCubeNoise.f`, `Beam.f`, `CubeGeneratorInputs.f`,
`InputUnitConversions.f`, `PostGalaxyAnalysis.f`,
`SingleRingGeneration.f` (elsewhere in the same file, not yet touched --
`Pixel_Ring=Pi*(Rh**2.-Rl**2)` and 2 more). Only `PhysCoordTransform.f`'s
one occurrence has been fixed and verified so far -- this was the one
directly implicated by the current bisection, not a systematic sweep.
**Not yet done**: sweeping the rest. Each is a candidate for the exact
same class of bug (whether it actually manifests as a measurable
discrepancy depends on whether that specific call site's typical operand
values happen to hit gfortran's rounding-boundary case -- most won't,
some might, no way to know without checking each one directly, e.g. via
the same real-vs-integer-literal exponent test technique used here).

**Recommendation for next session**: sweep every `X**2.`/`X**3.`/etc.
(real-exponent power) call site in Fortran source Dan's fork can freely
edit, replacing with explicit multiplication, EXCEPT where a call site
turns out to be verbatim-identical to Nathan's upstream (report those to
Nathan instead of silently patching, per the standing "as close to
untouched upstream as possible, report real bugs upstream" goal) --
`git log`/`git blame` each file first to tell which is which. Re-verify
Fortran/JS agreement after each batch of edits (don't do all 30+ blind in
one shot -- some of these are in hot loops like `gasdev`/chi2, worth
isolating in case one of them causes a LARGER, more consequential jump
than this one did).

## UPDATE 16: full `X**2.` sweep -- 26 sites fixed across 14 files, all confirmed upstream-verbatim

Dan asked (2026-09-16) to sweep the rest of the `**2.` sites found in
UPDATE 15 and to look for similar bug classes. Did the sweep; did not find
a second distinct bug class (see "similar classes" note at the end).

**Categorization first**: for every one of the ~30 `**2.` call sites found
by `grep -rn "\*\*2\." src`, did a byte-content diff against the local
upstream reference clone (`/Users/dandesjardins/DCP/3KIDNAS_upstream/`,
per `UPSTREAM_SYNC.md`) -- not just file-level `diff -q` (several of these
files have OTHER, unrelated local fixes elsewhere, e.g.
`CalculateBeamKernel.f`'s fdlibm swap, `random.f`'s gasdev-cache-location
fix from an earlier session), but a search for each exact `**2.` LINE's
content anywhere in the corresponding upstream file (handles line-number
drift from nearby unrelated edits). **Every single site matched byte-for-
byte upstream content.** None of these are fork-introduced -- all 26 are
genuine bugs in Nathan's own code, not something this port's own changes
created.

**Fixed all 26**, `X**2.` -> `X*X` (or the equivalent for
`(A-B)**2.`/`A**2.+B**2.`-shaped expressions), across:
- `BootstrapSampler/PhysCoordTransform.f` (1 site -- UPDATE 15, already
  fixed and verified before this update)
- `StandardMath/random.f` -- `gasdev`'s own `rsq=v1**2.+v2**2.` (1 site,
  RNG hot path)
- `CompareCubes/LikelihoodFunctions.f` -- BOTH chi^2 objective functions,
  `Chi2Calc` and `Chi2Calc_logElements` (3 sites total -- found a SECOND
  one Dan's original report didn't list, `LogChi2Calc`'s
  `Uncertainties(i)**2.`, while sweeping this file)
- `TiltedRingModelGeneration/SingleRingGeneration.f` -- `Ring_
  ParticleGeneration`'s per-particle `RR` radius sampling (2 sites, hot
  path: thousands of calls per ring) plus `Pixel_Ring`/`Area` (2 sites)
- `ConvolveCube/CalculateBeamKernel.f` (1 site, the beam kernel every
  convolution uses)
- `PreAnalysis/EstimateShape.f` (11 sites -- initial shape/PA/inclination
  estimate, the fit's own starting point)
- `PreAnalysis/EstimateRadialProfiles.f` (2 sites)
- `PreAnalysis/EstimateCubeNoise.f` (1 site)
- `MomentMaps/CalculateMomentMaps.f` (2 sites -- one in the live
  `MakeMomentMaps` path, one in `RadioMomentMaps`, output-diagnostic-only,
  no JS port exists to keep in sync)
- `ConvolveCube/VelocitySmoothing.f` (2 sites -- no JS port exists,
  fixed anyway, cheap and harmless)
- `ObjectDefinitions/Beam.f`, `Inputs/CubeGeneratorInputs.f`,
  `Inputs/InputUnitConversions.f`, `PostAnalysis/PostGalaxyAnalysis.f`
  (1 site each)

For every site with a live JS port, checked the port FIRST -- every single
one already used direct multiplication (`x*x`, never `Math.pow(x,2)`), so
none needed a JS-side change; this really is a one-sided Fortran bug
class, not a porting gap. Confirmed one incorrect assumption from an
EARLIER session's own comment in `GetMomentMaps.js` claiming
`moments[1]*moments[1]` was safe to leave as a plain multiply because
"Fortran's `**2.` constant-folds to a plain multiply for a literal
exponent" -- UPDATE 15's own reproduction disproves this (it does NOT
reliably constant-fold or reduce to `X*X`); the JS code itself was
already fine (already using multiply), just the comment's reasoning was
wrong -- left uncorrected in the JS file itself since fixing prose-only
comments across every affected file wasn't the priority, but flagging it
here so it isn't cited as evidence for a future "is `**2.` safe" question.

**Rebuilt and re-verified after the full sweep** (not verified after each
individual file -- given every site provably JS-matched already and the
change is mechanically identical everywhere, batched the verification
instead of the originally-planned incremental approach): anchor fit
`FINALVEC` still bit-exact (`41ABAAD3`/.../`3AB6767D`, unchanged).
Resampled-cube diff (seed=42, nBootstraps=1): unchanged at
`1.862645149230957e-09` from UPDATE 15 -- expected, none of the other 25
sites feed the resampling step itself (they're all in the FIT's own
moment-map/shape-estimate/chi2/particle-generation path, not
`BuildPhysCoordsArray`). 5-bootstrap comparison, both seeds (42, 99):
**zero change from UPDATE 15's post-single-fix numbers** -- realization
0 (seed 42) and realization 1 (seed 99, convergence failure) are BOTH
unchanged. Not a red flag: these are the same already-documented
Nelder-Mead chaos-sensitivity and SoFiA-threshold cases, and this sweep's
value is closing LATENT risk (other galaxies/seeds/cube sizes whose
particular operand values might hit gfortran's rounding-boundary case)
rather than fixing these two specific already-known hard cases.

**"Similar classes of bug"**: looked for, did not find, a second distinct
class this session. Specifically checked: (1) `**3.`/`**0.5`/other real-
exponent literals besides `2.` -- none found anywhere in `src` (grep for
`\*\*[0-9]+\.[0-9]*|\*\*\.[0-9]+` after the sweep returns zero live-code
hits). (2) `LOG`/`EXP`/`SQRT` intrinsics without the fdlibm substitution --
already systematically swept in an earlier session (`fd_log`/`fd_exp`/
`fd_atanh` etc. already applied wherever needed; `CalculateBeamKernel.f`'s
diff above shows this exact substitution already in place for `cos`/`sin`/
`exp` there). (3) Real-valued INTEGER-exponent forms like `**2` (no
decimal point) -- these ARE algebraically simplified by gfortran (verified
via the same `transfer()`-based standalone-program technique used for
`**2.`: `X**2` bit-identical to `X*X` for the same rounding-boundary value
that broke `X**2.`) and are NOT part of this bug class; left untouched
where they occur (e.g. `EstimateShape.f`'s `bb2**2`/`bb1**2`, both now
also converted to explicit multiplication anyway during the sweep, purely
for consistency/readability, not because they needed the fix).

**Recommendation, worth raising with Nathan**: this is a clean, well-
evidenced, reproducible upstream bug report -- `X**REAL_LITERAL` (a real,
not integer, exponent) is not guaranteed bit-identical to repeated
multiplication in gfortran, confirmed with a minimal standalone repro
(`transfer()`'d bit pattern, `-O0 -ffp-contract=off`, same flags this
project's own native build uses). Affects 26 sites across his own
codebase, all still present in `NateDeg/3KIDNAS` Dev branch as of this
session. Doesn't affect single-platform Fortran-only reproducibility (the
same binary always gives the same answer), but is a real portability/
numerical-hygiene issue (different gfortran versions, or a different
compiler entirely, could plausibly get a different -- still "valid" --
answer at any of these 26 sites) worth knowing about independent of this
project's own cross-language parity motivation for finding it.

## UPDATE 17: pathological JS-optimizer bisection, seed=1000/realization=4

Dan reported a NEW symptom from a large live batch: Fortran converges fast
on every single bootstrap realization, but JS occasionally burns ~60-70x
longer and fails to converge (`converged=false`, ITMAX=5000 hit) on
specific realizations, while Fortran sails through the exact same
realization normally. Reproducible instance found: `seed=1000`,
`nBootstraps=5`, `WALLABY_J100336-262923` (PA=81.271, Inc=31.49,
cdens=20) -- **realization 4** (0-indexed, 5th of 5). Fortran: fast,
pass-2 `iter=39`. JS: `converged=false`, ~100-125s.

**Methodology (Dan's proposal)**: dump Fortran's EXACT pass-2 starting
simplex (all 14 vertices x 13 params, full float32 hex) to a file via
`FORTRAN_SIMPLEX_DUMP_PATH` (new env var, `GalaxyFit.f`), force-feed the
identical simplex into JS's `amoeba` via `JS_SIMPLEX_OVERRIDE_PATH` (new
env var, `GalaxyFit.js`). If JS then reproduces Fortran's
iteration/convergence behavior, `amoeba` itself is exonerated and the gap
is upstream (objective-function-level noise accumulating over the
trajectory). If JS STILL gets stuck, the bug is in `amoeba` or in
`funk`/`tiltedRingModelComparison` itself.

**Result so far**: JS given Fortran's EXACT starting simplex still shows
`converged=false` (~102s, materially unchanged from the unforced run) --
ruling out "just a different starting point plus Nelder-Mead chaos" as
the full explanation. This directly matches Dan's own framing: "if we
STILL get something different, then the problem lies in the JS
implementation... maybe something is lingering between computations,
maybe something isn't updating, maybe something is computing with too
much or not enough accuracy."

**A costly false trail, now resolved**: bisecting call-by-call against
Fortran's own `FittingLog.txt` (`FULLVEC call=N` / `FULLVECPARAM` lines,
already-existing `TraceSwitch`-gated instrumentation, enabled here via
`WRKP_TRACE_DEBUG=1`, which `FitDriverScripts/RunWRKP.py` wires to
`TraceSwitch`) initially looked like a Fortran-internal inconsistency:
the dumped vertex-2 (`41A46ED2`, 20.55) didn't match what call=48 in the
log showed (`41850F70`, 16.63) -- and the dumped value only showed up
later, at call=63. This was NOT a bug. `TraceCallCounter` is a single
running counter shared across BOTH optimizer passes (pass 1's own initial
14-vertex loop + all of its amoeba iterations, THEN pass 2's), so naively
assuming "call 47 = pass-2 vertex 1, call 48 = pass-2 vertex 2" was wrong
-- call 47 just happened to coincide with pass 1's own converged best
point (chi2=113802.797, unchanged since pass 1's iter=31) being
re-evaluated, purely by coincidence of value, not position. The genuine
pass-2 vertex loop only starts once pass 1 fully finishes (confirmed via
the `ITER_F`/`CONVERGED_VECTOR`/`Param Guess Array Creation` markers
around it) -- for this run that's **call=62 through call=75** (14
vertices), not 47-60. Re-verified end to end with a clean, from-scratch
rerun (`node js/tools/run_both.js --seed 1000 --nBootstraps 5
--skip-js-dcp --cloudDensity 20 --objName WALLABY_J100336-262923 --cube
... --mask ... --pa 81.271 --inc 31.49`, `FORTRAN_SIMPLEX_DUMP_PATH` +
`WRKP_TRACE_DEBUG=1` set): dump and log are fully self-consistent once
calls are correctly attributed -- call=62's params == dump vertex 1
byte-for-byte, call=63 == dump vertex 2, etc. No Fortran bug here; this
was a bisection-methodology error on Dan's assistant's part, now
corrected. Fortran's 14 initial pass-2 vertex chi2 values for this
override (calls 62-75): 113979.195, 114029.508, 113835.164, 113839.617,
113793.039, 114036.594, 113944.852, 114076.828, 114091.867, 113794.438,
114152.555, 113875.203, 113916.289, 114386.992 -- this is the ground
truth the JS side's own 14 override-vertex chi2 values need to be diffed
against next, to see whether the FIRST divergence is already inside
`funk`'s very first (non-iterative) evaluation of these vertices, or only
appears once `amoeba` starts actually iterating.

**Incidental discovery while chasing this**: `run_both.js`'s default
galaxy (no `--objName`) is `WALLABY_J103538-484832`, NOT
`WALLABY_J100336-262923` -- an invocation with no args (e.g. testing
`--help`) silently runs a full, real, unseeded end-to-end fit rather than
printing usage and exiting. Worth a cheap follow-up (print usage and
exit when `--seed` is entirely absent AND no recognized flag is present,
or at least require an explicit `--seed`/`--help` before doing real work)
but not fixed yet -- out of scope for this investigation, noted here so
it isn't re-discovered by surprise later.

**Status**: in progress. Next step is comparing JS's own 14 override-
vertex chi2 values (same forced simplex, `JS_SIMPLEX_OVERRIDE_PATH` +
`TRACE_DEBUG=1` + per-realization `TRACE_DEBUG_FINALVEC_FILE` suffixing)
against the Fortran ground-truth list above, vertex by vertex.

### UPDATE 17 continued: root cause found -- RNG desync, not an amoeba/funk bug

**First pass (params forced, idum NOT forced)**: JS's 14 override-vertex
chi2 values came back wildly different from Fortran's -- ~113979-114387
(Fortran) vs ~121030-121692 (JS), a uniform ~6.6% gap across EVERY
vertex. That's too large and too uniform to be float32 rounding noise.
The reason: `funk()`/`tiltedRingModelComparison` has a SECOND hidden
input besides the parameter vector -- the RNG's mutable `idum` state
(particle placement is Monte Carlo). The override only forced the
parameters; `idum` was left to whatever pass 1's own natural, unforced
run happened to leave it at, and pass 1's own trajectory is NOT
externally forced at all. Checked directly: JS's idum at vertex 1
(266674843) bore no resemblance to Fortran's (1885505620) -- confirmed
desynced.

**Second pass (params AND idum forced, per-vertex)**: built a proper
controlled experiment. Extracted Fortran's own idum stream from its
`TRACE` lines (idum is one continuously-advancing stream shared across
BOTH passes' vertex loops AND iterations -- the printed value is the
POST-call state, i.e. call N's OUTPUT idum is call N+1's INPUT idum).
For pass 2's 14 vertices (Fortran calls 62-75), that gives the exact
INPUT idum each vertex evaluation should start from: 310281553,
1885505620, 1101243192, 1701357900, 1717828047, 1806624060, 989084484,
1372689030, 1099151444, 1704241316, 352627044, 960420574, 1362287062,
1355376704 (vertex k's input = this list's k-th value). Added
`JS_IDUM_OVERRIDE_SEQUENCE_PATH` (new env var, `FullModelComparison.js`)
to inject the matching value into `state.rng.state.ran2State.idum`
before each of the 14 vertex evaluations.

**Two methodology bugs found and fixed along the way** (both real,
both worth remembering if this apparatus gets reused):
1. Both overrides initially applied unconditionally to EVERY
   `galaxyFit_Simple` call -- including the INITIAL FIT, which has no
   meaningful realization index and got corrupted ("No best fit model
   made"), cascading into every bootstrap realization derived from it.
   Fixed with `JS_OVERRIDE_REALIZATION_INDEX`, checked in both
   `GalaxyFit.js` (simplex override) and `FullModelComparison.js` (idum
   override).
2. Realization-index scoping alone wasn't enough for the idum override:
   pass 1 shares the same `realizationIndex` as pass 2, and pass 1's own
   natural (unforced) `evalCount` can pass through the SAME numeric range
   later used for pass 2's vertices (this actually happened -- adding
   instrumentation shifted pass 1's own length between reruns), silently
   corrupting pass 1 instead of pass 2 and shifting where pass 2 actually
   starts. Fixed with a `state._simplexOverrideActive` flag, set only at
   the exact point the simplex override applies (pass 2, post-guess-
   array), checked by the idum override alongside realization scoping.
   Also replaced the fixed `JS_IDUM_OVERRIDE_START` call-number guess with
   a dynamic anchor (`_idumSeqStartCall`, set to whatever `evalCount` is
   the first time the flag is seen active) so this no longer silently
   misaligns if pass 1's natural length changes again.

**Result, with both bugs fixed**: JS's 14 pass-2 vertex chi2 values
(calls 60-73 this run): 113886.992, 113981.664, 113782.273, 113755.758,
113874.578, 114009.344, 113816.156, 113943.297, 114149.938, 113806.500,
114047.359, 114011.656, 113954.797, 114291.086 -- every one within
0.01%-0.12% of Fortran's corresponding value (113979.195, 114029.508,
113835.164, 113839.617, 113793.039, 114036.594, 113944.852, 114076.828,
114091.867, 113794.438, 114152.555, 113875.203, 113916.289, 114386.992).
That's the SAME residual float32-rounding-level noise this entire
session has been chasing down (X**2., PA, fdlibm, ...), not a new bug.
And critically: ZERO `badModelCheck` rejections this run, vs. TWO out of
14 in the first (unfixed-idum) pass -- both of those were a `Sigma<0`
rejection on a surface-density value that was only ever a few times
10^-4 from zero, i.e. the SAME rounding noise pushed a near-zero value
across the sign boundary on one platform but not the other (verified
Fortran has the identical `Sigma.lt.0.` check -- not a missing-check
mismatch, a genuine boundary-sensitivity artifact of the noise itself).

**Conclusion**: `amoeba` is not buggy, `funk`/`tiltedRingModelComparison`
is not buggy (beyond the already-known, already-minimized residual
float32 noise), and `badModelCheck` is not buggy (Fortran has the same
check and is equally exposed to the same boundary sensitivity in
principle). The seed=1000/realization=4 symptom Dan asked about --
"Fortran converges fast on every bootstrap, JS is pathologically slow on
some" -- is NOT a JS implementation defect. It's an emergent consequence
of two things compounding: (1) Nelder-Mead's well-known sensitivity to
its objective function's exact values, and (2) the objective function
being a Monte Carlo particle simulation whose randomness is a SINGLE,
continuously-advancing stream shared across the ENTIRE fit (both passes,
thousands of evaluations). A ~0.05%-level chi2 difference on any ONE
early evaluation (the residual noise floor already documented) can
change which iteration path pass 1's own amoeba takes, which changes how
many `ran2()` draws pass 1 consumes in total, which fully desyncs
Fortran's and JS's idum streams for everything from that point forward
-- not a small perturbation once that happens, but a completely
different random particle sample for every subsequent evaluation.
Nelder-Mead run on a different underlying random sample can legitimately
converge fast on one platform and stall/fail to converge on the other;
that's expected stochastic-optimizer behavior, not a cross-platform bug.
This also directly explains why it's realization-specific and
unpredictable rather than a systematic, always-reproducible-percentage
gap: it depends on whether pass 1's own iteration count happens to drift
apart between platforms for that specific realization's specific
resampled data, which is itself sensitive to the same ~0.05% noise floor
present in every realization.

**Practical takeaway**: this class of "pathological slow/non-converging
realization" cannot be fixed by further bit-parity work -- the remaining
noise is already near the floor of what's achievable given float32
arithmetic and two different platforms' transcendental-function
implementations. It's an inherent property of pairing a chaotic
deterministic optimizer with a stochastic objective function, present in
some form on EITHER platform in isolation too (Fortran's own fit would
be similarly sensitive to a hypothetical tiny perturbation). Not
recommended as an upstream bug report -- there's no bug to report.

**Cleanup status**: `JS_OVERRIDE_REALIZATION_INDEX`,
`JS_IDUM_OVERRIDE_SEQUENCE_PATH`, `state._simplexOverrideActive`, the
`badModelCheck._lastReason` tagging, and the `EVALCHI2`/idum file-log
lines in `FullModelComparison.js`/`GalaxyFit.js` are all still in place,
env-var-gated (zero cost when unset), following this project's established
pattern of leaving diagnostic tooling in place rather than ripping it out
-- not yet explicitly decided with Dan whether to keep permanently or
strip before the next commit.

### UPDATE 18: Dan's pushback -- asymmetry is real, hunted further, found a fix

Dan correctly rejected UPDATE 17's framing: "ALL of fortran's 200
bootstraps complete very fast whereas JS has SEVERAL blowups... performance
is clearly wrong... figure that out." Right call -- symmetric chaos would
predict Fortran occasionally failing too, and it doesn't. Kept hunting.

**Checked and ruled out, each with direct evidence, not assumption:**
- `amoeba`'s algorithm itself: read both implementations side by side,
  line for line. Faithful port, including the already-shared/fixed
  shrink-step bug. No divergence in the mechanics.
- `ftol`: identical on both platforms (0.005 pass 1, /5=0.001 pass 2) --
  confirmed by direct source read, not just comment-trusting.
- `badModelCheck`: found and fixed a REAL gap -- `vSinI` used native
  `Math.sin()` instead of `fdSin` (every other trig call in the hot path
  already uses the fdlibm substitutes; this one boundary check, evaluated
  on every single optimizer call, was missed). Fixed in
  `FullModelComparison.js`. Tested: zero effect on realization 4's own
  outcome (identical final chi2 before/after) -- inclination values in
  this trajectory never hit a value where native Math.sin and fdSin
  actually disagree. A real, worthwhile fix, just not the answer for this
  specific case.
- `nParticles` truncation-boundary sensitivity (the same class of bug as
  the `Sigma<0` badmodel boundary flip): added file-based NPTRACE logging
  keyed to `evalCount`/`realizationIndex` (new cross-module
  `global.__TRACE_REALIZATION_INDEX`/`__TRACE_EVAL_COUNT` tags, since
  `ring_CalcNumParticles` has no direct `state` access) and compared
  against Fortran's own NPTRACE output for the identical forced vertex:
  EXACT match, all 4 rings, both nParticles counts AND sigma hex values
  bit-identical (396/1187/1977/2768, same sigma hex on both sides). Not
  the source.
- FFT rank-2 (row-then-column vs FFTW's native 2D plan) accumulation-order
  difference: already investigated and documented in an EARLIER session
  (`FFTW3WasmRank2.js`'s own header) -- checked against a real compiled
  FFTW ground-truth harness, found ~1-96 ULP differences in DOUBLE
  PRECISION (~1e-14 relative), 1000-100000x smaller than a single float32
  ULP. Ruled out as too small to explain a ~0.4% noise floor.
- Ring/particle iteration order in `FillDataCubeByTiltedRing`: read both
  implementations, confirmed identical ring-major/particle-minor nested
  order and per-cell accumulation sequence on both platforms. Not a
  source of non-associative-summation divergence.
- Broader sweep of `CubeComparison.js`/`CubeKernelConvolution.js`/
  `CalculateBeamKernel.js`/`ParameterToTiltedRingInterface/*.js` for any
  remaining native `Math.sin/cos/atan/log/exp/pow` calls: clean, none
  found beyond the one already fixed.

**New, decisive measurement**: extracted Fortran's own `Current tolerance`
(rtol) trajectory for THIS SAME realization's pass 2, straight from its
FittingLog:
```
iter  rtol         y_hi        y_lo
0     0.00521      114386.99   113793.04
10    0.00291      113979.20   113648.02
14    0.00200      113875.20   113648.02
29    0.00261      113944.71   113648.02   <- got WORSE, not monotonic
31    0.00165      113828.73   113641.32
39    0.000996     113754.58   113641.32   <- converged, barely under 0.001
```
Fortran is NOT cruising to a clean answer here -- its own rtol is
non-monotonic and hovers in the same 0.1%-0.3% range JS gets stuck in,
before happening to dip a hair under ftol=0.001 at iter=39. It won this
specific case by a near-miss, not a comfortable margin. Combined with
JS's own measured chi2 spread (flat at 0.385%-0.403% across five
1000-call windows near the end of a stuck run, NOT narrowing -- direct
measurement, not inferred), the conclusion: `ftol=0.001` sits at or below
the intrinsic Monte-Carlo-plus-rounding noise floor of this objective
function for hard resampled datasets, on EITHER platform. Nelder-Mead
cannot satisfy a tolerance tighter than its own measurement noise, no
matter how many more iterations run once truly stalled (not slowly
converging -- confirmed flat, not shrinking, over 1000+ evals).

**Honest conclusion on root cause**: exhausted every discrete,
find-and-fix mechanism checked above. The residual ~0.05%-0.4%
cross-platform noise most likely isn't ONE bug but the compounding tail
of many already-individually-negligible rounding differences across
thousands of arithmetic operations per evaluation (particle placement +
convolution + chi2) -- the same long-tail pattern this whole session's
`X**2.`/PA/fdlibm work has been closing, now at a point where remaining
gaps are individually too small to find by code review and would need a
much more exhaustive automated per-operation trace-diff to fully close.
JS's noise floor being slightly higher than Fortran's (enough to lose
this near-miss coin-flip more often across many bootstrap realizations)
is the most likely explanation for the asymmetry Dan observed, but this
is not fully proven to the same standard as the earlier X**2./PA fixes.

**Practical fix applied** (`amoeba`, `GalaxyFit.js`): a stall detector,
NOT a numerical-parity fix -- added on top of, not instead of, the
investigation above. Tracks the simplex's best value (`y[ilo]`); if it
hasn't improved by more than `ftol/10` relative over a
`max(300, 20*ndim)`-iteration window, exits early (same `noConvergence`
signal as hitting ITMAX) instead of grinding to ITMAX=5000 once genuine
progress has demonstrably stopped. Chosen conservatively (window wide,
epsilon tight) specifically so it cannot fire during real, if slow,
convergence -- only once measured progress has flatlined. Verified on
seed=1000/nBootstraps=5: realizations 0-3 (all normally-converging)
produce BIT-IDENTICAL chi2/convergence results before and after (detector
never fires for them, as designed) -- realization 4 drops from ~101s to
~10s (10x), still honestly reported as `converged=false` (does not
fake success), landing at essentially the same chi2 (113338.680 vs
113607.586, both inside the already-measured-flat ~0.4% noise band).

**Not yet done**: Fortran's own `DownhillSimplex.f` has no equivalent
stall detector -- intentionally left untouched (Fortran doesn't appear to
need it, per Dan's own 200/200 success report, and this session found no
evidence Fortran's ITMAX=5000 cap is ever actually hit). Worth
reconsidering only if Fortran itself is later found to stall on some
harder dataset.
