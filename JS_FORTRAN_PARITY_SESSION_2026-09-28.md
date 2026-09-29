# JS/Fortran bit-exactness session notes -- 2026-09-28

Continuation of `JS_FORTRAN_PARITY_SESSION_2026-09-15.md`'s tracked work.
Started from a user question ("is FFTW the single source of errors?") that
turned into a from-scratch re-verification: FFTW was fully exonerated this
session (forward r2c and inverse c2r both proven bit-exact), and the real
bugs found were all elsewhere. Working galaxy: `WALLABY_J103538-484832`
(the default test galaxy -- has correctly-named `WALLABY_..._mask.fits`,
unlike `WALLABY_J100336-262923` which does not and causes an early,
confusing "No best fit model made" failure if used by mistake).
`cloudDensity=20`, `BootstrapSeed=42` throughout (user's explicit
preference -- density=500 was tried briefly to dodge density=20's inherent
optimizer noise, but the user correctly redirected back to 20: density=500
runs are much slower and don't actually help verify bitwise parity, since
parity should be checked via targeted intermediate-value dumps, not by
waiting for full noisy convergence).

## Central methodology (same as 2026-09-15, re-proven effective)

Bisect by dumping intermediate values from both platforms to matching
`TRACE_DUMP_PRECONV`-gated files (NOT `print*`/`console.log` -- Fortran's
stdout gets redirected/buffered unreliably for bootstrap sub-processes, a
recurring dead end this session before switching to files), then diff. This
found every real bug below. Column-major-vs-row-major array traversal
order is a recurring, sneaky bug class distinct from ordinary rounding --
see bugs #2/#3.

## Confirmed real bugs found + fixed this session

1. **`js/src/StandardMath/fdlibm.js`: stale FMA emulation.** `__kernel_sin`/
   `__kernel_cos`/`fdLog`/`fdLog1p`/`fdExp`/`fdAtan` emulated hardware FMA
   (via `fma.js`) based on a comment claiming the project's C build used
   FMA. Disassembled the actual current `fdlibm_k_sin.o`/`fdlibm_k_cos.o`
   objects: zero FMA instructions (the build now uses
   `-ffp-contract=off` in `src/makeflags`, added after that comment was
   written). Rewrote all six functions to plain, unfused arithmetic
   matching the C source's literal operation order. Verified 0/20000, then
   0/50000 domain-realistic mismatches against the real compiled C fdlibm
   (was ~0.3-0.6% before). `fma.js` itself is now dead code (kept, unused).

2. **`js/src/ConvolveCube/CalculateBeamKernel.js`: kernel-sum traversal
   order.** JS summed the beam kernel with a flat sequential scan; Fortran's
   `sum(B%Kernel)` on its column-major 2D array walks in `j`-outer/`i`-inner
   memory order. Non-associative float addition made this a real, if tiny,
   bug. Fixed the JS sum to `for j { for i { ... } }` explicitly. Verified:
   `rawsum`/`center`/`corner`/`r1c1`/`r2c3` all now match Fortran bit-exact.

3. **`js/src/PreAnalysis/InitialAnalysis.js`: `SInt`/`nCells` traversal
   order (3D version of bug #2, bigger blast radius -- full observed
   cube).** First fix attempt used `i`-outer/`j`-middle/`k`-inner, which is
   WRONG for a 3D array and (confirmed empirically) numerically identical
   to a flat scan, not real column-major order. Corrected to `k`-outer/
   `j`-middle/`i`-inner (true Fortran column-major: first index fastest).
   Verified: `SInt` now matches Fortran bit-exact
   (`6.9436478614807128906E-01` both sides).

4. **`js/bootstrap-realization-launcher.js`: `FITAchieved` gated on
   `noConvergence`.** Was `report.FITAchieved = !noConvergence`, meaning
   JS aborted the ENTIRE pipeline (no bootstrap realizations at all) any
   time `amoeba` hit `ITMAX` without formally converging. Fortran has no
   such gate -- Python's own `FITAchieved` (`ReadWRKPFit.LoadBestFitModelFile`)
   is just `os.path.isfile(...)`, and Fortran always writes its output
   using whatever best point `amoeba` found. This is a REAL, user-caught
   behavioral bug (user: "if it works in fortran but not in js, there's
   already a problem") -- not a precision issue. Fixed to
   `report.FITAchieved = true` unconditionally at that point;
   `report.converged` kept as a separate diagnostic field. This is what
   makes `cloudDensity=20` usable for testing at all now (JS no longer
   spuriously aborts on the same noise floor Fortran tolerates).

5. **`src/Outputs/FitOutput.f`: `F16.6` (6 decimal places) insufficient for
   `X_kin`/`Y_kin`/`RA`/`DEC`/`Inc_kin`/`PA_kin`/`VSys_kin`/`VDisp_kin`
   text output.** `VSys_kin` wrote `5745.641602` for true value
   `5745.6416015625` (~4.4e-7 truncation) -- small enough to round-trip to
   the same float32 bit pattern, but NOT reabsorbed by
   `RunBootstrapsDCP.ComputeBsCent`'s `VCenter=DeltaV/dV+RefChan`, which
   operates in float64 directly on the parsed text. First fix attempt
   (switching to bare list-directed `write(ValStr,*)`) was ALSO
   insufficient -- gfortran's list-directed `real(4)` formatting caps at
   ~9 significant digits regardless of magnitude, which happened to be
   WORSE for `VSys` specifically (`5745.64160`, only 5 decimals). Final fix:
   wrap every write in `DBLE(...)` to force full double-precision
   list-directed formatting (~17 sig figs) -- `VSys_kin` now writes
   `5745.6416015625000`, exact. `character(16) ValStr,ErrStr` widened to
   `character(30)` to hold the longer output. **This fix is real and
   correct (verified: the written text is now exact) but turned out NOT to
   explain the `CentV`/`PA` bug below -- keep it anyway, it's a genuine
   precision improvement independent of that bug.**

6. **`js/bootstrap-realization-launcher.js`: `bsCent` missing float32
   rounding -- the actual explanation for the `CentV`/`PA` mismatch bug #5
   didn't explain.** `BootstrapGlobals.f`'s `Type BootstrapCenter` declares
   `CentX`/`CentY`/`CentV`/`PA`/`Inc` as `real` (single precision) --
   Fortran's `read(10,*) BS_Cent%CentX,...` forces an implicit
   double-to-float32 rounding the moment Python's full-double-precision
   `ComputeBsCent` output is parsed back in. JS's `bsCent` object carried
   the full, un-rounded JS double straight through into
   `genFlipBootstrapSample`/`axisFlip`/`getCubeCoords`, silently using more
   precision than Fortran's own data type ever had. **Confirmed
   empirically, decisively**: `Math.fround(js_centV) === fortran_centV`
   and `Math.fround(js_pa) === fortran_pa`, exactly, for values that had
   differed by ~1 ULP. Fixed: added a `bsCentF32` object
   (`{...bsCent, centX: f32(...), centY: f32(...), centV: f32(...),
   pa: f32(...), inc: f32(...)}`) right after `bsCent` is destructured from
   `payload`, and pass `bsCentF32` (not `bsCent`) to
   `genFlipBootstrapSample`. Verified: `CentV` now matches Fortran
   bit-exact (`5.0531715393066406250E+01` both sides, was
   `...5393066406250` vs `...3807234346802` before).

## Independently proven bit-exact this session (not bugs, just verification)

- Forward r2c and inverse c2r FFTW transforms: bit-exact (2112/2112 and
  4096/4096 bins), using the project's own patched FFTW (routes twiddle
  trig through the shared fdlibm C source on both native and wasm builds --
  `third_party/fftw-3.3.8/kernel/trig.c`'s own header comment documents
  this, predates this session).
- The ENTIRE anchor-fit `amoeba` trajectory, both passes, to full
  completion: Pass 1 (65/65 iterations, converges) and Pass 2 (1715/1715
  iterations, hits `ITMAX=5007`) -- `rtol`/`y[ihi]`/`y[ilo]` match Fortran
  exactly at EVERY iteration, not just spot-checked. This is a much
  stronger proof than "same final chi2" -- it's iteration-by-iteration.
- The resulting best-fit parameter vector: all 15 params bit-exact
  (`PvBestTraceF.txt`/`PvBestTraceJS.txt`).
- `BS_Cent` geometry (`CentX`/`CentY`/`CentV`/`PA`) after bug #6's fix:
  bit-exact.

## RESOLVED: model-cube resynthesis divergence (bug #7)

**Full bisection chain, in order, each step proven before moving to the
next:**
1. `amoeba` trajectory: bit-exact, every iteration, both passes (see above).
2. Best-fit parameter vector: bit-exact (see above).
3. `BS_Cent` geometry: bit-exact after bug #6.
4. JS's *natural* `ran2` RNG state at the model resynthesis point: bit-exact
   against Fortran's captured state (`idum`/`idum2`/`iy`/`iv[32]`, all 32
   values) -- ruled out RNG drift entirely. (A prior RNG-transplant test via
   `TRACE_OVERRIDE_RAN2_STATE` gave a confusing literal `0.0`, later found
   to be a red herring from tracing the wrong array in a
   `--skip-fortran --local` anchor-only invocation -- abandoned once the
   natural-state comparison gave a clean answer instead.)
5. Pre-convolution model cube (`PreConvModel.fits`): **0/179520 diffs**,
   bit-exact.
6. Post-convolution-pre-scale model cube (`PostConvPreScaleModel.fits`):
   0 meaningful diffs (144228 raw diffs, all at the ~1e-22 numerical-zero
   noise floor -- FFT/convolution "ringing" far below any real flux level,
   inherent to non-associative floating point summation, not a bug).
7. **Final (post-`beamArea`-rescale) model cube: root cause found here.**
   `js/bootstrap-realization-launcher.js`'s `runInitialFit` reused
   `fitBeam.beamAreaPixels` (the "official" field, computed once at Beam
   allocation via `2.*Pi/(2.355*2.355)*Major*Minor` -- matches
   `Beam.f:62`'s `B%BeamAreaPixels` exactly) for the Jy/pixel->Jy/beam
   rescale (`modelFlux[i] *= beamAreaPixels`) before writing the model
   cube. But `FitOutput.f`'s `OutputCube` does NOT use that field for this
   specific rescale -- it computes a separate, local `BeamArea` variable
   inline via a different (mathematically equivalent, numerically NOT
   identical) formula: `2.*Pi*abs(BeamSigmaVector(0)*BeamSigmaVector(1))`,
   built from the already-divided sigma values (`Major/2.355`,
   `Minor/2.355`, each independently rounded) rather than one precomputed
   `2*Pi/2.355^2` constant times `Major*Minor`. A uniform scalar
   multiplier differing by ~1 ULP produces exactly the observed signature:
   26492 voxels differing by a uniform ~1.4e-7 relative amount (right at
   float32 machine epsilon) across the ENTIRE cube, roughly evenly split
   inside/outside the mask (ruling out a mask-related explanation tried
   first). **Fixed**: added a freshly-computed `outputBeamArea` in
   `runInitialFit`, replicating `FitOutput.f`'s exact local formula and
   operator grouping (using `fitBeam.beamSigmaVector`, already correctly
   populated), and multiply by THAT instead of `beamAreaPixels` for this
   one rescale. `beamAreaPixels` itself is untouched and still correct for
   its own (different) purpose elsewhere.

**Verified end-to-end**: after the fix, the final model cube shows
**0 meaningful diffs** (down from 26492), and the full bootstrap-resampled
cube (`ResampleFullF.txt`/`ResampleFullJS.txt`) shows **0/179520 diffs,
including zero of the previously-present numerical-noise-floor
artifacts** -- genuinely, completely bit-identical, not just
"bit-identical above some magnitude threshold." This closes the entire
investigation: every stage from the optimizer's own internal evaluations
through to the final bootstrap-resampled cube used for error analysis is
now proven bit-exact between Fortran and JS at `cloudDensity=20`,
`BootstrapSeed=42`, for `WALLABY_J103538-484832`.

**Bug count for this session: 7 real, confirmed bugs found and fixed**
(6 in the earlier list + this one), all independently verified by direct
before/after comparison, not just code review.

## FOUND + FIXED 2026-09-29: bugs #8 and #9 (harder test case: WALLABY_J100336-262923)

Discovered while re-running the exact reproducible test case from
`HANDOFF_JS_FORTRAN_NOISE_FLOOR_2026-09-16.md` (`WALLABY_J100336-262923`,
`pa≈81.271`, `inc≈31.49`, `cloudDensity=20`, `nBootstraps=10`) -- a harder
case than the `WALLABY_J103538-484832` anchor test bug #7 was proven
against. `run_both.js --local` showed `Inc_model` off by up to ~4e-6
absolute and `RHI_AS` off by up to **4.45%** despite `X_model`/`PA_model`/
`Vsys_model` often matching exactly, and despite everything checked (binary
build timestamps, `--local`'s direct `./src/*` requires, output-file
timestamps) confirming nothing was stale.

**Bug #8 -- `RAD2DEG` double-precision ratio collapses Fortran's two
sequential float32 roundings into one.**
`js/bootstrap-realization-launcher.js` (both `runBootstrapRealization` and
`runInitialFit` copies) computed degrees-from-radians via a precomputed
`const RAD2DEG = 180.0 / Math.PI` (native double `Math.PI`, not the
already-fixed float32 `Pi` from `BasicConstants.js`), then did
`f32(rad * RAD2DEG)` -- ONE double-precision multiply against a precise
ratio, rounded to float32 once at the end. Fortran's actual formula
(`FitOutput.f`, `Inc_kin`/`PA_kin`): `Inclination*180./Pi` -- `*` and `/`
are equal precedence, left-to-right, so this is
`(Inclination*180.)/Pi` -- TWO SEPARATE float32 roundings (multiply, round;
divide, round), using the float32 `Pi` constant, not one double-rounded
constant multiply. Same class of bug as #7 (mismatched rounding sequence
for a mathematically-equivalent formula). Explains why PA happened to
match exactly for some realizations (double-rounding doesn't always cross
a boundary) while Inc, for this specific test case, did.
**Fixed**: imported `Pi` from `BasicConstants` into both blocks, replaced
`RAD2DEG`-based conversion with a `radToKinDeg(rad)` helper doing
`f32(f32(rad * f32(180.0)) / Pi)`, matching Fortran's exact two-step
sequence. Used by both `INCLINATION`'s conversion and `toKinematicPA`'s
first step.
**Verified**: `Inc_model`/`PA_model` now show **0.000000 max diff** across
all 10 realizations (were up to ~4e-6 before).

**Bug #9 -- radial-profile text output (`Rad`/`VRot_model`/`SD_model`) was
lossy, amplified into a 4.45% `RHI_AS` error.**
`FitOutput.f`'s `StandardModelOutput` wrote the per-ring radial profile
table with fixed, narrow formats: `Rad`/`VRot_kin` as `F8.2` (2 decimal
places) and `SD_kin` as `G9.2` (**TWO SIGNIFICANT FIGURES** -- e.g. "2.9"
for a true value of 2.9300459...). `ExtractScalingParams.py`'s RHI
extraction (shared Python code, used identically for both the
Fortran-local and JS-local legs) interpolates the model's surface-density
profile against a threshold to find the HI radius -- feeding it a
2-significant-figure-quantized `SD_model` (Fortran side) against a full
double-precision one (JS side, never quantized) was enough to shift the
interpolated radius by up to 4.45%. Confirmed directly: Fortran's own
`BootstrapFits.csv` showed `SD_model` values like `"2.9, 2.0, 0.31, 1.2"`
for the exact same realization JS reported
`"2.9300459036646607, 2.0146256636770046, 0.31207922029800184,
1.1935189127727546"`.
**Fixed**: same pattern as bug #5 -- widened `RadialProfStr` from
`character(20)` to `character(30)` and switched all six per-ring writes to
list-directed (`write(...,*)DBLE(...)`) instead of `F8.2`/`F5.2`/`G9.2`.
**Verified**: `RHI_AS` max diff dropped from **1.18 absolute / 4.45%** to
**0.000001 absolute / 0.00%** -- essentially bit-exact.

**Build-system fixes needed to compile/link bug #9's fix (see "Gotchas" below
for full detail):** `fftw_wisdom_helper.o` was permanently added to
`src/ObjectLists`' `StandardMathObj` (closes a long-standing Makefile gap),
and `FitOutput.o` had to be recompiled by hand before relinking, since this
Makefile has no per-file dependency tracking (`make` alone silently
relinks with a stale `.o` after a source edit).

**`RA_model`/`DEC_model` diffs (~1e-6) are NOT a bug**: traced to two
independently-correct astrometry paths (Fortran's `ArcSecToDegrees` gets
overwritten by `GeometryFix.py`'s own `astropy.wcs.pixel_to_world` +
explicit `round(...,7)`; JS's goes through a separate `RunBootstrapsDCP.py`
`all_pix2world` call) -- different formulas by design, not a parity defect.

## FOUND + FIXED 2026-09-29: bug #10 (X_model/Y_model, isolated to ONE realization)

After bugs #8/#9, one last residual: `X_model`/`Y_model` differed by
~1.5e-5/2.9e-5 absolute (0.00% relative -- easy to dismiss as noise) on the
same `WALLABY_J100336-262923` 10-bootstrap run. Bisected by checking every
realization individually: **9 of 10 were already exactly 0.0000000** --
only realization index 8 (`BS_8`) differed at all. A uniform "just noise"
explanation doesn't fit a signature that's exactly zero on 9/10 draws and
nonzero on one -- that's the same boundary-crossing signature bugs #7/#8
had (a constant, systematic bias that only flips a near-tied decision for
specific inputs), so kept looking rather than writing it off.

**Root cause**: `runInitialFit`'s beam construction --
`fitBeam.beamPositionAngle = f32(f32(bpaDeg * Math.PI) / f32(180.0))` --
used native double `Math.PI` instead of the module's own float32 `Pi`
constant (`BasicConstants.js`, `f32(4.0*Math.atan(1.0))`). Fortran's
matching routine (`UnitConversions.f`'s `DegreesToRadians`:
`L_Rad=L_Deg*Pi/180.`) multiplies by its OWN real4-rounded `Pi`, not a more
precise value -- `bpaDeg(f32) * Pi(f32)` (Fortran's real4*real4 multiply,
using an already-quantized constant) is not guaranteed to round to the
same float32 result as `bpaDeg * Math.PI` (JS's double-precision multiply
by a much more precise constant, rounded to float32 only once at the end).
This `fitBeam` (and its `beamPositionAngle`) is built ONCE in
`runInitialFit` and reused for every bootstrap realization's convolution
kernel (per that function's own comment) -- so the bias is constant across
all 10 realizations, but only realization 8's optimizer landscape happened
to have a decision sitting close enough to a tie for that bias to flip it.
**Fixed**: `f32(f32(bpaDeg * Pi) / f32(180.0))`, using the already-fixed
`Pi` import (already in scope from bug #8's fix, same function).

**Verified**: re-ran the exact 10-bootstrap repro -- `X_model`/`Y_model`
now **0.000000 max diff across all 10 realizations** (were 1.5e-5/2.9e-5 on
realization 8 alone). Every field this run is now 0.00% except
`RA_model`/`DEC_model`, which are a known, by-design non-bug (above).

**Bug count for this session: 10 real, confirmed bugs found and fixed.**

## Also fixed 2026-09-29: run_both_report.html corner-plot NaN handling

Separate from the numerical parity bugs above -- a display bug in
`js/tools/run_both_report.html`'s corner plot, noticed on a VHI panel
where extraction had failed for every bootstrap realization in the batch
(a real, legitimate data gap, not a bug: faint/noisy source, not enough
S/N to reach the extrapolation region). Two issues in the same file:

- `gbpDrawScatter` never filtered non-finite points before drawing (unlike
  `gbpDrawHist`, which already did). `gbpScaleFor` degrades a field with
  zero finite values to a placeholder `[0,1]` range, but the actual `(NaN,
  NaN)` points still got passed to SVG's `cx`/`cy`; browsers silently
  coerce an invalid numeric SVG attribute to `0`, which rendered as a
  dead-straight line of dots pinned to the top of every panel in that
  row/column -- looking exactly like real, tightly-clustered data instead
  of "no data at all". Fixed: skip drawing (mirrors `gbpDrawHist`'s
  existing `isFinite` guard).
- Added an explicit "no data" label (both `gbpDrawScatter` and
  `gbpDrawHist`) when a field has zero finite values across every row, so
  an empty panel reads as "nothing to plot" rather than looking like a
  rendering bug.

## Historical section below, left as originally written (context for how
   bug #7 was found -- steps 1-6 above supersede the "in progress" framing)

Despite everything above being bit-exact, the **full resampled bootstrap
cube still shows ~5299/179520 voxels (~3%) differing by ~1 ULP each**
(`ResampleFullF.txt`/`ResampleFullJS.txt`, column-major-ordered flat dumps,
one value per line -- diff with a small Python script mapping flat index
back to `(i,j,k)` via `i=idx%nx; j=(idx//nx)%ny; k=idx//(nx*ny)`).

**Traced past `bsCent` entirely** (bug #6's fix had ZERO effect on this --
identical indices, identical values, bit-for-bit, before and after): the
**model cube itself** (the one-time resynthesis from the converged
parameter vector, used to build every bootstrap realization's difference
cube -- `FitOutput.f`'s `OutputCube` / JS's `runInitialFit`'s
`tiltedRingModelComparison(synthParams, state)` call) **already differs
at the affected voxels**, before any flip/interpolation code runs.
Confirmed via `VoxelTraceF.txt`/`VoxelTraceJS.txt` at `(i=29,j=11,k=2)`:
observed cube matches exactly; model cube does not
(`1.5355691402874072082E-06` Fortran vs `1.5355689129137317650e-6` JS,
~1.5e-4 relative -- much bigger than a simple ULP, this voxel is near the
noise floor).

This is surprising given the optimizer's own `chi2` evaluations (which
also resynthesize a model cube internally, every single call) were proven
bit-exact for 1715+ iterations. The working theory: the final output
resynthesis is a SEPARATE, one-off call
(not a replay of amoeba's last accepted evaluation), and if it consumes
RNG state (Monte Carlo particle placement) starting from a state that
naturally drifted even slightly, the resulting cube would differ despite
every prior optimizer step matching. **This exact question was
investigated in an earlier session too** -- extensive pre-built diagnostic
tooling already exists in `js/bootstrap-realization-launcher.js` for
exactly this: `TRACE_OVERRIDE_PARAMS`, `TRACE_OVERRIDE_IDUM`,
`TRACE_OVERRIDE_RAN2_STATE` (the last one takes Fortran's FULL `ran2`
state as JSON: `{"idum":N,"idum2":N,"iy":N,"iv":[32 ints]}`, matching
`FitOutput.f`'s `MaybeOverrideIdum`/`GetRan2State`).

**In progress when this doc was written**: a direct RNG-transplant test --
dump Fortran's exact `idum`/`idum2`/`iy`/`iv[32]` state right before its
own `BuildTiltedRingModel` call in `OutputCube` (added a new
`TRACE_DUMP_PRECONV`-gated trace, `OutputIdumTraceF.txt`, for this -- not
present before this session), then transplant that EXACT state into JS's
`runInitialFit` via `TRACE_OVERRIDE_RAN2_STATE` and re-check
`model_29_11_2`. First attempt at this test had a bug (traced the wrong
array -- `modelDC.flux`, which is `runBootstrapRealization`-specific and
never populated in an anchor-only `--skip-fortran --local` run, giving a
misleading literal `0.0`). Fixed to trace `fitModelDC.flux` (the actual
`runInitialFit` resynthesis output) right after its
`*= beamAreaPixels` conversion, matching Fortran's own
`ModelDC%Flux=ModelDC%Flux*BeamArea` unit convention. Re-running now via:

```bash
export TRACE_OVERRIDE_RAN2_STATE='{"idum":44935331,"idum2":1409304502,"iy":1161451613,"iv":[578099345,208640823,499590767,1148022595,1760857751,406777564,835806095,502605656,1985149339,386786030,2106073053,1263296582,789384200,280446953,1283079306,1644091015,1440078607,946610642,2142293930,1239805679,1745222282,169506799,44935331,663111705,563368577,1899430641,52157200,1577048331,1313558234,51483699,447954393,56650931]}'
TRACE_DUMP_PRECONV=1 node js/tools/run_both.js --local --skip-fortran --seed 42 --nBootstraps 1 --cloudDensity 20
# check js/DCPjobData/AnchorModelVoxelTraceJS.txt's model_29_11_2
```

**Next steps, in order**:
1. Get this RNG-transplant result. Compare `model_29_11_2` against
   Fortran's `1.5355691402874072082E-06`.
   - If it now matches exactly: RNG-state drift between the two platforms'
     accounting (something consumes draws differently between "amoeba
     returns" and "this specific resynthesis call") is the cause. Next:
     find WHERE the drift is introduced -- likely candidates are anything
     that calls the RNG between `amoeba`'s return and `BuildTiltedRingModel`
     (check `MaybeOverrideIdum`'s no-op path, `ParamToTiltedRing`, or
     whether `OutputBestFit_Simple` calls `TiltedRingModelComparison`/
     `StandardModelOutput` an extra, uncounted time with `FitNum=0` for the
     "IniEstimate" output BEFORE the `FitNum=2` call, consuming RNG state
     Fortran's own bookkeeping doesn't expect JS to replicate the same way).
   - If it still differs given IDENTICAL RNG state AND identical
     parameters: RNG state was never the cause, and the bug is in the
     resynthesis code proper (`SingleRingGeneration.f`/`.js`'s particle
     generation, or `FillDataCubeWithTiltedRing`/`CubeBeamConvolution` for
     this specific call path) despite that same code being proven correct
     when called as part of `chi2Calc` during optimization -- look for a
     divergent branch specific to the "final/output" invocation (e.g. a
     different `StrictEstimate` flag, particle count, or convolution
     padding decision that only triggers for this one call).
2. Once the model-cube divergence is fully explained and fixed, redo the
   full `ResampleFullF.txt`/`ResampleFullJS.txt` voxel diff (methodology
   above) to confirm 0/179520 diffs -- that's the actual end-to-end
   success criterion for "bootstrap resampling is bit-exact."
3. Only after that: consider whether the same rigor should be extended to
   the bootstrap realization's OWN fit (downstream of resampling, never
   checked this session) -- likely fine if everything upstream is finally
   exact and the same `amoeba`/`chi2Calc` code paths are reused, but not
   verified.

## Gotchas hit again this session (already documented 2026-09-15 but worth
   re-flagging since they cost real time)

- **`js/package/` (bravojs bundle) must be rebuilt after ANY `js/src/*.js`
  change**, via `cd js/package && node build-bravojs-package.js`, before
  `--local`/real DCP runs will reflect it. Found this bundle stale since
  Aug 18 early this session -- every fix before the rebuild silently never
  took effect in test runs. `js/bootstrap-realization-launcher.js` itself
  is NOT part of this bundle (top-level orchestrator file, edits apply
  immediately, no rebuild needed).
- **`fftw_wisdom_helper.o` wasn't wired into the Fortran Makefile** --
  `make all` failed the final link every time with an undefined-symbol
  error. FIXED PERMANENTLY (2026-09-29): added `fftw_wisdom_helper.o` to
  `StandardMathObj` in `src/ObjectLists`, so it's now part of `$(AllObj)`
  for every target. No more manual relinking needed.
- **This Makefile has NO per-file dependency tracking** -- `BootStrapSampler`/
  `SingleGalaxyFitter`'s targets only depend on their own main-program `.o`
  (`BootStrapGenerator.o`/`SingleGalaxyFitTests.o`), not on `$(AllObj)`, so
  editing e.g. `FitOutput.f` and running `make` does NOT recompile
  `FitOutput.o` -- it silently relinks with the STALE `.o`. Confirmed this
  bit a real fix (2026-09-29's radial-profile precision fix, below): `make`
  reported "Nothing to be done" and produced binaries with the pre-fix
  `FitOutput.o` still baked in. Workaround: manually recompile the specific
  changed `.f` file's `.o` (same flags as `makeflags`' `FLAGS`, with
  `-I../mods -I../../third_party/fftw-3.3.8/api`), THEN `rm` the target
  binaries before `make` (since existing binaries also short-circuit the
  link step) so it's forced to relink with the fresh object.
- **`WALLABY_J100336-262923` lacks a correctly-named mask file** (only has
  `SoFiA_J100336-262923_mask.fits`, not `WALLABY_J100336-262923_mask.fits`)
  -- using `--objName` to switch to it causes an early, confusing "No best
  fit model made" with no useful diagnostic. Stick to the default galaxy
  (`WALLABY_J103538-484832`) unless this is fixed.
- **Bash tool cwd does not persist reliably across calls in this
  environment** -- always use absolute paths or a single chained `&&`
  command; `cd X` in one call does not guarantee cwd=X in the next.
- **Disk space**: hit a genuine, machine-wide `ENOSPC` mid-session (98%
  full, 268MB free on a 228GB disk) -- unrelated to this project (freed
  ~140MB from session scratchpad files, but the real consumer was outside
  `/` visibility, likely Time Machine local snapshots or another APFS
  volume in the same container; resolved itself/was resolved externally
  during the session). If this recurs, check
  `tmutil listlocalsnapshots /` and Disk Utility, not this repo.

## OPEN, narrowed significantly 2026-09-29: cloudDensity=100 divergence isolated to the amoeba fit itself

Continuing the bug hunt on `WALLABY_J100336-262923`, 11 bootstraps,
`cloudDensity=100` (bugs #8/#9/#10 fixed the same test at `cloudDensity=20`;
this harder setting -- 5x the particles per ring -- exposed a NEW,
much larger divergence: realization 10 alone (of 11) differs by
X≈0.35/Inc≈0.91 absolute, ~1-2% relative, while the other 10 realizations
still match exactly).

**Ruled out, with direct proof, not just code review:**
1. `GeometryEstimates.js`'s `Math.acos` vs Python's `np.arccos` (initial
   inclination estimate): empirically swept 99,999 test values through
   both at float32 precision -- **zero differences**. Not the cause.
2. `FlipBootstrap.js`'s `thetaNew` rounding: initially suspected a missing
   intermediate float32 rounding (`f32(f32(2.0)*Pi - theta)`), but doubling
   a float32 value is always exact in IEEE-754, so the "fix" was a
   mathematical no-op (confirmed: adding it changed nothing; reverted).
3. `gasdev`'s rejection-sampling retry loop (`v1²+v2²<1`) -- the strongest
   structural candidate (a retry that COULD consume a different number of
   `ran2()` draws per platform) -- already matches Fortran's exact
   per-operation rounding sequence in `random.js`. No divergence possible.
4. `nParticles = int(CloudSurfDens * Sigma^cmode * PixelRingArea) + 1` --
   `cmode=0` for this test config (`Inputs/SingleGalaxyTestFittingOptions_
   Base.txt`), so `Sigma^0 = 1.0` exactly on both platforms (every `pow`
   implementation special-cases exponent 0). The rest of the formula
   already uses correctly-sequenced rounding matching Fortran op-for-op.
5. **Bootstrap resampling itself, proven bit-exact by direct reconstruction**
   (`js/tools/repro_resample.js` + a hand-built `BootstrapRuntimeInputs.f`
   input file replicating `MakeBootstrapSample.WriteBootstrapFile`'s exact
   geometry for realization 10, idum=-(42+10+1)=-53): ran JS's
   `genFlipBootstrapSample` and Fortran's real `Programs/BootStrapSampler`
   on identical inputs, dumped both flux arrays, sorted and compared --
   **zero differences** (the two arrays contain the identical multiset of
   values; a naive index-aligned diff showed nonzero only because the two
   dumps use different axis orderings, not because any value differs).
   Bonus: both sides also print an identical hardcoded `PXTRACE` debug line
   at voxel (14,15,29) -- matched to 8 decimals.
6. **SoFiA extraction itself, proven bit-exact between wasm and native**
   (`js/tools/repro_sofia.js`, run on the SAME `ReproBS10.fits` from #5):
   compared the WASM `sofia-wasm.js` build against the real native
   `third_party/SoFiA-2-master_2_5_1/sofia` binary -- **identical catalog
   row** (x, y, ell_maj, ell_min, ell_pa, kin_pa all matched to the
   catalog's own 6-decimal text precision). This was a previously-
   unexamined cross-platform boundary (wasm vs native compiled SoFiA-2 C
   code, not JS-port-vs-Fortran) -- worth remembering as a candidate for
   any FUTURE divergence too.

**Conclusion**: #5+#6 together prove the ENTIRE resample+SoFiA INPUT to
realization 10's fit -- the resampled cube AND the catalog-derived starting
guess -- is bit-identical between Fortran and JS. The divergence is
isolated to the amoeba fit itself (chi2 evaluation / particle generation
during optimization / convolution / simplex decisions), not anything
upstream. Given the fit's own RNG (`fitIdum`) is STATIC and NOT
realization-dependent (`bootstrap-realization-launcher.js`'s own comment:
"matches Fortran, RunWRKP never varies idum across realizations"), and
every particle-generation formula already checked out in earlier sessions'
extensive 2026-08-18 gasdev-desync bisection, the next step is a genuine
per-iteration amoeba trajectory trace for realization 10 specifically
(reusing `TRACE_OVERRIDE_IDUM`/`PARTTRACE`/`BISECTTRACE`, forcing the
EXACT verified starting guess into both a standalone `SingleGalaxyFitter`
run and a JS single-realization run) -- not yet done; this is a
substantially larger undertaking than anything fixed so far this session
and was not completed.

## FOUND + FIXED 2026-09-29 (continued): 3 real bugs in the pre-analysis stage, root cause still not fully closed

Picked the trace back up by loading the REAL payload the actual `--local`
run built (`js/DCPjobData/realization_payload.json`, copied out before
cleanup) and calling `runBootstrapRealization(10, payload)` directly in the
main thread (no worker pool) with `TRACE_DEBUG=1` -- a clean, non-
interleaved trace using the exact real inputs, no hand-transcription risk.
`js/tools/repro_realization10_trace.js` is this repro, kept for reuse.

**The trace immediately found the real divergence point, much earlier than
expected**: `FULLVECPARAM` call 1 (the FIRST simplex vertex, before ANY
amoeba iteration) already differs between platforms -- X=22.1069355
(Fortran) vs X=22.10097885 (JS), despite `idum` matching exactly at that
same call. So the divergence isn't in the fit's search trajectory at all --
it's in the STARTING GUESS itself, specifically the flux-weighted centre
estimate (`Iter_EstimateCenter`/`EstimateCenter`, `EstimateShape.f`) that
feeds it.

Traced one level deeper (added a temporary `TRACE_DEBUG`-gated print in
`EstimateShape.js`'s `estimateCenter`, and used the pre-existing
`PARITY_DEBUG` hooks already in `InitialAnalysis.js`): the flux-weighted
sum going into that centroid calc already differs -- JS's masked pre-
analysis cube sums to 0.8225861794126104, Fortran's to 0.208438516 (**NOT**
just a rounding difference: JS's own mask has **1041** nonzero pixels where
Fortran's SoFiA gives **1038**, on cubes and beams already independently
proven bit-identical). This is despite:
- The raw resampled cube itself matching exactly (still 14.072564761127083,
  same as `repro_resample.js`'s earlier proof).
- Native SoFiA and wasm SoFiA producing **pixel-for-pixel identical masks**
  when run on the SAME externally-written FITS file (`ReproBS10.fits`,
  Fortran's own resample output) -- confirmed by a direct numpy diff, zero
  differing pixels.

So wasm-vs-native SoFiA itself is NOT the problem (re-confirming the
earlier finding) -- the problem is specifically in the FITS file
`bootstrap-realization-launcher.js` itself hands to SoFiA. Diffed that
file's header against Fortran's own resampled-cube FITS output
(`ReproBS10.fits`) directly and found three real, confirmed bugs, all in
`js/src/BootstrapSampler/DataCubeFits.js`'s `dataCubeToFitsBytes` (which
builds every FITS file this pipeline hands to SoFiA) and its one caller in
`runBootstrapRealization`:

1. **`resampleBeam.beamPositionAngle` was never set** (`bootstrap-
   realization-launcher.js`, `runBootstrapRealization`'s own beam
   construction) -- stayed at `Beam2D`'s default of `0`, so every
   bootstrap realization's SoFiA call ran against a cube whose BPA FITS
   keyword claimed 0 deg instead of this galaxy's real ~12.3 deg. Fixed by
   reusing `observedBeam.beamSigma2` (already carries this value in
   radians, per `Beam.js`'s `beamSigmaVector[2]=beamPositionAngle`
   convention, and already serialized into the payload from the initial
   fit's own correctly-computed beam).
2. **`dataCubeToFitsBytes` wrote `beamPositionAngle` straight into the FITS
   `BPA` keyword without converting radians to degrees** -- a real,
   independent bug, MASKED for a long time by bug #1 above (0 rad
   coincidentally equals 0 deg, so this only became visible once bug #1
   started passing a real nonzero value through). Fixed:
   `BPA: (beam.beamPositionAngle || 0) * 180 / Math.PI`.
3. **`dataCubeToFitsBytes` hardcoded `CTYPE3: 'VELO-LSR', CUNIT3: 'km/s'`**
   regardless of the real cube's own convention -- confirmed directly
   against `WALLABY_J100336-262923_VelCube.fits`'s actual header:
   `CTYPE3=VOPT, CUNIT3=m/s`. Every FITS file this pipeline has EVER handed
   to SoFiA carried the wrong spectral-axis type/unit. Fixed by writing
   `CRVAL3`/`CDELT3` scaled to m/s (`dh`'s internal working unit is km/s)
   with `CTYPE3: 'VOPT', CUNIT3: 'm/s'`, matching Fortran's real output.

**Honest result**: all three are real, confirmed, independently-verified
bugs (each checked directly against Fortran's actual FITS header, not
guessed) -- but fixing them did **not** change the mask pixel count (still
1041, not 1038) or the final fit outcome for realization 10 (chi2, X/Y/Inc
all identical to before the fixes). So they were real correctness bugs
worth fixing regardless (every SoFiA call in this entire pipeline was
running against cubes with wrong beam angle AND wrong spectral-axis
convention), but they are NOT (by themselves) the explanation for this
specific mask-count divergence.

**Next lead, not yet chased**: a header diff also showed **`BITPIX -32`
(Fortran, real4) vs `BITPIX -64`** (JS, `cfitsio.writeImageDoubleWithHeader`
always writes double-precision) -- every FITS file this pipeline hands to
SoFiA is double-precision where Fortran's own output is single-precision.
If SoFiA's own noise/threshold statistics behave even slightly differently
reading float64 vs float32 pixel data, that could explain a genuinely
different segmentation despite identical underlying values. Checked
`cfitsio-wasm.js`'s exposed API: only `writeImageDouble`/
`writeImageDoubleWithHeader` exist, no float32 write path -- adding one
means touching the C driver source and rebuilding the wasm module, a
meaningfully bigger job than anything else fixed this session, and was
not attempted.

**Kept as reusable diagnostic tooling**: `js/tools/repro_realization10_
trace.js` (loads a saved `realization_payload.json` and calls
`runBootstrapRealization` directly, no worker pool, clean trace) and the
`PARITY_DEBUG`-gated trace prints added to `bootstrap-realization-
launcher.js` (mask nonzero count, resampleBeam.beamPositionAngle) --
these made this whole investigation tractable and should make the next
session's BITPIX lead much faster to check.

## CHASED 2026-09-29 (same day): BITPIX lead ruled out too -- root cause still open

Added a genuine float32 (BITPIX=-32) write path to the cfitsio wasm build,
to directly test the BITPIX=-32-vs-64 hypothesis above rather than leave it
as a guess:

- `third_party/cfitsio-4.6.3/wasm/cfitsio-driver.c`: added
  `cfits_create_image_float_wasm` (`FLOAT_IMG` instead of `DOUBLE_IMG`) and
  `cfits_write_image_data_float_wasm` (`TFLOAT` instead of `TDOUBLE`),
  mirroring the existing double-precision pair exactly. Header keywords
  still go through the same `TDOUBLE`/`TSTRING` write-key calls either way.
- `build.sh`: exported the two new driver functions, added `HEAPF32` to
  `EXPORTED_RUNTIME_METHODS`.
- Rebuilt via `emcc` (emsdk already installed at `~/DCP/Emscripten/emsdk`;
  `source emsdk_env.sh` first) -- built cleanly.
- `cfitsio-wasm.js`: added `writeImageFloatWithHeader`, mirroring
  `writeImageDoubleWithHeader` exactly (`HEAPF32.set` instead of
  `HEAPF64.set`, 4 bytes/element instead of 8).
- `DataCubeFits.js`'s `dataCubeToFitsBytes`: switched from
  `writeImageDoubleWithHeader` to `writeImageFloatWithHeader` -- confirmed
  via a fresh header read that the resulting FITS file now genuinely
  reports `BITPIX=-32`, matching Fortran.
- Also rebuilt `third_party/cfitsio-4.6.3/wasm/package/cfitsio-wasm.js`
  (the published-package bravojs bundle) via its own `build-bravojs-
  bundle.js`, so the real DCP-dispatch path picks this up too, not just
  `--local`.

**Result: no change.** The mask is still 1041 pixels (Fortran's SoFiA:
1038) on the exact same cube. Re-diffed the full FITS header against
Fortran's own output afterward: BPA/CTYPE3/CUNIT3/BITPIX now all match
exactly; the only remaining differences (BMAJ, BMIN, CDELT1/2, CRVAL1/2/3)
are sub-ULP text-formatting artifacts of the two FITS writers' own ASCII
keyword serialization (e.g. `CRVAL3 889520.562` vs `889520.568847656` --
the same float32 value, written with different numbers of decimal digits),
not real value differences. Confirmed the two catalogs' detected source
shares the EXACT SAME bounding box (`x_min..x_max, y_min..y_max,
z_min..z_max` = `16-30, 19-28, 10-31` in both) and the same f_min/f_max --
only n_pix (1038 vs 1041) and the resulting weighted centroid/ellipse
differ, meaning a handful of specific pixels right at the segmentation
threshold, inside an otherwise-identical detection, flip status between
Fortran's native SoFiA and the wasm build.

Also ruled out wasm-module state leakage as a red herring: ran the wasm
SoFiA module twice in a row, in the same process, against the same file --
identical n_pix both times, so it's not run-to-run nondeterminism.

**Where this leaves it**: with the FITS file now essentially byte-
equivalent to Fortran's own (same pixel values, same precision, same
header semantics), the remaining difference has to be inside SoFiA-2's
own C algorithm itself -- its internal noise estimation / smoothing /
thresholding arithmetic behaving differently when compiled to wasm
(Emscripten's musl libm) vs native (macOS's system libm), for pixels right
at a threshold boundary. That's a fundamentally different, much deeper
kind of investigation than anything else in this whole session: debugging
SoFiA-2's own upstream algorithm/build, not this project's JS-port-vs-
Fortran code. Not pursued further this session -- would need instrumenting
SoFiA-2's own source (a third-party dependency, not this project's code)
to find the exact internal computation responsible, most likely one of its
noise/RMS estimation or Gaussian-smoothing-kernel routines disagreeing at
the ULP level between musl and native libm for specific input values.

**Kept from this round**: the float32 write path itself
(`writeImageFloatWithHeader`) is a real, independent correctness fix
regardless of whether it solved this specific divergence -- every FITS
file this pipeline had ever handed to SoFiA was needlessly double-
precision where Fortran's real output is single-precision. Regression-
tested against the previously-fully-bit-exact `cloudDensity=20` case
(still 0.000000 max diff on every geometry field) and re-confirmed the
`cloudDensity=100` case is unchanged (still isolated to realization 10,
same magnitude) -- no regressions introduced.

## RESOLVED (root cause found, not fixable in this codebase) 2026-09-29: SoFiA itself definitively ruled out; the real cause is residual model-synthesis noise amplified by SoFiA's own threshold statistic

Pushed the bisection one level further per Dan's explicit instruction to
keep going. Two things closed this out:

**1. JS's own regenerated model cube compared directly against Fortran's**
(not a substitute, the ACTUAL cube each platform's own pipeline uses):
59,762/102,168 voxels bit-identical; the other 42,406 differ, but EVERY
one is in the noise floor (max value among them 2.3e-8, vs a real signal
peak of 0.0127 -- five orders of magnitude below anything meaningful).
This is the same FFT/convolution ringing already documented and accepted
as harmless in bug #7 earlier this session -- present independently on
BOTH platforms from non-associative summation, not new, not a bug. The
resampled cube inherits this same (even smaller, ~1e-19) residual.

**2. The decisive test Dan asked for**: fed the EXACT SAME cube file into
both native SoFiA and wasm SoFiA, for BOTH cubes (Fortran's own resample
output and JS's own resample output) -- a full 2x2 matrix:

|                 | Fortran's cube | JS's cube |
|-----------------|:--------------:|:---------:|
| **native SoFiA**|      1038      |   **1041**|
| **wasm SoFiA**  |      1038      |   **1041**|

The result depends ONLY on which cube goes in, never on which binary runs
it -- native SoFiA on JS's own generated cube gives byte-for-byte the SAME
catalog row as wasm SoFiA on that same file. **This definitively rules out
SoFiA's wasm build as the cause** -- it is fully deterministic and
behaviorally identical to the native build, confirmed directly, not
inferred from "everything else matched so it must be this."

**Actual root cause**: Fortran's and JS's independently-synthesized model
cubes agree to an already-accepted noise floor (~1e-15 to 1e-19 absolute,
utterly negligible next to real flux) but are not byte-identical -- and
SoFiA-2's own MAD-based noise/threshold statistic (`scfind.statistic=mad`,
`fluxRange=negative`, computed over a strided sample of the WHOLE cube,
i.e. including all those near-zero ringing voxels) is sensitive enough to
exactly which of those noise-floor values happen to be negative vs
positive to shift its own computed RMS by a tiny amount. Confirmed via the
full catalog diff: `std`/`skew`/`kurt`/`ell_maj`/`ell_pa`/etc. all shift
by small-but-nonzero amounts between the two cubes, not just n_pix --
consistent with a shifted threshold, not a segmentation-algorithm
difference. For this one faint, marginal source (the reproducible test
case is intentionally hard -- see the Sep-16 handoff), that tiny threshold
shift is enough to flip a handful of pixels right at the boundary, which
cascades into a different starting guess for the optimizer, which (once in
a harder cloudDensity=100 landscape) is enough to converge somewhere
measurably different.

**Where this leaves it**: this is not a bug in this project's JS port, the
Fortran code, or the DataCubeFits.js FITS-writing fixes made earlier today
(all of which were real and are kept) -- it's an inherent sensitivity of
SoFiA-2's own noise-estimation algorithm to floating-point noise this
project has already decided, elsewhere, is acceptable (bug #7's
resolution explicitly treated this exact magnitude of ringing as fine).
Closing this gap for real would mean either (a) making the model-cube
synthesis genuinely bit-identical at the ~1e-19 level too (a much harder
bar than anything targeted this session, and arguably chasing noise below
any physically meaningful threshold), or (b) patching SoFiA-2's own source
to make its noise statistic less sensitive to sub-threshold ringing (a
third-party dependency, out of scope for a JS/Fortran parity project).
Not pursued further -- this is the honest end of this particular thread.
