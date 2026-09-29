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
- **`fftw_wisdom_helper.o` isn't wired into the Fortran Makefile** -- `make
  all` fails the final link every time with an undefined-symbol error;
  manually relink both `BootStrapSampler` and `SingleGalaxyFitter` with it
  added to the object list (copy the exact failed link line from `make`'s
  own output, add `fftw_wisdom_helper.o`, then `mv` both binaries to
  `Programs/`).
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
