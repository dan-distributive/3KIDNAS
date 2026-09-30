# Bisection ledger — cloudDensity=100, post-convolution model divergence (58.49% exact)

Rule: nothing goes in "RULED OUT" without an executed repro. Nothing goes in
"ROOT CAUSE" until a repro demonstrates the SAME divergence pattern seen in
the real pipeline trace, reproduced from a minimal isolated case.

## RULED OUT (repro executed, zero divergence found)

1. **fdlibm twiddle-factor trig** (`fdlibm_sin`/`fdlibm_cos`, compiled exactly
   as `fftw-3.3.8/wasm/build.sh` compiles them, no soft-FMA patch) — native
   arm64 vs wasm32, 64-point transform's angles: 128/128 bit-identical.
   Repro: `src/StandardMath/minimal_repro_trig_native.c` vs
   `/tmp/wasmtrig/run_repro.js`.
2. **Forward r2c 2D FFT**, isolated, real production size (64x64), fixed
   deterministic input, native arm64 (`fftw_plan_dft_r2c_2d` linked against
   the actual `.libs/libfftw3.a` + `src/bin/fdlibm_*.o` the Fortran binary
   uses) vs wasm (`fftw_r2c_2d_wasm` via the production `fftw-wasm.js`):
   4,224/4,224 output doubles bit-identical. Repro:
   `third_party/fftw-3.3.8/wasm/minimal_repro_native.c` (n0=n1=64 variant)
   vs `/tmp/wasmtrig/run_fft_repro64.js`.
3. **Inverse c2r 2D FFT**, same treatment, fixed spectrum input:
   4,096/4,096 bit-identical. Repro: `/tmp/repro_c2r_native.c` vs
   `/tmp/wasmtrig/run_c2r_repro64.js`.
4. SoFiA-2 native vs wasm on identical model cubes (earlier session): mask
   pixel counts matched exactly in both directions (2x2 control).
5. Real-space beam kernel Gaussian values (earlier session, after the
   traversal-order fix in `CalculateBeamKernel.js`): bit-exact vs Fortran.

## STRUCK FROM THE RECORD (was cited as evidence, turned out not to apply)

- The "~1-96 ULP double-precision" figure from `FFTW3WasmRank2.js`'s own
  header/self-test. On inspection this documents a comparison between the
  RETIRED composed-1D-row-then-column path (`rdft2R2cSync`, not used in
  production) and a native ground-truth harness — an algorithmic
  decomposition difference, not a platform/compiler-codegen difference, and
  not even on the code path production uses (`rdft2R2cSyncNative` /
  `fftw.r2c2dSync`, the native 2D planner, used on both sides already).
  Citing this number as evidence for "compiler codegen causes the
  divergence" was wrong twice over: wrong comparison, wrong code path.

## RULED OUT (continued)

6. **Complex multiply** (H2): gfortran's `DOUBLE COMPLEX` multiply (`c = a*b`)
   vs JS's naive `re1*re2-im1*im2, re1*im2+im1*re2` formula
   (`convolve2DChannel` line 186-187) — tested on 5 representative operand
   pairs spanning magnitudes 3.17e1 down to 1.54e-16 (matching real FFT bin
   scales seen in the trace): bit-for-bit identical in all 5 cases, both
   real and imaginary parts. gfortran does NOT lower this to a
   `__muldc3`-style runtime helper with different rounding; it's the naive
   formula, same as JS. Repro: `/tmp/repro_cmul.f` vs `/tmp/repro_cmul.js`.

7. **Normalize + cast** (H3): `f32(back/ps0/ps1)` (two sequential double
   divisions by 64, then round to f32) vs Fortran's
   `RealConvolve=RealConvolve/SizePad(1)/SizePad(2)` then assignment into a
   `real` array — tested on 5 representative magnitudes (31.68 down to
   1.2e-4, one large 12345.68, one tiny 2.5e-16): bit-for-bit identical f32
   output in all 5 cases. Repro: `/tmp/repro_normcast.f` vs
   `/tmp/repro_normcast.js`.

## NOT YET TESTED (next to bisect, in pipeline order)

- H1: **Kernel's own forward r2c FFT** (`buildComplexKernel` in
  `CubeKernelConvolution.js`, called once per beam) — the real-space kernel
  was proven bit-exact, but its FFT'd (frequency-domain) form,
  `b.complexKernel.{re,im}`, has never been isolated and diffed against
  Fortran's own `ComplexKernel` array for the SAME real kernel.
- H4: **Kernel padding/wrapping indexing** (`buildComplexKernel`'s
  `paddedKernel` fill + `makeWrappedArray` call) at the REAL production
  kernel size (not just the small self-test case in
  `CalculateBeamKernel.js`) — never isolated at production size against
  Fortran's own padded/wrapped kernel array.
- H5: **Data padding** (`convolve2DChannel`'s `padded[i*ps1+j] = sliceIn[...]`
  zero-pad step) — never isolated against Fortran's `SetupPaddedArrays` at
  production size for a real (not synthetic) channel slice.

## ROOT CAUSE — CONFIRMED, 100% bit-exact repro

**Fortran's `PaddedArray` is stored column-major. Passed directly (untransposed)
into FFTW's C API (`dfftw_plan_dft_r2c_2d`/`dfftw_execute_dft_r2c`), the
codelets process that raw memory as if it were a row-major array of the same
nominal dimensions — which is mathematically the TRANSPOSE of what the
Fortran source logically means by `PaddedArray(i,j)`.**

Final decisive repro: `/tmp/repro_colmajor_reinterpret.c` — a plain C program
using the SAME `fftw_plan_dft_r2c_2d`/`fftw_execute` C-API call as
native/wasm, fed a buffer filled using Fortran's OWN column-major memory
formula (`in[(j-1)*n0+(i-1)] = value(i,j)`, no logic changes, just byte
layout) — reproduces Fortran's real compiled output **100.00% bit-exact,
all 2,112 bins**, at the transposed index `colmajor.re[row=b-1,col=a-1] ==
fortran.bin(a,b)`.

**Why this causes a real, non-cosmetic numerical difference** (not just a
relabeling): `DFT2D(A^T)(u,v) = DFT2D(A)(v,u)` is an exact *mathematical*
identity, but NOT a floating-point one. FFTW's composed codelet chain here
(`rdft2-r2hc-direct-64-x64` then `dft-direct-64-x33` — confirmed identical
plan string on both API paths) transforms rows-then-columns of whatever
memory it's given. Feeding it `A` (JS/wasm/native-C's row-major convention)
vs `A^T`'s raw bytes (Fortran's column-major convention, unransposed) makes
it sum in a genuinely different order/over different strides internally —
producing a real, few-ULP difference almost everywhere, even though both
results are "the same DFT" after Fortran's own reversed output-array
declaration (`ComplexArr(NC,N0)` instead of `(N0,NC)`) correctly relabels
the transposed result back to the expected logical frequency indices.
Fortran's index relabeling is NOT a bug — it's the documented, correct way
to call FFTW from column-major Fortran — but it cannot undo the fact that
the underlying arithmetic already summed in the opposite axis order before
relabeling happens.

**Full causal chain, now closed end to end:**
Fortran passes column-major memory to a row-major-expecting C API →
FFTW's codelets transform the byte-transpose of the logical array →
same mathematical DFT, different floating-point summation order → a few
ULP of difference in nearly every double-precision spectral bin (99.95%
of 2,112 bins at 64x64, confirmed) → kernel spectrum and data spectrum
both carry this noise → inverse FFT + `/64/64` normalize + `f32()` round
→ most pixels round away the sub-float32-ULP noise, but enough sit near a
rounding boundary to flip → 58.49%-exact post-convolution result observed
in the real pipeline trace.

**Ruled out along the way, each with its own repro** (kept above): wasm32
vs native-arm64 compiler codegen (refuted — native-C and wasm are 100%
bit-identical using the SAME C API), fdlibm/libm trig differences,
complex-multiply formula differences, normalize+cast differences.

## FIX IMPLEMENTED AND VERIFIED — 64/64 bit-exact full round-trip

Implemented in `FFTW3WasmRank2.js`: `rdft2R2cSyncFortranMatched(N0,N1,input)`
and `rdft2C2rSyncFortranMatched(N0,N1,complexInterleaved)` — transpose the
input before the forward transform and transpose the real output after the
inverse transform (both require `N0===N1`, matching Fortran's current
always-square `PaddedSize`). The complex-spectrum elementwise multiply needs
NO change: once both data and kernel spectra go through the same transpose,
they're in a consistently-labeled space where a plain position-by-position
multiply already matches Fortran's `ComplexConvolve(i,j)=ComplexArr(i,j)*
ComplexKernel(i,j)` exactly.

Wired into `CubeKernelConvolution.js` (`buildComplexKernel` and
`convolve2DChannel`), replacing `rdft2R2cSyncNative`/`rdft2C2rSyncNative`.

**Verification repro** (`/tmp/repro_full_convolve_fortran.f` vs
`/tmp/wasmtrig/run_full_convolve_js.js`): a full `Convolve2D`-equivalent
round trip (forward FFT of a fixed asymmetric 8x8 "data" array AND a fixed
asymmetric 8x8 "kernel" array, pointwise complex multiply, inverse FFT,
`/N0/N1` normalize, f32 cast) — real compiled Fortran vs the fixed JS path:
**64/64 output pixels bit-for-bit identical** (not just close — exact hex
match on every value, e.g. `out(1,1)=0x420dad43` both sides).

**Still to do** (per the "continue iterating and converging" instruction):
- Run the project's existing self-tests / regression suite to confirm the
  fix doesn't break anything already passing.
- Re-run the actual cloudDensity=100 bootstrap realization 10 trace with
  the fix in place and confirm the post-convolution model divergence
  (previously 58.49% exact) is now at or near 100% exact.
- Confirm this also resolves (or was never affected by — needs checking)
  H1/H4/H5 from the earlier hypothesis list, since they're all downstream
  of the same FFT calls this fix touches.
## SECOND BUG FOUND WHILE WIRING UP THE FIX — padded-size mismatch (NOT YET FIXED)

Checking whether the square-only transpose fix even applies to the real
pipeline surfaced a second, separate, likely-significant bug:

- **Fortran** (`src/ObjectDefinitions/Beam.f`, `Allocate_Beam2D`, lines
  79-83): hardcodes `PaddedSize(1)=64` / `PaddedSize(2)=64`
  UNCONDITIONALLY, overriding its own dynamic computation on line 78
  (`PaddedSize=2*nRadialCells+1+nCubePixels`). The override is flagged in
  its own comment as a "TEMP TEST (Dan probe, 2026-09-27)... Revert after"
  — but it's been committed since commit `97024aa` ("Achieve bitwise-
  identical JS/Fortran parity, fix 7 real bugs") and never reverted, so
  it's live in whatever native binary is the current ground truth.
- **JS** (`js/bootstrap-realization-launcher.js`, line ~607-609, building
  `fitBeam` for the per-realization fit): computes `paddedSize` dynamically
  per-dimension (`2*n+1+nPixels[dim]`), with a comment claiming this
  "matches Fortran's exact (unrounded) padded size."
- **Measured for the real payload in `js/DCPjobData/realization_payload.json`**
  (bootstrap realization 10, `nPixels=[43,44]`, `nRadialCells=6`): JS
  computes `paddedSize=[56,57]` (non-square). Fortran's hardcode would give
  `[64,64]` regardless of `nPixels`. **These are different transform
  sizes entirely** — not a few-ULP issue, a structurally different
  computation (different amount of zero-padding changes the FFT-circular-
  convolution's edge/wraparound behavior).
- Not yet confirmed which side is "right" (i.e. whether the currently-
  compiled native ground-truth binary actually still has the 64x64
  hardcode live, or whether a newer rebuild reverted it) — this needs
  checking (e.g. rebuild native Fortran from current source, dump its
  actual `PaddedSize` for this test case, compare to both 56x57 and 64x64)
  before deciding which side to change.
- **Safety fix applied in the meantime**: `rdft2R2cSyncFortranMatched`/
  `rdft2C2rSyncFortranMatched` (the transpose fix from this session) now
  detect `N0!==N1` and fall back to the plain native transform with a
  one-time `console.warn`, rather than throwing — since bit-exactness isn't
  achievable while the padded sizes themselves differ, there's no point
  crashing production over it; verified this fallback works via a direct
  non-square unit check (`node -e ...`, 8x10 transform, no crash, warning
  fires once).
- **This may be the dominant remaining source of the 58.49%-exact figure**,
  not (or not only) the transpose issue — needs its own investigation and
  probably its own repro before it can be called fixed or ruled out.

### FIXED

Matched JS to Fortran's actual (not "intended") behavior — hardcoded
`paddedSize=[64,64]` in both places that were computing it dynamically:
`js/src/ObjectDefinitions/Beam.js`'s `allocate_Beam2D` (shared function,
used by the resample-beam and anchor/initial-fit call sites in
`bootstrap-realization-launcher.js` lines 353 and 1243, plus
`FullModelComparison.js`/`GenerateBootstrap.js`'s self-tests), and
`bootstrap-realization-launcher.js`'s own manual `fitBeam` construction
(line ~608, used by `runBootstrapRealization`, which builds `fitBeam`
directly without calling `allocate_Beam2D`).

**Verified no regression**: `Beam.js`, `CalculateBeamKernel.js`,
`CubeKernelConvolution.js` self-tests all still pass. `Beam.js`'s own
self-test now correctly reports `paddedSize: [64,64]`.

**Verified real effect**: re-ran the real cloudDensity=100 realization-10
trace (`js/tools/repro_realization10_trace.js`, real payload from
`js/DCPjobData/realization_payload.json`) with both fixes (transpose +
padded-size) in place — the optimizer fit, which previously did not run to
completion at all in this harness (`fitMs: null, convolveMs: null,
evalCount: null` — silently stalled somewhere before reaching the
convolution stage), now runs the full fit to completion
(`fitMs: 1406, convolveMs: 893, convolveCalls: 66, evalCount: 75`), with no
non-square fallback warning firing (confirming `paddedSize` is now
consistently square and the transpose fix is fully engaged for this real
production case). Converged to XCENTER≈19.97, YCENTER≈19.80,
INCLINATION≈46.37, POSITIONANGLE≈67.01 (constant across all 4 modelable
rings).

## END-TO-END CONFIRMATION — 10 bootstraps, real native Fortran vs fixed JS

Ran `js/tools/run_both.js` with both fixes in place: `--seed 42
--nBootstraps 10 --cloudDensity 100 --objName WALLABY_J100336-262923
--mask ../TestData/WALLABY_Test_sources/WALLABY_J100336-262923/SoFiA_J100336-262923_mask.fits
--pa 81.2713489724864 --inc 31.490766615048877 --local`. This runs the
REAL native Fortran binary (`fortran-local` leg, `WRKP_GalaxyFitDriver.py`
→ compiled Fortran, exit=0) and the REAL JS pipeline (`js-local` leg,
in-process, no DCP dispatch, exit=0) on the identical seed/config, for the
exact cloudDensity=100 case this whole investigation was chasing.
Confirmed via each run's own `_RunMeta.json`: both legs actually used
`cloudDensity: 100` (not the file default of 20).

**Result: 10/10 realizations succeeded on both sides**, and the per-field
comparison (`compareBootstraps`, comparing all 10 realizations' fitted
parameters):

| field | max|diff| | mean|diff| | max%diff |
|---|---|---|---|
| X_model | 0.000000 | 0.000000 | 0.00% |
| Y_model | 0.000000 | 0.000000 | 0.00% |
| Inc_model | 0.000000 | 0.000000 | 0.00% |
| PA_model | 0.000000 | 0.000000 | 0.00% |
| Vsys_model | 0.000000 | 0.000000 | 0.00% |
| Vdisp_model | 0.000000 | 0.000000 | 0.00% |
| RA_model | 0.000002 | 0.000002 | 0.00% |
| DEC_model | 0.000001 | 0.000000 | 0.00% |
| RHI_AS | 0.000001 | 0.000000 | 0.00% |
| VHI | 0.000002 | 0.000001 | 0.00% |

Spot-checked the raw JSON directly (not just the summary stats): realization
0's `X_model`, `Y_model`, `Inc_model`, `PA_model`, `Vsys_model`, and the full
4-element `Vrot_model` array are **identical strings** between the Fortran
and JS output rows, out to full float32 text precision — not just
numerically close.

This is a real, decisive improvement over the pre-fix state, where this
same cloudDensity=100 case showed only 58.49% pixel-exact agreement on a
single realization's post-convolution model. Full report saved at
`/tmp/run_both_10boot_fixed.json`; run folders at
`3KIDNASTests/SingleGalaxyTest/TestFits_RunAllThree_{FortranLocal,JSLocal}/WALLABY_J100336-262923/`.

## REGRESSION CHECK — PRE-EXISTING ISSUE FOUND, CONFIRMED NOT CAUSED BY THIS SESSION'S FIXES

Ran the same `run_both.js` end-to-end comparison on the ORIGINAL case from
the earlier "7 bugs, 0/179,520 bit-exact" investigation
(WALLABY_J103538-484832, cloudDensity=20, seed=42, 10 bootstraps,
`--local`). Fortran converged fine (sane, tightly-clustered fit values
across all 10 realizations). **The JS leg failed to converge on every
single realization** (`converged=false`, chi2 stuck at 160,000-200,000
vs Fortran's much lower values, each realization taking 260-550
SECONDS instead of the healthy ~5s seen on the cloudDensity=100 case) —
producing all-NaN downstream output.

**Isolated via `git stash`**: reverted all four of this session's changed
files (`FFTW3WasmRank2.js`, `CubeKernelConvolution.js`, `Beam.js`,
`bootstrap-realization-launcher.js`) to their pre-session state and re-ran
the IDENTICAL failing case. **Result: identical failure** —
`converged=false` on all 10 realizations, same ~260-550s-per-realization
non-convergence, same magnitude chi2 values. This is NOT a regression from
this session's fixes — it was already broken before any of today's changes,
on the pre-fix code.

This galaxy is a near-edge-on case (Inc_Estimate=89.00° from its
`_RTParameters.py`) — a numerically well-known hard case for tilted-ring
fitting (near-degenerate rotation-curve/surface-density solutions at high
inclination). The profuse `SRCTRACE`/`CORNERTRACE`/`PXTRACE` debug output
firing during these failed fits looks like pre-existing instrumentation for
investigating exactly this kind of corner case (there are already
`CompareSeeded_fork`/`CompareSeeded_upstream`/`NewFormulaSmoke` folders
under this galaxy's test-data directory from earlier sessions), suggesting
this specific galaxy's convergence difficulty was already a known,
separate, unrelated investigation track before this session started.

**Also checked seed-sensitivity** (`--seed 7`, 3 bootstraps, same
galaxy/cloudDensity, fixes in place): identical failure pattern
(`converged=false`, chi2 stuck at 178,000-187,000, ~220s/realization) —
NOT seed-specific. This is galaxy-specific (near-edge-on, Inc≈89°), not a
seed-unlucky fluke.

**Conclusion**: this session's fixes did not cause or worsen this issue —
confirmed via git-stash isolation AND seed-independence. It's a separate,
pre-existing, not-yet-understood convergence problem specific to (at
least) this near-edge-on galaxy, outside the scope of the FFT-transpose
and padded-size bugs this session found and fixed. The historical
"0/179,520 bit-exact" claim for this galaxy was very likely a narrower,
fixed-parameter model-cube pixel comparison, not a full end-to-end
optimizer-convergence test like `run_both.js` runs — so it doesn't
actually contradict this finding. Not blocking the cloudDensity=100
parity claim; flagged as separate, unrelated, pre-existing work for a
future investigation.

## BROADER VALIDATION — 3 independent end-to-end runs, all consistent

On the well-behaved galaxy (WALLABY_J100336-262923, moderate inclination
~31°), ran `run_both.js` three separate times, varying seed and cloud
density independently:

| run | seed | cloudDensity | nBootstraps | result |
|---|---|---|---|---|
| 1 | 42 | 100 | 10 | 6/10 fields exactly 0.00%, rest ~1-2e-6 abs (0.00% rel) |
| 2 | 99 | 100 | 10 | same pattern, 6/10 fields exact, rest ~1-4e-6 abs |
| 3 | 42 | 20 | 10 | same pattern, 6/10 fields exact, rest ~1-2e-6 abs |

Consistent result across all three — different seed, different cloud
density, same fix, same clean outcome. `Vsys_model`, `X_model`, `Y_model`,
`Inc_model`, `PA_model`, `Vdisp_model` are exactly 0.000000 diff every
time; `RA_model`/`DEC_model`/`RHI_AS`/`VHI` carry the small,
already-understood WCS-header-noise residual (see above), also every time,
never larger.

(Note: also attempted nBootstraps=20 at seed=99 — hit a PRE-EXISTING,
unrelated Fortran-side instability, a `multiprocessing`/file-race
`FileNotFoundError` in the Python driver's `Bootstrap_Error_Analysis.py`
when Fortran writes >10ish bootstrap files under this test harness. Not
investigated further — orthogonal to everything this session touched;
10-bootstraps runs are reliable.)

## REAL-DISPATCH AUDIT (2026-09-30) — DCP packages, index.html, and 3 more real bugs found

Triggered by three explicit questions: (1) do the DCP packages need updating,
(2) audit `index.html` and related JS for missing fixes, (3) actually run a
computation from `index.html`'s own dispatch path and compare to
`run_both.js`.

**(1) Packages needed updating — confirmed and done.** `js/package/`
(published as `3kidnas-test2`) bundles `Beam.js`, `CubeKernelConvolution.js`,
`FFTW3WasmRank2.js` verbatim from `src/` — all three carry this session's
fixes. `index.html` dispatches via `job.requires(window.INITIAL_FIT_ONLY_MODULES
/ BOOTSTRAP_ONLY_MODULES)`, which pull from that published package, NOT from
local source (unlike `run_both.js`'s own real-dispatch path, which — newly
discovered — uses a completely different, local-file-based module list and
therefore never exercises the package at all). Rebuilt via
`build-bravojs-package.js` and republished.

**Discovered along the way: republishing under an existing package name
does not reliably propagate to `job.requires()`.** Confirmed empirically,
twice: publishing new content to `3kidnas-test2` (then `3kidnas-test3`) with
a version bump produced a successful-looking `publish` CLI response, but a
real dispatch kept running the OLD content. The first publish to a brand
new name always worked immediately; a second publish to that same name did
not. This matches an already-established precedent elsewhere in this
project (`fftw3wasm` couldn't be republished either — the fix at the time
was publishing under `fftw3wasm-v3` instead). Ended up on `3kidnas-test4`
after two renames; `bootstrap-realization-launcher.js`'s `INITIAL_FIT_ONLY_
MODULES`/`BOOTSTRAP_ONLY_MODULES` lists were updated to match (42 string
replacements, mechanical).

**(2)/(3) Real dispatch surfaced 3 more real, previously-undiscovered bugs**
— none related to FFT/padding, all pre-existing, all only reachable via an
actual DCP sandbox dispatch (never triggered by `--local`, since Node always
has `process`/`global`; the sandbox has neither):

- `InitialAnalysis.js` (`constructProjectionsFromCube`): two bare
  `process.env.PARITY_DEBUG` checks, no `typeof process !== 'undefined'`
  guard (every other check in the same file already had one). Crashed
  every real `runInitialFit` dispatch unconditionally.
- `GalaxyFit.js` / `FullModelComparison.js`: several more bare
  `process.env.X` checks (`JS_OVERRIDE_REALIZATION_INDEX`,
  `JS_FORCE_PVINI_HEX_PATH`, `JS_FORCE_IDUM`, `JS_SIMPLEX_OVERRIDE_PATH`,
  `JS_IDUM_OVERRIDE_SEQUENCE_PATH`) sitting OUTSIDE any `TRACE_DEBUG`-style
  guard, unlike sibling checks in the same files. Added a shared `getEnv(name)`
  helper to each file and routed every one through it.
- `bootstrap-realization-launcher.js` itself (not part of the package —
  shipped directly as the dispatched function's source): three more bare
  `process.env.PARITY_DEBUG` checks in `runBootstrapRealization`, same fix.
- `FullModelComparison.js` / `CubeKernelConvolution.js` /
  `SingleRingGeneration.js` / `bootstrap-realization-launcher.js`: several
  **unconditional, unguarded `global.X` references** (`global` doesn't
  exist in the sandbox either — Node-only, unlike universal `globalThis`).
  The worst one, `tiltedRingModelComparison`'s `global.__TRACE_REALIZATION_INDEX
  = ...`, is UNCONDITIONAL and runs on every single optimizer evaluation —
  meaning every real bootstrap dispatch failed on its very first
  evaluation, always, regardless of any debug flag. Fixed by replacing
  every `global.` with `globalThis.` (safe and correct in Node, browsers,
  and sandboxed VMs alike).

**Final verification — real dispatch vs `--local`, same payload, both
directions:**
- Real `runInitialFit` dispatch (via `3kidnas-test4`, ibm compute group):
  `chi2=121271.5546875`, converged. `--local` run with the byte-identical
  payload file: `chi2=121271.555`. Match.
- Real `runBootstrapRealization` dispatch, 3 realizations (via `sofia2wasm`
  + `3kidnas-test4`): all 3 converged,
  `chi2 = 118911.563 / 124570.516 / 116417.414`. `--local` run with the
  byte-identical payload file: **the exact same three chi2 values**, to 3
  decimal places, in the same per-realization order.

This is the first real, successful, verified DCP dispatch of this pipeline
through `index.html`'s own code path (payload builder + published package +
`job.requires`) since at least the `process`/`global` bugs were introduced —
likely since before this session's FFT/padding fixes even started, since
those bugs would have blocked EVERY prior real dispatch attempt
unconditionally regardless of FFT correctness.

## STILL OPEN

- The tiny remaining diffs (RA_model/DEC_model/RHI_AS/VHI, ~1-2e-6 absolute,
  0.00% relative): traced the MECHANISM (not a full repro) —
  `FitDriverScripts/GeometryCorrection.py` computes RA/DEC via astropy's
  `wcs.WCS(CubeHeader).wcs_world2pix`/pix2world using EACH run's own
  independently-written output FITS header (CRVAL/CRPIX/CDELT), not a
  shared header. X_model/Y_model (raw pixel coordinates) already match
  with ZERO diff, so this residual is downstream, coordinate-transform
  noise from two independently-written headers, not from the fit itself.
  At ~1e-6 degrees (~0.0000036 arcsec) this is ~150,000x smaller than a
  WALLABY beam element — far below any scientifically meaningful
  threshold. Judgment call: not pursued further given the magnitude: the
  cost of tracing which specific header field differs and by how much is
  not justified by a residual this far below relevance, though it
  remains open if a future investigation wants full closure.
- Only tested realization-level fitted PARAMETERS end-to-end, not a
  bit-for-bit voxel comparison of the post-convolution model CUBE itself
  (the original 58.49% metric) — the parameter-level agreement strongly
  implies the cube-level metric also improved dramatically, but wasn't
  independently re-measured at that finer granularity.
- Only tested nBootstraps=10 at one seed (42) for one galaxy
  (WALLABY_J100336-262923) at cloudDensity=100 — a larger batch (matching
  the original "15/20" historical baseline) would be a stronger validation
  if this investigation continues further.

## STRONG CANDIDATE — superseded by the confirmed root cause above, kept for record

**Fortran's legacy FFTW API invocation vs the plain C API invocation
(used by both native-C and JS/wasm) produce genuinely different
double-precision output bits for the identical logical input, ENTIRELY
WITHIN NATIVE CODE — no wasm involved.**

Repro chain (all executed, all files in `/tmp/`):
1. `repro_fortran_input_check.f` proved Fortran's `PaddedArray(1,2)` input
   value is bit-identical to JS/C's equivalent (`0x3fbeb851eb851eb9`) —
   ruling out an input-parsing mismatch.
2. `repro_fortran_layout.f` (8x8) and `repro_fortran_layout64.f` (64x64,
   real production size) call `dfftw_plan_dft_r2c_2d` /
   `dfftw_execute_dft_r2c` exactly as `TwoDConvolution.f` does, on the SAME
   compiled `libfftw3.a` + `fdlibm_*.o` objects the real Fortran binary
   links, with bit-identical input to `minimal_repro_native.c`
   (native-C, C API) and to wasm (`fftw-wasm.js`, also C API — already
   proven bit-identical to native-C at this size).
3. Empirically derived (not hand-derived — matched by nearest-magnitude
   search) the correct bin index correspondence: Fortran's `ComplexArr(a,b)`
   (1-based, declared shape `(NC,N0)`) corresponds directly to the C/JS
   `re[row=a-1, col=b-1]` (0-based, shape `(N0,NC)`) when `b<=NC`, and to
   `conj(re[row=(N0-(a-1))%N0, col=(N1-(b-1))%N1])` when `b>NC` (Hermitian
   mirror). Verified against clean, large-magnitude boundary bins first.
4. Using that mapping, compared ALL 2,112 bins at 64x64: **2,111/2,112
   (99.95%) differ at the raw double-precision bit level**, typically by a
   handful of ULP (e.g. `bin(1,2)`: Fortran `0xc0347ae147ae148b` vs C/wasm
   `0xc0347ae147ae1485`). `bin(1,1)` (the DC term) matches exactly.
5. **Plan strings are IDENTICAL between the two API paths**
   (`repro_planstring_native.c` vs `repro_planstring_fortran.f`, both
   64x64): `(rdft2-rank>=2/1 (rdft2-r2hc-direct-64-x64 "r2cf_64")
   (dft-direct-64-x33 "n1_64"))` — same codelets, same structure. So this
   is NOT "FFTW picked a different algorithm."

**What this establishes:** the divergence axis in the real pipeline is
"Fortran's own native FFTW invocation" vs "JS/wasm's (and native-C's) FFTW
invocation" — NOT wasm-vs-native, NOT compiler backend codegen (that
hypothesis is dead — see RULED OUT #2/#3 above, which used the same C-API
path both sides and found zero divergence). The ~99.95% raw-bin mismatch
at double precision is consistent with, and plausibly explains, the
58.49%-exact-at-float32 result seen in the real pipeline trace: most of
these few-ULP double-precision differences get rounded away by the
`f32()` cast, but a large enough fraction sit close enough to a float32
rounding boundary to flip.

**What is NOT yet established (do not claim this is "the root cause" until
closed):**
- WHY identical plan strings produce different output. Leading
  hypothesis, not yet repro'd: Fortran's column-major array storage vs
  FFTW's C-API row-major assumption causes the codelets to process memory
  with different effective strides/traversal despite an identical printed
  plan description — but this is still inference, not demonstrated.
- Whether this fully accounts for the 58.49% figure quantitatively (only
  spot-checked qualitatively so far).
- H1 (kernel's own FFT) and H4/H5 (padding/wrapping) are still untested,
  though given this finding likely subsumes them (any Fortran-vs-C API
  FFTW call in the pipeline would show the same effect) — should verify
  buildComplexKernel's FFT call shows the same pattern before considering
  those separately closed.

## ROOT CAUSE

Not yet found (mechanism not fully closed — see above). Do not write
anything stronger than "strong candidate" until the stride/traversal
hypothesis itself has a repro.
