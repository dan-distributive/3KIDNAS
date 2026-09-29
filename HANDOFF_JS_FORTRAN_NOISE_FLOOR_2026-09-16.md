# Handoff: JS optimizer noise-floor asymmetry investigation

**Written for:** a fresh Claude Code session picking this up with no prior context.
**Date:** 2026-09-16. **Owner:** Dan.

## The one-sentence problem

Across large batches of bootstrap realizations on the same galaxy/seed, **Fortran's optimizer converges fast on every single realization** (Dan's own report: 200/200 fast), while **JS's optimizer occasionally burns all 5000 Nelder-Mead iterations and fails to converge** on a handful of them (the same batches show "several" JS blowups). This is a real asymmetry, not imagined — it needs a real explanation, not a reassurance that "optimizers are sometimes chaotic."

A large amount of investigation already happened (see "What's already been ruled out" below) without finding one single decisive bug that explains the asymmetry. A practical mitigation (a stall detector) was shipped that fixes the *symptom* (wasted compute time), but the *root cause* of why JS's noise floor is apparently higher than Fortran's is still open. **That's the job for this next session**: keep hunting for the actual numerical cause, using the leads at the bottom of this doc as a starting point.

## Background you need

This is one thread in a much larger, multi-session effort to get the JS/DCP port of a Fortran galaxy-fitting pipeline (`3KIDNAS`) to bit-exact numerical parity with the original Fortran. That larger effort has already found and fixed several real, confirmed cross-platform bugs this month:

- **`X**2.` (real-literal exponent) is not bit-identical to `X*X` in gfortran** — a genuine upstream bug, found via a `transfer()`-based standalone repro, swept across 26 sites in 14 files. Worth reporting to Nathan (the upstream author) but not yet done.
- **A PA (position angle) "kinematic convention" round-trip bug** — geometry was being re-derived from a lossy, display-rounded value on both platforms independently, instead of carrying the raw pre-conversion angle through. Fixed via a new `_RawGeom_v1.txt` companion output file.
- **Systematic `fdlibm`-substitute sweep** — JS's native `Math.sin/cos/atan/log/exp` use a different internal algorithm than gfortran's runtime, and can disagree by more than 1 ULP for some inputs. Ported fdlibm equivalents (`fdSin`, `fdCos`, `fdAtan`, `fdLog`, `fdExp`, `fdAtanh` in `js/src/StandardMath/fdlibm.js`) are supposed to be used at every relevant call site instead of the native `Math.*` functions. **This sweep was NOT 100% complete** — see "New leads" below, found via a fresh grep this session.
- Every individual float32 arithmetic op needs its own `Math.fround()` wrap to match Fortran's per-op rounding (`REAL(4)` semantics) — done throughout, but easy to get subtly wrong when porting a compound Fortran expression into more JS statements than the original.

Read `/Users/dandesjardins/DCP/3KIDNAS/JS_FORTRAN_PARITY_SESSION_2026-09-15.md` for the full history if you need it (it's long — UPDATE 17 and UPDATE 18 are the sections directly relevant to this handoff's topic). This handoff document is meant to let you skip re-reading all of that.

## The reproducible test case

- Galaxy: `WALLABY_J100336-262923` (real WALLABY test data already in this repo)
- `BootstrapSeed=1000`, `nBootstraps=5`, `PA_Estimate=81.271`, `Inc_Estimate=31.49`, `cloudDensity=20` (cdens)
- **Realization 4** (0-indexed, 5th of 5) is the reliably-reproducible pathological case. Confirmed reproducible byte-for-byte across many reruns given the same seed.

Command to reproduce the JS-only leg (fast, no Fortran needed for most of this work):
```bash
cd /Users/dandesjardins/DCP/3KIDNAS
CUBE=/Users/dandesjardins/DCP/3KIDNAS/3KIDNASTests/TestData/WALLABY_Test_sources/WALLABY_J100336-262923/WALLABY_J100336-262923_VelCube.fits
MASK=/Users/dandesjardins/DCP/3KIDNAS/3KIDNASTests/TestData/WALLABY_Test_sources/WALLABY_J100336-262923/SoFiA_J100336-262923_mask.fits
node js/tools/run_both.js \
  --seed 1000 --nBootstraps 5 --skip-fortran --local --cloudDensity 20 \
  --objName WALLABY_J100336-262923 \
  --cube "$CUBE" --mask "$MASK" \
  --pa 81.271 --inc 31.49 \
  --json /tmp/some_output.json
```
Look for the `realization 4:` line in stdout. As of this session it shows `converged=false` with a ~10s fit time (after the stall-detector fix landed — see below; it used to be ~100s before that).

To get the Fortran-only leg for direct comparison, drop `--skip-fortran` and add `--skip-js-dcp` instead; **also `unset` any `JS_*`/`TRACE_*` env vars first** (see "Careful: env var gotchas" below).

**IMPORTANT — do not run `node js/tools/run_both.js` with no arguments.** Its default galaxy (no `--objName`) is `WALLABY_J103538-484832`, not the one you want, and it will kick off a real, slow, unseeded end-to-end fit. Always pass `--seed`, `--objName`, `--cube`, `--mask`, `--pa`, `--inc` explicitly as shown above. (This surprised the previous session — see UPDATE 17 in the main session log for the story. Might be worth fixing `run_both.js` to require `--seed` or print usage when no args are given, as a separate small task.)

## The forced-simplex + forced-idum methodology (reusable tooling, already built)

This is the most valuable thing to carry forward: a working experimental apparatus to directly compare Fortran's and JS's objective function (`funk`/`tiltedRingModelComparison`) given **byte-identical inputs**, isolating "does the math agree" from "do the two platforms happen to explore the same territory."

### Why this was necessary
`funk()` isn't a pure function of the parameter vector alone — it also depends on the RNG's mutable `idum` state (particle placement inside the model is Monte Carlo). Two platforms' `idum` streams desync from each other over the course of a fit (pass 1 runs un-forced on both sides, and any tiny per-op difference can nudge which of several near-tied Nelder-Mead moves gets taken, which changes how many `ran2()` draws get consumed, which fully desyncs the two idum streams going forward). So forcing ONLY the parameter vector to match isn't enough — you have to force `idum` too, or you're comparing two different stochastic models, not the same one.

### The pieces (all currently in the codebase, env-var gated, zero cost when unset)

1. **`FORTRAN_SIMPLEX_DUMP_PATH`** (Fortran, `src/GalaxyAnalysis/GalaxyFit.f` around line 120-144): dumps the exact pass-2 starting simplex (all `nParams+1` vertices × `nParams` params, full float32 hex) to `<path>.<ObjName>` (e.g. `..._Bootstrap_4`) right after `MakeParamGuessArray` builds it, before `DownhillSimplexRun` consumes it.

2. **`JS_SIMPLEX_OVERRIDE_PATH`** (JS, `js/src/GalaxyAnalysis/GalaxyFit.js`): reads that same dump file and force-injects it into `paramGuesses` right before pass 2's `downhillSimplexRun` call, overwriting whatever JS's own `makeParamGuessArray` would have produced.

3. **`JS_OVERRIDE_REALIZATION_INDEX`** (JS, both `GalaxyFit.js` and `FullModelComparison.js`): **critical scoping fix**. Without this, the overrides above apply to EVERY `galaxyFit_Simple` call in the whole run — the initial fit AND all 5 bootstrap realizations — corrupting all of them. Set this to `4` (or whichever realization you're investigating) to scope the override to just that one.

4. **`state._simplexOverrideActive`** (JS, set in `GalaxyFit.js` right where the simplex override applies): **second critical scoping fix**. Realization-index scoping alone isn't enough — pass 1 shares the same `realizationIndex` as pass 2, and pass 1's own natural (unforced) `evalCount` can pass through the same numeric range later used to identify pass 2's vertices. This flag makes the idum override (next item) strictly pass-2-only.

5. **`JS_IDUM_OVERRIDE_SEQUENCE_PATH`** (JS, `FullModelComparison.js`): a file of newline-separated decimal `idum` values, one per pass-2 vertex, injected into `state.rng.state.ran2State.idum` right before each vertex evaluation. **How to build this file from Fortran's log**: `idum` is one continuously-advancing stream; the value Fortran's `TRACE` line prints is the POST-call state, so vertex k's INPUT idum is the value printed on the PREVIOUS call's TRACE line. Extract with:
   ```bash
   awk '/^FULLVEC call= /{call=$3} /^ TRACE /{ if(call>=X && call<=Y) print call, $3 }' FittingLog.txt
   ```
   where X..Y is the call-number range of the 14 (or however many params+1) vertex evaluations minus one (i.e. one earlier than the first vertex's own call number) through the second-to-last vertex's call number. **Do not assume you know the call-number range in advance** — see "gotcha" below.
   The code anchors dynamically now (`tiltedRingModelComparison._idumSeqStartCall`, set to whatever `evalCount` is the first time `state._simplexOverrideActive` is seen true) — you should NOT need a hardcoded `JS_IDUM_OVERRIDE_START` any more; that was an earlier, since-removed hack that silently misaligned when pass 1's natural length changed between reruns.

6. **`TRACE_DEBUG=1` + `TRACE_DEBUG_FINALVEC_FILE=/path`**: enables per-eval `FULLVEC`/`FULLVECPARAM` dumps (13-param vectors, hex + decimal) into a file suffixed `.r<realizationIndex>` (parallel workers otherwise interleave their output into one unreadable stream). Also enables `EVALCHI2 call=N chi2=... idum=...` lines (added this session) into the same suffixed file, and `EVALCHI2 call=N BADMODEL reason=...` when `badModelCheck` rejects a vertex (also added this session, with the specific reason — e.g. `SIGMA ring=3 of nRings=4 sigma=... sigUse=...` — logged via `badModelCheck._lastReason`).

### GOTCHA — Fortran's `TraceCallCounter` is shared across BOTH optimizer passes

Do not assume "pass 2's vertex 1 is call N" based on a naive read of the log. `TraceCallCounter` in `FullModelComparison.f` is a single running counter shared across pass 1's entire run (its own 14-vertex initial loop AND all its amoeba iterations) and pass 2. The previous session burned real time on a false lead here: it assumed pass 2's vertex loop started right after a call whose value happened to coincidentally match pass 2's vertex 1 (because pass 1's own converged best point was, by coincidence of VALUE not position, identical to what became pass 2's unperturbed base vertex). The correct way to find where pass 2's real vertex loop starts: look for the `ITER_F`/`CONVERGED_VECTOR`/`Param Guess Array Creation` markers that mark pass 1's actual completion, and take the FULLVEC call immediately after those. Or, more robustly: search the log for the exact hex value of a KNOWN vertex from your dump file (e.g. vertex 2's first param) and confirm it appears exactly once in the expected position, not earlier by coincidence.

### GOTCHA — Fortran's `TraceSwitch` is not an env var

Unlike JS's `TRACE_DEBUG`, Fortran's per-eval `FULLVEC` tracing is gated by a `TraceSwitch` value read from a fitting-options text file (`SingleFitRuntimeInputs.f`), NOT an environment variable directly. However, `FitDriverScripts/RunWRKP.py` already wires `WRKP_TRACE_DEBUG=1` (an env var) to set `TraceSwitch=1` in the generated options file for you — just set `WRKP_TRACE_DEBUG=1` in the environment before invoking `run_both.js`'s Fortran leg and it'll work. (A previous, much slower path was accidentally taken this session before this was found: manually patching `run_both.js`'s `writeCloudDensityOptionsFile` function to append the switch lines directly to the WRONG file — the fitting-options file, not the top-level `WRKP_Input.txt`-style file that `SingleFitRuntimeInputs.f` actually reads `TraceSwitch` from. Don't repeat that; use `WRKP_TRACE_DEBUG=1`.)

### Careful: env var gotchas

- Always `unset JS_SIMPLEX_OVERRIDE_PATH JS_IDUM_OVERRIDE_SEQUENCE_PATH JS_OVERRIDE_REALIZATION_INDEX TRACE_DEBUG TRACE_DEBUG_FINALVEC_FILE WRKP_TRACE_DEBUG` before a "natural"/unforced run, or leftover exports from a previous investigation will silently corrupt it.
- `--cloudDensity 20` must be passed explicitly to `run_both.js` to match the original test's cdens=20 (its own default is `400`, per `Inputs/SingleGalaxyTestFittingOptions_Base.txt`).
- Use `--skip-js-dcp` (not `--local`) to get a Fortran-only run; use `--skip-fortran --local` to get a JS-only (no dispatch) run.

## What's already been ruled out this session (with direct evidence, not assumption)

Each of these was checked by reading both implementations side by side or by direct instrumented measurement — not inferred:

1. **`amoeba`'s algorithm itself** (`src/GeneralMinimizationRoutines/DownhillSimplex.f` vs `js/src/GalaxyAnalysis/GalaxyFit.js`'s `amoeba` function) — faithful line-for-line port, including a previously-fixed shared bug (the shrink step's `y(i)` not being updated after a shrink, now fixed identically on both sides — see the "Historical bug fix" comment in the JS source).
2. **`ftol` mismatch** — identical on both platforms: `0.005` for pass 1, `/5` = `0.001` for pass 2. Confirmed by direct source read (`src/GalaxyAnalysis/GalaxyFit.f:95,113` and `js/src/GalaxyAnalysis/GalaxyFit.js:329,346`).
3. **`nParticles` truncation-boundary sensitivity** — compared Fortran's and JS's particle counts directly (via file-based `NPTRACE` logging, byte-identical forced inputs): exact match, all rings, both counts AND sigma hex values bit-identical.
4. **FFT rank-2 (row-then-column vs FFTW's native 2D plan) accumulation-order difference** — already investigated and documented in an EARLIER session, see `js/src/ConvolveCube/FFTW3WasmRank2.js`'s own header comment: checked against a real compiled FFTW ground-truth harness, found ~1-96 ULP differences in DOUBLE PRECISION (~1e-14 relative) — 1000-100000x smaller than a single float32 ULP, and the float32-rounding step immediately after swallows almost all of it. Ruled out as too small to explain the observed noise.
5. **Ring/particle iteration order** in `FillDataCubeByTiltedRing.f`/`.js` — confirmed identical ring-major/particle-minor nested order and per-cell accumulation sequence on both platforms (floating-point summation isn't associative, so a different order would matter, but the order is the same).
6. **`badModelCheck`'s `vSinI` used native `Math.sin()` instead of `fdSin`** — this WAS a real, previously-missed gap (every other trig call in the hot path already used `fdSin`/`fdCos`; this one boundary check was missed). **Fixed** this session (`js/src/CompareCubes/FullModelComparison.js`, look for the "BUG FIX (Dan, 2026-09-16)" comment). Tested: zero measurable effect on realization 4's specific outcome (identical final chi2 before/after the fix) — the inclination values encountered in this specific trajectory never happened to hit a value where native `Math.sin` and `fdSin` actually disagree. Real fix, kept, but not the (or at least not the whole) answer.

## The key measurement that reframed the question

Dan's original framing implied Fortran is somehow immune to this class of problem. It isn't — it's just been lucky (so far, on the runs checked). Extracted Fortran's own `Current tolerance` (rtol) trajectory for realization 4's pass 2, straight from its own `print` statement in `DownhillSimplex.f` (enable via `WRKP_TRACE_DEBUG=1`, then `grep "Current tolerance" FittingLog.txt`):

```
iter  rtol         y_hi        y_lo
0     0.00521      114386.99   113793.04
10    0.00291      113979.20   113648.02
14    0.00200      113875.20   113648.02
29    0.00261      113944.71   113648.02   <- got WORSE, not monotonic
31    0.00165      113828.73   113641.32
39    0.000996     113754.58   113641.32   <- converged, barely under ftol=0.001
```

Fortran's own rtol is non-monotonic and hovers in the same 0.1%-0.3% range JS gets stuck in, before happening to dip a hair under `ftol=0.001` at iteration 39. **This is a near-miss, not a comfortable margin.**

Separately, measured JS's own chi2 spread during a genuine (unforced) blowup on realization 4, over five 1000-call windows near the end of a 5000-iteration run:
```
calls 4079-4278: relrange=0.3860%
calls 4279-4478: relrange=0.4018%
calls 4479-4678: relrange=0.3907%
calls 4679-4878: relrange=0.4032%
calls 4879-5078: relrange=0.3850%
```
**Flat, not narrowing.** This is a genuine stuck-oscillation, not slow-but-real convergence — over 1000+ evaluations there is zero visible trend toward the ~0.1% ftol threshold.

**Conclusion so far**: `ftol=0.001` sits at or below the intrinsic Monte-Carlo-plus-rounding noise floor of this objective function for hard resampled datasets, on EITHER platform. The open question is specifically: **why does JS's noise floor sit measurably higher than Fortran's often enough to lose this near-miss more frequently across many bootstrap realizations, when every individual mechanism checked above matches almost exactly?**

## The stall-detector fix already shipped (a mitigation, NOT a root-cause fix)

`js/src/GalaxyAnalysis/GalaxyFit.js`'s `amoeba` function now has a stall detector: tracks the simplex's best value (`y[ilo]`); if it hasn't improved by more than `ftol/10` relative over a `max(300, 20*ndim)`-iteration window, exits early (same `noConvergence` signal as hitting `ITMAX`) instead of grinding to `ITMAX=5000`. Verified on the seed=1000/nBootstraps=5 test case: realizations 0-3 (normally converging) produce **bit-identical** chi2/convergence results before and after (detector never fires for them — it's conservative by design). Realization 4 drops from ~101s to ~10s (10x), still honestly reports `converged=false`.

**This does not need to be undone or distrusted** — it's a safe, well-tested, production-worthy change that directly addresses the wasted-compute-time symptom. But it doesn't explain WHY the noise floor differs between platforms, which is the actual open question. Look for the "Stall detector (Dan, 2026-09-16)" comment block in `amoeba` for the full rationale.

Fortran's `DownhillSimplex.f` has NO equivalent stall detector — intentionally left untouched, since Fortran doesn't appear to need one (per Dan's 200/200 report) and no evidence was found that Fortran's own `ITMAX=5000` is ever actually hit. Worth reconsidering only if that changes.

## New leads for this session, from a fresh repo-wide grep

Run this to reproduce (and to catch anything new introduced since this handoff was written):
```bash
cd /Users/dandesjardins/DCP/3KIDNAS
grep -rn "Math\.\(sin\|cos\|atan2\?\|log\|exp\|pow\|sinh\|cosh\|tanh\|asin\|acos\)\b" js/src/ js/bootstrap-realization-launcher.js js/tools/*.js \
  | grep -v "node_modules\|/verify/\|self-test\|require.main === module"
```

Most hits are comments, self-tests (guarded by `if (require.main === module)`), or fixture/test-data generators unrelated to the real fit pipeline (double-check anything you're unsure about by reading its enclosing function — several genuinely are just self-test helpers, e.g. `FFTW3JS/BluesteinSolver.js`, `DhtRaderSolver.js`, `RealEngine1D.js`'s native trig calls are all inside `if (require.main === module)` blocks, not production code). The following are **real, live-code, unverified candidates** worth checking first:

1. **`Math.acos` vs gfortran's `acos()` — CONFIRMED to disagree, but NOT YET FIXED. Read this whole item before touching acos again, it saves real time.**

   **Confirmed the gap is real** (checked same session, immediately after the handoff above was first written): built a standalone comparison harness — `gfortran -O0 -funroll-loops -cpp -fbounds-check -Wl,-ld_classic -ffp-contract=off` (this project's exact `FLAGS`, from `src/makeflags`) compiling a tiny Fortran program that reads a shared binary file of float32 inputs and writes `acos()` of each to another binary file, compared byte-for-byte against JS. Swept 200,003 float32 values across `[-1,1]` (uniform random, edge cases at ±1/0, dense sampling near the ±1 boundary, and oversampled `[0,1)` matching the two call sites' realistic range). **Result: 3158/200003 (1.579%) disagree, always by exactly 1 ULP.** Reproduction scripts (kept for reuse, not deleted): `acos_test.f`, `acos_gen.js`, `acos_compare.js`, `acos_candidates.js`, `acos_sf.js` in this session's scratchpad (path was `/private/tmp/claude-501/-Users-dandesjardins/234b575e-72da-4de3-abce-1a6aae9023b9/scratchpad/` — scratchpads are session-specific and may not survive into a new session, so **regenerate from the descriptions below rather than assuming those files still exist**).

   **Four fix hypotheses tried, ALL FAILED to close the gap — do not re-attempt these without a new idea:**
   - Native `Math.acos()`: 1.579% mismatch (the baseline, i.e. doing nothing).
   - A from-scratch double-precision port of the classic Sun fdlibm/FreeBSD-msun `e_acos.c` algorithm (rational Padé approximation, same family as `fdSin`/`fdCos`/`fdAtan` already in `fdlibm.js`), computed in JS doubles then rounded to float32: **identical** 1.579% mismatch, same exact set of failing inputs, same bit patterns. i.e. this hand-written port is behaviorally identical to native `Math.acos` for every tested input — V8's `Math.acos` already IS essentially this algorithm.
   - The same port with FMA fusion inserted at the natural Horner-chain multiply-add points (same technique that fixed `fdAtan` — see that function's "FMA NOTE" comment in `fdlibm.js` for the methodology): **identical** 1.579% mismatch. FMA contraction is not the cause here.
   - `acos(x) = atan2(sqrt(1-x*x), x)` (a structurally different reformulation, in case gfortran's acos is computed this way internally): **identical** 1.579% mismatch, same failing inputs.
   - A genuine **single-precision** fdlibm port (FreeBSD msun's `s_acosf.c`, real float32 constants and arithmetic throughout, not double-then-round — tested the hypothesis that gfortran's `REAL*4 acos()` calls `acosf()` directly): **made it WORSE**, 20281/200003 (10.14%) mismatch. So gfortran's `REAL*4 acos()` is NOT calling a naive single-precision routine either — if anything this rules out that theory rather than supporting it, and suggests double-precision-then-round (whatever exact double algorithm produces it) is closer to correct, just not fully there.

   **What this means**: whatever macOS's system `libm`'s `acos()` (which is what gfortran actually links against and calls) does internally, it is NOT exactly the standard fdlibm/msun algorithm in any of the four variants tried, in either precision. This is a **harder problem than sin/cos/atan were** (those all turned out to match one of the "standard" fdlibm approaches once found). Two remaining avenues, neither attempted yet:
   - **Binary reverse-engineering**: disassemble the actual `acos`/`acosf` symbols in `/usr/lib/libSystem.B.dylib` (or wherever they resolve on this machine — `otool -L` on the compiled Fortran test binary, or `dtruss`/`nm`, would find the exact path) to see what algorithm Apple actually ships. This is a materially bigger undertaking than anything else tried in this investigation so far — likely a session in its own right.
   - **Empirical curve-fitting**: since only 1.579% of inputs disagree and always by exactly 1 ULP, it may be tractable to find the specific decision boundary (e.g. is it correlated with a specific bit pattern in the rational approximation's rounding, like "cases where `p/q`'s exact mathematical value falls within `2^-53` of a float64 rounding boundary"?) without needing the exact algorithm — essentially reverse-engineering just the ROUNDING TIE-BREAK behavior rather than the whole algorithm. Not attempted; unclear if tractable.
   - Given the size of this specific gap (1 ULP, 1.6% of inputs, on a value used only ONCE per fit as an initial guess) relative to the effort needed to close it, **it may not be worth pursuing further as the explanation for the noise-floor asymmetry** — it's a real, confirmed, unfixed gap, but whether it's a MEANINGFUL contributor to the pathological-convergence problem (vs. e.g. mostly getting washed out by pass 1's own subsequent iterations regardless of the exact starting inclination) has not been tested. If you want to test that specific question before investing in a full fix: temporarily hack in a lookup-table override for the ~3158 known-mismatching float32 inputs (regenerate via the harness above) covering the actual range of `ellip`/`sqrt-ratio` values this pipeline produces, and see if realization 4's outcome changes at all. That would tell you whether this is worth the disassembly effort before doing it.

2. **`Math.pow(10.0, r.sigUse)` for log-mode surface density** (`js/src/CompareCubes/FullModelComparison.js` — search for `linearLogSDSwitch === 1`, and also duplicated in `js/bootstrap-realization-launcher.js:598,1270`). **Not active for the current test case** (this galaxy's fitting options use linear SD mode, switch=0) — but if this investigation moves to a different galaxy/config that uses log-mode SD fitting, check this. `10**x` for a general real `x` is a genuine transcendental (not the `X**2.` class of bug, which was specifically about integer-valued real literals) — worth verifying against Fortran's `10.**SigUse` the same way as `Math.acos` above.

3. **`constructModelBasedPV`'s `Math.cos(-angUse), Math.sin(-angUse)`** (`js/bootstrap-realization-launcher.js:1468`) — used only to build position-velocity diagram OUTPUTS (`pvMajorData`, `pvMinorModel`, etc., called around line 1524-1529), which looks like post-fit output/plotting generation, not part of the optimizer's hot loop. Lower priority, but unverified — confirm it's genuinely output-only (not fed back into anything the optimizer sees) before deciding whether it matters for this investigation.

## Where to look if the acos lead doesn't pan out

If a rigorous acos check comes back clean (bit-exact, like sin/cos), the next places to look, roughly in order of how directly they sit in the per-evaluation hot path (particle generation → binning → convolution → chi2, run thousands of times per fit):

- **`src/TiltedRingModelGeneration/SingleRingGeneration.f` vs `.js`** — the actual particle-placement loop. Already has fdlibm substitutes wired in (including a WASM-accelerated variant swapped in after warmup — see `fdlibmWasm.warmUp()` in that file). Worth double-checking the WASM variant (`fdSinWasm`, `fdCosWasm`, `fdAtanhWasm` in `../StandardMath/fdlibm-wasm.js`) actually produces bit-identical results to the plain JS fdlibm functions it replaces — an unverified WASM reimplementation could itself be a source of divergence nobody's checked from this specific angle.
- **`src/CompareCubes/CubeComparison.f` vs `js/src/CompareCubes/CubeComparison.js`** (the chi2/likelihood calculation itself) — already swept for native Math calls this session (clean), but not re-verified for f32-wrap-per-op discipline as rigorously as `amoeba`/`amotry` were.
- Consider whether the RNG itself (`ran2`/`gasdev` in `src/StandardMath/random.f` vs `js/src/StandardMath/random.js`) has any remaining subtle difference beyond the `gasdev` `v1**2.+v2**2.` fix already applied this month — the RNG is the literal source of every downstream difference once `idum` desyncs, so it's worth an extra-careful re-read even though it was already checked once.
- A more exhaustive, automated approach (not yet attempted): use the forced-simplex + forced-idum apparatus above, but instead of just comparing final chi2 per vertex, add fine-grained instrumentation that dumps the FULL model cube (or a meaningful checksum/hash of it) at each stage (post-particle-generation, post-binning, post-convolution) for ONE single forced-identical evaluation, on both platforms, and diff stage by stage. This is exactly the methodology that found the original `X**2.` bug (bisecting "call 99 -> call 1 -> call 3 -> resolved") — it hasn't yet been applied to a full single evaluation's every intermediate stage for this specific noise-floor question, only to the final chi2 output.

## Files touched this session (all still in place)

- `js/src/CompareCubes/FullModelComparison.js` — `fdSin` fix, `EVALCHI2`/`BADMODEL` file logging, `badModelCheck._lastReason`, idum-sequence override + scoping, `global.__TRACE_REALIZATION_INDEX`/`__TRACE_EVAL_COUNT` tags.
- `js/src/GalaxyAnalysis/GalaxyFit.js` — simplex override + scoping (`state._simplexOverrideActive`), **stall detector** (permanent, not diagnostic).
- `js/src/TiltedRingModelGeneration/SingleRingGeneration.js` — file-based `NPTRACE` particle-count logging.
- `src/GalaxyAnalysis/GalaxyFit.f` — `FORTRAN_SIMPLEX_DUMP_PATH` dump mechanism, `ITER_F` trace print.
- `JS_FORTRAN_PARITY_SESSION_2026-09-15.md` — UPDATE 17 and UPDATE 18 have the full blow-by-blow if you need more detail than this handoff gives.

All diagnostic instrumentation is env-var-gated (zero cost when the relevant env var is unset) — none of it needs to be stripped before continuing, but be aware it's there so you're not confused by unfamiliar code when reading these files.
