# Upstream sync tracker: NateDeg/3KIDNAS (Dev branch) vs. local fork

Tracks reconciliation between the local Fortran/Python pipeline and
[NateDeg/3KIDNAS](https://github.com/NateDeg/3KIDNAS) `Dev` branch (216 commits
as of this doc's creation), plus propagating any adopted changes into the JS
port (`js/`). Background/methodology: `~/.claude/plans/breezy-launching-nova.md`.

**Upstream reference clone**: `/Users/dandesjardins/DCP/3KIDNAS_upstream/`
(plain checkout, not a git remote of this repo). Refresh with
`git -C /Users/dandesjardins/DCP/3KIDNAS_upstream pull`.

**Status values**: `not started` / `no action needed` (upstream already
functionally equivalent, or diff is local-only) / `needs review` (real
unabsorbed upstream change) / `flagged -- needs a decision` (real upstream
change with a real risk/tradeoff, do not apply without explicit sign-off) /
`applied+verified`.

Of 64 local Fortran files, 43 are byte-identical to upstream and are not
listed. **Full triage pass complete as of 2026-08-17** -- every file below has
been read in full, not just classified by keyword.

## Fortran (`src/`) -- triage complete

| File | Status | Note |
|---|---|---|
| `TiltedRingModelGeneration/TiltedRingModelGeneration.f` | **applied+verified (plumbing only)** | `BuildTiltedRingModel` gained `Noise,DC,BUse` args + new `CalcAvgChanPerPix` subroutine, ported to JS (`js/src/TiltedRingModelGeneration/TiltedRingModelGeneration.js`) and confirmed bit-matching Fortran. Signature/plumbing is live on both platforms; see `SingleRingGeneration.f` below for why the formula itself isn't switched on yet. |
| `TiltedRingModelGeneration/SingleRingGeneration.f` | **deferred (formula), applied+verified (plumbing)** | `Ring_CalcNumParticles` formula changed upstream: `DensMultiplications=CloudSurfDens*((Sigma/Noise)**cmode)` (was `*(Sigma**cmode)`, no noise term) and `nParticles=int(DensMultiplications*Pixel_Ring*AvgChannelsPerPix)+1` (was `*Pixel_Ring` only). This is upstream's newest commit (`76ade48`, tip of Dev, "Changed how the number of particles in each ring is calculated"), touching the exact particle-count subsystem the earlier bit-exactness investigation spent a full session proving matches Fortran/JS to 2245/2245 particles across an entire optimizer trajectory. `Noise`/`AvgChannelsPerPix` are threaded through the signature (ported to JS too) but **deliberately not used in the calc yet** -- Nathan's own attached example config used `cdens=10`, but his email text states the actual new default is `cdens=100`; 10 looks like a typo in the example rather than the intended value (also consistent with `GalaxyFitParameters.py`'s dict-based default of 100, see below). Revisit once confirmed, and re-verify Fortran/JS agreement at the real default before switching the new term on. Separately, a real, unrelated Fortran/JS parity bug was found and fixed while investigating this: `R%Sigma` differs from the JS port by a few float32 ULPs (~5e-7 relative, ordinary unavoidable cross-platform FP noise), which was landing on different sides of the `int(X)+1` truncation boundary often enough to desync `idum` permanently partway through a fit. Fixed with a new `RoundForParticleStability` function (masks the low 12 mantissa bits before truncating, ~900x safety margin over the observed noise), applied on both platforms, confirmed to hold bit-exact `idum` lockstep through an entire fit to convergence. A second, independent bug was found and fixed alongside it: the historical `cdens=400`(Fortran)/`500`(JS) asymmetry (an old undocumented fudge) was silently producing different particle *counts* per ring from the very first optimizer evaluation -- removed in favor of matched `cdens=400` on both platforms. **UPDATE 2026-09-16**: `RoundForParticleStability` REMOVED -- a later, unrelated fix (bootstrap resampling geometry's PA "kinematic convention" round-trip, see `JS_FORTRAN_PARITY_SESSION_2026-09-15.md` UPDATE 13) apparently closed the actual `R%Sigma`-noise source; re-tested the masking's removal across 2 seeds x 5 bootstrap realizations each with it fully disabled and found zero behavior change vs masked, on both the deterministic anchor fit and bootstrap realizations. `Ring_CalcNumParticles` is back to plain `int(DensMultiplications*Pixel_Ring)+1`, matching upstream's un-masked style exactly. |
| `CompareCubes/FullModelComparison.f` | **applied+verified (plumbing only)** | Same `BuildTiltedRingModel` call-site update (computes `NoiseSpec=ObservedDC%DH%Uncertainty*abs(ObservedDC%DH%ChannelSize)`), ported to JS and verified. Rest of the 96-line diff is 100% local `TraceSwitch`/`TraceCallCounter`/`PrintStageChecksum` instrumentation (STAGE POSTFILL/POSTCONV, BADMODEL reason prints) -- all local-only, keep. |
| `Outputs/FitOutput.f` | **applied+verified** | 3-way, all resolved: (1) local-only `character(8)->(16)`/`F8.2->F16.6` precision fix, kept; (2) same `BuildTiltedRingModel` call-site update as above, applied; (3) the stale local comment referencing the removed `bootstrap-fit-launcher.js` fixed to say `bootstrap-realization-launcher.js`. |
| `GalaxyAnalysis/GalaxyFit.f` | **needs review** | Real, isolated tuning change: `IniGuessWidth=0.5` (second-pass simplex guess width) -> upstream `0.25`. Everything else in the 78-line diff is local-only (`RunObjectiveProbe`, `CONVERGED_VECTOR`/`TraceSwitch` dumps) -- keep. |
| `BootstrapSampler/FlippingBootstrap.f` | **declined** | Upstream swaps the seedable `call ran2(idum)` for `call RANDOM_NUMBER(RandVal)` in the block-resampling flip decision. Decision (Dan, 2026-08-17): keep `ran2` -- needed to keep checking numerical accuracy against the JS port, which bit-matches `ran2` specifically. Not adopted; revisit only if a JS-side equivalent of `RANDOM_NUMBER` is ever built. |
| `PreAnalysis/EstimateRadialProfiles.f` | needs review (trivial) | New safety cap: `if(nRings .gt. nRingsMax) nRings=nRingsMax` -- `nRingsMax` already exists locally (declared+set at the same point), trivial/safe to adopt. Also drops a large commented-out dead-code block (harmless, optional cleanup). Rest of the 65-line diff is local fdlibm (`fd_cos`/`fd_sin`, 2 sites) and the local `FoundFillIn` bug fix -- keep both. |
| `Inputs/InputUnitConversions.f` | needs review (trivial) | New accepted `FUnit` string: `'Jy Beam-1'` (alongside existing `'Jy/beam'`/`'Jy/Beam'`) -- purely additive, zero risk, trivial to adopt. Rest is the local BPA deg->rad bug fix -- keep. |
| `ConvolveCube/CalculateBeamKernel.f` | no action needed | 100% local fdlibm forcing (`fd_cos`/`fd_sin`/`fd_exp`) -- keep. |
| `PreAnalysis/EstimateShape.f` | no action needed | 100% local fdlibm forcing (`fd_cos`/`fd_sin`, 2 sites) -- keep. |
| `StandardMath/FullCircTrig.f` | no action needed | 100% local fdlibm forcing (`fd_atan`, well-documented bug fix) -- keep. |
| `StandardMath/random.f` | no action needed | Local forces `fd_log` in `gasdev`; upstream uses native `log`. Local-only fdlibm forcing, consistent with the rest of the codebase -- keep. |
| `Inputs/BootstrapRuntimeInputs.f` | no action needed | **Corrects the plan's earlier assumption.** Local already fixed the PA double-conversion bug -- by removing the *Fortran-side* re-conversion (opposite mechanism from upstream's fix, which removes the *Python-side* pre-conversion in `MakeBootstrapSample.py`). Both converge on doing the conversion exactly once. `js/ARCHITECTURE.md` §7's "reported to Nathan, not fixed" note is stale and should be updated. |
| `GeneralMinimizationRoutines/DownhillSimplex.f` | no action needed | Local shrink-step `y(i)` fix is comment-only vs. upstream -- upstream's code (commit `0839b49`) already matches ours exactly. |
| `MomentMaps/CalculateMomentMaps.f` | no action needed | Local `nChannels` off-by-one fix is comment-only vs. upstream -- upstream's code (commit `8e83f19`) already matches ours exactly. |
| `ProgramMains/SingleGalaxyFitTests.f` | no action needed | Entire 422-line diff is the local `DumpFittingFixture`/`DumpBestFitModelCube` JSON bridge -- no upstream equivalent. |
| `Globals/PipelineGlobals.f` | no action needed | 100% local diagnostic-switch declarations (`FixtureOnlySwitch`/`ProbeSwitch`/`TraceSwitch`/`DumpFixtureSwitch`) -- keep. |
| `Inputs/SingleFitRuntimeInputs.f` | no action needed | 100% local: backward-compatible optional reads for the 4 switches above -- keep. |

**Build-config only, confirmed no algorithmic content** (local customization for the JS/DCP port + fdlibm build):
- `src/StandardMath/Makefile`, `src/makeflags`, `src/ObjectLists` (adds `fdlibm_*.o` to the link list + a comment pointing at `js/src/StandardMath/fdlibm.js`).

**Local-only, no upstream counterpart at all**:
- `src/StandardMath/fdlibm_*.c`, `fdlibm_private.h` (11 files), `src/StandardMath/wasm/`.

## Python (`FitDriverScripts/`) -- triage complete

| File | Status | Note |
|---|---|---|
| `MakeBootstrapSample.py` | no action needed | Confirms the `BootstrapRuntimeInputs.f` finding above -- local already writes PA/Inc unconverted (functionally matches upstream's fix). Separately, local adds a `BootstrapSeed`-derived per-realization `idum` (documented as mirrored in `bootstrap-resample-launcher.js`) -- local-only, keep, this is what `run_both.js --seed` depends on. |
| `RunWRKP.py` + `GalaxyFitParameters.py` | needs review (low priority) | Upstream refactored default fitting-options loading from a static text file (`Inputs/SingleGalaxyTestFittingOptions_Base.txt`) to a generated-from-dict approach (`GalaxyFitParameters.DefaultRuntimeOptions()`/`RunWRKP.GenerateDefaultFittingOptionsFile()`). **Confirmed pure mechanism change** -- compared every default value against local's current `SingleGalaxyTestFittingOptions_Base.txt`, all match except `cdens` (local=400, upstream default=100), which is local's own intentional test-galaxy tuning, not a missed sync. Local's `WriteWRKPMainFile` also still appends the 4 local-only diagnostic-switch lines (`FixtureOnlySwitch` etc.) on top of whichever mechanism is used -- must be preserved either way. Low priority: no numerical-parity impact, purely a maintainability question of whether to adopt upstream's new mechanism. |
| `SetFileLocations.py` | no action needed | 100% local: DCP path definitions (`js/...`, `js/app/` dissolved 2026-08-17 -- see below) + `BootstrapSeed` default + `UseDCP` key registration -- keep. |
| `FullSingleGalaxyFit.py` | no action needed | 100% local: the `UseDCP` branch point (JS/DCP dispatch vs. fortran-local), a `chdir`-to-config-dir convenience feature (well-documented, real usability fix), `pool.close()/pool.join()` cleanup (minor resource-leak fix upstream lacks), and the DCP timing-checkpoint breakdown -- keep all. |
| `Bootstrap_Error_Analysis.py` | no action needed | 100% local: per-realization timing instrumentation for `run_both.js` -- keep. |
| `Bootstrap_Outputs.py` | no action needed | 100% local: defensive `os.makedirs` calls (real observed issue) + `StoreBootstrapTimings_JSON` (new function, needed for DCP/JS timing comparison) -- keep. |
| `SoFiA_Driver.py` | no action needed | 100% local: per-realization SoFiA log redirection + `WRKP_TRACE_DEBUG`-gated trace prints -- keep. |
| `BootstrapBoxPlot.py`, `BootstrapModelPlot.py` | no action needed | 100% local: defensive `os.makedirs` calls -- keep. |
| `CubeAnalysis.py`, `GeometryCorrection.py` | no action needed | 100% local: `FITSFixedWarning` suppression (cosmetic) -- keep. |
| `ReadWRKPFit.py` | no action needed | Upstream added 2 trivial print statements -- harmless, could adopt for parity of console output but zero functional stakes. Not worth a dedicated commit. |

**Local-only additions, no upstream counterpart**: `RunBootstrapsDCP.py`, `RunInitialFitDCP.py` (DCP dispatch bridge).

## `Inputs/`

| File | Status | Note |
|---|---|---|
| `SingleGalaxyTestFittingOptions_Base.txt` | no action needed | Only value that would differ vs. upstream's new dict-based defaults is `cdens` (400 vs 100) -- confirmed local test-specific tuning, not a missed sync (see `GalaxyFitParameters.py` above). |
| `SingleGalaxyTestFittingOptions.in` | local-only | No upstream counterpart. |
| upstream's `TiltedRingModel_ori.in`, `WRKP_GalaxyFitDriver_Parameters_10AS.py` | no action needed | Confirmed: Nathan's own personal test config referencing his own `/Users/nate/Dropbox/...` paths -- not portable, not applicable here. |

## Summary

- **2 real decision points surfaced; both resolved:**
  1. **Particle-count-per-ring formula change** (`SingleRingGeneration.f`, `TiltedRingModelGeneration.f`, `FullModelComparison.f`, `FitOutput.f`) -- upstream's newest commit, real physics change (adds a noise- and channel-spread-aware particle density term). Plumbing (new `BuildTiltedRingModel` args, `CalcAvgChanPerPix`) applied+verified on both Fortran and JS; the formula itself deferred pending `cdens=100` confirmation (see `SingleRingGeneration.f` row above). Investigating this surfaced and fixed two real, independent parity bugs (float32-truncation `idum` desync; asymmetric `cdens=400/500`) -- both fixed, both validated to hold bit-exact Fortran/JS `idum` lockstep through a full fit to convergence on the pre-76ade48 formula.
  2. **RNG swap in `FlippingBootstrap.f`** (`ran2`->`RANDOM_NUMBER`) -- declined (Dan, 2026-08-17): keep `ran2` for bit-exact JS parity; revisit only if a JS-side `RANDOM_NUMBER` equivalent is built.
- **3 trivial, safe, real changes ready to adopt**: `EstimateRadialProfiles.f`'s `nRingsMax` cap, `InputUnitConversions.f`'s `'Jy Beam-1'` string, `GalaxyFit.f`'s `IniGuessWidth` 0.5->0.25. Not yet applied.
- **1 low-priority refactor** (fitting-options loading mechanism) confirmed to have zero numerical impact -- can be deferred indefinitely.
- **1 stale doc note to fix**: `js/ARCHITECTURE.md` §7's "reported to Nathan, not fixed" -- local already fixed this (see `BootstrapRuntimeInputs.f`/`MakeBootstrapSample.py` above).
- ~~1 stale code comment to fix: `FitOutput.f`'s reference to the no-longer-existing `bootstrap-fit-launcher.js`~~ -- fixed.
- Everything else (bulk of both trees) is confirmed local-only or functionally-already-equivalent to upstream -- no action needed.

## Progress log

- 2026-09-15: Switched on upstream 76ade48's particle-count formula on both
  Fortran and JS (confirmed `cdens=100` is upstream's real default, resolving
  the open question below), applied the 3 other "ready to adopt" fixes, and
  fixed a real, separate bug: `cdens=` in a galaxy's `RTParameters.py` was a
  silent no-op on the fork (now wired up, mirroring `BootstrapSeed`). Isolated
  formula bit-exactness re-confirmed. Hit a NEW blocking bug re-verifying
  end-to-end on `WALLABY_J100336-262923`: the JS port crashes in pre-analysis
  (`NaN` flux-weighted center, zero-flux window) on a galaxy the earlier
  76ade48 investigation never tested against -- root-caused same day (see
  below) and fixed.
- 2026-09-15 (continued): root-caused the blocking bug -- JS's FITS reader
  never zeroed NaN/blanked pixels (unlike Fortran's `nullval` substitution)
  and never populated `nValid`/`flattendValidIndices`, silently breaking the
  "handled upstream" assumption documented in `FullModelComparison.js`. Fixed
  in `DataCubeFits.js`; verified bit-exact against Fortran's masked-cube/
  moment-0 stats. Found + fixed two more real bugs via full optimizer-
  trajectory trace diffing (`idum`/chi2/PA per call, numeric tolerance): a
  forgotten JS-side mirror of the `IniGuessWidth` change, and
  `RoundForParticleStability` flooring instead of rounding to nearest.
  Pushed verified bit-exact `idum` lockstep from an immediate crash to 99
  calls deep. Concluded further chasing is a structural dead end, not a bug
  list -- see the session doc for the exact hex-level evidence (a genuine
  round-to-nearest tie straddle) and a recommended alternative success
  criterion (statistical bootstrap agreement, not exact RNG-state lockstep).
- 2026-09-15 (RESOLVED): the "structural dead end" above was wrong. Per
  Dan's instruction to bisect the actual divergence instead of guessing at
  candidate functions, widened the FULLVEC hex dump to every call instead
  of just call 99 -- found the true origin was call 1 (the unperturbed
  starting vector), not something accumulated over 98 iterations. Two real
  bugs, both the "combine multiply+divide in one f32(), round once instead
  of twice" class: `getGalaxyShape`'s degrees->radians conversion
  (`InitialAnalysis.js`) and a SECOND gfortran compile-time constant-
  folding artifact in `JyAS_To_MsolPC` (`BasicConstants.js` -- a THIRD,
  already-correct copy of this same constant existed in `CommonConsts.js`
  and was reused instead of hardcoding a literal). After both fixes: all
  166 calls of Fortran's own initial-fit trajectory are bit-identical to
  JS's, every one of 13 parameters, hex value for hex value -- confirmed
  programmatically, not spot-checked. 11 real bugs found and fixed total
  across this investigation. Bootstrap realizations (vs. the one initial
  fit just proven) not yet re-verified -- natural next step. Full details,
  exact repro, and all hex-level evidence:
  `JS_FORTRAN_PARITY_SESSION_2026-09-15.md`.

- 2026-08-17: Full triage pass complete across all 18 differing Fortran files, all 13 differing Python files, and `Inputs/`. Every file read in full. Original plan's assumption that the PA double-conversion bug was still unfixed locally was wrong (corrected above) -- the actual first priority is the particle-count formula change and the two trivial/safe fixes, pending user direction on the RNG-swap question.
- 2026-08-17: `FlippingBootstrap.f` RNG swap declined. Particle-count formula plumbing (`BuildTiltedRingModel` signature, `CalcAvgChanPerPix`) applied and ported to JS. Attempted to switch on the new formula itself (`cmode=1, cdens=10`, matching Nathan's example config) -- caused catastrophic JS-side non-convergence. Root-caused via bit-level tracing to a pre-existing, unavoidable ~5e-7-relative float32 difference in `R%Sigma` between Fortran and the JS port landing on different sides of the `int(X)+1` truncation in `Ring_CalcNumParticles`, permanently desyncing `idum`. Fixed with `RoundForParticleStability`/`roundForParticleStability` (mantissa bit-masking, ~900x safety margin), added on both platforms. Separately found and fixed a second, independent bug: the historical asymmetric `cdens=400`(Fortran)/`500`(JS) calibration was silently causing different particle counts per ring from the first optimizer evaluation; matched to `cdens=400` on both platforms. With both fixes in place on the *pre-76ade48* formula, Fortran and JS ran in bit-exact `idum` lockstep through an entire fit to convergence (chi2 agreeing to ~8e-8 relative). Decision (Dan, 2026-08-17): preserve and ship this validated old-formula+fixes state now; defer switching on the new 76ade48 formula until `cdens` is confirmed as `100` (per Nathan's email text, not the `10` in his example config -- also consistent with `GalaxyFitParameters.py`'s dict-based default of 100) and the same tight-agreement bar is re-verified.
- 2026-08-17: Dissolved the `js/app/` layer -- everything that lived under it (`galaxy-fit.html`, `run-galaxy-fit-cli.js`, `bootstrap-realization-launcher.js`, `dcp-client-browser/`, `sample-data/`, `DCPjobData/`, etc.) now sits directly under `js/`, as a true sibling of `src/`/`package/`/`tools/`. Motivated by the people-page deploy URL: `galaxy-fit.html`'s `../src/...` script tags forced a nested `.../3KIDNAS/app/galaxy-fit.html` deploy path instead of the clean `.../3KIDNAS/` root the transcode demo uses. Fixed by making `src/` a direct child everywhere the page is served from (both local dev and remote deploy) -- `../src/` -> `src/` in `galaxy-fit.html`'s 2 script tags and the 2 JS files' `require()` calls, `SetFileLocations.py`'s `DCPDir` updated to `js/` (was `js/app`), plus doc/comment updates in `js/ARCHITECTURE.md`, `RunInitialFitDCP.py`, `RunBootstrapsDCP.py`, `FitOutput.f`. Deploy is now a single `rsync js/ people:~/public_html/3KIDNAS/` (excluding `package/`/`tools/`) instead of two separate commands for `app/`+`src/`.
- 2026-09-16: **New upstream bug class found, worth a report to Nathan**: `X**2.` (a REAL, not INTEGER, exponent literal -- e.g. `X**2.` vs `X**2`) is NOT guaranteed bit-identical to `X*X` in gfortran, confirmed with a minimal standalone repro (`transfer()`'d bit pattern, `-O0 -ffp-contract=off`, this project's own native build flags): for one specific value, `X**2.` gives `0x432B0014` while `X*X` gives `0x432B0015`. `X**2` (integer exponent, no decimal point) does NOT have this problem -- confirmed bit-identical to `X*X` for the SAME value; gfortran algebraically simplifies an integer exponent but apparently routes a real-literal exponent through a generic `pow()`-style runtime call instead. Swept the whole codebase (`grep -rn "\*\*2\."`): 26 call sites across 14 files (`StandardMath/random.f`'s `gasdev`, both `CompareCubes/LikelihoodFunctions.f` chi^2 functions, `TiltedRingModelGeneration/SingleRingGeneration.f`'s per-particle radius sampling, `ConvolveCube/CalculateBeamKernel.f`, `PreAnalysis/EstimateShape.f` x11, `PreAnalysis/EstimateRadialProfiles.f` x2, `PreAnalysis/EstimateCubeNoise.f`, `MomentMaps/CalculateMomentMaps.f` x2, `ConvolveCube/VelocitySmoothing.f` x2, `ObjectDefinitions/Beam.f`, `Inputs/CubeGeneratorInputs.f`, `Inputs/InputUnitConversions.f`, `PostAnalysis/PostGalaxyAnalysis.f`, `BootstrapSampler/PhysCoordTransform.f`) -- every single one confirmed byte-identical to `NateDeg/3KIDNAS` Dev branch (none fork-introduced), all still present upstream as of this date. Fixed all 26 locally (`X**2.` -> `X*X`, mathematically identical, removes the `pow()`-path) for Fortran/JS parity; the JS port never had this problem (`X*X` has no real-exponent ambiguity to begin with -- every corresponding JS site already used direct multiplication). Doesn't affect Fortran-only reproducibility (same binary, same answer every time), but IS a real portability/numerical-hygiene issue -- different gfortran versions or compilers could plausibly land on a different, still-valid, answer at any of these 26 sites. Full writeup, reproduction steps, and before/after measurements: `JS_FORTRAN_PARITY_SESSION_2026-09-15.md` UPDATE 15/16.
