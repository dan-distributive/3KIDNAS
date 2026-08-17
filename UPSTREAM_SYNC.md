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
| `TiltedRingModelGeneration/TiltedRingModelGeneration.f` | **needs review (priority)** | `BuildTiltedRingModel` gained `Noise,DC,BUse` args + new `CalcAvgChanPerPix` subroutine -- part of the particle-count-formula change, see below. |
| `TiltedRingModelGeneration/SingleRingGeneration.f` | **needs review (priority)** | `Ring_CalcNumParticles` formula changed: `DensMultiplications=CloudSurfDens*((Sigma/Noise)**cmode)` (was `*(Sigma**cmode)`, no noise term) and `nParticles=int(DensMultiplications*Pixel_Ring*AvgChannelsPerPix)+1` (was `*Pixel_Ring` only, no `AvgChannelsPerPix`). **This is upstream's newest commit (`76ade48`, tip of Dev, "Changed how the number of particles in each ring is calculated"), touching the exact particle-count subsystem the earlier bit-exactness investigation spent a full session proving matches Fortran/JS to 2245/2245 particles across an entire optimizer trajectory.** Local `TraceSwitch`/`NPTRACE`/`PARTTRACE` prints around this code are local-only, keep. |
| `CompareCubes/FullModelComparison.f` | **needs review (priority)** | Same `BuildTiltedRingModel` call-site update (computes `NoiseSpec=ObservedDC%DH%Uncertainty*abs(ObservedDC%DH%ChannelSize)`). Rest of the 96-line diff is 100% local `TraceSwitch`/`TraceCallCounter`/`PrintStageChecksum` instrumentation (STAGE POSTFILL/POSTCONV, BADMODEL reason prints) -- all local-only, keep. |
| `Outputs/FitOutput.f` | **needs review** | 3-way, already fully resolved: (1) local-only `character(8)->(16)`/`F8.2->F16.6` precision fix, keep; (2) same `BuildTiltedRingModel` call-site update as above, needs review; (3) a stale local comment referencing the removed `bootstrap-fit-launcher.js` -- should say `bootstrap-realization-launcher.js`, trivial fix while touching this file. |
| `GalaxyAnalysis/GalaxyFit.f` | **needs review** | Real, isolated tuning change: `IniGuessWidth=0.5` (second-pass simplex guess width) -> upstream `0.25`. Everything else in the 78-line diff is local-only (`RunObjectiveProbe`, `CONVERGED_VECTOR`/`TraceSwitch` dumps) -- keep. |
| `BootstrapSampler/FlippingBootstrap.f` | **flagged -- needs a decision** | Upstream swaps the seedable `call ran2(idum)` for `call RANDOM_NUMBER(RandVal)` in the block-resampling flip decision, and drops the `use BootstrapGlobals`/`use BasicRanNumGen` imports. `ran2` is still used elsewhere on both sides (`random.f`, `SingleRingGeneration.f`, `GalaxyFit.f` -- confirmed via grep), so this is localized to just this one call site, not a systemic RNG replacement. **Adopting it as-is would break bit-exact JS/Fortran parity for bootstrap resampling** -- `RANDOM_NUMBER`'s underlying algorithm isn't a portable spec the way `ran2` is (unlike everything else in this codebase that's been carefully bit-matched). Do not apply without an explicit decision on how (or whether) to preserve reproducibility. |
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
| `SetFileLocations.py` | no action needed | 100% local: DCP path definitions (`js/app/...`) + `BootstrapSeed` default + `UseDCP` key registration -- keep. |
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

- **2 real decision points surfaced, neither yet applied:**
  1. **Particle-count-per-ring formula change** (`SingleRingGeneration.f`, `TiltedRingModelGeneration.f`, `FullModelComparison.f`, `FitOutput.f`) -- upstream's newest commit, real physics change (adds a noise- and channel-spread-aware particle density term), touches the most heavily-verified subsystem in this whole project. High value, high care needed.
  2. **RNG swap in `FlippingBootstrap.f`** (`ran2`->`RANDOM_NUMBER`) -- real risk to bit-exact JS parity if adopted naively. Recommend NOT adopting without a specific plan for how to keep resampling reproducible.
- **3 trivial, safe, real changes ready to adopt**: `EstimateRadialProfiles.f`'s `nRingsMax` cap, `InputUnitConversions.f`'s `'Jy Beam-1'` string, `GalaxyFit.f`'s `IniGuessWidth` 0.5->0.25.
- **1 low-priority refactor** (fitting-options loading mechanism) confirmed to have zero numerical impact -- can be deferred indefinitely.
- **1 stale doc note to fix**: `js/ARCHITECTURE.md` §7's "reported to Nathan, not fixed" -- local already fixed this (see `BootstrapRuntimeInputs.f`/`MakeBootstrapSample.py` above).
- **1 stale code comment to fix**: `FitOutput.f`'s reference to the no-longer-existing `bootstrap-fit-launcher.js`.
- Everything else (bulk of both trees) is confirmed local-only or functionally-already-equivalent to upstream -- no action needed.

## Progress log

- 2026-08-17: Full triage pass complete across all 18 differing Fortran files, all 13 differing Python files, and `Inputs/`. Every file read in full. Original plan's assumption that the PA double-conversion bug was still unfixed locally was wrong (corrected above) -- the actual first priority is the particle-count formula change and the two trivial/safe fixes, pending user direction on the RNG-swap question.
