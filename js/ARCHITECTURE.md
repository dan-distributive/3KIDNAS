# Local (Fortran) vs. DCP bootstrap pipeline

This document explains how the two ways of running the bootstrap error
analysis — fully local/native (`UseDCP=0`) and DCP-distributed (`UseDCP=1`)
— actually work under the hood: which files run in what order, exactly
where the DCP path branches off from the native one, what gets written to
disk along the way, and how to run each side (and the comparison harness)
yourself.

Both paths are driven by the same entry point and share the same initial
fit; they only diverge for the bootstrap loop.

```
python3 WRKP_GalaxyFitDriver.py <config>.py
          |
          v
FitDriverScripts/FullSingleGalaxyFit.py :: GalaxyFit()
```

## 1. Shared initial fit (identical on both paths)

This part runs exactly the same way regardless of `UseDCP`.

1. `RunWRKP.RunWRKP(GeneralDict, GalaxyDict, BSSwitch=0)` shells out to the
   compiled **`Programs/SingleGalaxyFitter`** binary — the real optimizer,
   run once on the actual observed cube.
2. That binary always dumps two side artifacts to its working directory
   (bare filenames, hardcoded on the Fortran side):
   - `diskfit_fixture.json` — the fit's internal parameter/fixture state
     (`DumpFittingFixture`)
   - `model_cube_bestfit.json` — the dense best-fit model cube flux array
     (`DumpBestFitModelCube`)
3. `ReadWRKPFit.ReadWRKPOutputFile()` parses the fit results into
   `GalaxyDict['BestFitModel']`.

**Branch point** — `FullSingleGalaxyFit.py`, right after the initial fit:

```python
if GalaxyDict.get('UseDCP', 0) == 1:
    os.makedirs(GeneralDict['DCPJobDir'], exist_ok=True)
    for FName in ['diskfit_fixture.json', 'model_cube_bestfit.json']:
        os.replace(FName, os.path.join(GeneralDict['DCPJobDir'], FName))
...
if GalaxyDict.get('UseDCP', 0) == 1:
    from . import RunBootstrapsDCP as RBD
    BootstrapModels = RBD.RunBootstrapsDCP(GeneralDict, GalaxyDict)
else:
    # local path, below
```

`DCPJobDir` is `js/app/DCPjobData/` (see `SetFileLocations.py`).

## 2. Local/native path (`UseDCP=0`)

`FitDriverScripts/Bootstrap_Error_Analysis.py :: GetBootstrapModel()`, once
per realization, sequentially, on this machine:

| Step | Calls | What it does |
|---|---|---|
| 1 | `MakeBootstrapSample.MakeBootstrapSample()` | Writes a bootstrap input `.txt` (`WriteBootstrapFile`), shells out to **`Programs/BootStrapSampler`** (`FlippingBootstrap.f`'s `GenFlipBootstrapSample`), gets a resampled cube FITS file back |
| 2 | `SoFiA_Driver.RunSoFiA()` | Shells out to the native **`sofia`** binary on the resampled cube |
| 3 | `SoFiA_Driver.LoadSoFiAOutput()` | Parses the SoFiA catalogue into an Inc/PA estimate |
| 4 | `RunWRKP.RunWRKP(..., BSSwitch=1)` | Shells out to **`Programs/SingleGalaxyFitter`** again — the *full* optimizer, pinned to the original fit's ring count |
| 5 | `ReadWRKPFit.ReadWRKPOutputFile()` | Parses this realization's fit result |

Three native binary calls per realization, all sequential on one machine.

## 3. DCP-distributed path (`UseDCP=1`)

`FitDriverScripts/RunBootstrapsDCP.py :: RunBootstrapsDCP()` replaces the
loop above with two distributed "rounds" plus a slimmer local step:

```
RunBootstrapsDCP()
 |
 |-- BuildBaseFixture()        reads diskfit_fixture.json, strips the
 |                             observed cube's full flux (round 2 doesn't
 |                             need it -- each worker gets its own delta)
 |-- GetCubeWCS()              reads the real cube's WCS via astropy
 |
 |== ROUND 1 (DCP job, parallel across all realizations) ==============
 |-- BuildResamplePayload()    reads diskfit_fixture.json (UN-stripped) +
 |                             model_cube_bestfit.json, computes bsCent
 |-- LaunchDCPResampleJob()    spawns: node bootstrap-resample-launcher.js
 |                             --> dispatches resampleAndSoFiA() to workers
 |                                   FlipBootstrap.genFlipBootstrapSample   (JS port of FlippingBootstrap.f)
 |                                   DataCubeFits.dataCubeToFitsBytes      (write real FITS bytes)
 |                                   sofia2wasm                            (replaces native `sofia`)
 |                                   ParseSoFiACatalog + GeometryEstimates (replaces LoadSoFiAOutput)
 |                             <-- per realization: cube+mask FITS (base64),
 |                                 inc/PA, timings
 |                             writes: DCPjobData/resample_payload.json,
 |                                     DCPjobData/dcp_resample_results.json
 |
 |== LOCAL (sequential, but only ONE native call per realization) =====
 |-- PrepareBootstrapDeltaFromRoundOne()   per realization:
 |     - decode round-1's cube/mask FITS bytes to disk
 |     - SoFiA_Driver.AdjustMaskFile()          isolate the source in the mask
 |     - SoFiA_Driver.WriteSoFiACatFileForWRKP()  from round-1's own Inc/PA
 |     - RunFixtureOnlyFit()   shells out to Programs/SingleGalaxyFitter,
 |                             FixtureOnlySwitch=1 (stops BEFORE the
 |                             optimizer) -- supplies this realization's own
 |                             noise estimate + Jy/beam->Jy/pixel factor.
 |                             THE ONE PIECE THAT STAYS NATIVE: no JS/wasm
 |                             port exists for this.
 |                             dumps a per-realization fixture, moved to
 |                             BootstrapFolder/<Obj>_Bootstrap_N_fixture.json
 |     - make_bootstrap_delta.make_delta()   builds the small "delta" payload
 |                             round 2 needs (whatever varies per realization)
 |
 |== ROUND 2 (DCP job, parallel across all realizations, UNCHANGED) ====
 |-- LaunchDCPFitJob()         spawns: node bootstrap-fit-launcher.js
 |                             --> dispatches bootstrap(delta, fixtureStr)
 |                                   galaxyFit_Simple            (JS tilted-ring optimizer)
 |                                   cubeBeamConvolution          (FFTW3JS-based beam convolution,
 |                                                                 called once per objective-function eval)
 |                             <-- per realization: fitted params, chi2, timings
 |                             writes: DCPjobData/deltas.ndjson,
 |                                     DCPjobData/base_fixture.json,
 |                                     DCPjobData/dcp_results.json
 |
 |-- AdaptResults()            worker results -> BootstrapModels (same shape
                               GetBootstrapModel() produces on the local path)
```

One native call per realization instead of three; the two removed calls
(resample, SoFiA) run in parallel across the DCP network instead of
sequentially on one laptop.

## 4. What gets dumped, where

| Artifact | Written by | Location | Used for |
|---|---|---|---|
| `<Obj>_AverageModel_v1.fits` | Initial fit (both paths) | `TargFolder/<Obj>/` | Reference best-fit model |
| `diskfit_fixture.json` | Initial fit (`DumpFittingFixture`) | CWD → moved to `DCPjobData/` if `UseDCP=1` | Round 2's shared base fixture; round 1's observed-cube payload |
| `model_cube_bestfit.json` | Initial fit (`DumpBestFitModelCube`) | CWD → moved to `DCPjobData/` if `UseDCP=1` | Round 1's model-cube payload |
| `resample_payload.json` | `LaunchDCPResampleJob` | `DCPjobData/` | Round 1 job's static args |
| `dcp_resample_results.json` | `LaunchDCPResampleJob` | `DCPjobData/` | Round 1 job's per-realization results |
| `<Obj>_Bootstrap_N_fixture.json` | `RunFixtureOnlyFit` | `TargFolder/<Obj>/BootstrapCubes/` | Per-realization noise + flux-conversion factor |
| `deltas.ndjson`, `base_fixture.json` | `LaunchDCPFitJob` | `DCPjobData/` | Round 2 job's inputs |
| `dcp_results.json` | `LaunchDCPFitJob` | `DCPjobData/` | Round 2 job's per-realization results |
| `<Obj>_BootstrapFits.csv` | `Bootstrap_Outputs.StoreBootstrappedModels_CSV` | `TargFolder/<Obj>/` | Final per-realization results (both paths, same format) |
| `<Obj>_FitTimeCheck.csv` | `FullSingleGalaxyFit.py` | `TargFolder/<Obj>/` | Timing checkpoints (both paths) |

## 5. Serving `galaxy-fit.html` locally

`js/app/galaxy-fit.html` loads `../src/PipelineConfig/defaultFittingOptions.js`
and `../src/PayloadBuilder/buildFitPayloads.js` as plain `<script>` tags --
`src/` is a *sibling* of `app/`, not nested inside it, since the `js/`
reorg. A static file server rooted at `js/app/` (e.g. `python3 -m
http.server` run from inside `app/`) refuses to serve anything above its
own root, so those two script tags 404 silently and every global they
define (`DefaultFittingOptions`, `PayloadBuilder`) comes up `undefined`.

Root the server one level up, at `js/` itself, and load the page via
`/app/galaxy-fit.html`:

```bash
cd /Users/dandesjardins/DCP/3KIDNAS/js
python3 -m http.server 8000
# then open http://localhost:8000/app/galaxy-fit.html
```

## 6. Running each side yourself

Both need a config `.py` file (see `3KIDNASTests/SingleGalaxyTest/Sample_KIDNAS_SingleGalaxyInput*.py`
for examples) — the only difference that matters here is `UseDCP`.

### Local only (no DCP, no credentials needed)

```bash
cd /Users/dandesjardins/DCP/3KIDNAS/3KIDNASTests/SingleGalaxyTest
python3 ../../WRKP_GalaxyFitDriver.py Sample_KIDNAS_SingleGalaxyInput_Small.py   # UseDCP=0
```

### DCP only

**`DCP_API_KEY` is required — the job will not dispatch without it.**
`DCP_COMPUTE_GROUP` and `DCP_SLICE_PRICE` are optional but important: omit
them and you get the **public** compute group and DCP's own **market-rate**
pricing, which may not be what you want for a real run.

```bash
cd /Users/dandesjardins/DCP/3KIDNAS/3KIDNASTests/SingleGalaxyTest
export DCP_API_KEY=0x<your identity>
export DCP_COMPUTE_GROUP=<joinKey>[,<joinSecret>]   # optional -- omit for public compute
export DCP_SLICE_PRICE=<price>                       # optional -- omit for market rate
python3 ../../WRKP_GalaxyFitDriver.py Sample_KIDNAS_SingleGalaxyInput_DCP_Small.py   # UseDCP=1
```

These three env vars are read by `RunBootstrapsDCP.py`'s `LaunchDCPResampleJob`
and `LaunchDCPFitJob` for **both** rounds (resample+SoFiA, and the fit) —
one export covers both DCP jobs in the pipeline.

### Comparison harness (`run_both.js`)

Runs local then DCP **sequentially** (not in parallel — both hit native
Fortran on this machine for part of their work, so running them
concurrently would contend for CPU and corrupt the timing comparison),
then diffs their per-realization results.

```bash
cd /Users/dandesjardins/DCP/3KIDNAS/js/tools
export DCP_API_KEY=0x<your identity>
export DCP_COMPUTE_GROUP=<joinKey>[,<joinSecret>]   # optional
export DCP_SLICE_PRICE=<price>                       # optional
node run_both.js --seed 42 --nBootstraps 5
```

Flags:

| Flag | Purpose |
|---|---|
| `--seed <idum>` | Threads the same `BootstrapSeed` into both configs so round-1 resampling draws identical per-realization random numbers on both sides (bit-exact-to-float32, verified). Omit to run both sides with today's unseeded, time-based randomness — the comparison is then only a loose statistical sanity check. |
| `--nBootstraps <N>` | Number of realizations (default 5) |
| `--match-known-bug` | See below |
| `--json <path>` | Also write the full timing + comparison report to a JSON file |
| `--skip-wipe` | Don't clear previous `TestFits_RunBoth_*` output first |
| `--skip-local` / `--skip-dcp` | Run only one side |

## 7. Known open issue: PA double-conversion (local path only)

`MakeBootstrapSample.WriteBootstrapFile` writes the resampling PA already
converted to radians (its own comment in the file says "PA & INC in
radians"). `BootstrapRuntimeInputs.f`'s `BootstrapIn()` then converts it
**again** (`BS_Cent%PA=BS_Cent%PA*Pi/180.`), so the local/native path's
actual internal resampling PA has been landing at ~5.78° instead of the
intended ~331.3°. This is a real, pre-existing bug independent of anything
in this document — reported to Nathan for confirmation, not fixed, since
fixing it changes existing local-pipeline behavior.

The DCP path never had this bug (it sends PA straight over JSON, no
text-file round-trip), so **without extra handling, `run_both.js` compares
two sides that are resampling with genuinely different geometry** — a real
gap, not a porting artifact. Pass `--match-known-bug` to make the DCP side
apply the same extra conversion for comparison purposes only (never affects
normal, non-`run_both` DCP runs), so the rest of the pipeline's fidelity can
be validated while Nathan's answer is pending:

```bash
node run_both.js --seed 42 --nBootstraps 5 --match-known-bug --json results.json
```
