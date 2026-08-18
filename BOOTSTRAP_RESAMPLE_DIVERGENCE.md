# Bootstrap-resampling Fortran/JS divergence — investigation log

**Status as of 2026-08-17: unresolved.** Four distinct, well-evidenced hypotheses
tested and disproven. Root cause not yet identified. This doc exists so the
investigation can be picked up cold, by anyone (including a future session
with no memory of how we got here).

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

## Hypotheses tested and disproven (all four, with hard evidence)

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

### 5. `idum` divergence from differing Nelder-Mead search paths — TESTED AND DISPROVEN
**Theory:** even when two independent fits converge to nearly the same
final answer, the *sequence* of evaluations each platform's optimizer
takes to get there is not guaranteed identical (this project's own
well-documented Nelder-Mead instability). Since `idum` is a running total
of every random draw consumed since the fit began, and the model cube is
resynthesized *after* convergence using whatever `idum` happens to be at
that point, the two platforms could easily be drawing from completely
different, uncorrelated positions in their own (individually bit-exact)
`ran2`/`gasdev` streams.
**Direct test built and run:** added `TRACE_OVERRIDE_IDUM` on both
platforms — an env var read right before the final model-cube resynthesis
call (`OutputCube` in `src/Outputs/FitOutput.f` via a new
`MaybeOverrideIdum` subroutine; `runInitialFit` in
`js/bootstrap-realization-launcher.js`, replacing `state.rng` with a fresh
`makeRng(overrideIdum)` right before `tiltedRingModelComparison`). Ran both
platforms with **the same fixed idum (-12345) AND the same fixed 15-param
vector** (Fortran's own converged values).
**Result: DISPROVEN.** If this were the cause, forcing both idum and
params identical should have collapsed the divergence to near-zero. It did
not — max abs diff at the tracked pixels was actually slightly *larger*
(0.00115 vs 0.00049 in the params-only test), and the overall cube diff
did not meaningfully improve.
**Verified the override mechanism was actually wired correctly** before
trusting this result: `tiltedRingModelComparison` destructures `state.rng`
at its own top and passes it straight to `buildTiltedRingModel` — the
override was confirmed to reach particle generation, not silently
bypassed.

## Where this leaves us

With idum and parameters both forced identical, the model cubes should be
bit-identical (up to the accepted noise floor) if `Ring_ParticleGeneration`
is a faithful port and nothing else differs. They aren't. That means
**something else still differs between the two `BuildTiltedRingModel` call
environments that hasn't been identified yet.** Candidates not yet tested:

- **The model cube's own header fields** (`start`, `channelSize`,
  `refVal`/`refLocation`, etc.) — `FindParticleCellLocation`'s channel
  binning depends directly on `DC%DH%Start(2)`/`dh.start[2]` and
  `ChannelSize`. If these aren't actually identical between the two
  `state`/`GalaxyDict` setups (even for "the same" galaxy), every
  particle's channel assignment would shift systematically — which would
  produce exactly this kind of large, localized, non-noise-floor
  divergence. **Not yet checked directly.**
- **`generalizedParamVectorToTiltedRing`'s deserialization** — converts
  the raw 15-number vector into per-ring `Rmid`, `Rwidth`, `Inclination`,
  `PositionAngle`, `VSys`, etc. Even with the raw 15 numbers forced
  identical, this deserialization step itself does real arithmetic
  (unit conversions, trig) that hasn't been checked line-by-line the way
  `PhysCoordTransform` was. A bug or an unmatched-trig-function here would
  explain a real, continuous ring-geometry difference feeding directly
  into `Ring_ParticleGeneration`.
- Not yet ruled out: some other input to `Ring_ParticleGeneration` besides
  idum/Sigma/ring-geometry (e.g. `CloudBaseSurfDens`, `cmode` — though
  these are simple scalars, unlikely but not directly re-verified in this
  specific controlled test).

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
| `TRACE_OVERRIDE_IDUM` | forces `idum` used for that same resynthesis call, on **both** platforms | `src/Outputs/FitOutput.f` (`MaybeOverrideIdum`), `js/bootstrap-realization-launcher.js` (new, added during this investigation) |

**To retarget the hardcoded pixel/index values above to a different cell**,
just edit the literal `if` conditions at each trace site — they're
intentionally narrow (avoid flooding output across ~180k cells/tens of
thousands of particles) rather than parameterized via env var. If you need
a different pixel, grep for `PXTRACE`/`CORNERTRACE`/`SRCTRACE` and change
the hardcoded indices, matching a divergent cell from a fresh pixel diff.

## Fixes kept from this investigation (real, worth keeping, not the root cause of the above)

1. **`RoundForInterpStability`/`roundForInterpStability`** — protects
   `GetFluxAtPoint`'s corner-selection truncation against the general class
   of float32-boundary-flip bug (even though it wasn't the cause here, the
   underlying risk it protects against is real, same as the particle-count
   fix from earlier in this session).
2. **`PhysCoordTransform` fdlibm fix** — closes a real gap in the
   project's fdlibm-forcing coverage; `BootstrapSampler/` was simply missed
   by the earlier sweep.
3. **`RoundForBinStability`/`roundForBinStability`** — same protective
   technique applied to per-particle nearest-pixel binning.

All three are committed as legitimate improvements independent of whether
they explain this specific divergence.

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
