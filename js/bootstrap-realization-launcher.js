'use strict';
// =============================================================================
// bootstrap-realization-launcher.js
//
// ONE DCP job, ONE slice per bootstrap realization; each slice runs the
// entire per-realization pipeline (resample -> SoFiA -> InitialAnalysis ->
// fit) atomically, mirroring Bootstrap_Error_Analysis.GetBootstrapModel and
// pool.starmap(BootstrapRunStep, ...) on the native side. Replaces an older
// two-DCP-dispatch design (bootstrap-resample-launcher.js +
// bootstrap-fit-launcher.js) that handed a delta/NDJSON payload between them.
//
// `payload` (built once by RunBootstrapsDCP.py) carries everything identical
// across realizations: observedDC, modelFlux, beamMajorAxis, beamMinorAxis,
// bsCent, velBlockSize, bootstrapSeed, sofiaParTemplate, observedBeam
// {beamSigma0/1/2, nRadialCells, sigmaLengths, pixelSizeX/Y, beamMajorAxis,
// beamMinorAxis}, centerSource, nRingsPerBeam, nTargRings, radGridArcsec
// (pinned ring grid), sdSwitch, linSDLims, logSDLims, vRotLims, sizeLims,
// noiseSigmaLim, cmode, cloudBaseSurfDens, fitIdum (the STATIC per-galaxy
// fit seed -- see InitialAnalysis.js's header for why it doesn't vary per
// realization), ftol, likelihoodSwitch.
//
// `runBootstrapRealization` is module-scope and self-contained (no closures,
// all requires inside its own body) because DCP's compute.for ships only
// this function's source to remote sandboxes, which can't see this file's
// module scope. That's also what lets it run directly for local testing:
//   const { runBootstrapRealization } = require('./bootstrap-realization-launcher');
//   const result = await runBootstrapRealization(0, payload);
// =============================================================================

// Node-only globals, used by the CLI/local-pool tail at the bottom of this
// file. Guarded so this same file can also be loaded via a plain <script>
// tag in a browser (galaxy-fit.html): compute.for() never actually EXECUTES
// runBootstrapRealization/runInitialFit locally in either environment (it
// only ships their source to a remote DCP sandbox), so the only code here
// that needs to run in-process is this module-scope setup and the CLI tail
// -- both no-ops in a browser, where isNode is false.
const isNode = typeof process !== 'undefined' && !!process.versions && !!process.versions.node;
let fs, os, execSync, Worker, isMainThread, parentPort, workerData;
if (isNode) {
  fs = require('fs');
  os = require('os');
  ({ execSync } = require('child_process'));
  ({ Worker, isMainThread, parentPort, workerData } = require('worker_threads'));
}

// ---------------------------------------------------------------------------
// getPerformanceCoreCount
// On heterogeneous-core hardware (Apple Silicon: performance + efficiency
// cores), os.cpus().length counts both core types as equivalent, but
// efficiency cores are dramatically slower for sustained CPU-bound
// FFTW/fdlibm work like this pipeline's fit loop -- confirmed directly on an
// M2 (4 P-cores + 4 E-cores, hw.perflevel0/1.physicalcpu): sizing the local
// worker pool to logical-cores-1 (7 workers) forces at least 3 of them onto
// E-cores or into P-core time-sharing, and a realization run under that
// full 7-way contention took 6.6-7.7x longer per objective-function
// evaluation than the exact same realization run solo -- close to real DCP
// dispatch's own per-realization speed (each of ITS workers gets a
// separate physical machine, not a shared core, so there's no P/E split to
// worry about there). Returns null (caller falls back to the original
// logical-cores-1 sizing) on non-Darwin platforms or Intel Macs, where
// there's no perflevel0/1 split and cores are already homogeneous -- the
// original sizing is already correct there.
// ---------------------------------------------------------------------------
function getPerformanceCoreCount() {
  if (process.platform !== 'darwin') return null;
  try {
    const out = execSync('sysctl -n hw.perflevel0.physicalcpu', { stdio: ['ignore', 'pipe', 'ignore'] })
      .toString().trim();
    const n = parseInt(out, 10);
    return Number.isFinite(n) && n > 0 ? n : null;
  } catch (e) {
    return null;
  }
}

// -----------------------------------------------------------------------
// LOCAL WORKER-THREAD ENTRY POINT -- only reached when this file is loaded
// as a worker_threads Worker (the --local pool below); never part of a real
// DCP dispatch. Checked ahead of `require.main === module` below because a
// Worker's own entry point ALSO satisfies that check -- isMainThread is
// what actually tells "spawned local worker" apart from "CLI process".
//
// payload is sent once at Worker construction (workerData), not resent per
// task -- it can be several MB (raw FITS bytes, base64); only lightweight
// {realizationIndex} messages flow after that.
// -----------------------------------------------------------------------
if (isNode && !isMainThread) {
  const { payload, kind } = workerData;
  parentPort.on('message', async (msg) => {
    if (!msg || msg.type !== 'task') return;
    try {
      // kind is only set for the (single-job) initial-fit local pool below;
      // the bootstrap pool's Worker never sets it, so this preserves the
      // existing bootstrap dispatch's exact behavior.
      const result = kind === 'initialFit'
        ? await runInitialFit(msg.realizationIndex, payload)
        : await runBootstrapRealization(msg.realizationIndex, payload);
      parentPort.postMessage({ type: 'result', realizationIndex: msg.realizationIndex, ok: true, result });
    } catch (e) {
      parentPort.postMessage({ type: 'result', realizationIndex: msg.realizationIndex, ok: false, error: e.message });
    }
  });
}

// -----------------------------------------------------------------------
// WORK FUNCTION
// -----------------------------------------------------------------------
async function runBootstrapRealization(realizationIndex, payload) {
  const report = { realizationIndex };
  const tEntry = Date.now();
  const timings = {
    resampleMs: null, sofiaMs: null, fixtureFitMs: null,
    fitMs: null, convolveMs: null, convolveCalls: null, evalCount: null, totalMs: null,
  };

  try {
    // Liveness ping for DCP's sandbox supervisor, which kills a slice with
    // ENOPROGRESS if it never sees a progress() call in time -- same fix as
    // runInitialFit's own (see that function's identical comment for the
    // real-dispatch failure that motivated it originally). This function
    // never had one: real bootstrap dispatch always failed at module
    // resolution before today, so nothing ever ran long enough to expose
    // the gap until the module-resolution fix let it actually start
    // executing -- confirmed directly, same "Sandbox never emitted a
    // progress event" / ENOPROGRESS failure, real dispatch, before this.
    if (typeof progress === 'function') progress(0);

    // No `?? default` fallbacks on payload fields anywhere in this function:
    // a silent default (e.g. `payload.centerSource ?? 0`) is exactly what
    // let a real bug ship unnoticed once (JS quietly fit 65 free params
    // instead of Fortran's 16, with no error). Missing data should throw,
    // naming the field, not produce a plausible-looking wrong answer.
    //
    // Declared INSIDE this function, not module scope: compute.for ships
    // only this function's own source to remote sandboxes, which can't see
    // module scope -- a module-scope version threw "requirePayloadField is
    // not defined" on real workers while working fine locally. Same reason
    // isTraceDebug() guards `process` instead of touching it directly:
    // `process` isn't guaranteed to exist in a DCP worker sandbox.
    function isTraceDebug() {
      return typeof process !== 'undefined' && process.env && process.env.TRACE_DEBUG === '1';
    }

    function requirePayloadField(payload, name) {
      const v = payload[name];
      if (v === undefined || v === null) {
        throw new Error(`runBootstrapRealization: payload.${name} is required but missing `
          + `(no silent default -- see this function's header comment for why).`);
      }
      return v;
    }

    // ---- Every field this realization reads, pulled out once here (see
    //      the file header for what each means). observedDC -> aliased to
    //      observedDCHeader because a DataCube INSTANCE gets built under the
    //      name `observedDC` a few lines down -- destructuring the payload
    //      field under its own name would collide with it. ----
    const observedDCHeader = payload.observedDC;
    const {
      beamMajorAxis, beamMinorAxis, bsCent, velBlockSize, bootstrapSeed,
      sofiaParTemplate, observedBeam, nRingsPerBeam, nTargRings, radGridArcsec,
      linSDLims, logSDLims, vRotLims, sizeLims, noiseSigmaLim,
      constParams, fixedParams, fitIdum, ftol,
    } = payload;
    const observedCubeRawFitsB64 = requirePayloadField(payload, 'observedCubeRawFitsB64');
    const modelCubeRawFitsB64    = requirePayloadField(payload, 'modelCubeRawFitsB64');
    const centerSource           = requirePayloadField(payload, 'centerSource');
    const sdSwitch               = requirePayloadField(payload, 'sdSwitch');
    const cmode                  = requirePayloadField(payload, 'cmode');
    const cloudBaseSurfDens      = requirePayloadField(payload, 'cloudBaseSurfDens');
    const likelihoodSwitch       = requirePayloadField(payload, 'likelihoodSwitch');

    const f32 = Math.fround;
    let DataCube, allocateDataCube, Beam2D, allocate_Beam2D, genFlipBootstrapSample, makeRng,
        dataCubeToFitsBytes, fitsBytesToDataCube, parseSoFiACatalog, getGeometryEstimates,
        initialAnalysis, TiltedRingModel, tiltRing_Allocate, ParameterVector, allocateParamVector,
        calculate2DBeamKernel, resetConvolveStats, getConvolveStats, warmUp,
        generalizedParamVectorToTiltedRing, warmUpWasmTrig, galaxyFit_Simple, resetEvalStats,
        getEvalStats, JyAS_To_MsolPC;
    try {
      const entry = require('initialFitEntry.js');
      ({ DataCube, allocateDataCube } = entry.DataCube);
      ({ Beam2D, allocate_Beam2D } = entry.Beam);
      ({ genFlipBootstrapSample } = entry.FlipBootstrap);
      ({ makeRng } = entry.random);
      ({ dataCubeToFitsBytes, fitsBytesToDataCube } = entry.DataCubeFits);
      ({ parseSoFiACatalog } = entry.ParseSoFiACatalog);
      ({ getGeometryEstimates } = entry.GeometryEstimates);
      ({ initialAnalysis } = entry.InitialAnalysis);
      ({ TiltedRingModel, tiltRing_Allocate } = entry.TiltedRing);
      ({ ParameterVector, allocateParamVector } = entry.ParameterVector);
      ({ calculate2DBeamKernel } = entry.CalculateBeamKernel);
      ({ resetConvolveStats, getConvolveStats, warmUp } = entry.CubeKernelConvolution);
      ({ generalizedParamVectorToTiltedRing } = entry.ParameterToTiltedRingVector);
      ({ warmUpWasmTrig } = entry.TiltedRingModelGeneration);
      ({ galaxyFit_Simple } = entry.GalaxyFit);
      ({ resetEvalStats, getEvalStats } = entry.FullModelComparison);
      ({ JyAS_To_MsolPC } = entry.BasicConstants);
    } catch (publishedPackageError) {
      // See runInitialFit's identical catch block for why both errors get
      // surfaced on a double failure instead of just the fallback's own.
      try {
        ({ DataCube, allocateDataCube } = require('./src/ObjectDefinitions/DataCube'));
        ({ Beam2D, allocate_Beam2D } = require('./src/ObjectDefinitions/Beam'));
        ({ genFlipBootstrapSample } = require('./src/BootstrapSampler/FlipBootstrap'));
        ({ makeRng } = require('./src/StandardMath/random'));
        ({ dataCubeToFitsBytes, fitsBytesToDataCube } = require('./src/BootstrapSampler/DataCubeFits'));
        ({ parseSoFiACatalog } = require('./src/SoFiACatalog/ParseSoFiACatalog'));
        ({ getGeometryEstimates } = require('./src/GeometryEstimates/GeometryEstimates'));
        ({ initialAnalysis } = require('./src/PreAnalysis/InitialAnalysis'));
        ({ TiltedRingModel, tiltRing_Allocate } = require('./src/ObjectDefinitions/TiltedRing'));
        ({ ParameterVector, allocateParamVector } = require('./src/ObjectDefinitions/ParameterVector'));
        ({ calculate2DBeamKernel } = require('./src/ConvolveCube/CalculateBeamKernel'));
        ({ resetConvolveStats, getConvolveStats, warmUp } = require('./src/ConvolveCube/CubeKernelConvolution'));
        ({ generalizedParamVectorToTiltedRing } = require('./src/ParameterToTiltedRingInterface/ParameterToTiltedRingVector'));
        ({ warmUpWasmTrig } = require('./src/TiltedRingModelGeneration/TiltedRingModelGeneration'));
        ({ galaxyFit_Simple } = require('./src/GalaxyAnalysis/GalaxyFit'));
        ({ resetEvalStats, getEvalStats } = require('./src/CompareCubes/FullModelComparison'));
        ({ JyAS_To_MsolPC } = require('./src/StandardMath/BasicConstants'));
      } catch (localFallbackError) {
        throw new Error(`Could not load pipeline modules via the published package OR the local fallback. `
          + `Published-package error: ${publishedPackageError.message}. `
          + `Local-fallback error: ${localFallbackError.message}`);
      }
    }
    // Published-package path (real DCP worker) tried first, falls back to
    // the local relative path for standalone `node` runs -- same try/catch
    // pattern FFTW3WasmRank2.js uses for fftw-wasm.js.
    let cfitsio, sofia;
    try {
      cfitsio = require('cfitsio-wasm.js');
    } catch (e) {
      cfitsio = require('../third_party/cfitsio-4.6.3/wasm/cfitsio-wasm.js');
    }
    try {
      sofia = require('sofia-wasm.js');
    } catch (e) {
      sofia = require('../third_party/SoFiA-2-master_2_5_1/wasm/sofia-wasm.js');
    }

    // One-time async wasm instantiation, kicked off now to overlap with
    // resample/SoFiA/InitialAnalysis below; awaited just before the
    // optimizer loop starts. fdlibm falling back to its JS port is fine
    // (slower, not wrong -- see SingleRingGeneration.js's warmUpWasmTrig).
    const fftwWarmUpPromise = warmUp();
    const trigWarmUpPromise = warmUpWasmTrig();

    // ---- Build the observed cube header (same for every realization) ----
    // Flux is loaded separately below, from RAW (Jy/beam) FITS bytes, NOT
    // from d.validFlux -- see the raw-bytes loading block for why.
    const observedDC = new DataCube();
    const odh = observedDC.dh;
    const d = observedDCHeader;
    odh.nPixels[0] = d.nPixelsX; odh.nPixels[1] = d.nPixelsY; odh.nChannels = d.nChannels;
    odh.pixelSize[0] = d.pixelSizeX; odh.pixelSize[1] = d.pixelSizeY; odh.channelSize = d.channelSize;
    odh.refLocation[0] = d.refLocX; odh.refLocation[1] = d.refLocY; odh.refLocation[2] = d.refLocV;
    odh.refVal[0] = d.refValX; odh.refVal[1] = d.refValY; odh.refVal[2] = d.refValV;
    odh.uncertainty = d.uncertainty;
    allocateDataCube(observedDC);

    // ---- Build the model cube header (same as observed) ----
    const modelDC = new DataCube();
    const mdh = modelDC.dh;
    mdh.nPixels[0] = d.nPixelsX; mdh.nPixels[1] = d.nPixelsY; mdh.nChannels = d.nChannels;
    mdh.pixelSize[0] = d.pixelSizeX; mdh.pixelSize[1] = d.pixelSizeY; mdh.channelSize = d.channelSize;
    mdh.refLocation[0] = d.refLocX; mdh.refLocation[1] = d.refLocY; mdh.refLocation[2] = d.refLocV;
    mdh.refVal[0] = d.refValX; mdh.refVal[1] = d.refValY; mdh.refVal[2] = d.refValV;
    mdh.uncertainty = d.uncertainty;
    allocateDataCube(modelDC);

    // ---- Load RAW (Jy/beam, un-brightness-converted) flux for both cubes.
    //      Native resamples and runs SoFiA on the raw cube; brightness
    //      conversion (/beamAreaPixels) only happens AFTER SoFiA, on the
    //      fit's own ObservedDC (SingleGalaxyFitTests.f:45-51). Converting
    //      earlier (as this file used to) fed SoFiA data ~beamAreaPixels-
    //      times too small and caused a real ~0.1 deg catalogue PA
    //      divergence from native -- conversion happens after SoFiA below
    //      instead. ----
    // atob + charCodeAt, not Buffer.from -- Buffer isn't guaranteed to exist
    // in a DCP worker sandbox (matches test-cfitsio-worker.js).
    const b64ToBytes = (b64) => Uint8Array.from(atob(b64), (c) => c.charCodeAt(0));
    const rawObservedBytes = b64ToBytes(observedCubeRawFitsB64);
    const rawModelBytes    = b64ToBytes(modelCubeRawFitsB64);
    const rawObservedDC = await fitsBytesToDataCube(cfitsio, rawObservedBytes, observedDC);
    const rawModelDC    = await fitsBytesToDataCube(cfitsio, rawModelBytes, modelDC);
    observedDC.flux = rawObservedDC.flux;
    modelDC.flux    = rawModelDC.flux;

    const resampleBeam = new Beam2D();
    resampleBeam.beamMajorAxis = beamMajorAxis;
    resampleBeam.beamMinorAxis = beamMinorAxis;
    allocate_Beam2D(resampleBeam, [d.nPixelsX, d.nPixelsY]);

    // ---- Resample (replaces the native BootStrapSampler binary) ----
    // Per-realization idum, same formula as MakeBootstrapSample.py's
    // WriteBootstrapFile -- independent of the fit's own STATIC idum below
    // (see InitialAnalysis.js's header for why that one doesn't vary).
    const resampleIdum = bootstrapSeed
      ? -(Math.abs(bootstrapSeed) + realizationIndex + 1)
      : -(Math.trunc(Date.now() % 2000000000) + realizationIndex + 1);
    const resampleRng = makeRng(resampleIdum);

    const tResampleStart = Date.now();
    const bootstrapCube = genFlipBootstrapSample(
      observedDC, modelDC, bsCent, velBlockSize, resampleRng
    );
    timings.resampleMs = Date.now() - tResampleStart;

    if (isTraceDebug()) {
      const flux = bootstrapCube.flux;
      let sum = 0, mn = Infinity, mx = -Infinity;
      for (let i = 0; i < flux.length; i++) {
        const v = flux[i];
        if (v === v) { // skip NaN, matches np.nansum/nanmin/nanmax
          sum += v;
          if (v < mn) mn = v;
          if (v > mx) mx = v;
        }
      }
      console.error('TRACE resampled cube', 'shape', bootstrapCube.dh.nChannels, bootstrapCube.dh.nPixels[1], bootstrapCube.dh.nPixels[0],
        'sum', sum, 'min', mn, 'max', mx,
        'px[0,0,0]', flux[0], 'px[-1,-1,-1]', flux[flux.length - 1]);
    }

    // ---- Write the bootstrap cube as real FITS bytes (feeds SoFiA) ----
    const cubeFitsBytes = await dataCubeToFitsBytes(cfitsio, bootstrapCube, resampleBeam);

    // One-off diagnostic dump (Fortran-vs-JS resampled-cube divergence
    // investigation) -- writes the exact bytes just built above to disk when
    // set, so they can be numpy-diffed against Fortran's own
    // BootstrapCubes/*_Bootstrap_N.fits for the same realization/seed. No
    // effect unless TRACE_DUMP_RESAMPLE_PATH is set.
    if (typeof process !== 'undefined' && process.env && process.env.TRACE_DUMP_RESAMPLE_PATH) {
      require('fs').writeFileSync(process.env.TRACE_DUMP_RESAMPLE_PATH, Buffer.from(cubeFitsBytes));
    }

    // BUG FIX (2026, flagged by Dan, reported to Nathan): the geometry
    // estimate (paEst/incEst below) needs the beam size in pixels at the
    // SoFiA catalogue step, same as native's SoFiA_Driver.py.LoadSoFiAOutput,
    // which computes it as abs(BMAJ/CDELT1) read from
    // GalaxyDict['CubeName']'s header -- the ORIGINAL, un-resampled source
    // cube (Bootstrap_Error_Analysis.py:GetCubeHeader, called on CubeName,
    // NOT the per-realization resampled CubeNameU), read once and identical
    // for every bootstrap realization of a given galaxy. That is NOT the
    // `beamMajorAxis` payload value above (the product of Fortran's own
    // separate, multi-step float32 InputUnitConversions.f pixel-conversion
    // chain, computed once during the initial fit) -- traced directly, the
    // two disagree by ~8.8e-8 relative for a real realization, enough for
    // arccos to produce a measurably different inclination, propagate
    // through the initial parameter vector for ~130 evaluations, and
    // eventually flip a near-tied Nelder-Mead simplex comparison.
    //
    // A first attempt read BMAJ/CDELT1 back from cubeFitsBytes (this
    // realization's own resampled cube, just written above) instead --
    // WRONG: dataCubeToFitsBytes derives that cube's BMAJ from
    // resampleBeam.beamMajorAxis, i.e. the very same suspect payload value,
    // so reading it back just recovers ~beamMajorAxis plus FITS round-trip
    // noise (confirmed: gave 4.994680881499988, not native's
    // 4.9946813188520105). Reading BMAJ/CDELT1 from rawObservedBytes (the
    // untouched original cube already in memory, decoded above at
    // "Load RAW ... flux for both cubes") instead reproduces native's value
    // exactly (verified bit-for-bit: 4.9946813188520105).
    const bmajForSofia   = await cfitsio.readKeyDouble(rawObservedBytes, 'BMAJ');
    const cdelt1ForSofia = await cfitsio.readKeyDouble(rawObservedBytes, 'CDELT1');
    const beamPixForSofia = Math.abs(bmajForSofia / cdelt1ForSofia);

    // ---- Run SoFiA (replaces the native sofia binary) ----
    const lines = sofiaParTemplate.split('\n').map((l) => l + '\n');
    if (lines.length && lines[lines.length - 1] === '\n') lines.pop();
    lines[17]  = 'input.data                 = /work/cube.fits\n';
    lines[142] = 'output.directory           = /work/out\n';
    lines[143] = 'output.filename            = result\n';
    lines[149] = 'output.writeMask           =  true \n';
    const par = lines.join('');

    const tSofiaStart = Date.now();
    const { exitCode, files } = await sofia.run({ cube: cubeFitsBytes, par });
    timings.sofiaMs = Date.now() - tSofiaStart;

    if (exitCode !== 0 || !files.has('result_cat.txt')) {
      report.sofiaFailed = true;
      timings.totalMs = Date.now() - tEntry;
      report.timings = timings;
      return report;
    }

    const parsed = parseSoFiACatalog(files.get('result_cat.txt'));
    if (!parsed.sofiaSuccess) {
      report.sofiaFailed = true;
      timings.totalMs = Date.now() - tEntry;
      report.timings = timings;
      return report;
    }
    if (!files.has('result_mask.fits')) {
      report.sofiaFailed = true;
      report.error = 'SoFiA succeeded but produced no result_mask.fits';
      timings.totalMs = Date.now() - tEntry;
      report.timings = timings;
      return report;
    }

    const { paEst, incEst } = getGeometryEstimates(
      { ell_maj: parsed.ellMaj, ell_min: parsed.ellMin, kin_pa: parsed.pa, ell_pa: parsed.pa },
      beamPixForSofia, 'WALLABY_Like'
    );
    if (isTraceDebug()) {
      console.error('TRACE SoFiA catalog', 'ellMaj', parsed.ellMaj, 'ellMin', parsed.ellMin, 'kin_pa', parsed.pa, 'maskVal', parsed.maskVal, 'x', parsed.x, 'y', parsed.y, 'z', parsed.z);
      console.error('TRACE geometry estimate', 'paEst', paEst, 'incEst', incEst);
      console.error('TRACE beamPix', beamPixForSofia, 'BMAJ', bmajForSofia, 'CDELT1', cdelt1ForSofia);
    }

    // ---- Mask: read SoFiA's FITS bytes into a DataCube, then isolate the
    //      selected source ID (Fortran/SoFiA_Driver.py's AdjustMaskFile:
    //      MDataNew = (MData == TargVal).astype(int)) ----
    const maskDC = await fitsBytesToDataCube(cfitsio, files.get('result_mask.fits'), bootstrapCube);
    for (let i = 0; i < maskDC.flux.length; i++) {
      maskDC.flux[i] = (maskDC.flux[i] === parsed.maskVal) ? f32(1.0) : f32(0.0);
    }

    // ---- Brightness conversion to Jy/pixel, now that SoFiA has already
    //      seen the raw cube (matches Fortran's DCBrightnessConversion,
    //      InputUnitConversions.f:273-294 -- called right after SoFiA, before
    //      fitting). Everything downstream needs Jy/pixel. ----
    const beamAreaPixels = resampleBeam.beamAreaPixels;
    for (let i = 0; i < bootstrapCube.flux.length; i++) {
      bootstrapCube.flux[i] = f32(bootstrapCube.flux[i] / beamAreaPixels);
    }
    if (isTraceDebug()) {
      console.error('TRACE brightness conversion', 'beamAreaPixels', beamAreaPixels);
    }

    // ---- InitialAnalysis: moment maps, shape, radial profiles, noise,
    //      initial parameter vector + fitting options -- run fresh on THIS
    //      realization's own resampled cube, matching Fortran's RunWRKP
    //      (every bootstrap fit, not just the initial one). ----
    const tInitAnalysisStart = Date.now();
    const initResult = initialAnalysis(bootstrapCube, maskDC, {
      centerSource,
      catalogueEllipseIncDeg: incEst,
      catalogueEllipsePADeg:  paEst,
      nRingsPerBeam,
      nTargRings,
      radGridArcsec,
      sdSwitch,
      linSDLims,
      logSDLims,
      vRotLims,
      sizeLims,
      noiseSigmaLim,
      beamMajorAxis,
      // Required -- without these, every realization silently fit 65 free
      // params instead of Fortran's real 16 (found via a chi2 cross-check
      // against Fortran's ProbeSwitch output).
      constParams,
      fixedParams,
    });
    timings.fixtureFitMs = Date.now() - tInitAnalysisStart; // schema name kept for run_both.js compatibility
    if (isTraceDebug()) {
      console.error('TRACE noise', initResult.noise, 'uncertainty', bootstrapCube.dh.uncertainty);
    }

    // ---- Build fit state ----
    const trfo = initResult.trFittingOptions;
    const pvIni = initResult.pvIni;

    const pvModel    = new ParameterVector();
    const pvFirstFit = new ParameterVector();
    pvModel.nParams  = pvIni.nParams;
    allocateParamVector(pvModel);
    pvModel.param.set(pvIni.param);
    pvModel.paramLowerLims.set(pvIni.paramLowerLims);
    pvModel.paramUpperLims.set(pvIni.paramUpperLims);
    pvModel.cyclicSwitch.set(pvIni.cyclicSwitch);
    pvModel.paramRange.set(pvIni.paramRange);

    const modelTR = initResult.modelTiltedRing;
    // cmode/cloudBaseSurfDens: FittingOptionsInputs.f config values, not
    // part of InitialAnalysis's own output (see GalaxyFit.js/FullModelComparison.js
    // for how these drive the chi2 objective's Monte Carlo cloud density).
    modelTR.cmode = cmode;
    modelTR.cloudBaseSurfDens = f32(cloudBaseSurfDens);

    // Beam for the FIT (convolution) -- distinct object from resampleBeam
    // above (that one only needed beamMajorAxis/beamMinorAxis for FITS
    // header + resampling geometry; this one needs the full kernel).
    const fitBeam = new Beam2D();
    fitBeam.beamSigmaVector[0] = f32(observedBeam.beamSigma0);
    fitBeam.beamSigmaVector[1] = f32(observedBeam.beamSigma1);
    fitBeam.beamSigmaVector[2] = f32(observedBeam.beamSigma2);
    fitBeam.nRadialCells = observedBeam.nRadialCells;
    fitBeam.sigmaLengths = f32(observedBeam.sigmaLengths);
    fitBeam.pixelSize[0] = f32(observedBeam.pixelSizeX);
    fitBeam.pixelSize[1] = f32(observedBeam.pixelSizeY);
    fitBeam.beamMajorAxis = f32(observedBeam.beamMajorAxis);
    fitBeam.beamMinorAxis = f32(observedBeam.beamMinorAxis);
    const n = fitBeam.nRadialCells;
    const kSz = 2 * n + 1;
    fitBeam.kernel = new Float32Array(kSz * kSz);
    // Fortran's exact (unrounded) padded size -- real fftw3wasm handles any
    // size natively, no next-power-of-2 fallback needed.
    fitBeam.paddedSize[0] = 2 * n + 1 + bootstrapCube.dh.nPixels[0];
    fitBeam.paddedSize[1] = 2 * n + 1 + bootstrapCube.dh.nPixels[1];
    fitBeam.complexSize[0] = Math.trunc(fitBeam.paddedSize[0] / 2) + 1;
    fitBeam.complexSize[1] = fitBeam.paddedSize[1];
    fitBeam.complexKernelCreated = false;
    calculate2DBeamKernel(fitBeam, fitBeam.pixelSize);

    const fitModelDC = new DataCube();
    fitModelDC.dh.nPixels[0]     = bootstrapCube.dh.nPixels[0];
    fitModelDC.dh.nPixels[1]     = bootstrapCube.dh.nPixels[1];
    fitModelDC.dh.nChannels      = bootstrapCube.dh.nChannels;
    fitModelDC.dh.pixelSize[0]   = bootstrapCube.dh.pixelSize[0];
    fitModelDC.dh.pixelSize[1]   = bootstrapCube.dh.pixelSize[1];
    fitModelDC.dh.channelSize    = bootstrapCube.dh.channelSize;
    fitModelDC.dh.refLocation[0] = bootstrapCube.dh.refLocation[0];
    fitModelDC.dh.refLocation[1] = bootstrapCube.dh.refLocation[1];
    fitModelDC.dh.refLocation[2] = bootstrapCube.dh.refLocation[2];
    fitModelDC.dh.refVal[0]      = bootstrapCube.dh.refVal[0];
    fitModelDC.dh.refVal[1]      = bootstrapCube.dh.refVal[1];
    fitModelDC.dh.refVal[2]      = bootstrapCube.dh.refVal[2];
    fitModelDC.dh.uncertainty    = bootstrapCube.dh.uncertainty;
    fitModelDC.dh.nValid         = bootstrapCube.dh.nValid;
    allocateDataCube(fitModelDC);
    fitModelDC.dh.start[0] = bootstrapCube.dh.start[0];
    fitModelDC.dh.start[1] = bootstrapCube.dh.start[1];
    fitModelDC.dh.start[2] = bootstrapCube.dh.start[2];
    fitModelDC.flattendValidIndices.set(bootstrapCube.flattendValidIndices);
    fitModelDC.dh.nValid = bootstrapCube.dh.nValid;

    const state = {
      pvIni, pvModel, pvFirstFit,
      modelTiltedRing:  modelTR,
      modelDC:          fitModelDC,
      observedDC:       bootstrapCube,
      observedBeam:     fitBeam,
      trFittingOptions: trfo,
      // Static per-galaxy fit seed, NOT derived from realizationIndex --
      // matches Fortran (RunWRKP never varies idum across realizations).
      rng:               makeRng(fitIdum),
      linearLogSDSwitch: sdSwitch,
      likelihoodSwitch,
      ftol:              f32(ftol),
      iniGuessWidth:     f32(1.0),
      paramToTiltedRing: generalizedParamVectorToTiltedRing,
    };

    // ---- Run optimizer ----
    await fftwWarmUpPromise;
    await trigWarmUpPromise;
    resetConvolveStats();
    resetEvalStats();
    const tFitStart = Date.now();
    const { pvModel: pvBest, noConvergence } = galaxyFit_Simple(state);
    timings.fitMs = Date.now() - tFitStart;
    const { convolveMs, convolveCalls } = getConvolveStats();
    timings.convolveMs = convolveMs;
    timings.convolveCalls = convolveCalls;
    const { evalCount } = getEvalStats();
    timings.evalCount = evalCount;


    const RAD2DEG = 180.0 / Math.PI;

    const outTR = new TiltedRingModel();
    outTR.nRings = trfo.nRings;
    tiltRing_Allocate(outTR);
    generalizedParamVectorToTiltedRing(pvBest, outTR, trfo);

    const R = outTR.r;
    const col = (fn) => R.map(fn);

    // Fortran (FitOutput.f:469-479, WRKP's "PA_kin" writer):
    function toKinematicPA(positionAngleRad) {
      let paOut = f32(f32(positionAngleRad * RAD2DEG) - f32(90.0));
      while (paOut < 0) paOut = f32(paOut + f32(360.0));
      while (paOut > 360) paOut = f32(paOut - f32(360.0));
      return paOut;
    }

    // Ring radii aren't a fit parameter -- converted to arcsec here
    // (Fortran's FitOutput.f formula) since bootstrap fits never reach
    // Fortran's own output stage to do it there.
    const arcsecPerPixel = Math.abs(bootstrapCube.dh.pixelSize[0]);
    const RRmid = trfo.radialProfiles.slice(0, trfo.nRings).map((rp) => f32(rp.rmid * arcsecPerPixel));

    // sigUse -> Msol/pc^2, same conversion Fortran applies at final WRKP
    // output time (FitOutput.f:368-372,511-518; UnitConversions.f).
    const sdConv1 = Math.abs(bootstrapCube.dh.channelSize)
      / (arcsecPerPixel * Math.abs(bootstrapCube.dh.pixelSize[1]));
    function sigUseToMsolPc2(sigUse) {
      const sdJyPerPixel = (sdSwitch === 1) ? Math.pow(10, sigUse) : sigUse;
      const sdJyPerArcsec2 = sdJyPerPixel * sdConv1;
      return sdJyPerArcsec2 / JyAS_To_MsolPC;
    }

    report.sofiaFailed = false;
    report.converged = !noConvergence;
    report.chi2 = pvBest.bestLike;
    report.XCENTER       = col((r) => r.centPos[0]);
    report.YCENTER       = col((r) => r.centPos[1]);
    report.INCLINATION   = col((r) => f32(r.inclination * RAD2DEG));
    report.POSITIONANGLE = col((r) => toKinematicPA(r.positionAngle));
    report.VSYS  = col((r) => r.vSys);
    report.VROT  = col((r) => r.vRot);
    report.VDISP = col((r) => r.vDisp);
    report.SURFDENS = col((r) => sigUseToMsolPc2(r.sigUse));
    report.SURFDENS_FACEON = report.SURFDENS; // matches runInitialFit's own alias, below
    report.R    = RRmid;
    report.R_SD = RRmid;
    timings.totalMs = Date.now() - tEntry;
    report.timings = timings;
  } catch (e) {
    report.sofiaFailed = true;
    report.error = (e && e.message) || String(e);
    timings.totalMs = Date.now() - tEntry;
    report.timings = timings;
  }

  return report;
}

// -----------------------------------------------------------------------
// estimateGeometry
//
// Runs SoFiA source-finding once, directly on the raw observed cube (no
// resample, no fit) -- the same catalogue-derived shape estimate a WALLABY-
// style survey pipeline would normally hand this driver ahead of time
// (SoFiA_Driver.WriteSoFiACatFileForWRKP just echoes whatever PA/Inc numbers
// it's given -- see that function's own comment for why there's no
// automated derivation on the initial-fit path today). This gives
// galaxy-fit.html a real, in-pipeline way to produce that same estimate
// instead of requiring it be done by hand beforehand: identical
// SoFiA-run/parse/GeometryEstimates chain runBootstrapRealization already
// runs live for every bootstrap realization, just once, on the un-resampled
// cube, with no fit afterward.
//
// Payload: observedCubeRawFitsB64 (raw, un-brightness-converted FITS
// bytes), sofiaParTemplate (same .par template text every other SoFiA run
// in this file uses).
//
// Self-contained (no closures, all requires inside its own body) for the
// same reason as runBootstrapRealization/runInitialFit -- see this file's
// header comment.
// -----------------------------------------------------------------------
async function estimateGeometry(inputSetElement, payload) {
  const report = {};
  try {
    // Liveness ping -- see runInitialFit's own comment on why this exists.
    if (typeof progress === 'function') progress(0);

    function requirePayloadField(payload, name) {
      const v = payload[name];
      if (v === undefined || v === null) {
        throw new Error(`estimateGeometry: payload.${name} is required but missing.`);
      }
      return v;
    }
    const observedCubeRawFitsB64 = requirePayloadField(payload, 'observedCubeRawFitsB64');
    const sofiaParTemplate = requirePayloadField(payload, 'sofiaParTemplate');

    const { parseSoFiACatalog } = require('./src/SoFiACatalog/ParseSoFiACatalog');
    const { getGeometryEstimates } = require('./src/GeometryEstimates/GeometryEstimates');
    let cfitsio, sofia;
    try {
      cfitsio = require('cfitsio-wasm.js');
    } catch (e) {
      cfitsio = require('../third_party/cfitsio-4.6.3/wasm/cfitsio-wasm.js');
    }
    try {
      sofia = require('sofia-wasm.js');
    } catch (e) {
      sofia = require('../third_party/SoFiA-2-master_2_5_1/wasm/sofia-wasm.js');
    }

    const b64ToBytes = (b64) => Uint8Array.from(atob(b64), (c) => c.charCodeAt(0));
    const rawObservedBytes = b64ToBytes(observedCubeRawFitsB64);

    // Same beam-pixels formula runBootstrapRealization's own SoFiA step
    // uses, read off this SAME (raw, un-resampled) cube -- see that
    // function's "BUG FIX (2026...)" comment for why it must be read this
    // way, not derived from any already-converted value.
    const bmaj = await cfitsio.readKeyDouble(rawObservedBytes, 'BMAJ');
    const cdelt1 = await cfitsio.readKeyDouble(rawObservedBytes, 'CDELT1');
    const beamPix = Math.abs(bmaj / cdelt1);

    // Non-fatal here (SoFiA's source-finding doesn't need a velocity axis
    // the way runInitialFit's kinematic modelling does, and a PA/Inc
    // estimate is still produced), but surfaced as a warning so a user
    // picking the wrong cube variant (e.g. a raw "_cube.fits" instead of
    // "_VelCube.fits") finds out immediately, before spending a real
    // dispatch on a fit that's guaranteed to fail this same check for real
    // (runInitialFit) -- see that function's own matching, fatal check for
    // why (src/Inputs/InputUnitConversions.f's ChannelVelocityConversion).
    const ACCEPTABLE_VELOCITY_CTYPE3 = ['VELO-LSR', 'VELOHEL', 'VOPT', 'VELO'];
    const ctype3ForWarning = (await cfitsio.readKeyString(rawObservedBytes, 'CTYPE3')).trim();
    let velocityAxisWarning = null;
    if (!ACCEPTABLE_VELOCITY_CTYPE3.includes(ctype3ForWarning)) {
      velocityAxisWarning =
        `This cube's spectral axis (CTYPE3=${JSON.stringify(ctype3ForWarning)}) isn't a velocity axis `
        + `-- the initial fit will fail on it. If CTYPE3 is FREQ, use the velocity-calibrated variant `
        + `instead (e.g. a "_VelCube.fits" file, not "_cube.fits").`;
    }

    // OBJECT is optional -- not every cube header carries it, and this
    // estimate is still useful without it (PA/Inc don't depend on it).
    let objectName = null;
    try {
      const raw = (await cfitsio.readKeyString(rawObservedBytes, 'OBJECT')).trim();
      if (raw) objectName = raw;
    } catch (e) { /* no OBJECT keyword in this header -- fine */ }

    const lines = sofiaParTemplate.split('\n').map((l) => l + '\n');
    if (lines.length && lines[lines.length - 1] === '\n') lines.pop();
    lines[17]  = 'input.data                 = /work/cube.fits\n';
    lines[142] = 'output.directory           = /work/out\n';
    lines[143] = 'output.filename            = result\n';
    lines[149] = 'output.writeMask           =  true \n';
    const par = lines.join('');

    const { exitCode, files } = await sofia.run({ cube: rawObservedBytes, par });
    if (exitCode !== 0 || !files.has('result_cat.txt')) {
      report.sofiaFailed = true;
      report.error = 'SoFiA found no source in this cube';
      return report;
    }
    const parsed = parseSoFiACatalog(files.get('result_cat.txt'));
    if (!parsed.sofiaSuccess) {
      report.sofiaFailed = true;
      report.error = 'SoFiA ran but produced no usable catalogue entry';
      return report;
    }

    const { paEst, incEst } = getGeometryEstimates(
      { ell_maj: parsed.ellMaj, ell_min: parsed.ellMin, kin_pa: parsed.pa, ell_pa: parsed.pa },
      beamPix, 'WALLABY_Like'
    );

    report.sofiaFailed = false;
    report.paEstDeg = paEst;
    report.incEstDeg = incEst;
    report.objectName = objectName;
    report.velocityAxisWarning = velocityAxisWarning;
  } catch (e) {
    report.sofiaFailed = true;
    report.error = (e && e.message) || String(e);
  }
  return report;
}

// -----------------------------------------------------------------------
// runInitialFit
// The ONE anchor/initial fit that seeds a galaxy's whole run -- Fortran's
// equivalent is SingleGalaxyTests.f's main program run with BSSwitch=0
// (RunWRKP.py), NOT a bootstrap realization. Two structural differences
// from runBootstrapRealization, both confirmed by reading the Fortran main
// program directly (SingleGalaxyTests.f:44-110), not assumed:
//   - NO resample: fits the raw observed/model cubes directly (Fortran
//     loads ObservedDC and fits it as-is -- genFlipBootstrapSample is a
//     bootstrap-realization-only step).
//   - NO live SoFiA: geometry (PA/Inc) is STATIC, payload-supplied degrees
//     (payload.paEstDeg/incEstDeg), matching Fortran's MakeCatalogue()/
//     PFlags%CatFlag just reading two numbers off a 2-line text file
//     (SoFiA_Driver.WriteSoFiACatFileForWRKP) -- there is no in-pipeline
//     shape estimation for the initial fit in Fortran either. Live SoFiA
//     only ever runs per bootstrap realization (Bootstrap_Error_Analysis.
//     GetBootstrapModel), which this function has no part of.
//
// From initialAnalysis() onward, the fit-state-build/optimizer/report-field
// logic below is intentionally the SAME code as runBootstrapRealization's
// (not extracted into a shared helper -- see this file's header: DCP's
// compute.for ships only a single target function's own source to remote
// sandboxes, so a call out to a module-scope helper would work locally but
// break on a real dispatch; duplication here is the correct, DCP-safe
// choice, not an oversight).
//
// Unlike runBootstrapRealization, this function does NOT write any files
// itself (same reasoning: a real DCP worker has no access to the local
// host's filesystem paths) -- the model-cube FITS and the diskfit_fixture.
// json-equivalent fixture are both returned as data (modelCubeFitsB64,
// fixtureJson) for the CALLER (RunInitialFitDCP.py, or the local/DCP CLI
// driver below) to persist wherever GeneralDict['FixtureFile']/
// BestFitModel['ModelCube'] expect them.
//
// Payload (built by RunInitialFitDCP.py's own payload builder -- NOT
// RunBootstrapsDCP.BuildRealizationPayload, which is bootstrap-realization-
// specific): observedCubeRawFitsB64, observedMaskRawFitsB64 (raw, un-
// brightness-converted FITS bytes for the observed cube and its mask --
// the pre-existing catalogue mask, NOT a live SoFiA output), paEstDeg,
// incEstDeg (static geometry, degrees), sigmaLengths (fitting-options
// config value, determines the convolution kernel's half-width),
// centerSource, nRingsPerBeam, nTargRings, radGridArcsec, sdSwitch,
// constParams, fixedParams, vRotLims, sizeLims, noiseSigmaLim, cmode,
// cloudBaseSurfDens, fitIdum, ftol, likelihoodSwitch -- same fit-side
// fields runBootstrapRealization's payload carries, since from
// initialAnalysis() onward the two paths are identical. Deliberately NOT
// given a pre-built observedDC header, observedBeam object, or linSDLims/
// logSDLims the way runBootstrapRealization's payload is/has (those come
// from a prior Fortran run's fixture, which doesn't exist yet for the
// initial fit) -- all three are derived here directly from the raw FITS
// bytes instead (see the unit-
// conversion block below), which is what makes this function usable for a
// genuinely fresh galaxy with zero Fortran involvement.
// -----------------------------------------------------------------------
async function runInitialFit(realizationIndex, payload) {
  // realizationIndex is unused -- there's only ever one anchor fit, no per-
  // slice index to act on -- but compute.for(inputSet, workFn, extraArgs)
  // always calls workFn(element, ...extraArgs), so the parameter list has
  // to match that shape (same as runBootstrapRealization's) for the real-
  // DCP dispatch below to bind payload correctly instead of accidentally
  // binding the inputSet element to it.
  const report = {};
  const tEntry = Date.now();
  const timings = {
    resampleMs: null, sofiaMs: null, fixtureFitMs: null,
    fitMs: null, convolveMs: null, convolveCalls: null, evalCount: null, totalMs: null,
  };

  try {
    function isTraceDebug() {
      return typeof process !== 'undefined' && process.env && process.env.TRACE_DEBUG === '1';
    }

    // Liveness ping for DCP's sandbox supervisor, which kills a slice with
    // ENOPROGRESS if it never sees a progress() call in time -- confirmed
    // directly: a real dispatch of this function was repeatedly killed with
    // "Sandbox never emitted a progress event" before this existed. Matches
    // amoeba's own first onProgress(0) call (GalaxyFit.js), so no backward
    // jump once the optimizer's own reporting takes over.
    if (typeof progress === 'function') progress(0);

    function requirePayloadField(payload, name) {
      const v = payload[name];
      if (v === undefined || v === null) {
        throw new Error(`runInitialFit: payload.${name} is required but missing `
          + `(no silent default -- see runBootstrapRealization's header comment for why).`);
      }
      return v;
    }

    const {
      nRingsPerBeam, nTargRings, radGridArcsec,
      vRotLims, sizeLims, noiseSigmaLim, constParams, fixedParams, fitIdum, ftol,
    } = payload;
    const observedCubeRawFitsB64 = requirePayloadField(payload, 'observedCubeRawFitsB64');
    const observedMaskRawFitsB64 = requirePayloadField(payload, 'observedMaskRawFitsB64');
    const paEstDeg               = requirePayloadField(payload, 'paEstDeg');
    const incEstDeg              = requirePayloadField(payload, 'incEstDeg');
    const centerSource           = requirePayloadField(payload, 'centerSource');
    const sdSwitch               = requirePayloadField(payload, 'sdSwitch');
    const cmode                  = requirePayloadField(payload, 'cmode');
    const cloudBaseSurfDens      = requirePayloadField(payload, 'cloudBaseSurfDens');
    const likelihoodSwitch       = requirePayloadField(payload, 'likelihoodSwitch');
    // sigmaLengths: a fitting-options config value (FittingOptionsInputs.f),
    // NOT derived from the cube/beam itself -- determines the convolution
    // kernel's half-width (nRadialCells = trunc(|beamSigma0*sigmaLengths|)).
    const sigmaLengths           = requirePayloadField(payload, 'sigmaLengths');

    const f32 = Math.fround;
    // Published-package path routes through ONE flat require -- confirmed
    // directly on a real dispatch that job.requires() itself (used at
    // job-deployment time, before the sandbox starts) mangles nested paths
    // ("Package 3kidnas-test version 0.1.0 does not contain file
    // src,BootstrapSampler,CubeDifference.js" -- slashes replaced with
    // commas), independently of module.declare()'s own dependency array
    // and in-sandbox require() calls, which DO resolve nested relative
    // paths correctly (Phase 0's probe, ModuleB.js -> ./sub/ModuleA.js).
    // initialFitEntry.js (package/initialFitEntry.js) is the one flat file
    // that internally require()s all 14 of these by relative path and
    // re-exports them -- see its own header comment. Falls back to the
    // ordinary local relative requires for standalone `node --local` runs,
    // same try/catch shape this file already uses for cfitsio-wasm.js/
    // sofia-wasm.js below. See ~/.claude/plans/breezy-launching-nova.md.
    let DataCube, allocateDataCube, Beam2D, allocate_Beam2D, makeRng, dataCubeToFitsBytes,
        fitsBytesToDataCube, initialAnalysis, computeSNDiagnostics, TiltedRingModel,
        tiltRing_Allocate, ParameterVector, allocateParamVector, calculate2DBeamKernel,
        resetConvolveStats, getConvolveStats, warmUp, generalizedParamVectorToTiltedRing,
        warmUpWasmTrig, galaxyFit_Simple, resetEvalStats, getEvalStats,
        tiltedRingModelComparison, JyAS_To_MsolPC, constructMomentMaps, flatIndxCalc;
    try {
      const entry = require('initialFitEntry.js');
      ({ DataCube, allocateDataCube, flatIndxCalc } = entry.DataCube);
      ({ Beam2D, allocate_Beam2D } = entry.Beam);
      ({ makeRng } = entry.random);
      ({ dataCubeToFitsBytes, fitsBytesToDataCube } = entry.DataCubeFits);
      ({ initialAnalysis, computeSNDiagnostics } = entry.InitialAnalysis);
      ({ TiltedRingModel, tiltRing_Allocate } = entry.TiltedRing);
      ({ ParameterVector, allocateParamVector } = entry.ParameterVector);
      ({ calculate2DBeamKernel } = entry.CalculateBeamKernel);
      ({ resetConvolveStats, getConvolveStats, warmUp } = entry.CubeKernelConvolution);
      ({ generalizedParamVectorToTiltedRing } = entry.ParameterToTiltedRingVector);
      ({ warmUpWasmTrig } = entry.TiltedRingModelGeneration);
      ({ galaxyFit_Simple } = entry.GalaxyFit);
      ({ resetEvalStats, getEvalStats, tiltedRingModelComparison } = entry.FullModelComparison);
      ({ JyAS_To_MsolPC } = entry.BasicConstants);
      ({ constructMomentMaps } = entry.GetMomentMaps);
    } catch (publishedPackageError) {
      // If the LOCAL fallback also fails (expected/normal in a real DCP
      // sandbox -- there's no local filesystem there at all), surface the
      // ORIGINAL published-package error too, not just the fallback's own
      // "module not available" -- that message alone is a red herring (it
      // names whatever this catch block's own first require() call was,
      // not the real cause), and swallowing the real error here is exactly
      // what made a real published-package failure look identical to
      // "local file not found" the first time this happened.
      try {
        ({ DataCube, allocateDataCube, flatIndxCalc } = require('./src/ObjectDefinitions/DataCube'));
        ({ Beam2D, allocate_Beam2D } = require('./src/ObjectDefinitions/Beam'));
        ({ makeRng } = require('./src/StandardMath/random'));
        ({ dataCubeToFitsBytes, fitsBytesToDataCube } = require('./src/BootstrapSampler/DataCubeFits'));
        ({ initialAnalysis, computeSNDiagnostics } = require('./src/PreAnalysis/InitialAnalysis'));
        ({ TiltedRingModel, tiltRing_Allocate } = require('./src/ObjectDefinitions/TiltedRing'));
        ({ ParameterVector, allocateParamVector } = require('./src/ObjectDefinitions/ParameterVector'));
        ({ calculate2DBeamKernel } = require('./src/ConvolveCube/CalculateBeamKernel'));
        ({ resetConvolveStats, getConvolveStats, warmUp } = require('./src/ConvolveCube/CubeKernelConvolution'));
        ({ generalizedParamVectorToTiltedRing } = require('./src/ParameterToTiltedRingInterface/ParameterToTiltedRingVector'));
        ({ warmUpWasmTrig } = require('./src/TiltedRingModelGeneration/TiltedRingModelGeneration'));
        ({ galaxyFit_Simple } = require('./src/GalaxyAnalysis/GalaxyFit'));
        ({ resetEvalStats, getEvalStats, tiltedRingModelComparison } = require('./src/CompareCubes/FullModelComparison'));
        ({ JyAS_To_MsolPC } = require('./src/StandardMath/BasicConstants'));
        ({ constructMomentMaps } = require('./src/PreAnalysis/GetMomentMaps'));
      } catch (localFallbackError) {
        throw new Error(`Could not load pipeline modules via the published package OR the local fallback. `
          + `Published-package error: ${publishedPackageError.message}. `
          + `Local-fallback error: ${localFallbackError.message}`);
      }
    }
    let cfitsio;
    try {
      cfitsio = require('cfitsio-wasm.js');
    } catch (e) {
      cfitsio = require('../third_party/cfitsio-4.6.3/wasm/cfitsio-wasm.js');
    }

    const fftwWarmUpPromise = warmUp();
    const trigWarmUpPromise = warmUpWasmTrig();

    const b64ToBytes = (b64) => Uint8Array.from(atob(b64), (c) => c.charCodeAt(0));
    // btoa + String.fromCharCode, not Buffer.from(...).toString('base64') --
    // same reasoning as b64ToBytes above (Buffer isn't guaranteed to exist
    // in a DCP worker sandbox; confirmed directly on a real dispatch,
    // "Buffer is not defined"). Chunked (not one String.fromCharCode(...
    // bytes) spread) so a real cube's byte count doesn't blow the engine's
    // max-arguments limit on a single call.
    const bytesToB64 = (bytes) => {
      let bin = '';
      const chunkSize = 8192;
      for (let i = 0; i < bytes.length; i += chunkSize) {
        bin += String.fromCharCode.apply(null, bytes.subarray(i, i + chunkSize));
      }
      return btoa(bin);
    };
    const rawObservedBytes = b64ToBytes(observedCubeRawFitsB64);
    const rawMaskBytes     = b64ToBytes(observedMaskRawFitsB64);

    // ---- Raw-FITS unit conversion (Fortran: src/Inputs/InputUnitConversions.f
    // 's DataCubeUnitConversions/ChannelVelocityConversion, invoked once per
    // galaxy inside PreGalaxyAnalysis -- every OTHER JS code path (bootstrap
    // realizations) receives the already-converted result via the fixture/
    // payload a prior Fortran run produced. This is the initial fit, so
    // there's no prior run to inherit it from -- done here, from the raw
    // FITS bytes, instead.
    //
    // Fortran reads every one of these header keys via cfitsio's single-
    // precision accessors (FTGKYE for BMAJ/BMIN/BPA, FTGKNE for CDELT/CRPIX/
    // CRVAL -- confirmed directly in src/Inputs/DataCubeInput.f), so the
    // FIRST representation of each value is already float32-rounded, before
    // any further arithmetic -- Math.fround() on cfitsio-wasm's double read
    // reproduces that identically (float32 has far less precision than
    // double, so rounding the same stored value to float32 lands on the
    // same result via either path). Each conversion step below is then its
    // own separately-rounded f32 op, matching Fortran's real4 arithmetic
    // throughout this codebase, not a single double-precision expression
    // narrowed once at the end. ----
    const rk = async (key) => f32(await cfitsio.readKeyDouble(rawObservedBytes, key));
    const [cdelt1Deg, cdelt2Deg, cdelt3Raw, crval1Deg, crval2Deg, crval3Raw,
      crpix1, crpix2, crpix3, bmajDeg, bminDeg, bpaDeg] = await Promise.all([
      rk('CDELT1'), rk('CDELT2'), rk('CDELT3'), rk('CRVAL1'), rk('CRVAL2'), rk('CRVAL3'),
      rk('CRPIX1'), rk('CRPIX2'), rk('CRPIX3'), rk('BMAJ'), rk('BMIN'), rk('BPA'),
    ]);
    const cunit3 = (await cfitsio.readKeyString(rawObservedBytes, 'CUNIT3')).trim();
    const ctype3 = (await cfitsio.readKeyString(rawObservedBytes, 'CTYPE3')).trim();

    // Reject unsupported spectral axes the same way Fortran does --
    // ChannelVelocityConversion (src/Inputs/InputUnitConversions.f:311-404)
    // only ever accepts these 4 CTYPE3 values; its own PossibleChannelKeys
    // list also has a 5th entry, 'FREQ', but the acceptability loop only
    // runs i=1..4, so a frequency-axis cube is unconditionally rejected in
    // Fortran too, not just unimplemented here. Confirmed as the real cause
    // of a "FillInMisingRings: no ring with an acceptable surface density
    // found" failure Dan hit by picking a plain "_cube.fits" (CTYPE3=FREQ,
    // the raw un-calibrated observation) instead of "_VelCube.fits"
    // (CTYPE3=VOPT) -- that failure mode gave no indication the cube itself
    // was the problem, so this check exists to fail immediately and clearly
    // instead of several pipeline stages later.
    const ACCEPTABLE_VELOCITY_CTYPE3 = ['VELO-LSR', 'VELOHEL', 'VOPT', 'VELO'];
    const ACCEPTABLE_VELOCITY_CUNIT3 = ['m/s', 'km/s'];
    if (!ACCEPTABLE_VELOCITY_CTYPE3.includes(ctype3)) {
      throw new Error(
        `This cube's spectral axis (CTYPE3=${JSON.stringify(ctype3)}) is not a velocity axis, `
        + `so it can't be used for kinematic fitting. Acceptable CTYPE3 values: `
        + `${ACCEPTABLE_VELOCITY_CTYPE3.join(', ')}. If CTYPE3 is FREQ, this is likely the raw, `
        + `un-calibrated cube -- use the velocity-calibrated variant instead `
        + `(e.g. a "_VelCube.fits" file, not "_cube.fits").`
      );
    }
    if (!ACCEPTABLE_VELOCITY_CUNIT3.includes(cunit3)) {
      throw new Error(
        `This cube's spectral axis units (CUNIT3=${JSON.stringify(cunit3)}) aren't supported for `
        + `kinematic fitting. Acceptable CUNIT3 values: ${ACCEPTABLE_VELOCITY_CUNIT3.join(', ')}.`
      );
    }

    // Spatial: degrees -> arcsec (DataCubeUnitConversions/DegreesToArcSec).
    // FITS WCS angular axes (CTYPE1/2 RA---SIN/DEC--SIN) are always degrees
    // per the FITS standard -- Fortran's own DEGREE/DEGREES/deg unit switch
    // is a validation check on user-text-format cubes, not something a real
    // FITS file varies (matches dataCubeToFitsBytes's own CUNIT1/2='deg').
    const DEG_TO_AS = f32(3600.0);
    const pixelSize0 = f32(cdelt1Deg * DEG_TO_AS);
    const pixelSize1 = f32(cdelt2Deg * DEG_TO_AS);
    const refVal0 = f32(crval1Deg * DEG_TO_AS);
    const refVal1 = f32(crval2Deg * DEG_TO_AS);

    // Spectral: m/s -> km/s (ChannelVelocityConversion/MToKM: LKM=LM/1000.,
    // a DIRECT division -- NOT multiplication by a pre-computed 1/1000
    // reciprocal, which rounds differently: verified directly, dividing
    // this galaxy's real CDELT3 by 1000 gives -4.059688091278076 -- matching
    // Fortran's real fixture value, -4.0596880913, to float32 precision --
    // while multiplying by f32(1/1000) instead gives -4.059688568115234, off
    // by more than 1 ULP) -- only when the cube says m/s; already-km/s cubes
    // need no conversion. Other CTYPE3/CUNIT3 combinations (frequency, etc.)
    // aren't supported for kinematic fitting by Fortran either
    // (ChannelVelocityConversion's own AcceptableType/AcceptableUnits
    // checks) -- not handled here either.
    const channelSize = cunit3 === 'm/s' ? f32(cdelt3Raw / f32(1000.0)) : cdelt3Raw;
    const refVal2 = cunit3 === 'm/s' ? f32(crval3Raw / f32(1000.0)) : crval3Raw;

    // CRPIX -> 0-indexed RefLocation -- an indexing-convention shift, not a
    // unit conversion (FITS is 1-indexed).
    const refLocation0 = f32(crpix1 - 1);
    const refLocation1 = f32(crpix2 - 1);
    const refLocation2 = f32(crpix3 - 1);

    // ---- Build the observed cube header, load RAW flux ----
    // Fortran (SingleGalaxyTests.f:44-51): ReadFullDataCube, THEN brightness-
    // convert immediately -- no SoFiA run sits between them for the initial
    // fit, unlike a bootstrap realization.
    const observedDC = new DataCube();
    const odh = observedDC.dh;
    odh.pixelSize[0] = pixelSize0; odh.pixelSize[1] = pixelSize1; odh.channelSize = channelSize;
    odh.refLocation[0] = refLocation0; odh.refLocation[1] = refLocation1; odh.refLocation[2] = refLocation2;
    odh.refVal[0] = refVal0; odh.refVal[1] = refVal1; odh.refVal[2] = refVal2;
    odh.uncertainty = f32(0.0); // overwritten by initialAnalysis()'s own getNoise() below
    const rawObservedDC = await fitsBytesToDataCube(cfitsio, rawObservedBytes, observedDC);
    observedDC.dh = rawObservedDC.dh;       // nPixels/nChannels filled in from the FITS data itself
    observedDC.flux = rawObservedDC.flux;
    observedDC.flattendValidIndices = rawObservedDC.flattendValidIndices;
    observedDC.channels = rawObservedDC.channels;

    // Beam: degrees -> arcsec -> pixels (major/minor), degrees -> radians
    // (position angle) -- InputUnitConversions.f:219-222,249-250,259-262,273-274.
    // allocate_Beam2D (already used identically for runBootstrapRealization's
    // resampleBeam) then derives beamSigmaVector/nRadialCells/paddedSize/
    // complexSize from these, exactly like it would from a fresh Fortran
    // BMAJ/BMIN/BPA read -- this is genuinely the SAME "build a beam from
    // scratch" path, not a bootstrap-only shortcut.
    const beamMajorAxisAS = f32(bmajDeg * DEG_TO_AS);
    const beamMinorAxisAS = f32(bminDeg * DEG_TO_AS);
    const fitBeam = new Beam2D();
    fitBeam.pixelSize[0] = pixelSize0;
    fitBeam.pixelSize[1] = pixelSize1;
    fitBeam.beamMajorAxis = f32(beamMajorAxisAS / Math.abs(pixelSize0));
    fitBeam.beamMinorAxis = f32(beamMinorAxisAS / Math.abs(pixelSize0));
    fitBeam.beamPositionAngle = f32(f32(bpaDeg * Math.PI) / f32(180.0));
    fitBeam.sigmaLengths = f32(sigmaLengths);
    allocate_Beam2D(fitBeam, [observedDC.dh.nPixels[0], observedDC.dh.nPixels[1]]);
    calculate2DBeamKernel(fitBeam, fitBeam.pixelSize);
    const beamAreaPixels = fitBeam.beamAreaPixels;

    // Brightness conversion to Jy/pixel -- SingleGalaxyTests.f:51,
    // DCBrightnessConversion, right after load (matches Fortran's order
    // exactly here; a bootstrap realization does this AFTER its own SoFiA
    // call instead -- see runBootstrapRealization's comment on why).
    for (let i = 0; i < observedDC.flux.length; i++) {
      observedDC.flux[i] = f32(observedDC.flux[i] / beamAreaPixels);
    }
    if (isTraceDebug()) {
      console.error('TRACE initial-fit brightness conversion', 'beamAreaPixels', beamAreaPixels);
      console.error('TRACE initial-fit beam pixels', fitBeam.beamMajorAxis, fitBeam.beamMinorAxis,
        'PA rad', fitBeam.beamPositionAngle);
    }

    // ---- Mask: the pre-existing catalogue mask (NOT a live SoFiA output).
    //      Fortran (SingleGalaxyTests.f:57-74, MaskSwitch==1): ANY
    //      |value|>=1e-7 becomes 1., else 0. -- deliberately NOT the
    //      per-source-ID equality test runBootstrapRealization uses on a
    //      freshly-SoFiA'd mask (there's no SoFiA source ID here). ----
    const maskDC = await fitsBytesToDataCube(cfitsio, rawMaskBytes, observedDC);
    for (let i = 0; i < maskDC.flux.length; i++) {
      maskDC.flux[i] = (Math.abs(maskDC.flux[i]) >= 1e-7) ? f32(1.0) : f32(0.0);
    }

    // ---- SD (surface density) limits -- PreAnalysis.f:26-39's
    //      PreGalaxyAnalysis preamble, recomputed fresh from THIS cube's own
    //      pixel/channel geometry every time it runs (Fortran doesn't cache
    //      this in the fixture either -- RunBootstrapsDCP.py's ComputeSDLims
    //      re-derives it per bootstrap realization the same way, already
    //      proven not precision-sensitive in practice: LIN_SD_LIMS_MSOL is a
    //      fixed PipelineGlobals.f constant, and the MSol/pc^2->Jy/px chain
    //      below only ever feeds optimizer parameter BOUNDS, not the
    //      objective function itself). Computed here (not passed via
    //      payload) so RunInitialFitDCP.py's payload builder never needs
    //      this cube's pixel/channel size at all -- only JS, which already
    //      derived them bit-exactly above, needs to know them. ----
    const LIN_SD_LIMS_MSOL = [f32(0.1), f32(50.0)];
    const sdConv1ForLims = f32(1.0 / Math.abs(f32(pixelSize0 * pixelSize1)));
    const linSDLims = LIN_SD_LIMS_MSOL.map((vMsol) => {
      const sdJyAS2 = f32(vMsol * JyAS_To_MsolPC);      // MSolPc2_To_JyAS2
      const sd1 = f32(sdJyAS2 / sdConv1ForLims);
      return f32(sd1 / Math.abs(channelSize));
    });
    const logSDLims = linSDLims.map((sd) => f32(Math.log10(sd)));
    if (isTraceDebug()) {
      console.error('TRACE initial-fit SD limits', linSDLims, logSDLims);
    }

    // ---- InitialAnalysis: same call shape as runBootstrapRealization's,
    //      static paEst/incEst instead of a live SoFiA catalogue. ----
    const tInitAnalysisStart = Date.now();
    const initResult = initialAnalysis(observedDC, maskDC, {
      centerSource,
      catalogueEllipseIncDeg: incEstDeg,
      catalogueEllipsePADeg:  paEstDeg,
      nRingsPerBeam,
      nTargRings,
      radGridArcsec,
      sdSwitch,
      linSDLims,
      logSDLims,
      vRotLims,
      sizeLims,
      noiseSigmaLim,
      beamMajorAxis: fitBeam.beamMajorAxis,
      constParams,
      fixedParams,
    });
    timings.fixtureFitMs = Date.now() - tInitAnalysisStart;
    if (isTraceDebug()) {
      console.error('TRACE noise', initResult.noise, 'uncertainty', observedDC.dh.uncertainty);
    }

    // ---- Build fit state (identical shape to runBootstrapRealization's) ----
    const trfo = initResult.trFittingOptions;
    const pvIni = initResult.pvIni;

    const pvModel    = new ParameterVector();
    const pvFirstFit = new ParameterVector();
    pvModel.nParams  = pvIni.nParams;
    allocateParamVector(pvModel);
    pvModel.param.set(pvIni.param);
    pvModel.paramLowerLims.set(pvIni.paramLowerLims);
    pvModel.paramUpperLims.set(pvIni.paramUpperLims);
    pvModel.cyclicSwitch.set(pvIni.cyclicSwitch);
    pvModel.paramRange.set(pvIni.paramRange);

    const modelTR = initResult.modelTiltedRing;
    modelTR.cmode = cmode;
    modelTR.cloudBaseSurfDens = f32(cloudBaseSurfDens);

    // fitBeam (beamSigmaVector/nRadialCells/kernel/paddedSize/complexSize)
    // was already fully built above, before brightness conversion -- see
    // that block's comment for why building it early (via allocate_Beam2D
    // from the raw FITS BMAJ/BMIN/BPA) replaces the copy-from-payload
    // construction runBootstrapRealization uses (it has no fixture to copy
    // from).

    const fitModelDC = new DataCube();
    fitModelDC.dh.nPixels[0]     = observedDC.dh.nPixels[0];
    fitModelDC.dh.nPixels[1]     = observedDC.dh.nPixels[1];
    fitModelDC.dh.nChannels      = observedDC.dh.nChannels;
    fitModelDC.dh.pixelSize[0]   = observedDC.dh.pixelSize[0];
    fitModelDC.dh.pixelSize[1]   = observedDC.dh.pixelSize[1];
    fitModelDC.dh.channelSize    = observedDC.dh.channelSize;
    fitModelDC.dh.refLocation[0] = observedDC.dh.refLocation[0];
    fitModelDC.dh.refLocation[1] = observedDC.dh.refLocation[1];
    fitModelDC.dh.refLocation[2] = observedDC.dh.refLocation[2];
    fitModelDC.dh.refVal[0]      = observedDC.dh.refVal[0];
    fitModelDC.dh.refVal[1]      = observedDC.dh.refVal[1];
    fitModelDC.dh.refVal[2]      = observedDC.dh.refVal[2];
    fitModelDC.dh.uncertainty    = observedDC.dh.uncertainty;
    fitModelDC.dh.nValid         = observedDC.dh.nValid;
    allocateDataCube(fitModelDC);
    fitModelDC.dh.start[0] = observedDC.dh.start[0];
    fitModelDC.dh.start[1] = observedDC.dh.start[1];
    fitModelDC.dh.start[2] = observedDC.dh.start[2];
    fitModelDC.flattendValidIndices.set(observedDC.flattendValidIndices);
    fitModelDC.dh.nValid = observedDC.dh.nValid;

    const state = {
      pvIni, pvModel, pvFirstFit,
      modelTiltedRing:  modelTR,
      modelDC:          fitModelDC,
      observedDC:       observedDC,
      observedBeam:     fitBeam,
      trFittingOptions: trfo,
      rng:               makeRng(fitIdum),
      linearLogSDSwitch: sdSwitch,
      likelihoodSwitch,
      ftol:              f32(ftol),
      iniGuessWidth:     f32(1.0),
      paramToTiltedRing: generalizedParamVectorToTiltedRing,
    };

    // ---- Run optimizer ----
    await fftwWarmUpPromise;
    await trigWarmUpPromise;
    resetConvolveStats();
    resetEvalStats();
    const tFitStart = Date.now();
    const { pvModel: pvBest, noConvergence } = galaxyFit_Simple(state);
    timings.fitMs = Date.now() - tFitStart;
    const { convolveMs, convolveCalls } = getConvolveStats();
    timings.convolveMs = convolveMs;
    timings.convolveCalls = convolveCalls;
    const { evalCount } = getEvalStats();
    timings.evalCount = evalCount;

    const RAD2DEG = 180.0 / Math.PI;

    const outTR = new TiltedRingModel();
    outTR.nRings = trfo.nRings;
    tiltRing_Allocate(outTR);
    generalizedParamVectorToTiltedRing(pvBest, outTR, trfo);

    const R = outTR.r;
    const col = (fn) => R.map(fn);

    function toKinematicPA(positionAngleRad) {
      let paOut = f32(f32(positionAngleRad * RAD2DEG) - f32(90.0));
      while (paOut < 0) paOut = f32(paOut + f32(360.0));
      while (paOut > 360) paOut = f32(paOut - f32(360.0));
      return paOut;
    }

    const arcsecPerPixel = Math.abs(observedDC.dh.pixelSize[0]);
    const RRmid = trfo.radialProfiles.slice(0, trfo.nRings).map((rp) => f32(rp.rmid * arcsecPerPixel));

    const sdConv1 = Math.abs(observedDC.dh.channelSize)
      / (arcsecPerPixel * Math.abs(observedDC.dh.pixelSize[1]));
    function sigUseToMsolPc2(sigUse) {
      const sdJyPerPixel = (sdSwitch === 1) ? Math.pow(10, sigUse) : sigUse;
      const sdJyPerArcsec2 = sdJyPerPixel * sdConv1;
      return sdJyPerArcsec2 / JyAS_To_MsolPC;
    }

    report.sofiaFailed = false;
    report.FITAchieved = !noConvergence;
    report.converged = !noConvergence;
    report.CHI2 = pvBest.bestLike;
    report.chi2 = pvBest.bestLike;
    report.XCENTER       = col((r) => r.centPos[0]);
    report.YCENTER       = col((r) => r.centPos[1]);
    report.INCLINATION   = col((r) => f32(r.inclination * RAD2DEG));
    report.POSITIONANGLE = col((r) => toKinematicPA(r.positionAngle));
    report.VSYS  = col((r) => r.vSys);
    report.VROT  = col((r) => r.vRot);
    report.VDISP = col((r) => r.vDisp);
    report.SURFDENS = col((r) => sigUseToMsolPc2(r.sigUse));
    report.SURFDENS_FACEON = report.SURFDENS;

    // Raw (un-unit-converted) spectral header values, needed only by a
    // SUBSEQUENT bootstrap-loop payload builder's ComputeBsCent-equivalent
    // (resample-geometry velocity-center calc, RunBootstrapsDCP.py's
    // ComputeBsCent: VCenter = (VSYS - CRVAL3/1000)/(CDELT3/1000) + CRPIX3,
    // deliberately the RAW header values, not this function's own
    // f32-rounded/unit-converted channelSize/refVal2 above). Surfaced here
    // so a payload-builder running OUTSIDE a DCP sandbox (a Node CLI, or a
    // browser page building the next job's payload) never needs its own
    // cfitsio instance just to re-read these three numbers off the cube a
    // second time -- this function already read them once, above.
    report.rawSpectralHeader = { crval3: crval3Raw, crpix3, cdelt3: cdelt3Raw };
    // Raw spatial WCS reference, analogous to rawSpectralHeader above --
    // needed client-side for GeometryCorrection.py:GetGlobalPositionAngle's
    // pixel->sky deprojection (arcsec-scaled CRVAL1/2 and CRPIX1/2-1, from
    // which crval/cdelt in degrees are recovered by /3600 -- see
    // wcsHeader.js for the actual transform).
    report.wcsHeader = {
      refValArcsec: [refVal0, refVal1],
      refPix0Indexed: [refLocation0, refLocation1],
      pixelSizeArcsec: [pixelSize0, pixelSize1],
    };
    report.R    = RRmid;
    report.R_SD = RRmid;

    // ---- Extras a bootstrap realization never needs ----

    // Re-synthesize the model cube AT THE BEST-FIT PARAMS -- state.modelDC
    // otherwise holds whatever the LAST amoeba evaluation computed, which is
    // not necessarily the best vertex. Matches Fortran's OutputBestFit_Simple/
    // FitOutput.f, which does the same re-synthesis (ParamToTiltedRing ->
    // BuildTiltedRingModel -> FillDataCubeWithTiltedRing -> CubeBeamConvolution)
    // before writing the model cube. Deliberately AFTER getEvalStats() above,
    // so timings.evalCount stays exactly what a bootstrap realization reports
    // (pure optimizer evaluations), not polluted by this one extra call.
    // One-off diagnostic (Fortran-vs-JS model-cube divergence investigation):
    // substitute Fortran's own converged parameter vector in place of this
    // fit's own (independently converged, <0.03% different) vector, so the
    // resynthesized model cube can be compared apples-to-apples against
    // Fortran's AverageModel_v1.fits using IDENTICAL inputs. No effect
    // unless TRACE_OVERRIDE_PARAMS is set (a JSON array of 15 numbers).
    let synthParams = pvBest.param;
    if (typeof process !== 'undefined' && process.env && process.env.TRACE_OVERRIDE_PARAMS) {
      synthParams = Float64Array.from(JSON.parse(process.env.TRACE_OVERRIDE_PARAMS));
      console.error('TRACE using overridden params for resynthesis:', Array.from(synthParams));
    }
    tiltedRingModelComparison(synthParams, state);
    // BUG FIX (2026, flagged by Dan): tiltedRingModelComparison leaves
    // fitModelDC.flux in Jy/pixel (the units the fit's own chi2 comparison
    // needs, to match brightness-converted observedDC) -- but
    // RunBootstrapsDCP.BuildRealizationPayload's own header comment is
    // explicit that BestFitModel['ModelCube'] must be RAW, un-brightness-
    // converted Jy/beam, same convention as the observed cube, since the
    // bootstrap loop's resample+SoFiA step needs both cubes in the same raw
    // units. Fortran's own OutputCube (FitOutput.f) does exactly this
    // conversion before writing (`ModelDC%Flux=ModelDC%Flux*BeamArea`,
    // its comment: "Because the cube is in units of Jy/pixel, convert back
    // to Jy/beam") -- missing it here left every js-dcp bootstrap
    // realization resampling against a model cube ~beamAreaPixels times too
    // faint (confirmed directly: a real run's fortranLocal vs jsDcp model
    // cube FITS differed in both sum and max by ~28.2x, matching
    // beamAreaPixels~28.196 almost exactly), which fed a systematically
    // wrong residual into every subsequent flip-resample.
    const modelFlux = fitModelDC.flux;
    for (let i = 0; i < modelFlux.length; i++) {
      modelFlux[i] = f32(modelFlux[i] * beamAreaPixels);
    }
    const modelCubeFitsB64 = bytesToB64(
      await dataCubeToFitsBytes(cfitsio, fitModelDC, fitBeam)
    );

    // Moment maps for the browser's moment-map visualization (see
    // ~/.claude/plans/breezy-launching-nova.md, Phase 1). Observed map
    // (masked) already comes out of initialAnalysis() as
    // initResult.observedMaps -- computed there for shape/center/PA/inc
    // estimation, just never surfaced until now. Model map computed here,
    // deliberately UNmasked (matches MomentMapPlotFncs.py's own convention:
    // the model's Mom0 contour is meant to extend past the data's masked
    // edge in the reference plots, not stop at it) -- fitModelDC at this
    // point is already in Jy/beam (post the beamAreaPixels conversion just
    // above), the same unit convention constructMomentMaps expects.
    const modelMaps = constructMomentMaps(fitModelDC);
    // Flattened to plain [nPixX*nPixY] arrays in simple row-major (x +
    // y*nPixX) order -- NOT DataCube's own channel-innermost flatIndxCalc
    // layout -- so the browser-side canvas code never needs to know about
    // flatIndxCalc at all, just nPixX/nPixY and 3 flat arrays per cube.
    function extractMoments(maps) {
      const { nPixels } = maps.dh;
      const nx = nPixels[0], ny = nPixels[1];
      const mom0 = new Array(nx * ny);
      const mom1 = new Array(nx * ny);
      const mom2 = new Array(nx * ny);
      for (let i = 0; i < nx; i++) {
        for (let j = 0; j < ny; j++) {
          const base = flatIndxCalc(i, j, 0, maps.dh);
          const out = i + j * nx;
          mom0[out] = maps.flux[base];
          mom1[out] = maps.flux[base + 1];
          mom2[out] = maps.flux[base + 2];
        }
      }
      return { mom0, mom1, mom2 };
    }
    report.momentMaps = {
      nPixX: observedDC.dh.nPixels[0],
      nPixY: observedDC.dh.nPixels[1],
      pixelSizeX: observedDC.dh.pixelSize[0],
      pixelSizeY: observedDC.dh.pixelSize[1],
      beamMajPix: Math.abs(fitBeam.beamMajorAxis),
      observed: extractMoments(initResult.observedMaps),
      model: extractMoments(modelMaps),
    };

    // Position-velocity diagrams (see ~/.claude/plans/breezy-launching-nova.md
    // Phase 3). Port of CubeAnalysis.py's ConstructModelBasedPVDiagram --
    // NOT the differently-named-but-unused ConstructPVDiagram (that one
    // takes an externally-supplied center/sizing and isn't part of the
    // BootstrapModelPlot call chain at all; every PV panel there, data
    // AND model, goes through ConstructModelBasedPVDiagram, sized to the
    // MODEL's own radial extent regardless of which cube's flux is fed
    // in). A beam-width pseudo-slit reprojection, not a 1-pixel slice or
    // an interpolated cut -- rotate each cube spaxel into the slit frame,
    // round to an integer output-pixel index, accept only spaxels within
    // +-beamSize/2 of the slit axis, sum every channel of every accepted
    // spaxel.
    //
    // Units: observedDC.flux is Jy/PIXEL here (divided by beamAreaPixels
    // above, for the fit's own chi2 comparison); fitModelDC.flux is
    // Jy/BEAM (multiplied by beamAreaPixels above, to match
    // BestFitModel['ModelCube']'s convention). A residual (data-model)
    // needs both in the SAME units -- converts a temporary Jy/beam copy of
    // the observed flux for PV purposes only, doesn't mutate observedDC.flux
    // itself (still needed in Jy/pixel form by nothing after this point,
    // but safer not to assume that stays true).
    const pixelSizeXAbs = Math.abs(observedDC.dh.pixelSize[0]);
    const dRForPV = report.R_SD.length === 1 ? report.R_SD[0] : report.R_SD[1] - report.R_SD[0];
    const rTestForPV = report.R_SD[report.R_SD.length - 1] + 2 * dRForPV;
    const nRPixForPV = Math.trunc(rTestForPV / pixelSizeXAbs);
    const nSpatialPixForPV = 2 * nRPixForPV + 1;
    const nxPV = observedDC.dh.nPixels[0], nyPV = observedDC.dh.nPixels[1], nChanPV = observedDC.dh.nChannels;

    const observedFluxJyBeam = new Float32Array(observedDC.flux.length);
    for (let i = 0; i < observedFluxJyBeam.length; i++) {
      observedFluxJyBeam[i] = f32(observedDC.flux[i] * beamAreaPixels);
    }
    const diffFluxJyBeam = new Float32Array(observedFluxJyBeam.length);
    for (let i = 0; i < diffFluxJyBeam.length; i++) {
      diffFluxJyBeam[i] = f32(observedFluxJyBeam[i] - fitModelDC.flux[i]);
    }

    function constructModelBasedPV(cubeFlux, dh, angleDeg, centerX, centerY, beamSizePix) {
      let angUse = angleDeg + 90;
      if (angUse > 360) angUse -= 360;
      if (angUse < 0) angUse += 360;
      angUse = (angUse * Math.PI) / 180;
      const cosA = Math.cos(-angUse), sinA = Math.sin(-angUse);
      const pv = new Float64Array(nSpatialPixForPV * nChanPV);
      const half = nSpatialPixForPV / 2; // deliberately NOT Math.trunc -- matches Python's own float division (nSpatialPix/2 in Python 3)
      for (let i = 0; i < nxPV; i++) {
        const x = i - centerX;
        for (let j = 0; j < nyPV; j++) {
          const y = j - centerY;
          const xp = x * cosA - y * sinA;
          const yp = x * sinA + y * cosA;
          const k = Math.round(xp + half);
          if (k < 0 || k >= nSpatialPixForPV) continue;
          if (Math.abs(yp) > beamSizePix / 2) continue;
          const base = flatIndxCalc(i, j, 0, dh);
          const pvBase = k * nChanPV;
          for (let m = 0; m < nChanPV; m++) {
            const v = cubeFlux[base + m];
            if (Number.isFinite(v)) pv[pvBase + m] += v;
          }
        }
      }
      return pv;
    }

    // 4-corner-RMS noise estimator (GetPVNoise/SetCornerLims/CornerSum) --
    // a specific, bespoke estimator (RMS of nonzero pixels in the 4 corners
    // of the PV array), not a global stddev. Zero-valued pixels excluded
    // (treated as "off the edge of the diagram", not real noise) -- exact
    // corner-index math replicated, not approximated.
    function getPVNoise(pv, nSpatial, nChan) {
      const cornerSize = 0.25;
      function lims(shapeDim, step) {
        if (step === 0) return [0, Math.trunc(cornerSize * shapeDim)];
        const high = Math.trunc(shapeDim - 1);
        return [Math.trunc(high - cornerSize * shapeDim), high];
      }
      let nPix = 0, rmsTot = 0;
      for (let si = 0; si < 2; si++) {
        const [xLo, xHi] = lims(nSpatial, si);
        for (let sj = 0; sj < 2; sj++) {
          const [yLo, yHi] = lims(nChan, sj);
          for (let i = xLo; i < xHi; i++) {
            for (let j = yLo; j < yHi; j++) {
              const val = pv[i * nChan + j];
              if (val !== 0) { nPix++; rmsTot += val * val; }
            }
          }
        }
      }
      if (nPix === 0) nPix = 1;
      return Math.sqrt(rmsTot / nPix);
    }

    const centerX = report.XCENTER[0], centerY = report.YCENTER[0];
    const paMajor = report.POSITIONANGLE[0];
    const paMinor = paMajor + 90;

    const pvMajorData = constructModelBasedPV(observedFluxJyBeam, observedDC.dh, paMajor, centerX, centerY, fitBeam.beamMajorAxis);
    const pvMajorModel = constructModelBasedPV(fitModelDC.flux, fitModelDC.dh, paMajor, centerX, centerY, fitBeam.beamMajorAxis);
    const pvMajorDiff = constructModelBasedPV(diffFluxJyBeam, observedDC.dh, paMajor, centerX, centerY, fitBeam.beamMajorAxis);
    const pvMinorData = constructModelBasedPV(observedFluxJyBeam, observedDC.dh, paMinor, centerX, centerY, fitBeam.beamMajorAxis);
    const pvMinorModel = constructModelBasedPV(fitModelDC.flux, fitModelDC.dh, paMinor, centerX, centerY, fitBeam.beamMajorAxis);
    const pvMinorDiff = constructModelBasedPV(diffFluxJyBeam, observedDC.dh, paMinor, centerX, centerY, fitBeam.beamMajorAxis);

    report.pvDiagrams = {
      nSpatialPix: nSpatialPixForPV,
      nChan: nChanPV,
      pixelSizeX: pixelSizeXAbs,
      channelVelsKmS: Array.from({ length: nChanPV }, (_, m) => observedDC.channels[m]),
      major: {
        data: Array.from(pvMajorData), model: Array.from(pvMajorModel), diff: Array.from(pvMajorDiff),
        noiseData: getPVNoise(pvMajorData, nSpatialPixForPV, nChanPV),
        noiseDiff: getPVNoise(pvMajorDiff, nSpatialPixForPV, nChanPV),
      },
      minor: {
        data: Array.from(pvMinorData), model: Array.from(pvMinorModel), diff: Array.from(pvMinorDiff),
        noiseData: getPVNoise(pvMinorData, nSpatialPixForPV, nChanPV),
        noiseDiff: getPVNoise(pvMinorDiff, nSpatialPixForPV, nChanPV),
      },
    };

    // SN_*/RMS diagnostics -- see InitialAnalysis.js's computeSNDiagnostics
    // header for why these are computed here and not inside initialAnalysis()
    // itself (dead weight on the bootstrap-realization hot path).
    const diag = computeSNDiagnostics(initResult.maskedObservedDC, maskDC, initResult.noise, beamAreaPixels);
    // mJy/beam, matching FitOutput.f:424-425's own conversion exactly
    // (ObservedDC%DH%Uncertainty*ObservedBeam%BeamAreaPixels*1000.) --
    // initResult.noise itself is the raw per-pixel internal-unit RMS
    // (EstimateCubeNoise.js), never converted before this point. Missing
    // this conversion previously left report.RMS (and everything
    // downstream reading it, e.g. RunInitialFitDCP.py's BestFitModel)
    // ~28,000x too small versus the Fortran-written _BSModel.txt/
    // _AvgModel_v1.txt for this galaxy (beamAreaPixels~28.2, so
    // ~28.2*1000) -- confirmed by comparing this session's matched-seed
    // fortran-local vs js-dcp run_both.js output.
    report.RMS = f32(f32(initResult.noise * beamAreaPixels) * 1000);
    report.SN_Integrated = diag.snInt;
    report.SN_Peak = diag.snPeak;
    report.SN_Avg = diag.snAvg;
    report.SN_Median = diag.snMedian;
    report.modelCubeFitsB64 = modelCubeFitsB64;

    // diskfit_fixture.json-equivalent -- field names match Fortran's
    // DumpFittingFixture (SingleGalaxyTests.f:174-509) exactly, since
    // RunBootstrapsDCP.BuildRealizationPayload reads this file by those
    // field names regardless of which side wrote it.
    report.fixtureJson = {
      nParams: pvIni.nParams,
      param: Array.from(pvIni.param),
      paramLowerLims: Array.from(pvIni.paramLowerLims),
      paramUpperLims: Array.from(pvIni.paramUpperLims),
      paramRange: Array.from(pvIni.paramRange),
      cyclicSwitch: Array.from(pvIni.cyclicSwitch),
      nRings: trfo.nRings,
      nFittedParamsTotal: trfo.nFittedParamsTotal,
      constParams: Array.from({ length: 13 }, (_, i) => !!trfo.constParams[i]),
      fixedParams: Array.from({ length: 13 }, (_, i) => !!trfo.fixedParams[i]),
      paramLowerLims13: Array.from(trfo.paramLowerLims).slice(0, 13),
      paramUpperLims13: Array.from(trfo.paramUpperLims).slice(0, 13),
      paramRange13: Array.from(trfo.paramRange).slice(0, 13),
      cyclicSwitch13: Array.from(trfo.cyclicSwitch).slice(0, 13),
      radialProfiles: trfo.radialProfiles.slice(0, trfo.nRings).map((rp) => ({
        rmid: rp.rmid, rwidth: rp.rwidth, centPos0: rp.centPos[0], centPos1: rp.centPos[1],
        inclination: rp.inclination, positionAngle: rp.positionAngle, vSys: rp.vSys,
        vRot: rp.vRot, vRad: rp.vRad, vDisp: rp.vDisp, vvert: rp.vvert, dvdz: rp.dvdz,
        sigUse: rp.sigUse, z0: rp.z0, zGradiantStart: rp.zGradiantStart,
      })),
      observedDC: {
        nPixelsX: observedDC.dh.nPixels[0], nPixelsY: observedDC.dh.nPixels[1],
        nChannels: observedDC.dh.nChannels,
        pixelSizeX: observedDC.dh.pixelSize[0], pixelSizeY: observedDC.dh.pixelSize[1],
        channelSize: observedDC.dh.channelSize,
        startX: observedDC.dh.start[0], startY: observedDC.dh.start[1], startV: observedDC.dh.start[2],
        refLocX: observedDC.dh.refLocation[0], refLocY: observedDC.dh.refLocation[1], refLocV: observedDC.dh.refLocation[2],
        refValX: observedDC.dh.refVal[0], refValY: observedDC.dh.refVal[1], refValV: observedDC.dh.refVal[2],
        uncertainty: observedDC.dh.uncertainty,
        nValid: observedDC.dh.nValid,
        flattendValidIndices: Array.from(observedDC.flattendValidIndices).slice(0, observedDC.dh.nValid),
        validFlux: Array.from(observedDC.flattendValidIndices).slice(0, observedDC.dh.nValid)
          .map((idx) => observedDC.flux[idx]),
      },
      observedBeam: {
        beamSigma0: fitBeam.beamSigmaVector[0], beamSigma1: fitBeam.beamSigmaVector[1], beamSigma2: fitBeam.beamSigmaVector[2],
        nRadialCells: fitBeam.nRadialCells, sigmaLengths: fitBeam.sigmaLengths,
        pixelSizeX: fitBeam.pixelSize[0], pixelSizeY: fitBeam.pixelSize[1],
        beamMajorAxis: fitBeam.beamMajorAxis, beamMinorAxis: fitBeam.beamMinorAxis,
        kernel: Array.from(fitBeam.kernel),
      },
      ftol, idum: fitIdum, linearLogSDSwitch: sdSwitch, likelihoodSwitch,
      cmode, cloudBaseSurfDens,
      nRingsPerBeam,
    };

    timings.totalMs = Date.now() - tEntry;
    report.timings = timings;
  } catch (e) {
    report.sofiaFailed = true;
    report.FITAchieved = false;
    report.error = (e && e.message) || String(e);
    if (typeof process !== 'undefined' && process.env && process.env.TRACE_DEBUG === '1') {
      console.error('TRACE runInitialFit stack', e && e.stack);
    }
    timings.totalMs = Date.now() - tEntry;
    report.timings = timings;
  }

  return report;
}

// module.exports for Node (require()'d by the CLI tail below and by
// RunInitialFitDCP.py/RunBootstrapsDCP.py's subprocess dispatch); globals
// for a browser <script> include (galaxy-fit.html), where `module` doesn't
// exist at all.
if (typeof module !== 'undefined' && module.exports) {
  module.exports = { runBootstrapRealization, runInitialFit, estimateGeometry, parseComputeGroups };
} else if (typeof globalThis !== 'undefined') {
  globalThis.runBootstrapRealization = runBootstrapRealization;
  globalThis.runInitialFit = runInitialFit;
  globalThis.estimateGeometry = estimateGeometry;
  globalThis.parseComputeGroups = parseComputeGroups;
}

// =============================================================================
// CLI / DCP job driver -- only runs when this file is executed directly
// (`node bootstrap-realization-launcher.js ...`), not when required as a
// module (RunBootstrapsDCP.py invokes it as a subprocess).
// =============================================================================
function formatResultLine(r) {
  const t = r.timings
    ? ` [total=${r.timings.totalMs}ms resample=${r.timings.resampleMs}ms sofia=${r.timings.sofiaMs}ms initAnalysis=${r.timings.fixtureFitMs}ms fit=${r.timings.fitMs}ms]`
    : '';
  return r.sofiaFailed
    ? `  realization ${r.realizationIndex}: FAILED${r.error ? ' -- ' + r.error : ''}${t}`
    : `  realization ${r.realizationIndex}: converged=${r.converged} chi2=${r.chi2.toFixed(3)}${t}`;
}

function printSummary(results, nBootstraps, dispatchMs, label) {
  const nOk = results.filter((r) => !r.sofiaFailed).length;
  console.log(`\n${nOk}/${nBootstraps} realizations succeeded.`);

  const withTimings = results.filter((r) => r.timings);
  if (withTimings.length) {
    const execTimes = withTimings.map((r) => r.timings.totalMs);
    const sumExec = execTimes.reduce((a, v) => a + v, 0);
    console.log(`\n${label} wall time: ${(dispatchMs / 1000).toFixed(1)}s`);
    console.log(`Worker execution time: min=${Math.min(...execTimes)}ms max=${Math.max(...execTimes)}ms `
      + `mean=${(sumExec / execTimes.length).toFixed(0)}ms sum=${(sumExec / 1000).toFixed(1)}s`);
    console.log(`${label === 'Dispatch' ? 'DCP overhead' : 'Overhead'}: ${((dispatchMs - sumExec) / 1000).toFixed(1)}s`);
  }
}

// Module list for real-DCP dispatch of EITHER work function -- runInitialFit
// is a subset of runBootstrapRealization's own call graph (no
// GeometryEstimates/ParseSoFiACatalog/sofia2wasm, since there's no live
// SoFiA run), so one shared superset list is simpler to keep in sync than
// two nearly-identical ones. Pure data (module-scope constant reference is
// fine here -- job.requires() just reads this array locally, before
// dispatch; it doesn't ship with the work function's own source the way a
// helper FUNCTION call from inside runBootstrapRealization/runInitialFit
// would).
const INITIAL_FIT_AND_BOOTSTRAP_MODULES = [
  './src/ObjectDefinitions/DataCube',
  './src/ObjectDefinitions/Beam',
  './src/ObjectDefinitions/TiltedRing',
  './src/ObjectDefinitions/ParameterVector',
  './src/ObjectDefinitions/Particle',
  './src/BootstrapSampler/FlipBootstrap',
  './src/BootstrapSampler/GenerateBootstrap',
  './src/BootstrapSampler/CubeDifference',
  './src/BootstrapSampler/PhysCoordTransform',
  './src/BootstrapSampler/DataCubeFits',
  './src/StandardMath/random',
  './src/StandardMath/Interpolation',
  './src/StandardMath/BasicConstants',
  './src/StandardMath/fdlibm',
  './src/StandardMath/FullCircTrig',
  './src/SoFiACatalog/ParseSoFiACatalog',
  './src/GeometryEstimates/GeometryEstimates',
  './src/PreAnalysis/InitialAnalysis',
  './src/PreAnalysis/EstimateCubeNoise',
  './src/PreAnalysis/GetMomentMaps',
  './src/PreAnalysis/EstimateShape',
  './src/PreAnalysis/VelProfileAnalysis',
  './src/PreAnalysis/EstimateRadialProfiles',
  './src/PreAnalysis/ModellingInitializations',
  './src/CompareCubes/MaskCube',
  './src/UnitConversions/UnitConversions',
  './src/ConvolveCube/CubeKernelConvolution',
  './src/ConvolveCube/CalculateBeamKernel',
  './src/ConvolveCube/FFTW3WasmRank2',
  './src/ParameterToTiltedRingInterface/ParameterToTiltedRingVector',
  './src/TiltedRingModelGeneration/SingleRingGeneration',
  './src/TiltedRingModelGeneration/TiltedRingModelGeneration',
  './src/GalaxyAnalysis/GalaxyFit',
  './src/CompareCubes/FullModelComparison',
  // Published DCP packages
  'fftw3wasm-v3/fftw-wasm.js',
  'sofia2wasm/sofia-wasm.js',
  'cfitsio4wasm/cfitsio-wasm.js',
  // NOT YET PUBLISHED as a DCP package -- bundled locally in-tree instead
  // (fdlibm-wasm.js + the compiled fdlibm-module.js it requires, both
  // copied into src/StandardMath/ alongside the canonical fdlibm.js).
  // Publishing this as its own DCP package later (mirroring fftw3wasm's
  // publish steps above) would let these two local copies be replaced by a
  // bare package require.
  './src/StandardMath/fdlibm-wasm',
  './src/StandardMath/fdlibm-module',
];
// Exposed for galaxy-fit.html, which needs the exact same list for its own
// job.requires() call -- kept as this ONE array (not a second, hand-copied
// list) so the two dispatch entry points (this file's CLI tail, the
// browser page) can never drift apart. No-op in Node (module.exports above
// already covers that consumer).
if (typeof module === 'undefined' && typeof globalThis !== 'undefined') {
  globalThis.INITIAL_FIT_AND_BOOTSTRAP_MODULES = INITIAL_FIT_AND_BOOTSTRAP_MODULES;
}

// Real-dispatch module list for runInitialFit ONLY. Lists just ONE file
// from 3kidnas-test -- its own flat entry point (package/
// initialFitEntry.js) -- plus the 2 already-published wasm packages it
// needs. NOT a per-file list of all 38 published pipeline files: confirmed
// directly on a real dispatch that job.requires() itself mangles nested
// paths ("Package 3kidnas-test version 0.1.0 does not contain file
// src,BootstrapSampler,CubeDifference.js"), so job.requires() can only
// ever reference this one flat file -- initialFitEntry.js's own
// module.declare() dependency array (not job.requires()) is what pulls in
// the other 38, via the already-proven-working intra-package relative
// require() mechanism. See initialFitEntry.js's own header comment and
// ~/.claude/plans/breezy-launching-nova.md for the full story.
//
// Deliberately separate from INITIAL_FIT_AND_BOOTSTRAP_MODULES above, not
// a replacement for it -- that list is still what the bootstrap loop's
// real dispatch and any Node-CLI real dispatch use, both unaffected by
// today's browser-only module-resolution bug (Node's own local-module
// bundling papers over it there), so leaving them on the old list carries
// no regression risk.
//
// Every file in the 3kidnas-test1 package is listed individually below,
// not just initialFitEntry.js -- confirmed via isolated probes
// (package-probe3 through package-probe6) that the sandbox only
// auto-discovers module.declare()-declared dependencies transitively up
// to 2 deep; beyond that, each file must appear in job.requires() itself
// or it fails with "Module '<path>' is not available."
const INITIAL_FIT_ONLY_MODULES = [
  '3kidnas-test1/initialFitEntry.js',
  '3kidnas-test1/CubeDifference.js',
  '3kidnas-test1/DataCubeFits.js',
  '3kidnas-test1/FlipBootstrap.js',
  '3kidnas-test1/GenerateBootstrap.js',
  '3kidnas-test1/PhysCoordTransform.js',
  '3kidnas-test1/CubeComparison.js',
  '3kidnas-test1/FullModelComparison.js',
  '3kidnas-test1/LikelihoodFunctions.js',
  '3kidnas-test1/MaskCube.js',
  '3kidnas-test1/CalculateBeamKernel.js',
  '3kidnas-test1/CubeKernelConvolution.js',
  '3kidnas-test1/FFTW3WasmRank2.js',
  '3kidnas-test1/GalaxyFit.js',
  '3kidnas-test1/Beam.js',
  '3kidnas-test1/DataCube.js',
  '3kidnas-test1/ParameterVector.js',
  '3kidnas-test1/Particle.js',
  '3kidnas-test1/TiltedRing.js',
  '3kidnas-test1/ParameterToTiltedRingVector.js',
  '3kidnas-test1/EstimateCubeNoise.js',
  '3kidnas-test1/EstimateRadialProfiles.js',
  '3kidnas-test1/EstimateShape.js',
  '3kidnas-test1/GetMomentMaps.js',
  '3kidnas-test1/InitialAnalysis.js',
  '3kidnas-test1/ModellingInitializations.js',
  '3kidnas-test1/VelProfileAnalysis.js',
  '3kidnas-test1/BasicConstants.js',
  '3kidnas-test1/FullCircTrig.js',
  '3kidnas-test1/Interpolation.js',
  '3kidnas-test1/fdlibm-module.js',
  '3kidnas-test1/fdlibm-wasm.js',
  '3kidnas-test1/fdlibm.js',
  '3kidnas-test1/fma.js',
  '3kidnas-test1/random.js',
  '3kidnas-test1/SingleRingGeneration.js',
  '3kidnas-test1/TiltedRingModelGeneration.js',
  '3kidnas-test1/FillDataCubeByTiltedRing.js',
  '3kidnas-test1/UnitConversions.js',
  // runInitialFit itself never calls these two -- they're only here because
  // initialFitEntry.js's own module.declare() dependency array unconditionally
  // requires them (added for runBootstrapRealization's benefit, since both
  // functions share one entry point), so the sandbox needs them listed
  // regardless of which function is doing the requiring. Omitting them
  // reproduces the exact "must list every file" bug this whole list exists
  // to work around, just for these 2 files specifically -- confirmed
  // directly: real initial-fit dispatch failed with "Module
  // './ParseSoFiACatalog.js' is not available." until they were added here.
  '3kidnas-test1/ParseSoFiACatalog.js',
  '3kidnas-test1/GeometryEstimates.js',
  // Published DCP packages
  'fftw3wasm-v3/fftw-wasm.js',
  'cfitsio4wasm/cfitsio-wasm.js',
];
if (typeof module === 'undefined' && typeof globalThis !== 'undefined') {
  globalThis.INITIAL_FIT_ONLY_MODULES = INITIAL_FIT_ONLY_MODULES;
}

// Real-dispatch module list for runBootstrapRealization -- INITIAL_FIT_ONLY_
// MODULES already covers every package file either function needs (see its
// own comment on why ParseSoFiACatalog.js/GeometryEstimates.js are there
// despite runInitialFit never calling them); this just adds the
// already-published sofia2wasm package runBootstrapRealization's own live
// per-realization SoFiA run needs. INITIAL_FIT_AND_BOOTSTRAP_MODULES above
// is untouched and still what Node-CLI real dispatch uses (unaffected by
// the browser-only module-resolution bug), so this is additive, not a
// replacement.
const BOOTSTRAP_ONLY_MODULES = [
  ...INITIAL_FIT_ONLY_MODULES,
  'sofia2wasm/sofia-wasm.js',
];
if (typeof module === 'undefined' && typeof globalThis !== 'undefined') {
  globalThis.BOOTSTRAP_ONLY_MODULES = BOOTSTRAP_ONLY_MODULES;
}

// Parses --computeGroups=joinKey[,joinSecret][:joinKey[,joinSecret]...] into
// the array of {joinKey, joinSecret} objects job.computeGroups expects --
// ':' between groups, ',' between a group's own joinKey/joinSecret. ':' was
// chosen over ';' specifically because ';' is a shell command separator --
// unquoted, it silently splits `--computeGroups a,b;c,d` into multiple shell
// commands instead of one argument (hit directly: a real dispatch ran with
// only the first group, no --slicePrice/--json, and no error visible in the
// scrollback). ':' isn't a shell metacharacter, so this works unquoted too.
// Falls back to the public group when no flag is given.
function parseComputeGroups(flagValue) {
  return (flagValue || 'public').split(':').map((g) => g.trim()).filter(Boolean).map((group) => {
    const [joinKey, joinSecret] = group.split(',').map((s) => s.trim());
    return joinSecret ? { joinKey, joinSecret } : { joinKey };
  });
}

// isMainThread guard needed here too: `require.main === module` is ALSO
// true inside a worker started via `new Worker(__filename)`, so without it
// every pool worker re-entered the CLI dispatcher on top of its own task
// handling -- observed directly as a run with real converged results mixed
// with a mid-run "--apiKey required" crash and a missing output file.
// Wrapped in an IIFE (cliMain), not just the outer `if`: Node's CommonJS
// loader implicitly wraps every module's whole body in a function
// (function(exports, require, module, __filename, __dirname) {...}), which
// is what makes the bare `return;` statements below legal there -- but a
// browser <script> tag has no such wrapper, so without this IIFE those same
// `return`s are a SyntaxError at PARSE time (not just unreached code) the
// moment this file is loaded via <script src="..."> in galaxy-fit.html,
// even though isNode being false means cliMain() itself is never called.
if (isNode && isMainThread && require.main === module) { (function cliMain() {
  function arg(name, def) {
    const eq = process.argv.find((a) => a.startsWith('--' + name + '='));
    if (eq) return eq.slice(name.length + 3);
    const i = process.argv.indexOf('--' + name);
    return i >= 0 ? process.argv[i + 1] : def;
  }
  const payloadPath = arg('payload');
  const outPath = arg('out', 'dcp_realization_results.json');
  const localFlag = arg('local');
  const isLocal = localFlag !== undefined;
  const isInitialFit = arg('initialFit') !== undefined;

  // -----------------------------------------------------------------------
  // --initial-fit: the ONE anchor fit (runInitialFit), not a bootstrap
  // realization pool -- see runInitialFit's own header for what makes it
  // different. Just one job, so no worker_threads pool for --local (no
  // parallelism to gain from one job); real DCP dispatch below still goes
  // through the exact same compute.for/job.requires/job.exec shape
  // runBootstrapRealization's dispatch uses (just inputSet=[0], workFn=
  // runInitialFit) -- same dual-mode design, ready for a future caller that
  // wants many galaxies' initial fits dispatched over the network in
  // parallel, not a bespoke local-only path.
  // -----------------------------------------------------------------------
  if (isInitialFit) {
    const payload = JSON.parse(fs.readFileSync(payloadPath, 'utf8'));

    if (isLocal) {
      (async () => {
        console.log('[local] running initial fit...');
        const dispatchStart = Date.now();
        const result = await runInitialFit(0, payload);
        const dispatchMs = Date.now() - dispatchStart;
        fs.writeFileSync(outPath, JSON.stringify(result));
        console.log(result.sofiaFailed
          ? `FAILED${result.error ? ' -- ' + result.error : ''}`
          : `converged=${result.converged} chi2=${result.chi2.toFixed(3)} [total=${result.timings.totalMs}ms]`);
        console.log(`Wall time: ${(dispatchMs / 1000).toFixed(1)}s`);
      })().catch((e) => { console.error(e); process.exit(1); });
      return;
    }

    require('dcp-client').init().then(async () => {
      const identity = require('dcp/identity');
      const compute = require('dcp/compute');

      const apiKey = arg('apiKey');
      const computeGroupsFlag = arg('computeGroups');
      const slicePriceFlag = arg('slicePrice');

      if (!apiKey) { console.error('ERROR: --apiKey=0x... is required to launch a DCP job'); process.exit(1); }
      await identity.set(apiKey);

      // compute.for(inputSet, workFn, extraArgs) calls workFn(element,
      // ...extraArgs) per element -- runInitialFit's signature now matches
      // that shape (realizationIndex, payload) explicitly, same as
      // runBootstrapRealization's, rather than relying on payload alone
      // being the inputSet element. (Earlier version of this line passed
      // [0] as inputSet with runInitialFit(payload) -- a single-parameter
      // signature -- which silently bound payload to 0 and left the real
      // payload as an unused second argument; confirmed directly on a real
      // dispatch: "payload.observedCubeRawFitsB64 is required but missing",
      // since (0).observedCubeRawFitsB64 is just undefined, not a throw.)
      const job = compute.for([0], runInitialFit, [payload]);
      // The 3kidnas-test published package, not INITIAL_FIT_AND_BOOTSTRAP_
      // MODULES's local paths -- those can never resolve from a real
      // dispatch at all (Node-only bundling machinery; see
      // INITIAL_FIT_ONLY_MODULES's own comment and
      // ~/.claude/plans/breezy-launching-nova.md).
      job.requires(INITIAL_FIT_ONLY_MODULES);
      job.computeGroups = parseComputeGroups(computeGroupsFlag);
      job.public = {
        name: '🌌 ' + (arg('jobName', '')) + ' (initial fit)',
        description: '3KIDNAS initial/anchor fit',
        link: 'https://www.candiapl.ca/',
      };
      job.on('readystatechange', (ev) => console.log(`Ready state: ${ev}`));
      job.on('accepted', () => console.log(`  Job id: ${job.id}\n  Awaiting result...`));
      job.on('error', (error) => console.error('  Job error:', error));
      job.on('nofunds', (ev) => console.log(ev));

      const dispatchStart = Date.now();
      const slicePrice = slicePriceFlag !== undefined ? parseFloat(slicePriceFlag) : compute.marketValue(1.0);
      const results = await job.exec(slicePrice);
      const dispatchMs = Date.now() - dispatchStart;

      fs.writeFileSync(outPath, JSON.stringify(results[0]));
      console.log(`Dispatch wall time: ${(dispatchMs / 1000).toFixed(1)}s`);
    }).catch((e) => { console.error(e); process.exit(1); });
    return;
  }

  // -----------------------------------------------------------------------
  // --estimateGeometry: a one-off SoFiA run for a PA/Inc starting-point
  // estimate. Local-only, deliberately no real-DCP dispatch branch --
  // measured well under a second on a real test cube (SoFiA source-finding
  // on one un-resampled cube, no fit), so dispatching it over the network
  // would only add DCP's own scheduling/sandbox overhead for zero benefit.
  // run-galaxy-fit-cli.js's estimate-geometry subcommand calls
  // estimateGeometry directly in-process instead of going through this CLI
  // at all; this branch exists for ad-hoc/manual use from a terminal.
  // -----------------------------------------------------------------------
  const isEstimateGeometry = arg('estimateGeometry') !== undefined;
  if (isEstimateGeometry) {
    const payload = JSON.parse(fs.readFileSync(payloadPath, 'utf8'));
    (async () => {
      console.log('estimating geometry...');
      const dispatchStart = Date.now();
      const result = await estimateGeometry(0, payload);
      const dispatchMs = Date.now() - dispatchStart;
      fs.writeFileSync(outPath, JSON.stringify(result));
      console.log(result.sofiaFailed
        ? `FAILED${result.error ? ' -- ' + result.error : ''}`
        : `PA=${result.paEstDeg.toFixed(3)} deg  Inc=${result.incEstDeg.toFixed(3)} deg${result.objectName ? '  object=' + result.objectName : ''}`);
      console.log(`Wall time: ${(dispatchMs / 1000).toFixed(1)}s`);
    })().catch((e) => { console.error(e); process.exit(1); });
    return;
  }

  const nBootstraps = parseInt(arg('nBootstraps'), 10);

  if (isLocal) {
    // Direct in-process execution via a worker_threads pool -- no
    // dcp-client, no job/Supervisor, no network. Each worker just calls the
    // exported work function directly, like a real DCP worker does
    // internally, minus the dispatch/sandboxing. Deliberately NOT
    // dcp-client's job.localExec(): it still calls identity.get() (throws
    // with no identity set, even though it's local-only) and still runs
    // _publishLocalModules, which names the bundled module package from a
    // content hash and requires owning that exact name -- already claimed
    // by another identity in this environment, so every identity hit
    // "Module name conflicts with an existing DCP Module."
    //
    // Parallelism mirrors fortran-local's multiprocessing.Pool (nProcessors
    // == nBootstraps, one OS process per realization) via a worker_threads
    // pool instead, capped at nBootstraps. Each worker is reused across
    // multiple realizations via task messages rather than one Worker per
    // realization -- Worker startup (fresh V8 isolate, re-requiring every
    // dependency, re-instantiating every wasm module) isn't free.
    //
    // Pool size: originally sized to logical-cores-1 (7 on this 8-logical-
    // core M2), which measured at only ~1.8x real speedup, not 7x. NOT
    // whole-chip thermal throttling as first suspected -- directly measured
    // (ps -M per-thread %CPU during a controlled 4-worker run: each thread
    // held ~87-90% CPU simultaneously, ~354% total, genuinely parallel
    // across cores, not serialized onto one). The actual cause: this is a
    // 4-performance + 4-efficiency core chip (hw.perflevel0/1.physicalcpu),
    // and logical-cores-1 (7) oversubscribes the 4 fast P-cores -- at least
    // 3 workers land on much-slower E-cores or time-share a P-core.
    // Confirmed with per-evaluation timing: realizations running under full
    // 7-way contention took 6.6-7.7x longer per objective-function
    // evaluation than the same realization run solo; realizations that
    // started once earlier workers had already freed a P-core ran at only
    // 1.4-1.9x. Sizing to the P-core count (getPerformanceCoreCount above)
    // keeps every worker on a fast core instead.
    (async () => {
      const payload = JSON.parse(fs.readFileSync(payloadPath, 'utf8'));
      const nCores = os.cpus().length;
      const pCores = getPerformanceCoreCount();
      const nWorkers = Math.max(1, Math.min(nBootstraps, pCores || (nCores - 1)));
      console.log(`[local] ${nBootstraps} realization(s) across ${nWorkers} worker thread(s) `
        + `(${nCores} logical cores available${pCores ? `, ${pCores} performance cores` : ''})`);

      const results = new Array(nBootstraps);
      let nextIndex = 0;
      const dispatchStart = Date.now();

      function runPoolWorker() {
        return new Promise((resolveWorker, rejectWorker) => {
          const worker = new Worker(__filename, { workerData: { payload } });
          function assignNext() {
            if (nextIndex >= nBootstraps) {
              worker.terminate().then(() => resolveWorker(), rejectWorker);
              return;
            }
            const i = nextIndex++;
            worker.postMessage({ type: 'task', realizationIndex: i });
          }
          worker.on('message', (msg) => {
            if (!msg || msg.type !== 'result') return;
            if (msg.ok) {
              results[msg.realizationIndex] = msg.result;
              console.log(formatResultLine(msg.result));
            } else {
              results[msg.realizationIndex] = { realizationIndex: msg.realizationIndex, sofiaFailed: true, error: msg.error };
              console.log(formatResultLine(results[msg.realizationIndex]));
            }
            assignNext();
          });
          worker.on('error', (err) => rejectWorker(err));
          assignNext();
        });
      }

      await Promise.all(Array.from({ length: nWorkers }, runPoolWorker));

      const dispatchMs = Date.now() - dispatchStart;
      fs.writeFileSync(outPath, JSON.stringify(results));
      printSummary(results, nBootstraps, dispatchMs, 'Local');
    })().catch((e) => { console.error(e); process.exit(1); });
    return;
  }

  require('dcp-client').init().then(async () => {
    const identity = require('dcp/identity');
    const compute = require('dcp/compute');

    const apiKey = arg('apiKey');
    const computeGroupsFlag = arg('computeGroups');
    const slicePriceFlag = arg('slicePrice');

    if (!apiKey) { console.error('ERROR: --apiKey=0x... is required to launch a DCP job'); process.exit(1); }
    await identity.set(apiKey);

    const inputSet = Array.from({ length: nBootstraps }, (_, i) => i);
    const payload = JSON.parse(fs.readFileSync(payloadPath, 'utf8'));
    const job = compute.for(inputSet, runBootstrapRealization, [payload]);

    // Published-package equivalent of runInitialFit's own dispatch above --
    // see BOOTSTRAP_ONLY_MODULES's own comment for what it adds on top of
    // INITIAL_FIT_ONLY_MODULES (ParseSoFiACatalog.js/GeometryEstimates.js/
    // sofia2wasm, needed for this loop's live per-realization SoFiA run).
    job.requires(BOOTSTRAP_ONLY_MODULES);

    job.computeGroups = parseComputeGroups(computeGroupsFlag);

    job.public = {
      name: '🌌 ' + (arg('jobName', '')),
      description: 'Full 3KIDNAS bootstrap pipeline',
      link: 'https://www.candiapl.ca/',
    };

    job.on('readystatechange', (ev) => console.log(`Ready state: ${ev}`));
    job.on('accepted', () => console.log(`  Job id: ${job.id}\n  Awaiting results...`));
    job.on('error', (error) => console.error('  Job error:', error));
    job.on('nofunds', (ev) => console.log(ev));
    job.on('result', (ev) => console.log(formatResultLine(ev.result)));

    const dispatchStart = Date.now();
    const slicePrice = slicePriceFlag !== undefined ? parseFloat(slicePriceFlag) : compute.marketValue(1.0);
    job.greedyEstimation = true;
    job.estimationSlices = 1000;
    const results = await job.exec(slicePrice);
    const dispatchMs = Date.now() - dispatchStart;

    fs.writeFileSync(outPath, JSON.stringify(results));
    printSummary(results, nBootstraps, dispatchMs, 'Dispatch');
  
  }).catch((e) => { console.error(e); process.exit(1); });
})(); }
