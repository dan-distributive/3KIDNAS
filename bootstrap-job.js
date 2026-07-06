/**
 * @file        bootstrap-job.js
 * @description DCP job driver for 3KIDNAS bootstrap fitting.
 *              Generates N bootstrap realizations of a galaxy kinematic fit
 *              in parallel across a DCP worker pool. Each worker sandbox receives
 *              one realization index and runs a full GalaxyFit_Simple on the
 *              corresponding bootstrap cube.
 *
 * @authors     Dan Desjardins <dan@distributive.network>
 * @date        July 2026
 * @copyright   2026 Distributive Corp.
 *
 * @usage       node bootstrap-job.js --apiKey=0x<identity> [options]
 *
 *   Flags:
 *     --apiKey=0x...              Identity / API key to run the job under (required)
 *     --computeGroup=key,secret   Compute group as joinKey,joinSecret
 *                                 (secret optional; e.g. --computeGroup=demo,dcp
 *                                 or --computeGroup=public). Defaults to public group.
 *     --bootstrapts=N             Number of bootstrap realizations to perform (default 1000)
 *     --slicePrice=N              Per-slice price in DCC (default 5.24)
 *
 *   Examples:
 *     node bootstrap-job.js --apiKey=0x45d7... --computeGroup=demo,dcp --bootstraps=50 --slicePrice=5.24
 *     node bootstrap-job.js --apiKey=0x45d7...
 *     node bootstrap-job.js --apiKey=0x45d7... --computeGroup=public
 *
 * @requires    node.js
 * @requires    dcp-client
 */

const fs = require('fs');
const path = require('path');

async function main() {
  const identity = require('dcp/identity');
  const compute  = require('dcp/compute');
  const RangeObject = require('dcp/range-object').RangeObject

  // ---- CLI FLAGS ----
  // --apiKey=0x...              identity / API key (required)
  // --computeGroup=key,secret   join key + secret (secret optional)
  // --bootstraps=1000           number of bootstrap realizations to compute (default 1000)
  // --slicePrice=5.24           per-slice price (default 5.24)
  function getFlag(name) {
    const pfx = `--${name}=`;
    const hit = process.argv.find(a => a.startsWith(pfx));
    return hit ? hit.slice(pfx.length) : undefined;
  }

  const apiKey     = getFlag('apiKey');
  const cg         = getFlag('computeGroup');
  const bootstraps = parseInt(getFlag('bootstraps') ?? '1000', 10);
  const slicePrice = parseFloat(getFlag('slicePrice') ?? '5.24');

  if (!apiKey) {
    console.error('ERROR: --apiKey=0x... is required');
    process.exit(1);
  }
  if (!Number.isInteger(bootstraps) || bootstraps < 1) {
    console.error('ERROR: --bootstraps must be a positive integer');
    process.exit(1);
  }
  if (Number.isNaN(slicePrice) || slicePrice <= 0) {
    console.error('ERROR: --slicePrice must be a positive number');
    process.exit(1);
  }

  // Parse computeGroup as "joinKey,joinSecret" (secret optional)
  let computeGroup;
  if (cg) {
    const [joinKey, joinSecret] = cg.split(',');
    if (!joinKey) {
      console.error('ERROR: --computeGroup must be joinKey or joinKey,joinSecret');
      process.exit(1);
    }
    computeGroup = joinSecret ? { joinKey, joinSecret } : { joinKey };
  }

  // API KEY
  await identity.set(apiKey);

  // STATIC INPUTS — load fixture and best-fit model cube from disk
  // These are passed as arguments to every worker so they don't need
  // to be fetched from a URL. Each is a plain JSON string.
  const fixtureStr   = fs.readFileSync(path.join(__dirname, 'inputs/diskfit_fixture.json'),    'utf8');
  const modelCubeStr = fs.readFileSync(path.join(__dirname, 'inputs/model_cube_bestfit.json'), 'utf8');

  // INPUT SET — one integer per bootstrap realization (used as RNG seed offset)
  const realizations = new RangeObject({ start: 1, end: bootstraps });

  // WORK FUNCTION — runs in each DCP worker sandbox
  // realizationIndex: integer 1..N (unique per worker sandbox)
  // fixtureStr:       JSON string of diskfit_fixture.json
  // modelCubeStr:     JSON string of model_cube_bestfit.json
  async function bootstrap(realizationIndex, fixtureStr, modelCubeStr) {
    progress(0.0);

    const f32 = Math.fround;

    // Load all required modules (must be listed in job.requires below)
    const { makeRng }            = require('./src/StandardMath/random');
    const { DataCube, allocateDataCube, flatIndxCalc }
                                 = require('./src/ObjectDefinitions/DataCube');
    const { Beam2D }             = require('./src/ObjectDefinitions/Beam');
    const { TiltedRingModel, TiltedRingFittingOptions,
            tiltRing_Allocate, tiltRingFittingOptions_Allocate,
            logicalTiltedRingIndexing }
                                 = require('./src/ObjectDefinitions/TiltedRing');
    const { ParameterVector, allocateParamVector }
                                 = require('./src/ObjectDefinitions/ParameterVector');
    const { calculate2DBeamKernel }
                                 = require('./src/ConvolveCube/CalculateBeamKernel');
    const { generalizedParamVectorToTiltedRing }
                                 = require('./src/ParameterToTiltedRingInterface/ParameterToTiltedRingVector');
    const { galaxyFit_Simple }   = require('./src/GalaxyAnalysis/GalaxyFit');
    const { genBootstrapSample } = require('./src/BootstrapSampler/GenerateBootstrap');

    // ---- Parse fixture ----
    const d = JSON.parse(fixtureStr);

    const pvIni = new ParameterVector();
    pvIni.nParams = d.nParams;
    allocateParamVector(pvIni);
    for (let i = 0; i < d.nParams; i++) {
      pvIni.param[i]          = f32(d.param[i]);
      pvIni.paramLowerLims[i] = f32(d.paramLowerLims[i]);
      pvIni.paramUpperLims[i] = f32(d.paramUpperLims[i]);
      pvIni.paramRange[i]     = f32(d.paramRange[i]);
      pvIni.cyclicSwitch[i]   = d.cyclicSwitch[i];
    }

    const trfo = new TiltedRingFittingOptions();
    trfo.nRings = d.nRings;
    trfo.nRingsPerBeam = 1;
    trfo.nTargRings = -1;
    for (let i = 0; i <= 12; i++) {
      trfo.constParams[i]    = d.constParams[i];
      trfo.fixedParams[i]    = d.fixedParams[i];
      trfo.paramLowerLims[i] = f32(d.paramLowerLims13[i]);
      trfo.paramUpperLims[i] = f32(d.paramUpperLims13[i]);
      trfo.paramRange[i]     = f32(d.paramRange13[i]);
      trfo.cyclicSwitch[i]   = d.cyclicSwitch13[i];
    }
    tiltRingFittingOptions_Allocate(trfo);
    logicalTiltedRingIndexing(trfo);
    for (let r = 0; r < d.nRings; r++) {
      const rp = trfo.radialProfiles[r];
      const dp = d.radialProfiles[r];
      rp.rmid = f32(dp.rmid); rp.rwidth = f32(dp.rwidth);
      rp.centPos[0] = f32(dp.centPos0); rp.centPos[1] = f32(dp.centPos1);
      rp.inclination = f32(dp.inclination); rp.positionAngle = f32(dp.positionAngle);
      rp.vSys = f32(dp.vSys); rp.vRot = f32(dp.vRot); rp.vRad = f32(dp.vRad);
      rp.vDisp = f32(dp.vDisp); rp.vvert = f32(dp.vvert); rp.dvdz = f32(dp.dvdz);
      rp.sigUse = f32(dp.sigUse); rp.z0 = f32(dp.z0);
      rp.zGradiantStart = f32(dp.zGradiantStart);
    }

    // Observed cube
    const obs = new DataCube();
    const dh  = obs.dh;
    dh.nPixels[0] = d.observedDC.nPixelsX;
    dh.nPixels[1] = d.observedDC.nPixelsY;
    dh.nChannels  = d.observedDC.nChannels;
    dh.pixelSize[0] = f32(d.observedDC.pixelSizeX);
    dh.pixelSize[1] = f32(d.observedDC.pixelSizeY);
    dh.channelSize  = f32(d.observedDC.channelSize);
    dh.refLocation[0] = f32(d.observedDC.refLocX);
    dh.refLocation[1] = f32(d.observedDC.refLocY);
    dh.refLocation[2] = f32(d.observedDC.refLocV);
    dh.refVal[0] = f32(d.observedDC.refValX);
    dh.refVal[1] = f32(d.observedDC.refValY);
    dh.refVal[2] = f32(d.observedDC.refValV);
    dh.uncertainty = f32(d.observedDC.uncertainty);
    dh.nValid = d.observedDC.nValid;
    const nCells = dh.nPixels[0] * dh.nPixels[1] * dh.nChannels;
    obs.pixels   = new Float32Array(2 * Math.max(dh.nPixels[0], dh.nPixels[1]));
    obs.channels = new Float32Array(dh.nChannels);
    obs.flux     = new Float32Array(nCells);
    obs.flattendValidIndices = new Int32Array(dh.nValid);
    for (let l = 0; l < dh.nValid; l++)
      obs.flattendValidIndices[l] = d.observedDC.flattendValidIndices[l];
    for (let l = 0; l < dh.nValid; l++)
      obs.flux[obs.flattendValidIndices[l]] = f32(d.observedDC.validFlux[l]);
    obs.dh.start[0] = f32(d.observedDC.startX);
    obs.dh.start[1] = f32(d.observedDC.startY);
    obs.dh.start[2] = f32(d.observedDC.startV);

    // Beam
    const beam = new Beam2D();
    beam.beamSigmaVector[0] = f32(d.observedBeam.beamSigma0);
    beam.beamSigmaVector[1] = f32(d.observedBeam.beamSigma1);
    beam.beamSigmaVector[2] = f32(d.observedBeam.beamSigma2);
    beam.nRadialCells  = d.observedBeam.nRadialCells;
    beam.sigmaLengths  = f32(d.observedBeam.sigmaLengths);
    beam.pixelSize[0]  = f32(d.observedBeam.pixelSizeX);
    beam.pixelSize[1]  = f32(d.observedBeam.pixelSizeY);
    beam.beamMajorAxis = f32(d.observedBeam.beamMajorAxis);
    beam.beamMinorAxis = f32(d.observedBeam.beamMinorAxis);
    const n = beam.nRadialCells;
    const kSz = 2 * n + 1;
    beam.kernel = new Float32Array(kSz * kSz);
    beam.paddedSize[0] = 2*n+1 + dh.nPixels[0];
    beam.paddedSize[1] = 2*n+1 + dh.nPixels[1];
    beam.complexSize[0] = Math.trunc(beam.paddedSize[0]/2)+1;
    beam.complexSize[1] = beam.paddedSize[1];
    beam.complexKernelCreated = false;
    calculate2DBeamKernel(beam, beam.pixelSize);

    // ---- Parse best-fit model cube ----
    const m = JSON.parse(modelCubeStr);
    const modelCube = new DataCube();
    modelCube.dh.nPixels[0] = m.nPixelsX;
    modelCube.dh.nPixels[1] = m.nPixelsY;
    modelCube.dh.nChannels  = m.nChannels;
    modelCube.dh.pixelSize[0]  = dh.pixelSize[0];
    modelCube.dh.pixelSize[1]  = dh.pixelSize[1];
    modelCube.dh.channelSize   = dh.channelSize;
    modelCube.dh.refLocation[0]= dh.refLocation[0];
    modelCube.dh.refLocation[1]= dh.refLocation[1];
    modelCube.dh.refLocation[2]= dh.refLocation[2];
    modelCube.dh.refVal[0] = dh.refVal[0];
    modelCube.dh.refVal[1] = dh.refVal[1];
    modelCube.dh.refVal[2] = dh.refVal[2];
    modelCube.dh.uncertainty = dh.uncertainty;
    modelCube.dh.nValid = dh.nValid;
    allocateDataCube(modelCube);
    modelCube.dh.start[0] = obs.dh.start[0];
    modelCube.dh.start[1] = obs.dh.start[1];
    modelCube.dh.start[2] = obs.dh.start[2];
    modelCube.flattendValidIndices.set(obs.flattendValidIndices);
    let idx = 0;
    for (let i = 0; i < m.nPixelsX; i++)
      for (let j = 0; j < m.nPixelsY; j++)
        for (let k = 0; k < m.nChannels; k++)
          modelCube.flux[flatIndxCalc(i, j, k, modelCube.dh)] = f32(m.flux[idx++]);

    // ---- Bootstrap geometry from best-fit radial profiles ----
    const bsCent = {
      centX: f32(trfo.radialProfiles[0].centPos[0]),
      centY: f32(trfo.radialProfiles[0].centPos[1]),
      centV: f32(trfo.radialProfiles[0].vSys),
      pa:    f32(trfo.radialProfiles[0].positionAngle),
      inc:   f32(trfo.radialProfiles[0].inclination),
    };

    // ---- Generate bootstrap realization ----
    const bootstrapCube = genBootstrapSample(
      obs, modelCube, beam, bsCent,
      /* spatialBlockSize */ Math.fround(1.0),
      /* velBlockSize     */ Math.fround(1.0)
    );

    // ---- Build fit state — unique seed per realization ----
    const bsIdum = d.idum - realizationIndex;

    const pvModel    = new ParameterVector();
    const pvFirstFit = new ParameterVector();
    pvModel.nParams  = pvIni.nParams;
    allocateParamVector(pvModel);
    pvModel.param.set(pvIni.param);
    pvModel.paramLowerLims.set(pvIni.paramLowerLims);
    pvModel.paramUpperLims.set(pvIni.paramUpperLims);
    pvModel.cyclicSwitch.set(pvIni.cyclicSwitch);
    pvModel.paramRange.set(pvIni.paramRange);

    const modelTR  = new TiltedRingModel();
    modelTR.nRings = trfo.nRings;
    modelTR.cmode  = 0;
    modelTR.cloudBaseSurfDens = Math.fround(100.0);
    tiltRing_Allocate(modelTR);
    for (let r = 0; r < trfo.nRings; r++) {
      modelTR.r[r].rmid   = f32(trfo.radialProfiles[r].rmid);
      modelTR.r[r].rwidth = f32(trfo.radialProfiles[r].rwidth);
    }

    const modelDC = new DataCube();
    modelDC.dh.nPixels[0]    = bootstrapCube.dh.nPixels[0];
    modelDC.dh.nPixels[1]    = bootstrapCube.dh.nPixels[1];
    modelDC.dh.nChannels     = bootstrapCube.dh.nChannels;
    modelDC.dh.pixelSize[0]  = bootstrapCube.dh.pixelSize[0];
    modelDC.dh.pixelSize[1]  = bootstrapCube.dh.pixelSize[1];
    modelDC.dh.channelSize   = bootstrapCube.dh.channelSize;
    modelDC.dh.refLocation[0]= bootstrapCube.dh.refLocation[0];
    modelDC.dh.refLocation[1]= bootstrapCube.dh.refLocation[1];
    modelDC.dh.refLocation[2]= bootstrapCube.dh.refLocation[2];
    modelDC.dh.refVal[0]     = bootstrapCube.dh.refVal[0];
    modelDC.dh.refVal[1]     = bootstrapCube.dh.refVal[1];
    modelDC.dh.refVal[2]     = bootstrapCube.dh.refVal[2];
    modelDC.dh.uncertainty   = bootstrapCube.dh.uncertainty;
    modelDC.dh.nValid        = bootstrapCube.dh.nValid;
    allocateDataCube(modelDC);
    modelDC.dh.start[0] = bootstrapCube.dh.start[0];
    modelDC.dh.start[1] = bootstrapCube.dh.start[1];
    modelDC.dh.start[2] = bootstrapCube.dh.start[2];
    modelDC.flattendValidIndices.set(bootstrapCube.flattendValidIndices);

    const state = {
      pvIni, pvModel, pvFirstFit,
      modelTiltedRing:  modelTR,
      modelDC,
      observedDC:       bootstrapCube,
      observedBeam:     beam,
      trFittingOptions: trfo,
      rng:              makeRng(bsIdum),
      linearLogSDSwitch: d.linearLogSDSwitch,
      ftol:             Math.fround(d.ftol),
      iniGuessWidth:    Math.fround(1.0),
      paramToTiltedRing: generalizedParamVectorToTiltedRing
    };

    // ---- Run optimizer ----
    const pvBest = galaxyFit_Simple(state);

    progress(1.0);

    return {
      realizationIndex,
      params: Array.from(pvBest.param),
      chi2:   pvBest.bestLike
    };
  }


  // JOB
  const job = compute.for(realizations, bootstrap, [fixtureStr, modelCubeStr]);


  // JOB REQUIRED PACKAGES — only top-level entry points needed
  job.requires([
    './src/StandardMath/random',
    './src/StandardMath/Interpolation',
    './src/ObjectDefinitions/DataCube',
    './src/ObjectDefinitions/Beam',
    './src/ObjectDefinitions/TiltedRing',
    './src/ObjectDefinitions/ParameterVector',
    './src/ObjectDefinitions/Particle',
    './src/ConvolveCube/CalculateBeamKernel',
    './src/ConvolveCube/CubeKernelConvolution',
    './src/TiltedModelGeneration/TiltedRingModelGeneration',
    './src/TiltedModelGeneration/SingleRingGeneration',
    './src/TiltedRingToDataCube/FillDataCubeByTiltedRing',
    './src/CompareCubes/LikelihoodFunctions',
    './src/CompareCubes/CubeComparison',
    './src/CompareCubes/FullModelComparison',
    './src/ParameterToTiltedRingInterface/ParameterToTiltedRingVector',
    './src/GalaxyAnalysis/GalaxyFit',
    './src/BootstrapSampler/CubeDifference',
    './src/BootstrapSampler/PhysCoordTransform',
    './src/BootstrapSampler/GenerateBootstrap',
  ]);


  // JOB COMPUTE GROUPS
  job.computeGroups = computeGroup
  ? [computeGroup]
  : [{ joinKey: 'public' }];  // default


  // JOB PUBLIC INFO
  job.public = {
    name:        '🌌 3KIDNAS Bootstrap',
    description: 'Galaxy kinematic bootstrap fitting — WALLABY HI survey',
    link:        'https://www.candiapl.ca/'
  };


  // JOB EVENTS
  job.on('readystatechange', (ev) => console.log(`Ready state: ${ev}`));
  job.on('accepted', ()   => console.log(`  Job id: ${job.id}\n  Awaiting results...`));
  job.on('error', (error) => console.error('  Job error:', error));
  job.on('nofunds', (ev)  => console.log(ev));
  job.on('result', (ev) => {
    if (typeof ev?.result?.chi2 !== 'number') return;
    console.log(`  Realization ${ev.result.realizationIndex}: chi²=${ev.result.chi2.toExponential(4)}`);
  });


  // JOB DEPLOYMENT
  const results = await job.exec(slicePrice);


  // COLLATE RESULTS
  const resultArr = Array.from(results).sort((a,b) => a.realizationIndex - b.realizationIndex);
  const nParams   = resultArr[0].params.length;

  console.log('\n=== Bootstrap Results ===');
  console.log(`N realizations: ${resultArr.length}`);
  console.log(`\n${'Param'.padEnd(6)} ${'Mean'.padEnd(12)} ${'Std'.padEnd(12)} ${'p16'.padEnd(12)} ${'p50'.padEnd(12)} ${'p84'.padEnd(12)}`);
  console.log('-'.repeat(66));

  const stats = [];
  for (let p = 0; p < nParams; p++) {
    const vals = resultArr.map(r => r.params[p]).sort((a,b) => a-b);
    const mean = vals.reduce((a,v) => a+v, 0) / vals.length;
    const std  = Math.sqrt(vals.reduce((a,v) => a+(v-mean)**2, 0) / vals.length);
    const p16  = vals[Math.floor(0.16 * vals.length)];
    const p50  = vals[Math.floor(0.50 * vals.length)];
    const p84  = vals[Math.floor(0.84 * vals.length)];
    stats.push({ mean, std, p16, p50, p84 });
    console.log(`[${String(p).padEnd(3)}]  ${mean.toFixed(4).padEnd(12)} ${std.toFixed(4).padEnd(12)} ${p16.toFixed(4).padEnd(12)} ${p50.toFixed(4).padEnd(12)} ${p84.toFixed(4).padEnd(12)}`);
  }

  const chi2vals = resultArr.map(r => r.chi2).sort((a,b) => a-b);
  const chi2mean = chi2vals.reduce((a,v) => a+v, 0) / chi2vals.length;
  console.log(`\nchi² mean: ${chi2mean.toExponential(4)}`);
  console.log(`chi² range: [${chi2vals[0].toExponential(4)}, ${chi2vals[chi2vals.length-1].toExponential(4)}]`);

  // Save results — timestamped so runs don't overwrite each other
  const outDir = path.join(__dirname, 'outputs');
  fs.mkdirSync(outDir, { recursive: true });   // no-op if it already exists

  // Local time, filename-safe: YYYY-MM-DD_HH-MM-SS
  const now = new Date();
  const pad = (n) => String(n).padStart(2, '0');
  const stamp = `${now.getFullYear()}-${pad(now.getMonth() + 1)}-${pad(now.getDate())}`
              + `_${pad(now.getHours())}-${pad(now.getMinutes())}-${pad(now.getSeconds())}`;

  const outPath = path.join(outDir, `bootstrap_results_${stamp}.json`);
  fs.writeFileSync(outPath, JSON.stringify({
    nRealizations: resultArr.length,
    nParams,
    results: resultArr,
    stats
  }, null, 2));
  console.log(`\nResults saved to ${outPath}`);
}

require('dcp-client')
  .init()
  .then(main)
  .catch(console.error);
