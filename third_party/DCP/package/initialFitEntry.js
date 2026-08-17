/**
 * @file        initialFitEntry.js (published bundle, hand-written)
 * @description Flat, package-ROOT entry point into 3kidnas-test -- the ONLY
 *              file this package's own job.requires() ever names directly.
 *
 *              Real dispatch revealed TWO separate mechanisms this needs to
 *              work around, not one:
 *              1. job.requires() itself (used at job-deployment time,
 *                 before the sandbox even starts) mangles nested paths --
 *                 "Package 3kidnas-test version 0.1.0 does not contain file
 *                 src,BootstrapSampler,CubeDifference.js" (slashes replaced
 *                 with commas). Fix: job.requires() lists only this ONE
 *                 flat file.
 *              2. Even module.declare()'s own intra-package relative
 *                 require() -- which Phase 0's probe DID prove works for
 *                 ONE level of nesting (ModuleB.js's own require('./sub/
 *                 ModuleA.js')) -- broke on a real dispatch once the target
 *                 was TWO levels deep: "Module './src/ObjectDefinitions/
 *                 DataCube.js' is not available." Fix: every published file
 *                 (this one included) sits FLAT at the package root, zero
 *                 directory levels, and every require() references a
 *                 sibling by bare './<basename>.js' -- the one pattern
 *                 actually proven to work end-to-end.
 *
 * @usage       const entry = require('initialFitEntry.js');
 *              const { DataCube, allocateDataCube } = entry.DataCube;
 */
module.declare([
  './DataCube.js',
  './Beam.js',
  './random.js',
  './DataCubeFits.js',
  './InitialAnalysis.js',
  './TiltedRing.js',
  './ParameterVector.js',
  './CalculateBeamKernel.js',
  './CubeKernelConvolution.js',
  './ParameterToTiltedRingVector.js',
  './TiltedRingModelGeneration.js',
  './GalaxyFit.js',
  './FullModelComparison.js',
  './BasicConstants.js',
  './FlipBootstrap.js',
  './ParseSoFiACatalog.js',
  './GeometryEstimates.js',
  './GetMomentMaps.js',
], function (require, exports, module) {
  exports.DataCube = require('./DataCube.js');
  exports.Beam = require('./Beam.js');
  exports.random = require('./random.js');
  exports.DataCubeFits = require('./DataCubeFits.js');
  exports.InitialAnalysis = require('./InitialAnalysis.js');
  exports.TiltedRing = require('./TiltedRing.js');
  exports.ParameterVector = require('./ParameterVector.js');
  exports.CalculateBeamKernel = require('./CalculateBeamKernel.js');
  exports.CubeKernelConvolution = require('./CubeKernelConvolution.js');
  exports.ParameterToTiltedRingVector = require('./ParameterToTiltedRingVector.js');
  exports.TiltedRingModelGeneration = require('./TiltedRingModelGeneration.js');
  exports.GalaxyFit = require('./GalaxyFit.js');
  exports.FullModelComparison = require('./FullModelComparison.js');
  exports.BasicConstants = require('./BasicConstants.js');
  // Bootstrap-realization-only exports -- runInitialFit doesn't need these
  // (no live SoFiA/resampling for the anchor fit), runBootstrapRealization
  // does (see that function's own require() block for how it's consumed).
  exports.FlipBootstrap = require('./FlipBootstrap.js');
  exports.ParseSoFiACatalog = require('./ParseSoFiACatalog.js');
  exports.GeometryEstimates = require('./GeometryEstimates.js');
  // runInitialFit needs this one directly (to build the model cube's own
  // moment maps for the browser's moment-map visualization) -- unlike the
  // 3 above, this genuinely is a runInitialFit dependency, not just carried
  // along for runBootstrapRealization's sake.
  exports.GetMomentMaps = require('./GetMomentMaps.js');
});
