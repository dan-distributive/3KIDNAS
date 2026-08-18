/**
 * @file        build-bravojs-package.js
 * @description Generalized version of cfitsio-4.6.3/wasm/package/
 *              build-bravojs-bundle.js's per-file wrapping, scaled from 2
 *              hand-flattened files to N separately-declared, FLAT files
 *              that cross-require each other by basename -- see
 *              ~/.claude/plans/breezy-launching-nova.md for full context.
 *
 *              Every published file sits directly at the package root, no
 *              subdirectories at all -- NOT preserving the source tree's
 *              directory structure (an earlier version of this script did;
 *              real dispatch proved that wrong). Phase 0's probe only ever
 *              proved ONE level of nesting resolves (sub/ModuleA.js); a
 *              real dispatch of the nested-directory version failed with
 *              "Module './src/ObjectDefinitions/DataCube.js' is not
 *              available" -- two levels deep breaks, independently of
 *              job.requires()'s own separate nested-path bug. Flattening
 *              to zero levels removes the question entirely by only ever
 *              using the one pattern actually proven to work.
 *
 *              For each source file: reads it, extracts its own relative
 *              require('./...')/require('../...') calls, REWRITES each one
 *              to './<basename>' (all 38 basenames are confirmed unique --
 *              see the collision check below), and wraps the rewritten
 *              content in:
 *
 *                module.declare(['./Dep1.js', ...], function (require, exports, module) {
 *                  <rewritten file content>
 *                });
 *
 *              Non-relative requires (e.g. require('cfitsio-wasm.js')) and
 *              relative requires that resolve OUTSIDE this package's own
 *              source root are left completely untouched (see
 *              isInPackageRequire()'s own comment).
 *
 * @usage       node build-bravojs-package.js
 *              (run from js/package/; reads from ../src/,
 *              writes flat .js files + package.dcp directly into ./)
 */
'use strict';

const fs = require('fs');
const path = require('path');

const SRC_ROOT = path.join(__dirname, '..'); // js/
const OUT_ROOT = __dirname;                   // js/package/

// The verified, self-test-guard-aware transitive closure of runInitialFit's
// own dependencies (see the plan file for how this was derived and
// cross-checked against the existing, under-specified
// INITIAL_FIT_AND_BOOTSTRAP_MODULES list). Paths relative to SRC_ROOT.
const FILES = [
  'src/BootstrapSampler/CubeDifference.js',
  'src/BootstrapSampler/DataCubeFits.js',
  'src/BootstrapSampler/FlipBootstrap.js',
  'src/BootstrapSampler/GenerateBootstrap.js',
  'src/BootstrapSampler/PhysCoordTransform.js',
  'src/CompareCubes/CubeComparison.js',
  'src/CompareCubes/FullModelComparison.js',
  'src/CompareCubes/LikelihoodFunctions.js',
  'src/CompareCubes/MaskCube.js',
  'src/ConvolveCube/CalculateBeamKernel.js',
  'src/ConvolveCube/CubeKernelConvolution.js',
  'src/ConvolveCube/FFTW3WasmRank2.js',
  'src/GalaxyAnalysis/GalaxyFit.js',
  'src/GeometryEstimates/GeometryEstimates.js',
  'src/ObjectDefinitions/Beam.js',
  'src/ObjectDefinitions/DataCube.js',
  'src/ObjectDefinitions/ParameterVector.js',
  'src/ObjectDefinitions/Particle.js',
  'src/ObjectDefinitions/TiltedRing.js',
  'src/ParameterToTiltedRingInterface/ParameterToTiltedRingVector.js',
  'src/PreAnalysis/EstimateCubeNoise.js',
  'src/PreAnalysis/EstimateRadialProfiles.js',
  'src/PreAnalysis/EstimateShape.js',
  'src/PreAnalysis/GetMomentMaps.js',
  'src/PreAnalysis/InitialAnalysis.js',
  'src/PreAnalysis/ModellingInitializations.js',
  'src/PreAnalysis/VelProfileAnalysis.js',
  'src/SoFiACatalog/ParseSoFiACatalog.js',
  'src/StandardMath/BasicConstants.js',
  'src/StandardMath/FullCircTrig.js',
  'src/StandardMath/Interpolation.js',
  'src/StandardMath/fdlibm-module.js',
  'src/StandardMath/fdlibm-wasm.js',
  'src/StandardMath/fdlibm.js',
  'src/StandardMath/fma.js',
  'src/StandardMath/random.js',
  'src/TiltedRingModelGeneration/SingleRingGeneration.js',
  'src/TiltedRingModelGeneration/TiltedRingModelGeneration.js',
  'src/TiltedRingToDataCube/FillDataCubeByTiltedRing.js',
  'src/UnitConversions/UnitConversions.js',
];

// Collision check -- flattening only works if every file's basename is
// unique across the whole set.
{
  const seen = new Map();
  for (const relPath of FILES) {
    const base = path.basename(relPath);
    if (seen.has(base)) {
      throw new Error(`Basename collision: ${relPath} and ${seen.get(base)} both flatten to ${base}`);
    }
    seen.set(base, relPath);
  }
}

const REQUIRE_RE = /require\(\s*'([^']+)'\s*\)/g;
const GUARD_RE = /if\s*\(\s*require\.main\s*===\s*module\s*\)/;

// True if `req` (a relative require string found inside the file at
// relPath) resolves to another file that IS part of this package (i.e.
// stays under SRC_ROOT) -- false for e.g. FFTW3WasmRank2.js's
// `require('../../../fftw-3.3.8/wasm/fftw-wasm.js')` local-fallback path
// (in a try/catch after the bare require('fftw-wasm.js') already succeeds
// against the published fftw3wasm-v3 package instead) -- this package has
// no file at that path, so it's left as its original string, unrewritten
// (harmless: that branch never actually executes in a real dispatch).
function isInPackageRequire(relPath, req) {
  const resolvedTarget = path.resolve(path.join(SRC_ROOT, path.dirname(relPath)), req);
  return resolvedTarget.startsWith(SRC_ROOT + path.sep);
}

function flattenRequiresInContent(content, relPath) {
  return content.replace(REQUIRE_RE, (fullMatch, req) => {
    if (!req.startsWith('.')) return fullMatch; // bare/published package -- untouched
    if (!isInPackageRequire(relPath, req)) return fullMatch; // out-of-package -- untouched, see above
    const resolvedTarget = path.resolve(path.join(SRC_ROOT, path.dirname(relPath)), req);
    const flatName = './' + path.basename(resolvedTarget);
    return `require('${flatName}')`;
  });
}

function getFlatDeps(content, relPath) {
  // Only requires BEFORE a require.main===module self-test guard (if any)
  // count -- matches the same convention the dependency tracer used to
  // derive FILES in the first place (see the plan file).
  const guardMatch = GUARD_RE.exec(content);
  const topLevelContent = guardMatch ? content.slice(0, guardMatch.index) : content;
  const deps = [];
  let m;
  REQUIRE_RE.lastIndex = 0;
  while ((m = REQUIRE_RE.exec(topLevelContent))) {
    const req = m[1];
    if (!req.startsWith('.')) continue;
    if (!isInPackageRequire(relPath, req)) {
      console.log(`  (excluding out-of-package dep from ${relPath}: '${req}' -- resolves outside the package root)`);
      continue;
    }
    const resolvedTarget = path.resolve(path.join(SRC_ROOT, path.dirname(relPath)), req);
    deps.push('./' + path.basename(resolvedTarget));
  }
  return [...new Set(deps)];
}

const manifestFiles = {};
let wrapped = 0;

for (const relPath of FILES) {
  const srcPath = path.join(SRC_ROOT, relPath);
  if (!fs.existsSync(srcPath)) {
    throw new Error(`FILES lists ${relPath} but it doesn't exist at ${srcPath}`);
  }
  const content = fs.readFileSync(srcPath, 'utf8');
  const deps = getFlatDeps(content, relPath);
  const rewrittenContent = flattenRequiresInContent(content, relPath);

  const wrappedContent = `module.declare(${JSON.stringify(deps)}, function (require, exports, module) {\n${rewrittenContent}\n});\n`;

  const flatName = path.basename(relPath);
  fs.writeFileSync(path.join(OUT_ROOT, flatName), wrappedContent);
  manifestFiles[flatName] = flatName;
  wrapped++;
  console.log(`wrapped ${relPath} -> ${flatName} (${deps.length} local dep${deps.length === 1 ? '' : 's'})`);
}

manifestFiles['initialFitEntry.js'] = 'initialFitEntry.js';

const packageDcp = {
  name: '3kidnas-test2',
  version: '0.8.0',
  files: manifestFiles,
};
fs.writeFileSync(path.join(OUT_ROOT, 'package.dcp'), JSON.stringify(packageDcp, null, 4) + '\n');

console.log(`\nWrapped ${wrapped}/${FILES.length} files (flat) + 1 hand-written entry point. Wrote package.dcp (version ${packageDcp.version}).`);
