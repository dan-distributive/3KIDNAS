/**
 * @file        sofia-wasm.js
 * @description Single entry point for running SoFiA-2 (the HI source
 *              finding pipeline, ../sofia.c + ../src/*.c) inside a DCP
 *              worker. job.requires(['./sofia-wasm']) is enough -- this
 *              file require()s sofia-module.js, the Emscripten build of
 *              SoFiA-2 statically linked against a wasm build of wcslib
 *              (see build.sh), with the wasm binary embedded directly in
 *              the JS (SINGLE_FILE=1) so there's no separate .wasm fetch
 *              for the sandbox to reject.
 *
 *              SoFiA-2 is a CLI tool: it takes one argument (a parameter
 *              file) and does the rest via file I/O -- reading the input
 *              cube, writing catalogues/masks/moment maps to an output
 *              directory. run() below does the sandbox-specific plumbing
 *              that implies: writing the cube and parameter file into
 *              Emscripten's in-memory filesystem (MEMFS), invoking main()
 *              via callMain() exactly as if from a shell, and reading the
 *              output directory back out into an in-memory result object.
 *              None of this touches a real filesystem or network; it's
 *              the same buffer-in/buffer-out shape as duckdb-wasm.js's
 *              registerFileBuffer(), just for a full CLI pipeline instead
 *              of a SQL engine.
 *
 * @usage       const sofia = require('./sofia-wasm');
 *              const { exitCode, log, files } = await sofia.run({
 *                cube: cubeBytes,           // Uint8Array, a FITS cube
 *                par: 'scfind.enable = true\n...',  // parameter file text
 *              });
 *              // files is a Map of output filename -> Uint8Array
 *              // (test_cat.txt, test_mask.fits, test_mom0.fits, etc.,
 *              // named per whatever output.filename the par text sets).
 */
'use strict';

// require() exists in Node and inside a DCP sandbox (via job.requires()'s
// own module bundle) but not in a raw browser tab -- sofia-module.js itself
// is already Emscripten-UMD (a bare `var SofiaModule = ...` at module
// scope, so loading it via a plain <script> tag already puts it on
// globalThis for free; only THIS hand-written wrapper needs the guard.
const createSofiaModule = (typeof require === 'function')
  ? require('./sofia-module')
  : globalThis.SofiaModule;

/**
 * Run one SoFiA-2 invocation to completion. Each call gets a fresh
 * Emscripten instance (fresh MEMFS, fresh heap), so concurrent or
 * repeated calls in the same worker never see each other's files --
 * the same isolation you'd get from separate worker processes, just
 * without the overhead of actually spawning them.
 *
 * @param {Object} opts
 * @param {Uint8Array} opts.cube - FITS data cube bytes.
 * @param {string} opts.par - Parameter file text. Must set
 *   `input.data = /work/cube.fits` (the fixed path this function writes
 *   the cube to) and should set `output.directory = /work/out` (created
 *   for you) plus `output.filename` and whichever `output.write*` flags
 *   you want products for.
 * @param {Object<string,Uint8Array|string>} [opts.extraFiles] - Extra
 *   input files to preload (weights/mask/noise cubes, catalogues for
 *   flagging, etc.), keyed by the absolute MEMFS path your par text
 *   references, e.g. `{'/work/weights.fits': weightsBytes}`.
 * @returns {Promise<{exitCode: number, log: string, files: Map<string, Uint8Array>}>}
 */
async function run({ cube, par, extraFiles = {} }) {
  const log = [];
  const Module = await createSofiaModule({
    print: (line) => log.push(line),
    printErr: (line) => log.push(line),
  });

  Module.FS.mkdir('/work');
  Module.FS.mkdir('/work/out');
  Module.FS.writeFile('/work/cube.fits', cube);
  Module.FS.writeFile('/work/sofia.par', par);
  for (const [path, contents] of Object.entries(extraFiles)) {
    Module.FS.writeFile(path, contents);
  }

  let exitCode = 0;
  try {
    exitCode = Module.callMain(['/work/sofia.par']);
  } catch (e) {
    // Emscripten's callMain throws an ExitStatus (not a real error) for
    // any exit() call, including SoFiA's own ensure()-triggered failures
    // (see ../src/common.h's ERR_* codes) -- surface the code, don't throw.
    if (e && e.name === 'ExitStatus') {
      exitCode = e.status;
    } else {
      throw e;
    }
  }

  const files = new Map();
  for (const name of Module.FS.readdir('/work/out')) {
    if (name === '.' || name === '..') continue;
    files.set(name, Module.FS.readFile('/work/out/' + name));
  }

  return { exitCode, log: log.join('\n'), files };
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { run };
} else if (typeof globalThis !== 'undefined') {
  globalThis.SofiaWasm = { run };
}
