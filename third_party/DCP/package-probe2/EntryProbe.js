/**
 * @file EntryProbe.js -- v2: tests whether ONE LEVEL of subdirectory
 * nesting (matching Phase 0's actual working case, ModuleB.js -> ./sub/
 * ModuleA.js) fixes what a fully-flat 14-entry dependency array (v1,
 * published as 3kidnas-probe2@0.0.1) could not: ALL 14 flat siblings
 * failed uniformly with "Module './DepNN.js' is not available." This
 * version is identical except every dependency now lives one directory
 * down (./deps/DepNN.js instead of ./DepNN.js).
 */
module.declare([
  './deps/Dep01.js', './deps/Dep02.js', './deps/Dep03.js', './deps/Dep04.js', './deps/Dep05.js',
  './deps/Dep06.js', './deps/Dep07.js', './deps/Dep08.js', './deps/Dep09.js', './deps/Dep10.js',
  './deps/Dep11.js', './deps/Dep12.js', './deps/Dep13.js', './deps/Dep14.js',
], function (require, exports, module) {
  exports.tryAll = function () {
    const results = {};
    const deps = [
      './deps/Dep01.js', './deps/Dep02.js', './deps/Dep03.js', './deps/Dep04.js', './deps/Dep05.js',
      './deps/Dep06.js', './deps/Dep07.js', './deps/Dep08.js', './deps/Dep09.js', './deps/Dep10.js',
      './deps/Dep11.js', './deps/Dep12.js', './deps/Dep13.js', './deps/Dep14.js',
    ];
    for (const dep of deps) {
      try {
        results[dep] = { ok: true, val: require(dep).val };
      } catch (e) {
        results[dep] = { ok: false, error: e.message };
      }
    }
    return results;
  };
});
