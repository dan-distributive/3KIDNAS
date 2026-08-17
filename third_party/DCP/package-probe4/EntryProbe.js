/**
 * Bisecting the boundary found between package-probe3 (2 deps: WORKS)
 * and package-probe2 (14 deps: FAILS uniformly). This tests 8 flat
 * deps, all required synchronously at the factory's own top level
 * (same shape as probe3's proven-working 2-dep case).
 */
module.declare([
  './Dep01.js', './Dep02.js', './Dep03.js', './Dep04.js',
  './Dep05.js', './Dep06.js', './Dep07.js', './Dep08.js',
], function (require, exports, module) {
  const d01 = require('./Dep01.js');
  const d02 = require('./Dep02.js');
  const d03 = require('./Dep03.js');
  const d04 = require('./Dep04.js');
  const d05 = require('./Dep05.js');
  const d06 = require('./Dep06.js');
  const d07 = require('./Dep07.js');
  const d08 = require('./Dep08.js');
  exports.vals = [d01.val, d02.val, d03.val, d04.val, d05.val, d06.val, d07.val, d08.val];
});
