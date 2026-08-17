/**
 * Bisecting between probe3 (2 deps: WORKS) and probe5 (5 deps: FAILS).
 * Testing 3 flat deps, same shape (synchronous top-level requires).
 */
module.declare([
  './Dep01.js', './Dep02.js', './Dep03.js',
], function (require, exports, module) {
  const d01 = require('./Dep01.js');
  const d02 = require('./Dep02.js');
  const d03 = require('./Dep03.js');
  exports.vals = [d01.val, d02.val, d03.val];
});
