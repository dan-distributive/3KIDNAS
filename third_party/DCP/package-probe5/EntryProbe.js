/**
 * Bisecting between probe3 (2 deps: WORKS) and probe4 (8 deps: FAILS).
 * Testing 5 flat deps, same shape (synchronous top-level requires).
 */
module.declare([
  './Dep01.js', './Dep02.js', './Dep03.js', './Dep04.js', './Dep05.js',
], function (require, exports, module) {
  const d01 = require('./Dep01.js');
  const d02 = require('./Dep02.js');
  const d03 = require('./Dep03.js');
  const d04 = require('./Dep04.js');
  const d05 = require('./Dep05.js');
  exports.vals = [d01.val, d02.val, d03.val, d04.val, d05.val];
});
