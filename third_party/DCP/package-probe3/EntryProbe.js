/**
 * Tests the tightest possible increment from Phase 0's proven case: TWO
 * flat siblings instead of one, both required directly at the factory's
 * own top level (synchronously, during initial module load -- not
 * deferred into a nested function, matching Phase 0's ModuleB.js exactly,
 * unlike package-probe2's tryAll()-deferred version).
 */
module.declare(['./DepA.js', './DepB.js'], function (require, exports, module) {
  const a = require('./DepA.js');
  const b = require('./DepB.js');
  exports.aVal = a.val;
  exports.bVal = b.val;
});
