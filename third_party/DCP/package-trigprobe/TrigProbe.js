/**
 * Tests whether native Math.sin/Math.cos/Math.atan differ between local
 * Node (Dan's Mac) and a real DCP worker's V8 build, for the same
 * representative angle values used in the tilted-ring pipeline. Also
 * cross-checks the wasm-backed fdlibm sin/cos/atanh for the same inputs --
 * those SHOULD already be bit-identical cross-machine (Wasm float ops are
 * IEEE-754-strict by spec), so this doubles as a control.
 *
 * Returns exact bit patterns (via Float64Array/DataView) rather than
 * decimal strings, since a decimal round-trip could mask a genuine last-bit
 * difference.
 */
module.declare(["./fdlibm-wasm.js"], function (require, exports, module) {

function toHexBits(x) {
  const buf = new ArrayBuffer(8);
  new Float64Array(buf)[0] = x;
  const view = new DataView(buf);
  let hex = '';
  for (let i = 0; i < 8; i++) {
    hex += view.getUint8(i).toString(16).padStart(2, '0');
  }
  return hex;
}

// Representative angles: some "nice" pipeline-realistic values (ring
// position angles / inclinations in radians, roughly 0..2*Pi and 0..Pi/2),
// plus a few edge-case-y values (near Pi/2 multiples, small, negative,
// moderately large) since transcendental-function implementations tend to
// diverge most in large-argument range reduction.
const TEST_ANGLES = [
  0.0, 0.1, 0.5, 1.0,
  Math.PI / 6, Math.PI / 4, Math.PI / 3, Math.PI / 2,
  Math.PI, 3 * Math.PI / 2, 2 * Math.PI,
  -1.234567, 4.276057, 12.9834,
  1.5707963267948966, // Pi/2 exactly, as a double
  100.123456789,
];

exports.runTrigProbe = async function runTrigProbe() {
  const fdlibm = require('./fdlibm-wasm.js');
  await fdlibm.warmUp();

  const rows = TEST_ANGLES.map((x) => ({
    x,
    xBits: toHexBits(x),
    nativeSin: Math.sin(x),
    nativeSinBits: toHexBits(Math.sin(x)),
    nativeCos: Math.cos(x),
    nativeCosBits: toHexBits(Math.cos(x)),
    wasmSin: fdlibm.fdSinWasm(x),
    wasmSinBits: toHexBits(fdlibm.fdSinWasm(x)),
    wasmCos: fdlibm.fdCosWasm(x),
    wasmCosBits: toHexBits(fdlibm.fdCosWasm(x)),
  }));

  // atanh has domain (-1, 1) -- separate, in-domain test values.
  const ATANH_TEST_VALUES = [0.0, 0.1, -0.1, 0.5, -0.5, 0.9, -0.9, 0.999];
  const atanhRows = ATANH_TEST_VALUES.map((x) => ({
    x,
    xBits: toHexBits(x),
    nativeAtanh: Math.atanh(x),
    nativeAtanhBits: toHexBits(Math.atanh(x)),
    wasmAtanh: fdlibm.fdAtanhWasm(x),
    wasmAtanhBits: toHexBits(fdlibm.fdAtanhWasm(x)),
  }));

  return { rows, atanhRows };
};

});
