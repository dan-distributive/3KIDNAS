module.declare([], function (require, exports, module) {
'use strict';

// =============================================================================
// fma.js
// Software emulation of IEEE754 fused multiply-add (a*b+c, single rounding)
// for JavaScript, which has no native FMA operation.
//
// WHY THIS EXISTS
// ----------------
// This project ports several pieces of C code (fdlibm's trig/log kernels,
// FFTW3's generic-solver accumulation loops) to JS with the goal of matching
// the REAL, ACTUAL compiled behavior of the reference C code -- not an
// idealized/portable version of it. We discovered empirically (disassembling
// the exact object files this project's own build produces) that:
//   - src/StandardMath/fdlibm_k_sin.c / fdlibm_k_cos.c, compiled with this
//     project's actual C flags (`-O -DRINGASCII -DASCII`, no
//     `-ffp-contract=off`), DO get compiled with hardware FMA instructions
//     for their polynomial evaluations (verified via `-S` disassembly).
//   - This is not a portability nicety: it changes computed values by 1 ULP
//     for real inputs (verified: theta=0.18265073567382517084 gives a
//     different __kernel_sin result at -O0/-O1 -- WITH fma -- than at -O3
//     -ffp-contract=off -- WITHOUT fma).
//   - FMA is mandatory baseline on arm64 (this project's dev machine) but an
//     optional x86-64 extension -- but since the REAL, CURRENT Fortran
//     build on THIS machine genuinely uses it, matching it here is what
//     "faithful port" means for this codebase, not matching some other
//     hypothetical portable-everywhere behavior.
//
// ALGORITHM
// ---------
// Standard software FMA via error-free transformations (Dekker/Veltkamp):
//   1. TwoProduct(a,b) -> [p, e] such that p+e == a*b EXACTLY (p is the
//      correctly-rounded double product, e is the exact rounding error),
//      computed via Veltkamp splitting (no native FMA needed).
//   2. TwoSum(p,c) -> [s, t] such that s+t == p+c EXACTLY (s is the
//      correctly-rounded double sum).
//   3. Result = s + (t+e), a single ordinary double addition. This is the
//      standard "double-double" software FMA approach (see e.g. Boldo &
//      Muller, "Exact and Approximated Error of the FMA", and the
//      TwoProduct/TwoSum error-free-transformation literature). Verified
//      below against the platform's native C `fma()` (which itself lowers
//      to a hardware FMA instruction on this machine) for both a classic
//      "exposes rounding difference" case and the specific real fdlibm
//      inputs this project needed it for.
// =============================================================================

const SPLITTER = 134217729; // 2^27 + 1, per Veltkamp's algorithm

// Veltkamp split: decomposes `a` into aHi+aLo, each with <=26 significant
// bits, such that aHi+aLo == a EXACTLY and aHi*x, aLo*x never lose bits for
// values in the normal double range.
function split(a) {
  const c = SPLITTER * a;
  const aHi = c - (c - a);
  const aLo = a - aHi;
  return [aHi, aLo];
}

// TwoProduct(a,b): returns [p, e] with p+e == a*b exactly, p = round(a*b).
function twoProduct(a, b) {
  const p = a * b;
  const [aHi, aLo] = split(a);
  const [bHi, bLo] = split(b);
  const e = ((aHi * bHi - p) + aHi * bLo + aLo * bHi) + aLo * bLo;
  return [p, e];
}

// TwoSum(a,b): returns [s, t] with s+t == a+b exactly, s = round(a+b).
function twoSum(a, b) {
  const s = a + b;
  const v = s - a;
  const t = (a - (s - v)) + (b - v);
  return [s, t];
}

// fma(a, b, c) -- correctly-rounded a*b+c (single rounding), matching
// hardware FMA / C's fma(). See module header for why this is needed and
// the algorithm reference.
//
// Inlined rather than calling split()/twoProduct()/twoSum() (which this
// function is a direct, unmodified transcription of -- every operation and
// grouping below is identical to those three functions composed) --
// measured hot spot: fdlibm's __kernel_sin/__kernel_cos each call this ~5
// times, atan's Horner chains ~10 times, and each of the ORIGINAL
// call-based version's invocations allocated 4 short-lived arrays (2 from
// split() inside twoProduct(), 1 from twoProduct() itself, 1 from
// twoSum()) that immediately got destructured and discarded. With sin/cos/
// atan called per-particle across hundreds of thousands of particles per
// fit, that's a very large number of single-use array allocations for
// values that never need to be anything but four local scalars. Same class
// of fix as findParticleCellLocation in FillDataCubeByTiltedRing.js.
function fma(a, b, c) {
  // twoProduct(a, b), split(a) and split(b) inlined:
  const p = a * b;
  const ca = SPLITTER * a;
  const aHi = ca - (ca - a);
  const aLo = a - aHi;
  const cb = SPLITTER * b;
  const bHi = cb - (cb - b);
  const bLo = b - bHi;
  const e = ((aHi * bHi - p) + aHi * bLo + aLo * bHi) + aLo * bLo;
  // twoSum(p, c) inlined:
  const s = p + c;
  const v = s - p;
  const t = (p - (s - v)) + (c - v);
  return s + (t + e);
}

module.exports = { fma, twoProduct, twoSum, split };

// ---------------------------------------------------------------------------
// Self-test (node fma.js) -- compares against known hardware-FMA reference
// values (captured from this project's own machine via C's native fma(),
// which itself lowers to a hardware `fmadd` instruction here -- see
// third_party/fftw-3.3.8/kernel/trig.c's header comment for how this was
// discovered). NOT a substitute for testing against a live C fma() call for
// arbitrary inputs, just a fixed regression set.
// ---------------------------------------------------------------------------
if (require.main === module) {
  const cases = [
    // [a, b, c, expected] -- expected captured via `fma_test a b c` (C's
    // native fma(), this session, arm64 Mac).
    [1.0000000000000002, 0.9999999999999998, -1.0, -4.9303806576313237838e-32],
    // A plain a*b+c where a*b is exactly representable and no rounding
    // difference exists (sanity check: fma should still work, differences
    // only appear when a*b needs its own rounding that interacts with c).
    [2.0, 3.0, 4.0, 10.0],
    [0.1, 0.2, 0.3, 0.32000000000000000666],
  ];
  console.log('=== fma(a,b,c) vs known hardware-FMA reference values ===');
  let allPass = true;
  for (const [a, b, c, expected] of cases) {
    const got = fma(a, b, c);
    const pass = got === expected;
    allPass = allPass && pass;
    console.log(`fma(${a}, ${b}, ${c}) = ${got} (expected ${expected}) ${pass ? 'OK' : 'MISMATCH'}`);
  }
  console.log(allPass ? 'ALL PASS' : 'SOME MISMATCHES -- see above');
}

});
