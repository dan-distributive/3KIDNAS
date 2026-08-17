'use strict';

// =============================================================================
// CompositeSolver1D.js
// Faithful JS port of FFTW3's radix-3 Cooley-Tukey decimation-in-time
// combine step for the full-COMPLEX 1D DFT (dft/ct.c's apply_dit + dft/
// scalar/codelets/t1_3.c), for composite sizes n = 3*m where m is an odd
// prime already covered by GenericSolver1D.js's dftGeneric.
//
// Confirmed via fftw_print_plan (verify/ground_truth_harness.c) that for the
// COLUMN (full complex) pass of a 2D r2c/c2r transform, N0=57 (a real,
// production padded size: 2*6+1+44 for one of this pipeline's actual test
// galaxies) decomposes as:
//   (dft-ct-dit/3 (dftw-direct-3/4 "t1_3") (dft-vrank>=1-x3/1 (dft-generic-19)))
// wrapped in a "dft-buffered" cache-optimization layer that is PURELY a
// memory-layout reorganization (dft/buffered.c's apply(): buffer, run the
// SAME cld plan, copy out -- confirmed by reading its apply(), no arithmetic
// difference), so it's correctly ignored here -- the per-column math is
// exactly the unbuffered dft-ct-dit/3 decomposition below.
//
// Unlike RaderSolver.js's size-36 sub-transform (which itself needed a
// radix-6x6 decomposition because 36 has no useful prime factor), the m=19
// "direct" stage here needs NO new codelet at all: 19 is prime, so FFTW3
// (and this port) just calls the already-verified dftGeneric(19) directly.
// Only the twiddle/combine stage (t1_3) is new.
//
// SCOPE: this generalizes to any n = 3*m where m is an odd prime (e.g. a
// future 3*23=69), since the driver below is parametric in m -- not just a
// hardcoded n=57 special case -- but has only been verified against real
// FFTW3 output for n=57 specifically (see verify/run_phase1_tests.js).
//
// FMA note: same as RaderSolver.js -- t1_3.c's FMA/FNMS macros are plain C
// expressions, and this project's FFTW3 build has -ffp-contract=off, so
// plain unfused JS arithmetic is the byte-exact match (no software-FMA
// emulation needed, unlike fdlibm.js).
//
// Source: third_party/fftw-3.3.8/{dft/ct.c, dft/buffered.c, dft/scalar/
// codelets/t1_3.c}, FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const { realCexp } = require('./Trig.js');
const { dftGeneric } = require('./GenericSolver1D.js');

const KP866025403 = 0.866025403784438646763723170752936183471402627;
const KP500000000 = 0.5;

// ---------------------------------------------------------------------------
// dft/scalar/codelets/t1_3.c: t1_3 -- twiddle radix-3 combine, one "column"
// (fixed m, here called `col` to avoid clashing with the outer sub-DFT size
// m) at a time. br/bi are the 3 block values (one per phase p=0,1,2); Wc/Ws
// [p] = cos/sin(2*pi*col*p/n) for p=1,2 (p=0 is the untwiddled DC block, per
// kernel/twiddle.c's TW_FULL loop starting at i=1). Transcribed from the
// #else (non ARCH_PREFERS_FMA) branch -- this project's build compiles that
// branch (config.h: `#undef ARCH_PREFERS_FMA`), confirmed via -S disassembly
// showing zero fmadd/fmsub/fnmadd/fnmsub for these codelets under this
// project's actual -ffp-contract=off build flags.
// ---------------------------------------------------------------------------
function radix3Twiddle(br, bi, Wc, Ws) {
  const T1 = br[0], Ti = bi[0];

  const T3 = br[1], T5 = bi[1];
  const T6 = Wc[1] * T3 + Ws[1] * T5;
  const Te = Wc[1] * T5 - Ws[1] * T3;

  const T8 = br[2], Ta = bi[2];
  const Tb = Wc[2] * T8 + Ws[2] * Ta;
  const Tf = Wc[2] * Ta - Ws[2] * T8;

  const Tc = T6 + Tb;
  const Th = Te + Tf;

  const outR = new Float64Array(3), outI = new Float64Array(3);
  outR[0] = T1 + Tc;
  outI[0] = Th + Ti;

  const Td = T1 - KP500000000 * Tc;
  const Tg = KP866025403 * (Te - Tf);
  outR[2] = Td - Tg;
  outR[1] = Td + Tg;
  const Tj = KP866025403 * (Tb - T6);
  const Tk = Ti - KP500000000 * Th;
  outI[1] = Tj + Tk;
  outI[2] = Tk - Tj;

  return [outR, outI];
}

// ---------------------------------------------------------------------------
// dftRadix3CT -- radix-3 decimation-in-time complex DFT of size n=3*m, m an
// odd prime. Same "physical position == output index" in-place property
// derived in RaderSolver.js's dftSize36 header applies here (r=3 instead of
// r=6): after the twiddle/combine stage, physical position [p*m+col] holds
// X[col+p*m], which equals the physical index itself, so no final
// permutation step is needed.
// ---------------------------------------------------------------------------
function dftRadix3CT(n, reIn, imIn) {
  const m = n / 3;
  const workRe = new Float64Array(n), workIm = new Float64Array(n);

  // Stage 1 (direct, "dft-vrank>=1-x3/1 (dft-generic-<m>)"): decimate by
  // phase p=0,1,2 (element p+k*3 for k=0..m-1), full size-m DFT per phase.
  const riP = new Float64Array(m), iiP = new Float64Array(m);
  const roP = new Float64Array(m), ioP = new Float64Array(m);
  for (let p = 0; p < 3; p++) {
    for (let k = 0; k < m; k++) {
      riP[k] = reIn[p + k * 3];
      iiP[k] = imIn[p + k * 3];
    }
    dftGeneric(m, riP, 0, iiP, 0, 1, roP, 0, ioP, 0, 1);
    for (let k = 0; k < m; k++) {
      workRe[p * m + k] = roP[k];
      workIm[p * m + k] = ioP[k];
    }
  }

  // Stage 2 (twiddle, "dftw-direct-3/4 t1_3"): for each column col=0..m-1,
  // gather the 3 block values, twiddle-multiply by exp(-2*pi*i*col*p/n) for
  // p=1,2, radix-3-combine, write back in place.
  const br = new Float64Array(3), bi = new Float64Array(3);
  const Wc = new Float64Array(3), Ws = new Float64Array(3);
  for (let col = 0; col < m; col++) {
    for (let p = 0; p < 3; p++) {
      br[p] = workRe[p * m + col];
      bi[p] = workIm[p * m + col];
    }
    for (let p = 1; p < 3; p++) {
      const [c, s] = realCexp((col * p) % n, n);
      Wc[p] = c; Ws[p] = s;
    }
    const [outR, outI] = radix3Twiddle(br, bi, Wc, Ws);
    for (let p = 0; p < 3; p++) {
      workRe[p * m + col] = outR[p];
      workIm[p * m + col] = outI[p];
    }
  }

  return [workRe, workIm];
}

module.exports = { dftRadix3CT, radix3Twiddle };

// ---------------------------------------------------------------------------
// Self-test (node CompositeSolver1D.js) -- smoke test only, NOT a substitute
// for the real cross-language verification (see verify/run_phase1_tests.js,
// which exercises this via Rank2Orchestration.js's 57-size fixtures).
// DC-term sanity check: X[0] must equal the plain sum of the input.
// ---------------------------------------------------------------------------
if (require.main === module) {
  console.log('=== CompositeSolver1D self-test (n=57=3x19) ===');
  const n = 57;
  const ri = new Float64Array(n), ii = new Float64Array(n);
  for (let i = 0; i < n; i++) ri[i] = i + 1;
  const [roArr, ioArr] = dftRadix3CT(n, ri, ii);
  const sum = Array.from(ri).reduce((a, b) => a + b, 0);
  console.log(`DC (out[0]): ${roArr[0]} (expect sum = ${sum})`);
}
