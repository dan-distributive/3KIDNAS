'use strict';

// =============================================================================
// GenericSolver1D.js
// Faithful JS port of FFTW3's "generic" 1D solvers -- the ones FFTW3's
// planner actually selects (confirmed via fftw_print_plan, see
// verify/fixtures/plan_*.txt) for ODD PRIME transform sizes with
// FFTW_ESTIMATE, which is exactly what this pipeline's padded array sizes
// frequently are (e.g. 47, 43, 37 -- all arise from real galaxy cube
// dimensions + beam kernel radius). NOT a naive O(n^2) DFT: a "Hartley
// fold" of the input (sum/difference pairs indexed by i and n-i) followed
// by a twiddle-weighted dot product, using precomputed twiddle factors
// (see Trig.js's realCexp).
//
// Ports, one per source function, matching upstream exactly:
//   dft/generic.c:  hartley, cdot            -- full complex DFT (both
//                                                 directions; direction is
//                                                 selected by the caller
//                                                 swapping which array is
//                                                 "real" vs "imaginary",
//                                                 see Rank2Orchestration.js)
//   rdft/generic.c: hartley_r2hc, cdot_r2hc  -- real -> half-complex (fwd)
//                   hartley_hc2r, cdot_hc2r  -- half-complex -> real (inv)
//
// Twiddle table: FFTW3 builds this via kernel/twiddle.c's compute(), which
// for the `half_tw = {TW_HALF,1,0},{TW_NEXT,1,0}` instruction these
// solvers' awake() uses, works out to: for outer index i (1..(n-1)/2) and
// inner index k (1..(n-1)/2), twiddle[i][k] = realCexp(MULMOD(k, i-1, n), n)
// -- traced by hand from compute()'s TW_HALF case (kernel/twiddle.c) plus
// generic.c's awake() (X(twiddle_awake)(wakefulness, &ego->td, half_tw, n,
// n, (n-1)/2)). MULMOD(x,y,p) is plain (x*y)%p for our tiny sizes
// (kernel/ifftw.h's fast-path condition x <= 92681-y always holds here).
// This port builds that same table directly with realCexp rather than
// replicating FFTW3's twiddle-cache/bytecode machinery (twlist hashing,
// tw_instr programs) -- that machinery exists to let FFTW3 reuse tables
// across many plans; it doesn't change the NUMBERS, only how they're
// cached, and reproducing the numbers is all that matters for this port.
//
// FFT_SIGN = -1 (kernel/ifftw.h) -- the forward-transform sign convention;
// r2hc/hc2r's `#if FFT_SIGN == -1` branches are taken as written below.
//
// Source: third_party/fftw-3.3.8/{dft,rdft}/generic.c, FFTW3 (c) 2003,
// 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const { realCexp } = require('./Trig.js');

// ---------------------------------------------------------------------------
// buildTwiddleTable(n) -- shared by all three solvers below (same awake()
// logic in dft/generic.c and rdft/generic.c).
//
// Derivation (kernel/twiddle.c's compute(), TW_HALF case): tw_instr is
// {unsigned char op; signed char v; short i;} (kernel/ifftw.h), so
// half_tw's `{TW_HALF, 1, 0}` sets v=1, NOT v=0 -- easy to misread as a
// 0-initialized struct literal, but the 2nd field is v=1 (confirmed against
// real compiled FFTW3 output for a size-7 case, whose actual bin-1 twiddle
// values only match MULMOD(k, j+1, n), not MULMOD(k, j, n)). compute()'s
// outer loop `j` (0..(n-1)/2-1) is consumed by apply()/apply_r2hc()/
// apply_hc2r()'s outer loop `i` (1..(n-1)/2) via j = i-1 (first apply() call
// reads the first-written twiddle block, etc.), so the final exponent is
// MULMOD(inner_k, (outer_i - 1) + 1, n) = MULMOD(inner_k, outer_i, n).
//
// Returns a flat array `W` such that block i (1-indexed, i=1..(n-1)/2)
// starts at W[(i-1)*(n-1)] and contains (n-1)/2 [c,s] pairs for k=1..(n-1)/2,
// i.e. W[(i-1)*(n-1) + (k-1)*2] = c, W[(i-1)*(n-1) + (k-1)*2 + 1] = s,
// where [c,s] = realCexp(MULMOD(k, i, n), n) -- matching the exact layout
// `cdot`/`cdot_r2hc`/`cdot_hc2r` read sequentially (w[0],w[1] per iteration,
// w += 2), and `apply`/`apply_r2hc`/`apply_hc2r` advancing W += n-1 per
// outer i.
// ---------------------------------------------------------------------------
function buildTwiddleTable(n) {
  const half = (n - 1) / 2;
  const W = new Float64Array(half * (n - 1));
  for (let i = 1; i <= half; i++) {
    const base = (i - 1) * (n - 1);
    for (let k = 1; k <= half; k++) {
      const m = (k * i) % n; // MULMOD(k, i, n)
      const [c, s] = realCexp(m, n);
      W[base + (k - 1) * 2] = c;
      W[base + (k - 1) * 2 + 1] = s;
    }
  }
  return W;
}

// Twiddle tables are pure functions of n -- cache like FFTW3's own twlist
// does (reused across many calls at the same size within one convolution).
const _twiddleCache = new Map();
function getTwiddleTable(n) {
  let W = _twiddleCache.get(n);
  if (!W) { W = buildTwiddleTable(n); _twiddleCache.set(n, W); }
  return W;
}

// ---------------------------------------------------------------------------
// dft/generic.c: hartley(n, xr, xi, xs, o, &pr, &pi)
// Folds a length-n complex sequence (xr/xi, stride xs) into a
// real "Hartley-style" working buffer o (length 2n), returning the DC sums
// in pr/pi. Ported with explicit array + offset + stride args (JS has no
// pointer arithmetic) rather than raw pointers -- purely a mechanical
// translation, same reads/writes/order.
// ---------------------------------------------------------------------------
function hartley(n, xr, xrOff, xi, xiOff, xs, o, oOff) {
  let oi = oOff;
  let sr, si;
  o[oi] = sr = xr[xrOff];
  o[oi + 1] = si = xi[xiOff];
  oi += 2;
  for (let i = 1; i + i < n; ++i) {
    const a = xr[xrOff + i * xs];
    const b = xr[xrOff + (n - i) * xs];
    sr += (o[oi] = a + b);
    const c = xi[xiOff + i * xs];
    const d = xi[xiOff + (n - i) * xs];
    si += (o[oi + 1] = c + d);
    o[oi + 2] = a - b;
    o[oi + 3] = c - d;
    oi += 4;
  }
  return [sr, si];
}

// dft/generic.c: cdot(n, x, w, &or0, &oi0, &or1, &oi1)
function cdot(n, x, xOff, w, wOff) {
  let xi = xOff;
  let rr = x[xi], ri = 0, ir = x[xi + 1], ii = 0;
  xi += 2;
  let wi = wOff;
  for (let i = 1; i + i < n; ++i) {
    rr += x[xi] * w[wi];
    ir += x[xi + 1] * w[wi];
    ri += x[xi + 2] * w[wi + 1];
    ii += x[xi + 3] * w[wi + 1];
    xi += 4;
    wi += 2;
  }
  return [rr + ii, ir - ri, rr - ii, ir + ri]; // or0, oi0, or1, oi1
}

// ---------------------------------------------------------------------------
// dftGeneric -- port of dft/generic.c's apply().
// Full complex DFT of length n (n odd prime), in place semantics matched by
// the caller (Rank2Orchestration.js): ri/ii = input real/imag arrays
// (stride `is`, offset riOff/iiOff), ro/io = output (stride `os`, offset
// roOff/ioOff). Direction (forward vs inverse) is NOT selected here -- per
// upstream, the caller achieves that by swapping which array is passed as
// "real" vs "imaginary" (see Rank2Orchestration.js's INVERSE column step).
// ---------------------------------------------------------------------------
function dftGeneric(n, ri, riOff, ii, iiOff, is, ro, roOff, io, ioOff, os) {
  const W = getTwiddleTable(n);
  const buf = new Float64Array(n * 2);
  const [pr, pi] = hartley(n, ri, riOff, ii, iiOff, is, buf, 0);
  ro[roOff] = pr;
  io[ioOff] = pi;

  let wBase = 0;
  for (let i = 1; i + i < n; ++i) {
    const [or0, oi0, or1, oi1] = cdot(n, buf, 0, W, wBase);
    ro[roOff + i * os] = or0;
    io[ioOff + i * os] = oi0;
    ro[roOff + (n - i) * os] = or1;
    io[ioOff + (n - i) * os] = oi1;
    wBase += n - 1;
  }
}

// ---------------------------------------------------------------------------
// rdft/generic.c: hartley_r2hc / cdot_r2hc -- forward real -> half-complex.
// ---------------------------------------------------------------------------
function hartleyR2hc(n, x, xOff, xs, o, oOff) {
  let oi = oOff;
  let sr = x[xOff];
  o[oi] = sr;
  oi += 1;
  for (let i = 1; i + i < n; ++i) {
    const a = x[xOff + i * xs];
    const b = x[xOff + (n - i) * xs];
    sr += (o[oi] = a + b);
    o[oi + 1] = b - a; // FFT_SIGN == -1 branch
    oi += 2;
  }
  return sr;
}

function cdotR2hc(n, x, xOff, w, wOff) {
  let xi = xOff;
  let rr = x[xi], ri = 0;
  xi += 1;
  let wi = wOff;
  for (let i = 1; i + i < n; ++i) {
    rr += x[xi] * w[wi];
    ri += x[xi + 1] * w[wi + 1];
    xi += 2;
    wi += 2;
  }
  return [rr, ri]; // or0, oi1
}

// rdftGenericR2hc -- port of rdft/generic.c's apply_r2hc(), FUSED with
// rdft2-rdft.c's hc2c() (the bridge FFTW3 itself uses to feed a plain-RDFT
// r2hc solver's packed halfcomplex output into an RDFT2 problem's separate
// real/imaginary arrays, which the 2D rank>=2 orchestration needs). hc2c
// is pure data movement (rio[i*os]=r[i]; iio[i*os]=r[n-i]; no arithmetic),
// so fusing it here changes nothing numerically -- it just skips writing
// to an intermediate packed buffer only to immediately unpack it, while
// still producing bit-for-bit the same values hc2c would.
//
// Output here is a single INTERLEAVED complex array (matching fftw_complex
// / verify/ground_truth_harness.c's layout: O[2k]=Re, O[2k+1]=Im for bin
// k), not FFTW's internal packed-real halfcomplex format -- i.e. this
// function's OUTPUT is what rdft/generic.c's apply_r2hc + rdft2-rdft.c's
// hc2c together produce, expressed directly.
function rdftGenericR2hc(n, x, xOff, xs, O, oOff, os) {
  const W = getTwiddleTable(n);
  const buf = new Float64Array(n);
  const pr = hartleyR2hc(n, x, xOff, xs, buf, 0);
  O[oOff] = pr;         // hc2c: rio[0] = r[0]
  O[oOff + 1] = 0;       // hc2c: iio[0] = 0

  let wBase = 0;
  for (let i = 1; i + i < n; ++i) {
    const [or0, oi1] = cdotR2hc(n, buf, 0, W, wBase);
    // hc2c: rio[i*os] = r[i] ; iio[i*os] = r[n-i]  (r[i]=or0, r[n-i]=oi1
    // per cdot_r2hc's own hartley-trick derivation -- see cdot_r2hc above)
    O[oOff + i * os] = or0;
    O[oOff + i * os + 1] = oi1;
    wBase += n - 1;
  }
}

// ---------------------------------------------------------------------------
// rdft/generic.c: hartley_hc2r / cdot_hc2r -- inverse half-complex -> real.
// Like rdftGenericR2hc, FUSED with rdft2-rdft.c's c2hc() (the exact inverse
// of hc2c -- also pure data movement: packedR[0]=rio[0]; packedR[i]=rio[i];
// packedR[n-i]=iio[i]; no arithmetic). Input here is a single INTERLEAVED
// complex array X (X[2k]=Re(bin k), X[2k+1]=Im(bin k)), matching
// rdftGenericR2hc's output format and verify/ground_truth_harness.c's
// fftw_complex layout -- not FFTW's internal packed-real format.
// ---------------------------------------------------------------------------
function hartleyHc2r(n, X, xOff, xs, o, oOff) {
  let oi = oOff;
  // c2hc: packedR[0] = rio[0] = X[0]  (X[1], bin 0's imaginary part, is
  // unused -- a real signal's DC bin is always real; matches c2hc exactly,
  // which never reads iio[0] either.)
  let sr = X[xOff];
  o[oi] = sr;
  oi += 1;
  for (let i = 1; i + i < n; ++i) {
    // c2hc: packedR[i] = rio[i*is] = Re(bin i) ; packedR[n-i] = iio[i*is] = Im(bin i)
    const a = X[xOff + i * xs];       // Re(bin i)
    o[oi] = a + a;
    sr += o[oi];
    const b = X[xOff + i * xs + 1];   // Im(bin i)
    o[oi + 1] = b + b;
    oi += 2;
  }
  return sr;
}

function cdotHc2r(n, x, xOff, w, wOff) {
  let xi = xOff;
  let rr = x[xi], ii = 0;
  xi += 1;
  let wi = wOff;
  for (let i = 1; i + i < n; ++i) {
    rr += x[xi] * w[wi];
    ii += x[xi + 1] * w[wi + 1];
    xi += 2;
    wi += 2;
  }
  return [rr - ii, rr + ii]; // FFT_SIGN==-1: or0 = rr-ii, or1 = rr+ii
}

// rdftGenericHc2r -- port of rdft/generic.c's apply_hc2r(), fused with c2hc
// (see hartleyHc2r above). Interleaved-complex input (n, complex-stride xs)
// -> real output (stride os).
function rdftGenericHc2r(n, X, xOff, xs, O, oOff, os) {
  const W = getTwiddleTable(n);
  const buf = new Float64Array(n);
  const pr = hartleyHc2r(n, X, xOff, xs, buf, 0);
  O[oOff] = pr;

  let wBase = 0;
  for (let i = 1; i + i < n; ++i) {
    const [or0, or1] = cdotHc2r(n, buf, 0, W, wBase);
    O[oOff + i * os] = or0;
    O[oOff + (n - i) * os] = or1;
    wBase += n - 1;
  }
}

module.exports = {
  buildTwiddleTable, getTwiddleTable,
  hartley, cdot, dftGeneric,
  hartleyR2hc, cdotR2hc, rdftGenericR2hc,
  hartleyHc2r, cdotHc2r, rdftGenericHc2r,
};

// ---------------------------------------------------------------------------
// Self-test (node GenericSolver1D.js) -- smoke test only. Round-trips a
// small odd-prime-length real signal through rdftGenericR2hc then
// rdftGenericHc2r (unnormalized, so expect n*original back) and checks
// dftGeneric's DC term against a plain sum. NOT a substitute for the real
// bit-exact verification against verify/ground_truth_harness.c fixtures.
// ---------------------------------------------------------------------------
if (require.main === module) {
  const n = 7;
  const x = new Float64Array([1, 2, 3, 4, 5, 6, 7]);

  console.log('=== rdftGenericR2hc / rdftGenericHc2r round-trip (n=7) ===');
  // n/2+1 complex bins, interleaved [re,im] -> 2*(n/2+1) reals. os=2 is the
  // complex stride (each bin occupies 2 consecutive reals).
  const nbins = Math.floor(n / 2) + 1;
  const complexOut = new Float64Array(nbins * 2);
  rdftGenericR2hc(n, x, 0, 1, complexOut, 0, 2);
  console.log('r2hc output (interleaved complex):', Array.from(complexOut));

  const back = new Float64Array(n);
  rdftGenericHc2r(n, complexOut, 0, 2, back, 0, 1);
  const scaled = Array.from(back).map((v) => v / n);
  console.log('hc2r inverse / n:', scaled, '(expect original input back)');

  console.log('\n=== dftGeneric DC term check (n=7) ===');
  const ri = new Float64Array([1, 2, 3, 4, 5, 6, 7]);
  const ii = new Float64Array(7);
  const ro = new Float64Array(7);
  const io = new Float64Array(7);
  dftGeneric(n, ri, 0, ii, 0, 1, ro, 0, io, 0, 1);
  console.log('DC (ro[0]):', ro[0], '(expect sum =', ri.reduce((a, b) => a + b, 0), ')');
}
