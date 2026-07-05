'use strict';

// =============================================================================
// random.js
// High-fidelity port of src/StandardMath/random.f (BasicRanNumGen module)
// Press et al. "Numerical Recipes in Fortran 77"
//
// PRECISION
// ---------
// Fortran `real` is 32-bit single precision. All float constants are
// hardcoded as their exact IEEE 754 single-precision values (computed via
// Python struct.pack/unpack) to match Fortran compile-time evaluation.
// Runtime float results are wrapped in Math.fround() to truncate to 32-bit.
// Integer arithmetic uses Math.imul() for signed 32-bit overflow.
//
// 1-ULP ROUNDING DIFFERENCE vs gfortran on Apple Silicon (ARM)
// ------------------------------------------------------------
// gfortran uses fused multiply-add (FMA) on ARM NEON, computing AM*iy in
// a single rounding step. JavaScript has no FMA access, so f32(AM * iy)
// rounds differently by exactly 1 ULP in roughly half of all outputs.
// The integer state (idum, idum2, iv, iy) is confirmed bit-identical to
// Fortran. The 1-ULP difference is the minimum possible error between
// FMA and non-FMA implementations and does not affect the statistical
// properties of the RNG or optimizer convergence.
// =============================================================================

const f32     = Math.fround;
const imul32  = Math.imul;
const int32   = x => x | 0;

// Exact 32-bit float constants (Fortran parameter values)
// Computed as: struct.unpack('f', struct.pack('f', 1.0/IM))[0]
const AM_IM   = 4.6566128730773926e-10;   // 1./2147483647  (ran0, ran1)
const AM_IM1  = 4.6566128730773926e-10;   // 1./2147483563  (ran2) — same f32
const FAC_ran3 = 9.999999717180685e-10;   // 1./1000000000  (ran3)


// ---------------------------------------------------------------------------
// ran0  —  Park-Miller minimal standard (no shuffle)
// ---------------------------------------------------------------------------
function ran0(state) {
  const IA = 16897, IM = 2147483647, IQ = 127773, IR = 2836, MASK = 123459876;
  let idum = int32(state.idum ^ MASK);
  const k  = Math.trunc(idum / IQ);
  idum     = int32(imul32(IA, int32(idum - imul32(k, IQ))) - imul32(IR, k));
  if (idum < 0) idum = int32(idum + IM);
  state.idum = int32(idum ^ MASK);
  return f32(AM_IM * idum);
}


// ---------------------------------------------------------------------------
// ran1  —  Park-Miller + Bays-Durham shuffle (NTAB=32)
// ---------------------------------------------------------------------------
function makeRan1State(idum) {
  return { idum: int32(idum), iv: new Int32Array(32), iy: 0, initialized: false };
}

function ran1(state) {
  const IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836;
  const NTAB = 32, NDIV = 67108864;
  const RNMX = f32(9.99999880790710449219e-01);

  if (state.idum <= 0 || !state.initialized) {
    state.idum = Math.max(-state.idum, 1);
    for (let j = NTAB + 8; j >= 1; j--) {
      const k    = Math.trunc(state.idum / IQ);
      state.idum = int32(imul32(IA, int32(state.idum - imul32(k, IQ))) - imul32(IR, k));
      if (state.idum < 0) state.idum = int32(state.idum + IM);
      if (j <= NTAB) state.iv[j - 1] = state.idum;
    }
    state.iy = state.iv[0];
    state.initialized = true;
  }

  const k    = Math.trunc(state.idum / IQ);
  state.idum = int32(imul32(IA, int32(state.idum - imul32(k, IQ))) - imul32(IR, k));
  if (state.idum < 0) state.idum = int32(state.idum + IM);
  const j         = int32(1 + Math.trunc(state.iy / NDIV));
  state.iy        = state.iv[j - 1];
  state.iv[j - 1] = state.idum;
  return Math.min(f32(AM_IM * state.iy), RNMX);
}


// ---------------------------------------------------------------------------
// ran2  —  L'Ecuyer + Bays-Durham shuffle. Primary generator in pipeline.
//
// NOTE: Fortran source typo on primary LCG advance uses IQ2 instead of IQ1.
// Reproduced exactly for bit-identical output.
// ---------------------------------------------------------------------------
function makeRan2State(idum) {
  return {
    idum:        int32(idum),
    idum2:       int32(123456789),
    iv:          new Int32Array(32),
    iy:          0,
    initialized: false
  };
}

function ran2(state) {
  const IM1 = 2147483563, IM2 = 2147483399, IMM1 = 2147483562;
  const IA1 = 40014,  IA2 = 40692;
  const IQ1 = 53668,  IQ2 = 52774;
  const IR1 = 12211,  IR2 = 3791;
  const NTAB = 32, NDIV = 67108862;
  const RNMX = f32(9.99999880790710449219e-01);

  if (state.idum <= 0) {
    state.idum  = Math.max(-state.idum, 1);
    state.idum2 = state.idum;
    for (let j = NTAB + 8; j >= 1; j--) {
      const k    = Math.trunc(state.idum / IQ1);
      state.idum = int32(imul32(IA1, int32(state.idum - imul32(k, IQ1))) - imul32(k, IR1));
      if (state.idum < 0) state.idum = int32(state.idum + IM1);
      if (j <= NTAB) state.iv[j - 1] = state.idum;
    }
    state.iy = state.iv[0];
    state.initialized = true;
  }

  // Primary LCG — Fortran typo: uses IQ2 not IQ1 in subtraction term
  let k      = Math.trunc(state.idum / IQ1);
  state.idum = int32(imul32(IA1, int32(state.idum - imul32(k, IQ2))) - imul32(k, IR1));
  if (state.idum < 0) state.idum = int32(state.idum + IM1);

  // Secondary LCG
  k           = Math.trunc(state.idum2 / IQ2);
  state.idum2 = int32(imul32(IA2, int32(state.idum2 - imul32(k, IQ2))) - imul32(k, IR2));
  if (state.idum2 < 0) state.idum2 = int32(state.idum2 + IM2);

  // Shuffle table
  const j         = int32(1 + Math.trunc(state.iy / NDIV));
  state.iy        = int32(state.iv[j - 1] - state.idum2);
  state.iv[j - 1] = state.idum;
  if (state.iy < 1) state.iy = int32(state.iy + IMM1);

  return Math.min(f32(AM_IM1 * state.iy), RNMX);
}


// ---------------------------------------------------------------------------
// ran3  —  Knuth subtractive generator
// ---------------------------------------------------------------------------
function makeRan3State(idum) {
  return { idum: int32(idum), inext: 0, inextp: 0, ma: new Int32Array(56), iff: 0 };
}

function ran3(state) {
  const MBIG = 1000000000, MSEED = 161803398, MZ = 0;

  if (state.idum < 0 || state.iff === 0) {
    state.iff   = 1;
    let mj      = Math.abs(MSEED - Math.abs(state.idum)) % MBIG;
    state.ma[54] = mj;
    let mk = 1;
    for (let i = 1; i <= 54; i++) {
      const ii         = (21 * i) % 56;
      state.ma[ii - 1] = mk;
      mk = mj - mk;
      if (mk < MZ) mk += MBIG;
      mj = state.ma[ii - 1];
    }
    for (let k = 1; k <= 4; k++) {
      for (let i = 1; i <= 56; i++) {
        state.ma[i - 1] -= state.ma[((i + 30) % 56)];
        if (state.ma[i - 1] < MZ) state.ma[i - 1] += MBIG;
      }
    }
    state.inext = 0; state.inextp = 31; state.idum = 1;
  }

  if (++state.inext  === 56) state.inext  = 1;
  if (++state.inextp === 56) state.inextp = 1;
  let mj = state.ma[state.inext - 1] - state.ma[state.inextp - 1];
  if (mj < MZ) mj += MBIG;
  state.ma[state.inext - 1] = mj;
  return f32(FAC_ran3 * mj);
}


// ---------------------------------------------------------------------------
// gasdev  —  Box-Muller Gaussian deviate (mean 0, sigma 1)
// Pairs with ran2. Caches second deviate across calls (Fortran SAVE iset,gset)
// ---------------------------------------------------------------------------
function makeGasdevState(ran2State) {
  return { ran2State, iset: 0, gset: f32(0.0) };
}

function gasdev(state) {
  if (state.ran2State.idum < 0) state.iset = 0;

  if (state.iset === 0) {
    let v1, v2, rsq;
    do {
      v1  = f32(f32(2.0) * ran2(state.ran2State) - f32(1.0));
      v2  = f32(f32(2.0) * ran2(state.ran2State) - f32(1.0));
      rsq = f32(f32(v1 * v1) + f32(v2 * v2));
    } while (rsq >= f32(1.0) || rsq === f32(0.0));

    const fac  = f32(Math.fround(Math.sqrt(f32(f32(-2.0) * f32(Math.log(rsq))) / rsq)));
    state.gset = f32(v1 * fac);
    state.iset = 1;
    return f32(v2 * fac);
  } else {
    state.iset = 0;
    return state.gset;
  }
}


// ---------------------------------------------------------------------------
// makeRng  —  convenience bundle for pipeline use
// Single idum threaded through ran2 + gasdev, matching Fortran.
// ---------------------------------------------------------------------------
function makeRng(idum) {
  const ran2State   = makeRan2State(idum);
  const gasdevState = makeGasdevState(ran2State);
  return {
    ran2:   () => ran2(ran2State),
    gasdev: () => gasdev(gasdevState),
    state:  { ran2State, gasdevState }
  };
}


module.exports = {
  ran0, ran1, ran2, ran3, gasdev,
  makeRan1State, makeRan2State, makeRan3State, makeGasdevState,
  makeRng
};


// ---------------------------------------------------------------------------
// Self-test (node random.js) — compare against Fortran reference output
// ---------------------------------------------------------------------------
if (require.main === module) {
  console.log('=== ran2 (seed -1, first 10 values) ===');
  const rng = makeRng(-1);
  for (let i = 0; i < 10; i++) console.log(i, rng.ran2().toFixed(10));
  console.log('\n=== gasdev (continuing same state, first 10 values) ===');
  for (let i = 0; i < 10; i++) console.log(i, rng.gasdev().toFixed(10));
}
