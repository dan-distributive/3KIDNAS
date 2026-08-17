'use strict';

// =============================================================================
// q1.js -- faithful JS port of the FULL "q1_r embedded cofactor" computation:
// dft/dftw-directsq.c's q1_r codelet (cldw, the DIF twiddle-transpose step)
// FOLLOWED BY dft/vrank-geq1.c's peeling wrapper around the recursive "cld"
// sub-transform (dft/ct.c's mkplan DECDIF+TRANSPOSE branch). Both steps
// together replace what Composite1D.js's normal per-phase Stage-1 loop does
// for a plain t1/t2 cofactor -- see Planner/chooseDecomposition.js's
// chooseComplexEmbedded() header for the full derivation of exactly when
// real FFTW selects this path (an IN-PLACE cofactor nested under an outer
// radix-r split where r == the ambient vector count -- i.e. rOuter ==
// parentRadix, always true whenever q1 applies at all, since a q1_R
// codelet's own applicability requires its radix to equal that ambient v).
//
// STRUCTURE (reverse-engineered by literal transcription of q1_3.c's
// non-FMA branch's arithmetic + a careful trace of dft/ct.c's mkplan()
// DECDIF+TRANSPOSE branch's tensor construction for the recursive "cld",
// then validated bit-exact end-to-end against real fftw_execute() for an
// in-place n=27 DFT (dft-ct-dit/3(t1_3, dft-ct-dif/3(q1_3, (vrank>=1-x3
// n1_3)))) and an in-place n=125 DFT (dft-ct-dit/5(t2_5, dft-ct-dif/5
// (q1_5, (vrank>=1-x5 n1_5)))) -- both using the SAME shared (not
// independently-regenerated) input/output fixture to rule out floating-
// point noise from a libm-vs-Math.sin mismatch, per this project's
// standing bit-exactness bar):
//
// Given outer radix `rOuter` (== parentRadix == q1's own registered
// radix) and further-recursion size `k` (so this q1 layer's own cofactor,
// m_outer = rOuter*k, is what chooseComplexEmbedded(m_outer, rOuter)
// predicts), q1 operates on the FULL rOuter^2*k-element buffer directly
// -- NOT a per-phase-gathered slice, since it genuinely reads/writes
// across what would otherwise be `rOuter` independent phases:
//
//   position(mi, a, b) = mi*rOuter + a*1 + b*(k*rOuter)
//     mi = 0..k-1 (this q1 layer's own recursion index)
//     a  = 0..rOuter-1 (the AMBIENT phase index, inherited from the outer
//          radix-rOuter split -- this is the dimension a plain t1/t2
//          cofactor would have gathered independently per phase)
//     b  = 0..rOuter-1 (q1's OWN radix index)
//
//   1. cldw (q1_rOuter itself): for each mi, for each a: read row b=0..
//      rOuter-1 at the position above, compute a standard rOuter-point
//      DFT (FFT_SIGN=-1) of that row -- confirmed BIT-IDENTICAL to this
//      port's own noTwiddle[r] (n1_r) codelet's arithmetic (same operand
//      groupings, same order -- genfft generates the same base DFT-r
//      algorithm whether embedding it in "notw" or "twidsq" context), so
//      that's reused directly. freq[a][j] = DFT bin j of row a, j=0..
//      rOuter-1. TRANSPOSE: freq[a][j] becomes output position (j,a)
//      (i.e. output[j][a] = freq[a][j] -- the input's "a" becomes the
//      output's "b" role). j=0 (DC) written untwiddled; j=1..rOuter-1
//      multiplied by realCexp(j*mi, rOuter*k) (this port's standard
//      conj-multiply convention), depending on (j,mi) only, NOT on a.
//   2. cld (vrank-geq1 wrapper around chooseComplexEmbedded(k, undefined)
//      -- for k=r_codelet-sized base cases this is just direct n1_k, but
//      in general it's k's own FULL recursive sub-plan): for each of the
//      rOuter*rOuter (j,a) pairs left untouched by step 1's own row/column
//      roles, gather the k elements at stride rOuter starting at offset
//      (j + a*k*rOuter), run them through the sub-plan, write back in
//      place at the same k positions.
//
// The combined result, read at position p*m_outer+col (p=0..rOuter-1,
// col=0..m_outer-1, m_outer=rOuter*k), is IDENTICAL to what Composite1D.js
// 's normal Stage-1 loop produces for a plain (non-q1) cofactor -- i.e.
// this is a drop-in replacement for that loop's output, not for its
// per-call signature (it needs the FULL un-gathered buffer, not one
// phase at a time). Confirmed via the position-identity p=a, col=j+rOuter*s
// (s=mi): p*m_outer+col = a*rOuter*k + j + rOuter*mi = mi*rOuter+a+j*k*rOuter
// with b re-labeled... this is exactly position(mi,a,j) from step 1's own
// output-write, so no separate reshuffle step is needed.
//
// Source: FFTW3 (c) 2003, 2007-14 Matteo Frigo and MIT, GPLv2+.
// =============================================================================

const { realCexp } = require('../../Trig');
const { noTwiddle } = require('./index');

// Radixes 2, 3, 5, 8 are independently bit-exact verified end-to-end
// against real fftw_execute() (n=24, n=27, n=125, n=128 respectively).
// Radixes 4 and 6 are DELIBERATELY EXCLUDED here: q1_6 was tried (same
// generic executeQ1Embedded code path, no radix-specific logic) and found
// to diverge from real FFTW by ~1e-14 for n=216 in-place -- confirmed via
// a direct dftComposeEmbedded(216,...) vs real fftw_execute() comparison
// after a broad safety sweep surfaced ~1871 mismatches once
// isFullyPortedComplexEmbedded started trusting q1-based structures. q1_4
// was never independently numerically verified either (no clean, non-
// bluestein-buried real N was found to exercise it standalone), so it's
// excluded on the same "don't trust what wasn't proven" basis, even
// though its row-DFT arithmetic was separately confirmed byte-identical
// to noTwiddle[4] via direct C source comparison. Do not re-add either
// without first getting a real, isolated ground-truth match to 0 maxDiff
// (see Codelets/complex/q1.js's own header for the verification method:
// an isolated fftw_plan_guru_dft matching dft/ct.c's DECDIF+TRANSPOSE
// tensor construction, or a natural top-level N where fftw_fprint_plan
// shows the radix in question -- e.g. via a small C harness sweeping
// fftw_plan_dft_1d(n,io,io,...) in-place and grepping for "q1_r").
const Q1_RADIXES = [2, 3, 5, 8];

// executeQ1Embedded(rOuter, k, subExecute, reFull, imFull) -- reFull/imFull
// are length rOuter*rOuter*k, native "physical position" layout (see
// header). subExecute(k, reK, imK) -> [roK, ioK] executes the k-sized
// cofactor's own sub-plan (Composite1D.js passes its own recursive
// executePlan here). Returns [workRe, workIm], same length, laid out so
// that workRe[p*mOuter+col] (mOuter=rOuter*k) is Stage-2's expected input
// -- a direct drop-in for what a per-phase Stage-1 loop would have
// produced for a plain (non-q1) cofactor.
function executeQ1Embedded(rOuter, k, subExecute, reFull, imFull) {
  const n = reFull.length;
  const irs = k * rOuter;
  const rowFn = noTwiddle[rOuter];
  const nTwiddle = rOuter * k;

  const stage1Re = reFull.slice();
  const stage1Im = imFull.slice();

  for (let mi = 0; mi < k; mi++) {
    const base = mi * rOuter;
    const freqRe = [];
    const freqIm = [];
    for (let a = 0; a < rOuter; a++) {
      const rowRe = new Float64Array(rOuter);
      const rowIm = new Float64Array(rOuter);
      for (let b = 0; b < rOuter; b++) {
        const pos = base + a + b * irs;
        rowRe[b] = reFull[pos];
        rowIm[b] = imFull[pos];
      }
      const [fr, fi] = rowFn(rowRe, rowIm);
      freqRe.push(fr);
      freqIm.push(fi);
    }
    for (let a = 0; a < rOuter; a++) {
      const pos0 = base + 0 + a * irs;
      stage1Re[pos0] = freqRe[a][0];
      stage1Im[pos0] = freqIm[a][0];
    }
    for (let j = 1; j < rOuter; j++) {
      const [Wc, Ws] = realCexp((j * mi) % nTwiddle, nTwiddle);
      for (let a = 0; a < rOuter; a++) {
        const pos = base + j + a * irs;
        const xr = freqRe[a][j];
        const xi = freqIm[a][j];
        stage1Re[pos] = xr * Wc + xi * Ws;
        stage1Im[pos] = xi * Wc - xr * Ws;
      }
    }
  }

  const outRe = stage1Re.slice();
  const outIm = stage1Im.slice();
  for (let j = 0; j < rOuter; j++) {
    for (let a = 0; a < rOuter; a++) {
      const off = j + a * irs;
      const reK = new Float64Array(k);
      const imK = new Float64Array(k);
      for (let s = 0; s < k; s++) {
        reK[s] = stage1Re[off + s * rOuter];
        imK[s] = stage1Im[off + s * rOuter];
      }
      const [roK, ioK] = subExecute(k, reK, imK);
      for (let s = 0; s < k; s++) {
        outRe[off + s * rOuter] = roK[s];
        outIm[off + s * rOuter] = ioK[s];
      }
    }
  }

  return [outRe, outIm];
}

module.exports = { executeQ1Embedded, Q1_RADIXES };
