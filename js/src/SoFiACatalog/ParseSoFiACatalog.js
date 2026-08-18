'use strict';

// =============================================================================
// ParseSoFiACatalog.js
// Faithful port of FitDriverScripts/SoFiA_Driver.py's LoadSoFiAOutput --
// parses sofia2wasm's result_cat.txt (bytes/text, not a file path) the same
// way Python's np.loadtxt(...,skiprows=13,usecols=(15,27,28,33,3,4,5)) does.
//
// COLUMN INDICES -- verified against a REAL run of this project's actual
// template_par_file.par (not re-derived from the catalogue's own labeled
// header, which would have been wrong here): SoFiA's auto-generated source
// name (parameter.prefix='SoFiA') is "SoFiA J103538.30-484833.2" --
// TWO whitespace tokens, not one, because of the space after "SoFiA". That
// shifts every column's 0-based token index by +1 relative to naively
// reading the catalogue's own "labeled column number - 1". Checked directly:
// token[15]=f_sum, token[27]=ell_maj, token[28]=ell_min, token[33]=kin_pa,
// token[3]=x, token[4]=y, token[5]=z -- exactly Python's own
// usecols=(15,27,28,33,3,4,5), confirmed against real values (plausible
// flux/size/PA/pixel-position numbers, not just "a number landed there").
//
// Python's kin_pa and ell_pa (label 29, a DIFFERENT column) are NOT the same
// column -- SoFiA_Driver.py only ever reads kin_pa (token 33) and copies it
// into BOTH the 'kin_pa' and 'ell_pa' keys of the dict it hands to
// GetGeometryEstimates, so PAEstimate's nan-fallback never actually reads
// the catalogue's real ell_pa column in this pipeline's current usage.
// Replicated exactly, not "fixed", since diverging from what Python already
// does would be the actual bug here.
//
// skiprows=13 in Python's call is NOT the real reason header lines get
// skipped -- checked directly: a real catalogue has 21 lines (comments +
// one blank line) before the first data row, not 13. numpy's own default
// `comments='#'` behavior silently drops every '#'-prefixed line regardless
// of skiprows, and loadtxt skips blank lines by default too -- that's what
// actually makes Python's call work. Replicated here as "drop blank lines
// and lines starting with #", which is robust to the header being a
// different length than 13 or 21 (e.g. if reliability.enable added rows),
// rather than hardcoding a line count that happened to work once.
// =============================================================================

const COL = { FSUM: 15, ELL_MAJ: 27, ELL_MIN: 28, KIN_PA: 33, X: 3, Y: 4, Z: 5 };

/**
 * @param {Uint8Array|string} catalogText - sofia2wasm's result_cat.txt bytes, or already-decoded text
 * @returns {{sofiaSuccess: boolean, ellMaj?: number, ellMin?: number, pa?: number,
 *            x?: number, y?: number, z?: number, maskVal?: number}}
 */
function parseSoFiACatalog(catalogText) {
  const text = typeof catalogText === 'string' ? catalogText : new TextDecoder('utf-8').decode(catalogText);

  const dataLines = text.split('\n').filter((line) => {
    const trimmed = line.trim();
    return trimmed.length > 0 && !trimmed.startsWith('#');
  });

  if (dataLines.length === 0) {
    return { sofiaSuccess: false };
  }

  const rows = dataLines.map((line) => line.trim().split(/\s+/).map(Number));

  // Multiple detections -> use the entry with the maximum total flux
  // (FSum), matching Python's `np.where(FSum == np.max(FSum))[0][0]`.
  let maxFluxIndx = 0;
  let maxFSum = rows[0][COL.FSUM];
  for (let i = 1; i < rows.length; i++) {
    if (rows[i][COL.FSUM] > maxFSum) {
      maxFSum = rows[i][COL.FSUM];
      maxFluxIndx = i;
    }
  }
  const row = rows[maxFluxIndx];

  return {
    sofiaSuccess: true,
    ellMaj: row[COL.ELL_MAJ],
    ellMin: row[COL.ELL_MIN],
    pa: row[COL.KIN_PA],
    x: row[COL.X],
    y: row[COL.Y],
    z: row[COL.Z],
    // Python: WorkingDict['MaxFluxIndx']+1 -- the row's position in the
    // (already flux-ordered-by-detection, not re-sorted here) catalogue,
    // not a literal read of the catalogue's own "id" column. Replicated
    // as-is; see this file's header comment for why.
    maskVal: maxFluxIndx + 1,
  };
}

// BOTH assignments, unconditionally -- not an if/else picking one or the
// other. Two real, independent consumers need two different exports:
//   - module.exports: every require()-based consumer (Node, the real DCP
//     sandbox's module.declare(), the published package, the local
//     cjs2-node.js test shim). Must be UNCONDITIONAL -- confirmed directly
//     that the real sandbox's `module.exports` starts out falsy (unlike
//     plain Node's pre-populated {}), so gating this behind
//     `module.exports` being already-truthy (the previous, broken version
//     of this code) silently skipped the assignment there, and real
//     bootstrap dispatch failed with "parseSoFiACatalog is not a function".
//   - globalThis.ParseSoFiACatalog: index.html's own
//     runEstimateGeometryLocally() loads this file as a raw <script> tag
//     (no module wrapper, `module` genuinely undefined there) and its
//     require() shim reads window.ParseSoFiACatalog directly -- removing
//     this broke that path with "Cannot destructure property
//     'parseSoFiACatalog' of 'require(...)' as it is undefined."
// `typeof module !== 'undefined'` (not `&& module.exports`) is still the
// right guard for the first one -- module is a real function parameter
// (module.declare's factory) in every context except the raw <script> tag.
if (typeof module !== 'undefined') {
  module.exports = { parseSoFiACatalog };
}
if (typeof globalThis !== 'undefined') {
  globalThis.ParseSoFiACatalog = { parseSoFiACatalog };
}
