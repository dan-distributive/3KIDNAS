'use strict';

// =============================================================================
// extract_opcosts.js
// Dev-time-only offline extraction: reads FFTW3's own vendored codelet source
// (third_party/fftw-3.3.8/{dft,rdft}/scalar/...) and produces
// Planner/registrationTables.js, a checked-in JS data module with one entry
// per registered codelet: its ops-cost tuple ({add,mul,fma,other}, straight
// out of the codelet's own `desc` struct literal) and its registration-order
// position (LIFO matters -- see chooseDecomposition.js).
//
// Run manually: `node extract_opcosts.js` from this directory. Not required
// at runtime -- registrationTables.js is checked into git as plain data so
// this port has no runtime dependency on third_party/fftw-3.3.8 (important
// for eventual standalone-repo extraction, see FFTW3JS/README.md).
//
// Struct formats and registration order verified by direct source reading
// this session (see Planner/registrationTables.js's own header for the
// summary table). Key structural facts this parser depends on:
//   - Every codelet file has TWO copies of its desc struct, guarded by
//     `#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)
//      ... #else ... #endif`. We want ONLY the #else (non-FMA) copy, matching
//     this project's -ffp-contract=off build convention.
//   - Registration order = literal order of the SOLVTAB(...) array in each
//     family's codlist.c (dft/scalar/codelets/codlist.c,
//     rdft/scalar/r2cf/codlist.c, rdft/scalar/r2cb/codlist.c). FFTW's
//     register_solver() PREPENDS to a linked list, so the planner tries
//     solvers in REVERSE of this array order (last-registered first) --
//     chooseDecomposition.js is responsible for that reversal, not this file;
//     here we just record the literal source order faithfully.
//   - okp() stride-matching is a structural no-op for every scalar codelet in
//     this codebase (is/os/ivs/ovs and rs/vs/ms are always hardcoded 0 in
//     every desc literal, which okp() treats as "don't care") -- confirmed
//     by grep across all families, so this extractor does not need to model
//     it at all.
//   - Out of scope (not extracted): r2cfII_*/r2cbIII_* (R2HCII/HC2RIII kinds,
//     serve REDFT/RODFT only, irrelevant to ordinary r2c/c2r -- see
//     FFTW3JS/README.md's scope boundary).
// =============================================================================

const fs = require('fs');
const path = require('path');

const FFTW_ROOT = path.resolve(__dirname, '../../../../third_party/fftw-3.3.8');
const DFT_CODELETS_DIR = path.join(FFTW_ROOT, 'dft/scalar/codelets');
const R2CF_DIR = path.join(FFTW_ROOT, 'rdft/scalar/r2cf');
const R2CB_DIR = path.join(FFTW_ROOT, 'rdft/scalar/r2cb');

function readSolvtabOrder(codlistPath, solvtabName) {
  const text = fs.readFileSync(codlistPath, 'utf8');
  const marker = `const solvtab X(${solvtabName})`;
  const start = text.indexOf(marker);
  if (start === -1) throw new Error(`solvtab ${solvtabName} not found in ${codlistPath}`);
  const end = text.indexOf('SOLVTAB_END', start);
  const block = text.slice(start, end);
  const names = [];
  const re = /SOLVTAB\(X\(codelet_(\w+)\)\)/g;
  let m;
  while ((m = re.exec(block))) names.push(m[1]);
  return names;
}

// Family classification purely from the codelet name prefix -- confirmed
// this fully determines struct type / register-function / GENUS.kind, no
// need to parse register-call source (see extract_opcosts.js header).
function classify(name) {
  if (/^n1_\d+$/.test(name)) return { structType: 'kdft_desc', kind: 'DFT', group: 'n1' };
  if (/^t1_\d+$/.test(name)) return { structType: 'ct_desc', kind: 'DFT', group: 't1' };
  if (/^t2_\d+$/.test(name)) return { structType: 'ct_desc', kind: 'DFT', group: 't2' };
  if (/^q1_\d+$/.test(name)) return { structType: 'ct_desc', kind: 'DFT', group: 'q1' };

  if (/^r2cf_\d+$/.test(name)) return { structType: 'kr2c_desc', kind: 'R2HC', group: 'r2cf' };
  if (/^hf_\d+$/.test(name)) return { structType: 'hc2hc_desc', kind: 'R2HC', group: 'hf' };
  if (/^hf2_\d+$/.test(name)) return { structType: 'hc2hc_desc', kind: 'R2HC', group: 'hf2' };
  if (/^hc2cf_\d+$/.test(name)) return { structType: 'hc2c_desc', kind: 'R2HC', group: 'hc2cf', via: 'VIA_RDFT' };
  if (/^hc2cf2_\d+$/.test(name)) return { structType: 'hc2c_desc', kind: 'R2HC', group: 'hc2cf2', via: 'VIA_RDFT' };
  if (/^hc2cfdft_\d+$/.test(name)) return { structType: 'hc2c_desc', kind: 'R2HC', group: 'hc2cfdft', via: 'VIA_DFT' };
  if (/^hc2cfdft2_\d+$/.test(name)) return { structType: 'hc2c_desc', kind: 'R2HC', group: 'hc2cfdft2', via: 'VIA_DFT' };

  if (/^r2cb_\d+$/.test(name)) return { structType: 'kr2c_desc', kind: 'HC2R', group: 'r2cb' };
  if (/^hb_\d+$/.test(name)) return { structType: 'hc2hc_desc', kind: 'HC2R', group: 'hb' };
  if (/^hb2_\d+$/.test(name)) return { structType: 'hc2hc_desc', kind: 'HC2R', group: 'hb2' };
  if (/^hc2cb_\d+$/.test(name)) return { structType: 'hc2c_desc', kind: 'HC2R', group: 'hc2cb', via: 'VIA_RDFT' };
  if (/^hc2cb2_\d+$/.test(name)) return { structType: 'hc2c_desc', kind: 'HC2R', group: 'hc2cb2', via: 'VIA_RDFT' };
  if (/^hc2cbdft_\d+$/.test(name)) return { structType: 'hc2c_desc', kind: 'HC2R', group: 'hc2cbdft', via: 'VIA_DFT' };
  if (/^hc2cbdft2_\d+$/.test(name)) return { structType: 'hc2c_desc', kind: 'HC2R', group: 'hc2cbdft2', via: 'VIA_DFT' };

  throw new Error(`unrecognized codelet name ${name}`);
}

// Extract the #else (non-FMA) desc statement's {add,mul,fma,other} tuple,
// the leading radix/n integer, and (for hc2c_desc only) the VIA_RDFT/
// VIA_DFT literal from the register call, by isolating the #else...#endif
// block first (there are exactly two #if/#else/#endif-guarded copies of
// each; we want the second).
function parseDescFile(filePath, structType) {
  const text = fs.readFileSync(filePath, 'utf8');
  const ifIdx = text.indexOf('#if defined(ARCH_PREFERS_FMA)');
  if (ifIdx === -1) throw new Error(`no ARCH_PREFERS_FMA guard in ${filePath}`);
  const elseIdx = text.indexOf('#else', ifIdx);
  const endifIdx = text.indexOf('#endif', elseIdx);
  if (elseIdx === -1 || endifIdx === -1) throw new Error(`malformed #if/#else/#endif in ${filePath}`);
  const nonFmaBlock = text.slice(elseIdx, endifIdx);

  const descRe = new RegExp(`${structType}\\s+desc\\s*=\\s*\\{([^;]*)\\}\\s*;`);
  const descMatch = nonFmaBlock.match(descRe);
  if (!descMatch) throw new Error(`no ${structType} desc statement found in #else block of ${filePath}`);
  const descBody = descMatch[1];

  const radixMatch = descBody.match(/^\s*(\d+)/);
  if (!radixMatch) throw new Error(`no leading radix/n integer in desc of ${filePath}`);
  const radix = parseInt(radixMatch[1], 10);

  const opsMatch = descBody.match(/\{\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\}/);
  if (!opsMatch) throw new Error(`no {add,mul,fma,other} ops tuple found in desc of ${filePath}`);
  const ops = {
    add: parseInt(opsMatch[1], 10),
    mul: parseInt(opsMatch[2], 10),
    fma: parseInt(opsMatch[3], 10),
    other: parseInt(opsMatch[4], 10),
  };

  let via = null;
  const viaMatch = nonFmaBlock.match(/X\(khc2c_register\)\s*\([^)]*?,\s*(HC2C_VIA_\w+)\s*\)/);
  if (viaMatch) via = viaMatch[1];

  return { radix, ops, via };
}

// Out of scope by construction: R2HCII/HC2RIII kinds (serve REDFT/RODFT
// only). classify() has no case for these on purpose -- skip them here
// before they'd hit classify()'s "unrecognized" throw.
function isOutOfScope(name) {
  return /^r2cfII_\d+$/.test(name) || /^r2cbIII_\d+$/.test(name);
}

function extractFamily(dirPath, codlistName, solvtabName) {
  const order = readSolvtabOrder(path.join(dirPath, codlistName), solvtabName).filter((n) => !isOutOfScope(n));
  const entries = [];
  for (let i = 0; i < order.length; i++) {
    const name = order[i];
    const info = classify(name);
    const filePath = path.join(dirPath, `${name}.c`);
    const { radix, ops, via } = parseDescFile(filePath, info.structType);
    entries.push({
      name,
      registrationIndex: i,
      structType: info.structType,
      group: info.group,
      kind: info.kind,
      radix,
      ops,
      via: via || info.via || null,
    });
  }
  return entries;
}

function main() {
  const complex = extractFamily(DFT_CODELETS_DIR, 'codlist.c', 'solvtab_dft_standard')
    .filter((e) => e.group !== 'q1'); // q1 (twiddle-squared DIF, N=r^2 in-place) -- low-priority narrow path, see Phase 6; excluded from the main table for now, not deleted from source.
  const q1 = extractFamily(DFT_CODELETS_DIR, 'codlist.c', 'solvtab_dft_standard')
    .filter((e) => e.group === 'q1');

  const realForwardAll = extractFamily(R2CF_DIR, 'codlist.c', 'solvtab_rdft_r2cf');
  const realBackwardAll = extractFamily(R2CB_DIR, 'codlist.c', 'solvtab_rdft_r2cb');

  // r2cfII_*/r2cbIII_* were never emitted by classify() as anything other
  // than a thrown error -- so extractFamily() would already have crashed on
  // them. Confirm scope exclusion explicitly instead: skip by name pattern.
  // (classify() intentionally has no case for them; adjust extractFamily's
  // solvtab reader to skip those names before calling parseDescFile.)

  const out = [];
  out.push("'use strict';");
  out.push('');
  out.push('// AUTO-GENERATED by verify/extract_opcosts.js -- do not hand-edit.');
  out.push('// Regenerate with: node verify/extract_opcosts.js');
  out.push('// See that file\'s header and FFTW3JS/README.md for what this data means');
  out.push('// and the scope boundary (ordinary rank-2 r2c/c2r only; R2HCII/HC2RIII and');
  out.push('// DHT/REDFT/RODFT-only families excluded).');
  out.push('//');
  out.push('// Each entry: { name, registrationIndex, structType, group, kind, radix,');
  out.push('//               ops: {add,mul,fma,other}, via }');
  out.push('// registrationIndex is literal SOLVTAB array order (source order), LOWEST');
  out.push('// first. FFTW\'s planner tries solvers in REVERSE of this (LIFO -- most');
  out.push('// recently registered first); chooseDecomposition.js does that reversal.');
  out.push('');
  out.push(`exports.complexCodelets = ${JSON.stringify(complex, null, 2)};`);
  out.push('');
  out.push(`exports.q1Codelets = ${JSON.stringify(q1, null, 2)};`);
  out.push('');
  out.push(`exports.realForwardCodelets = ${JSON.stringify(realForwardAll, null, 2)};`);
  out.push('');
  out.push(`exports.realBackwardCodelets = ${JSON.stringify(realBackwardAll, null, 2)};`);
  out.push('');

  const outPath = path.resolve(__dirname, '../Planner/registrationTables.js');
  fs.writeFileSync(outPath, out.join('\n'));
  console.log(`Wrote ${outPath}`);
  console.log(`  complex (non-q1): ${complex.length}, q1: ${q1.length}`);
  console.log(`  realForward: ${realForwardAll.length}, realBackward: ${realBackwardAll.length}`);
}

main();
