'use strict';
// One-off standalone repro (2026-09-29, Dan) -- companion to
// repro_resample.js; reads /tmp/resample_repro/ReproBS10.fits, so run that
// script first. Kept as a template, not wired into any test runner.
//
// Runs the WASM SoFiA build on the (proven bit-identical, see
// repro_resample.js) resampled cube, to check whether wasm SoFiA's catalog
// output diverges from NATIVE SoFiA's for this specific noisy input -- an
// entirely different cross-platform boundary (wasm vs native compiled
// SoFiA-2 C code) never examined before this session; everything else
// audited so far was JS-port-vs-Fortran, not wasm-vs-native-of-the-SAME-
// third-party-C-source.
//
// RESULT (2026-09-29): proved bit-exact (to the catalog text's own
// precision). Native SoFiA (run by hand: `third_party/SoFiA-2-master_2_5_1/
// sofia sofia_native_filled.par` against the same ReproBS10.fits) produced
// an IDENTICAL catalog row to this script's wasm output -- x, y, ell_maj,
// ell_min, ell_pa, kin_pa all matched to the catalog's 6-decimal precision.
// Combined with repro_resample.js's result, this proves the ENTIRE
// resample+SoFiA input to realization 10's fit is bit-identical between
// platforms -- the divergence is isolated to the amoeba fit itself, not
// anything upstream of it. See JS_FORTRAN_PARITY_SESSION_2026-09-28.md.
const fs = require('fs');

async function main() {
  const sofia = require('../../third_party/SoFiA-2-master_2_5_1/wasm/sofia-wasm.js');
  await sofia.ready;

  const parTemplate = fs.readFileSync(
    '/Users/dandesjardins/DCP/3KIDNAS/third_party/SoFiA-2-master_2_5_1/template_par_file.par', 'utf8');
  const lines = parTemplate.split('\n').map((l) => l + '\n');
  if (lines.length && lines[lines.length - 1] === '\n') lines.pop();
  lines[17]  = 'input.data                 = /work/cube.fits\n';
  lines[142] = 'output.directory           = /work/out\n';
  lines[143] = 'output.filename            = result\n';
  lines[149] = 'output.writeMask           =  true \n';
  const par = lines.join('');

  const cubeFitsBytes = fs.readFileSync('/tmp/resample_repro/ReproBS10.fits');
  const { exitCode, files } = await sofia.run({ cube: cubeFitsBytes, par });
  console.log('exitCode', exitCode, 'files', [...files.keys()]);
  if (files.has('result_cat.txt')) {
    fs.writeFileSync('/tmp/resample_repro/wasm_cat.txt', files.get('result_cat.txt'));
    console.log(Buffer.from(files.get('result_cat.txt')).toString('utf8'));
  }
}

main().catch((e) => { console.error(e); process.exit(1); });
