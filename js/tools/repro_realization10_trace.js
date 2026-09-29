'use strict';
// One-off standalone repro (2026-09-29, Dan) -- kept as a template, not
// wired into any test runner. Loads a saved realization_payload.json and
// calls runBootstrapRealization(realizationIndex, payload) DIRECTLY in the
// main thread (no worker pool), with full TRACE_DEBUG/PARITY_DEBUG tracing
// -- gives a clean, non-interleaved per-call trace for one specific
// bootstrap realization's own resample+SoFiA+fit, using the REAL payload a
// real --local run built (no hand-transcription risk, unlike reconstructing
// inputs by hand -- see repro_resample.js/repro_sofia.js for why that
// matters).
//
// To regenerate the payload for a different run: run_both.js's --local
// jsLocal leg always writes js/DCPjobData/realization_payload.json (one
// shared payload covering every realization for that run) -- copy it out
// to /tmp/resample_repro/ IMMEDIATELY after the run finishes, before the
// next run overwrites it.
//
// RESULT (2026-09-29): this is what actually found bugs #11/#12/#13 (see
// JS_FORTRAN_PARITY_SESSION_2026-09-28.md) -- the FULLVEC trace showed
// call 1's own X/Y/Inc already differing from Fortran's, before ANY amoeba
// iteration, which redirected the whole investigation from "somewhere in
// the fit" to "the pre-fit centre/shape estimate", which led to the FITS-
// header bugs in DataCubeFits.js.
process.env.TRACE_DEBUG = '1';
process.env.TRACE_DEBUG_FINALVEC_FILE = '/tmp/resample_repro/js_trace';

const fs = require('fs');
const { runBootstrapRealization } = require('../bootstrap-realization-launcher.js');

async function main() {
  const payload = JSON.parse(fs.readFileSync('/tmp/resample_repro/realization_payload.json', 'utf8'));
  const report = await runBootstrapRealization(10, payload);
  console.log('report', JSON.stringify({
    chi2: report.chi2, converged: report.converged,
    XCENTER: report.XCENTER, YCENTER: report.YCENTER,
    INCLINATION: report.INCLINATION, POSITIONANGLE: report.POSITIONANGLE,
    timings: report.timings,
  }, null, 2));
}

main().catch((e) => { console.error(e); process.exit(1); });
