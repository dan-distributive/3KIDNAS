Synced copies, used only by `galaxy-fit.html`'s "run locally (testing)" checkbox for the
geometry-estimate step, so that page works fully self-contained under `third_party/DCP/` regardless
of where it's served from -- no dependency on a specific server root one level up.

The real dispatch path (Node CLI, real DCP jobs) never reads these -- it always uses the canonical
originals via `require(...)`/`job.requires()`:

- `cfitsio-module.js`, `cfitsio-wasm.js` <- `third_party/cfitsio-4.6.3/wasm/`
- `sofia-module.js`, `sofia-wasm.js` <- `third_party/SoFiA-2-master_2_5_1/wasm/`

If cfitsio or SoFiA-2 ever get rebuilt, re-copy these four files from their canonical locations
above -- nothing else needs to change.
