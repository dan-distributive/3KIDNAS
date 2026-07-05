# bootstrap-job.js

DCP job driver for 3KIDNAS bootstrap fitting. Generates N bootstrap realizations of a galaxy kinematic fit in parallel across a DCP worker pool. Each worker sandbox receives one realization index and runs a full `GalaxyFit_Simple` on the corresponding bootstrap cube.

- **Authors:** Dan Desjardins \<dan@distributive.network\>
- **Date:** July 2026
- **Requires:** `Node.js`, `dcp-client`

## Usage

```bash
node bootstrap-job.js --apiKey=0x<identity> [options]
```

### Flags

| Flag | Description |
| --- | --- |
| `--apiKey=0x...` | Identity / API key to run the job under (**required**) |
| `--computeGroup=key,secret` | Compute group as `joinKey,joinSecret` (secret optional; e.g. `--computeGroup=demo,dcp` or `--computeGroup=public`). Defaults to public group. |
| `--bootstraps=N` | Number of bootstrap realizations to perform (default 1000) |
| `--slicePrice=N` | Per-slice price in DCC (default 5.24) |

### Examples

```bash
# Full invocation
node bootstrap-job.js \
  --apiKey=0x45d7... \
  --computeGroup=demo,dcp \
  --bootstraps=1000 \
  --slicePrice=5.24

# Minimal — just the required key, defaults for the rest
node bootstrap-job.js --apiKey=0x45d7...

# Public compute group, no secret
node bootstrap-job.js --apiKey=0x45d7... --computeGroup=public
```
