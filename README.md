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
| `--apiKey=0x...` | DCP Identity / API key to run the job under (**required**) see `Obtaining an API Key` below |
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

## Obtaining an API key

The `--apiKey` is a DCP identity that authorizes job creation and pays for compute. Each deployed job debits compute credits (DCC) from the bank account you associate with the key. Create one on the DCP portal:

1. Create an account at [dcp.cloud](https://dcp.cloud) using the standard email and password sign-up.
2. In the left nav, open **Keys**.
3. In the **API Keys** section, click **Add**. The *Add API Key* dialog opens.
4. Fill in the dialog:
   - **Purpose:** select `Job Creation`.
   - **Description:** a name to identify the key later (e.g. `3kidnas bootstrap deploy`).
   - **Payment Account:** the bank account to debit when jobs run under this key. Credits are drawn from here on every deployment, so pick the account you intend to fund this work.
   - **Expires in:** choose a lifetime (e.g. `one year`).
5. Click **Save**.
6. Back in the **API Keys** table, use the **download** or **copy** icon in the *Actions* column to retrieve the key's `0x...` value. That string is what you pass to `--apiKey`.

<!-- ![Add API Key dialog with Purpose, Description, Payment Account and Expires in filled out](docs/img/add-api-key.png) -->

> **Keep the key private.** It can spend the associated account's credits. Don't commit it to the repo or paste it into shared logs; pass it on the command line (or via an environment variable / secrets manager in automated runs). Keys expire, so if deployments start failing on auth, check whether the key has lapsed and issue a new one.
