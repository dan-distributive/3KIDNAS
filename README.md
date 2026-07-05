# bootstrap-job.js

DCP job driver for 3KIDNAS bootstrap fitting. Generates N bootstrap realizations of a galaxy kinematic fit in parallel across a DCP worker pool. Each worker sandbox receives one realization index and runs a full `GalaxyFit_Simple` on the corresponding bootstrap cube.

This is a JavaScript port of the 3KIDNAS kinematic-fitting pipeline, adapted to run on the Distributive Compute Protocol (DCP). The original pipeline is written in Fortran (see [Attribution](#attribution) below).

- **Port author:** Dan Desjardins \<dan@distributive.network\> — JavaScript / DCP port
- **Date:** July 2026
- **Requires:** `Node.js`, `dcp-client`

## Attribution

This project ports the **3KIDNAS** pipeline to JavaScript so it can run distributed on DCP. It does not replace or fork the upstream science code; all credit for the underlying algorithm and its scientific validity belongs to the original authors.

**3KIDNAS** (3D Kinematic Data aNalysis Algorithm for Surveys — pronounced like "echidnas") is an automated pipeline that fits kinematic models of rotating disk galaxies and produces statistically meaningful uncertainties on their rotation curves. It is being developed by WALLABY Technical Working Group 5, led by Kristine Spekkens and Nathan Deg, and is intended to model the kinematics of thousands of galaxies resolved by the WALLABY HI survey.

- Original pipeline (Fortran): <https://github.com/NateDeg/3KIDNAS>
- Pipeline description and diagnostic plots: <https://wallaby-survey.org/tools/>

The bootstrap uncertainty estimation implemented here mirrors the upstream method; this port's contribution is parallelizing those realizations across a DCP worker pool rather than running them serially.

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
6. Back in the **API Keys** table, click the **copy to clipboard** icon in the *Actions* column to copy the key's `0x...` value. That string is what you pass to `--apiKey`. (The download icon exports the full keystore, not what you want here.)

<!-- ![Add API Key dialog with Purpose, Description, Payment Account and Expires in filled out](docs/img/add-api-key.png) -->

> **Keep the key private.** It can spend the associated account's credits. Don't commit it to the repo or paste it into shared logs; pass it on the command line (or via an environment variable / secrets manager in automated runs). Keys expire, so if deployments start failing on auth, check whether the key has lapsed and issue a new one.
>
> **Watch the account balance.** Jobs are paid for in DCC from the key's payment account, and that account can run dry. If a job won't start or stalls with no progress, check the account balance under **Bank** first. If it's low, top it up (buy credits) or earn more by running a worker on the public network.
