# 3KIDNAS on DCP
An automated kinematic pipeline that produces kinematic models of rotating disk galaxies with robust and statistically meaningful uncertainties.

Original github source: https://github.com/NateDeg/3KIDNAS

minimal js port of fortran components of 3KIDNAS used in bootstrap computations. Deployed on DCP to exploit parallelism.

/**
 * @file        bootstrap-job.js
 * @description DCP job driver for 3KIDNAS bootstrap fitting.
 *              Generates N bootstrap realizations of a galaxy kinematic fit
 *              in parallel across a DCP worker pool. Each worker sandbox receives
 *              one realization index and runs a full GalaxyFit_Simple on the
 *              corresponding bootstrap cube.
 *
 * @authors     Dan Desjardins <dan@distributive.network>
 * @date        July 2026
 * @copyright   2026 Distributive Corp.
 *
 * @usage       node bootstrap-job.js --apiKey=0x<identity> [options]
 *
 *   Flags:
 *     --apiKey=0x...              Identity / API key to run the job under (required)
 *     --computeGroup=key,secret   Compute group as joinKey,joinSecret
 *                                 (secret optional; e.g. --computeGroup=demo,dcp
 *                                 or --computeGroup=public). Defaults to public group.
 *     --bootstrapts=N             Number of bootstrap realizations to perform (default 1000)
 *     --slicePrice=N              Per-slice price in DCC (default 5.24)
 *
 *   Examples:
 *     node bootstrap-job.js --apiKey=0x45d7... --computeGroup=demo,dcp --bootstraps=50 --slicePrice=5.24
 *     node bootstrap-job.js --apiKey=0x45d7...
 *     node bootstrap-job.js --apiKey=0x45d7... --computeGroup=public
node bootstrap-job.js \
  --apiKey=0x45d77fb82cf7b021580c98512c6103cf23f08dacbb1e5ea604a670b388fef0a0 \
  --computeGroup=google,95rhwgha \
  --bootstraps=1000 \
  --slicePrice=5.24
 *
 * @requires    node.js
 * @requires    dcp-client
 */
