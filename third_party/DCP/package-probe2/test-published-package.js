/**
 * @usage node test-published-package.js --apiKey=0x<identity> [--computeGroup=key,secret]
 */
'use strict';

async function testProbe() {
  const setProgress = (p) => { if (typeof progress === 'function') progress(p); };
  setProgress(0.0);
  try {
    const { tryAll } = require('EntryProbe.js');
    const results = tryAll();
    setProgress(1.0);
    return results;
  } catch (e) {
    setProgress(1.0);
    return { entryPointError: (e && e.message) || String(e) };
  }
}

async function main() {
  const identity = require('dcp/identity');
  const compute = require('dcp/compute');

  function getFlag(name) {
    const pfx = `--${name}=`;
    const hit = process.argv.find((a) => a.startsWith(pfx));
    return hit ? hit.slice(pfx.length) : undefined;
  }

  const apiKey = getFlag('apiKey');
  const cg = getFlag('computeGroup');
  if (!apiKey) {
    console.error('ERROR: --apiKey=0x... is required');
    process.exit(1);
  }
  let computeGroup;
  if (cg) {
    const [joinKey, joinSecret] = cg.split(',');
    computeGroup = joinSecret ? { joinKey, joinSecret } : { joinKey };
  }

  await identity.set(apiKey);

  const job = compute.for([0], testProbe, []);
  job.requires(['3kidnas-probe2/EntryProbe.js']);
  job.computeGroups = computeGroup ? [computeGroup] : [{ joinKey: 'public' }];
  job.public = {
    name: '3kidnas-probe2 (breadth test)',
    description: 'Isolates whether a 14-entry module.declare() dependency array behaves differently than a 1-entry one',
  };

  job.on('readystatechange', (ev) => console.log(`Ready state: ${ev}`));
  job.on('accepted', () => console.log(`  Job id: ${job.id}\n  Awaiting result...`));
  job.on('error', (error) => console.error('  Job error:', error));
  job.on('nofunds', (ev) => console.log(ev));

  const [result] = await job.exec();
  console.log('\n=== RESULT ===\n');
  console.log(JSON.stringify(result, null, 2));
}

require('dcp-client').init().then(main).catch((e) => { console.error(e); process.exit(1); });
