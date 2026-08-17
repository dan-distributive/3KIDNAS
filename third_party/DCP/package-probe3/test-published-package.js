'use strict';
async function testProbe() {
  const setProgress = (p) => { if (typeof progress === 'function') progress(p); };
  setProgress(0.0);
  try {
    const entry = require('EntryProbe.js');
    setProgress(1.0);
    return { ok: true, aVal: entry.aVal, bVal: entry.bVal };
  } catch (e) {
    setProgress(1.0);
    return { ok: false, error: e.message };
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
  if (!apiKey) { console.error('ERROR: --apiKey=0x... is required'); process.exit(1); }
  await identity.set(apiKey);
  const job = compute.for([0], testProbe, []);
  job.requires(['3kidnas-probe3/EntryProbe.js']);
  job.computeGroups = [{ joinKey: 'public' }];
  job.public = { name: '3kidnas-probe3 (2-dep boundary test)', description: 'Does exactly 2 deps break vs Phase 0s proven 1?' };
  job.on('readystatechange', (ev) => console.log(`Ready state: ${ev}`));
  job.on('accepted', () => console.log(`  Job id: ${job.id}\n  Awaiting result...`));
  job.on('error', (error) => console.error('  Job error:', error));
  job.on('nofunds', (ev) => console.log(ev));
  const [result] = await job.exec();
  console.log('\n=== RESULT ===\n');
  console.log(JSON.stringify(result, null, 2));
}
require('dcp-client').init().then(main).catch((e) => { console.error(e); process.exit(1); });
