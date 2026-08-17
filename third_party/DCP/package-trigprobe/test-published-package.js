'use strict';
async function trigProbeWork() {
  const setProgress = (p) => { if (typeof progress === 'function') progress(p); };
  setProgress(0.0);
  try {
    const { runTrigProbe } = require('TrigProbe.js');
    const result = await runTrigProbe();
    setProgress(1.0);
    return result;
  } catch (e) {
    setProgress(1.0);
    return { entryPointError: (e && e.stack) || String(e) };
  }
}

async function main() {
  const identity = require('dcp/identity');
  const compute = require('dcp/compute');
  const apiKey = '0xf1512793d2dcb94a0102d53e6ab55ac8b145982342eae999be826aed54533ec7';
  await identity.set(apiKey);
  const job = compute.for([0], trigProbeWork, []);
  job.requires([
    '3kidnas-trigprobe/TrigProbe.js',
    '3kidnas-trigprobe/fdlibm-wasm.js',
    '3kidnas-trigprobe/fdlibm-module.js',
  ]);
  job.computeGroups = [{ joinKey: 'public' }];
  job.public = { name: '3kidnas-trigprobe', description: 'Does native Math.sin/cos differ between local Node and a real DCP worker?' };
  job.on('readystatechange', (ev) => console.log(`Ready state: ${ev}`));
  job.on('accepted', () => console.log(`  Job id: ${job.id}\n  Awaiting result...`));
  job.on('error', (error) => console.error('  Job error:', error));
  job.on('nofunds', (ev) => console.log(ev));
  const [result] = await job.exec();
  require('fs').writeFileSync(__dirname + '/trigprobe_remote_result.json', JSON.stringify(result, null, 2));
  console.log('\nWrote trigprobe_remote_result.json');
}
require('dcp-client').init().then(main).catch((e) => { console.error(e); process.exit(1); });
