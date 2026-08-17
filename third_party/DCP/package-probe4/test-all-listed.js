'use strict';
async function testProbe() {
  const setProgress = (p) => { if (typeof progress === 'function') progress(p); };
  setProgress(0.0);
  try {
    const entry = require('EntryProbe.js');
    setProgress(1.0);
    return { ok: true, vals: entry.vals };
  } catch (e) {
    setProgress(1.0);
    return { ok: false, error: e.message };
  }
}

async function main() {
  const identity = require('dcp/identity');
  const compute = require('dcp/compute');
  const apiKey = '0xf1512793d2dcb94a0102d53e6ab55ac8b145982342eae999be826aed54533ec7';
  await identity.set(apiKey);
  const job = compute.for([0], testProbe, []);
  job.requires([
    "3kidnas-probe4/EntryProbe.js", "3kidnas-probe4/Dep01.js", "3kidnas-probe4/Dep02.js",
    "3kidnas-probe4/Dep03.js", "3kidnas-probe4/Dep04.js", "3kidnas-probe4/Dep05.js",
    "3kidnas-probe4/Dep06.js", "3kidnas-probe4/Dep07.js", "3kidnas-probe4/Dep08.js",
  ]);
  job.computeGroups = [{ joinKey: 'public' }];
  job.public = { name: '3kidnas-probe4 (all-files-listed, 8 dep)', description: 'Confirming the all-files-listed fix scales past 3' };
  job.on('readystatechange', (ev) => console.log(`Ready state: ${ev}`));
  job.on('accepted', () => console.log(`  Job id: ${job.id}\n  Awaiting result...`));
  job.on('error', (error) => console.error('  Job error:', error));
  job.on('nofunds', (ev) => console.log(ev));
  const [result] = await job.exec();
  console.log('\n=== RESULT ===\n');
  console.log(JSON.stringify(result, null, 2));
}
require('dcp-client').init().then(main).catch((e) => { console.error(e); process.exit(1); });
