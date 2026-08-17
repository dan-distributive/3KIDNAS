'use strict';
require('/Users/dandesjardins/DCP/node_modules/bravojs/utility/cjs2-node.js');
const { runTrigProbe } = require('./TrigProbe.js');
runTrigProbe().then((result) => {
  require('fs').writeFileSync(__dirname + '/trigprobe_local_result.json', JSON.stringify(result, null, 2));
  console.log('Wrote trigprobe_local_result.json');
}).catch((e) => { console.error('THREW:', e.stack); process.exit(1); });
