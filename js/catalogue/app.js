'use strict';
// =============================================================================
// app.js -- catalogue.html: browse the WALLABY catalogue, select galaxies, fit them on DCP.
//
// Run shape: ONE initial-fit job for the selection, one slice per galaxy. Each slice is the
// URL of that galaxy's JSON on this server (workers fetch it; runInitialFit merges it over the
// shared options argument, and payload.galaxies[id] supplies its PA/Inc). As each fit lands,
// that galaxy's bootstrap job is deployed: nBoot slices of runBootstrapRealization with the
// payload as the job argument, as in index.html. Bootstrap deploys go through a small queue,
// so a burst of fits does not mean a burst of simultaneous uploads.
//
// Result events are not trusted to be complete or unique: each galaxy is handled once, and
// when a job finishes its exec() results are reconciled against what was handled.
//
// Results come back to this tab. With ?store=<upload key>, workers instead POST them to
// tools/results-receiver.py on dcp.network and the job results are short acknowledgements; the
// page reads fits and realizations back from RESULTS_BASE. ?simulate swaps DCP for a local fake with the same interface.
// =============================================================================

const PEOPLE = 'https://people.distributive.network';
const query = new URLSearchParams(location.search);
const SIMULATE = query.has('simulate');
const onPeople = location.origin === PEOPLE;
const DATASET_BASE = new URL(query.get('dataset') || (onPeople ? '../wallaby/' : PEOPLE + '/dandesjardins/wallaby/'), location.href);
// Results live on dcp.network for now (tools/results-receiver.py behind nginx there).
const RESULTS_BASE = new URL(query.get('results') || 'https://dcp.network/3kidnas-results/');
const UPLOAD_BASE = new URL(query.get('upload') || 'https://dcp.network/3kidnas-upload/');
const TABLE_PAGE = 300;

const $ = (id) => document.getElementById(id);

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------
function el(tag, attrs, ...children) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === undefined || v === null || v === false) continue;
    if (k === 'class') node.className = v;
    else if (k === 'text') node.textContent = v;
    else if (k.startsWith('on')) node.addEventListener(k.slice(2), v);
    else node.setAttribute(k, v === true ? '' : v);
  }
  for (const c of children) if (c !== null && c !== undefined) node.append(c);
  return node;
}

function logLine(text, kind) {
  const log = $('runLog');
  log.hidden = false;
  log.append(el('div', { class: 'line' + (kind ? ' ' + kind : ''), text: `${new Date().toLocaleTimeString()}  ${text}` }));
  while (log.childElementCount > 500) log.firstChild.remove();
  log.scrollTop = log.scrollHeight;
}

function first(v) { return Array.isArray(v) ? v[0] : v; }
function fixed(v, d) { const n = Number(first(v)); return Number.isFinite(n) ? n.toFixed(d) : '–'; }

function parseResult(value) {
  if (typeof value === 'string') {
    try { return JSON.parse(value); } catch (e) { return null; }
  }
  return value && typeof value === 'object' ? value : null;
}

async function fetchJson(url) {
  const res = await fetch(url, { cache: 'no-store' });
  if (!res.ok) throw new Error(`${url} returned ${res.status}`);
  return res.json();
}

function b64ToBytes(b64) {
  const bin = atob(b64);
  const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

class Limiter {
  constructor(n) { this.n = n; this.active = 0; this.waiting = []; }
  async acquire() {
    if (this.active < this.n) { this.active++; return; }
    await new Promise((resolve) => this.waiting.push(resolve));
    this.active++;
  }
  release() {
    this.active--;
    const next = this.waiting.shift();
    if (next) next();
  }
}

// Median + 1.4825*MAD combined with the (median - best fit) offset, as index.html and
// Bootstrap_Error_Analysis.BootsrapErrors_ForParam do. The reported value stays the best fit.
// wrap: position angles are unwrapped to within 180 deg of the best fit first.
function medianMadUncertainty(values, bestFitVal, wrap) {
  if (wrap) values = values.map((v) => (v - bestFitVal > 180 ? v - 360 : v - bestFitVal < -180 ? v + 360 : v));
  const median = (arr) => { const n = arr.length; return n % 2 ? arr[(n - 1) / 2] : (arr[n / 2 - 1] + arr[n / 2]) / 2; };
  const med = median(values.slice().sort((a, b) => a - b));
  const mad = median(values.map((v) => Math.abs(v - med)).sort((a, b) => a - b));
  const diff = bestFitVal - med;
  return { median: med, err: Math.sqrt((1.4825 * mad) ** 2 + diff * diff) };
}
function fmtValErr(val, err) {
  if (!Number.isFinite(val)) return '–';
  if (!Number.isFinite(err) || err === 0) return val.toFixed(2);
  const decimals = Math.max(0, 2 - (Math.floor(Math.log10(Math.abs(err))) + 1));
  return `${val.toFixed(decimals)} ± ${err.toFixed(decimals)}`;
}

// ---------------------------------------------------------------------------
// DCP account (same storage keys as index.html, so settings carry over)
// ---------------------------------------------------------------------------
const API_KEY_STORAGE_KEY = 'galaxy-fit:apiKey';
const COMPUTE_GROUPS_STORAGE_KEY = 'galaxy-fit:computeGroups';

function stored(key, fallback) { try { return localStorage.getItem(key) ?? fallback; } catch (e) { return fallback; } }
function store(key, value) { try { localStorage.setItem(key, value); } catch (e) { /* private window */ } }

const apiKeyInput = $('apiKeyInput');
apiKeyInput.value = stored(API_KEY_STORAGE_KEY, '');
apiKeyInput.addEventListener('change', () => store(API_KEY_STORAGE_KEY, apiKeyInput.value.trim()));

let computeGroupRowEls = [];
function loadStoredComputeGroups() {
  try {
    const groups = JSON.parse(stored(COMPUTE_GROUPS_STORAGE_KEY, 'null'));
    if (Array.isArray(groups) && groups.length) return groups;
  } catch (e) { /* default below */ }
  return [{ joinKey: '', joinSecret: '' }];
}
function persistComputeGroups() {
  store(COMPUTE_GROUPS_STORAGE_KEY, JSON.stringify(computeGroupRowEls.map((r) => ({ joinKey: r.keyInput.value.trim(), joinSecret: r.secretInput.value.trim() }))));
}
function makeComputeGroupRow(joinKey, joinSecret) {
  const keyInput = el('input', { type: 'text', placeholder: 'joinKey (blank = public)' });
  keyInput.value = joinKey || '';
  const secretInput = el('input', { type: 'password', placeholder: 'joinSecret (optional)' });
  secretInput.value = joinSecret || '';
  const row = el('div', { class: 'compute-group-row' });
  const removeBtn = el('button', { type: 'button', class: 'remove-group-btn', title: 'Remove this group', text: '×' });
  keyInput.addEventListener('change', persistComputeGroups);
  secretInput.addEventListener('change', persistComputeGroups);
  removeBtn.addEventListener('click', () => {
    if (computeGroupRowEls.length > 1) { computeGroupRowEls = computeGroupRowEls.filter((r) => r.el !== row); row.remove(); }
    else { keyInput.value = ''; secretInput.value = ''; }
    persistComputeGroups();
  });
  row.append(keyInput, secretInput, removeBtn);
  return { el: row, keyInput, secretInput };
}
computeGroupRowEls = loadStoredComputeGroups().map((g) => makeComputeGroupRow(g.joinKey, g.joinSecret));
for (const r of computeGroupRowEls) $('computeGroupRows').append(r.el);
$('addComputeGroupBtn').addEventListener('click', () => {
  const r = makeComputeGroupRow('', '');
  computeGroupRowEls.push(r);
  $('computeGroupRows').append(r.el);
});
function getComputeGroups() {
  const groups = computeGroupRowEls.map((r) => ({ joinKey: r.keyInput.value.trim(), joinSecret: r.secretInput.value.trim() })).filter((g) => g.joinKey);
  if (!groups.length) return [{ joinKey: 'public' }];
  return groups.map((g) => (g.joinSecret ? g : { joinKey: g.joinKey }));
}
function getSlicePaymentOffer() {
  const v = parseFloat($('slicePaymentOfferInput').value);
  return Number.isFinite(v) && v >= 0 ? v : 1;
}

// Results: Local comes back through DCP to this tab; Remote has workers post to the
// receiver on dcp.network with an upload key minted by tools/new-upload-key.sh (ssh).
// Simulation can use Remote only with ?upload= and ?results= pointing at a local receiver.
const UPLOAD_KEY_STORAGE_KEY = 'galaxy-catalogue:uploadKey';
const RESULTS_MODE_STORAGE_KEY = 'galaxy-catalogue:resultsMode';
const REMOTE_ALLOWED = !SIMULATE || (query.has('upload') && query.has('results'));
const RESULTS_HELP = {
  local: 'Results come back to this tab. Fine for a few hundred galaxies; download them before closing.',
  remote: 'Workers post results to dcp.network, so they survive the tab and a run can resume. Needed for big runs.',
};
function remoteMode() { return $('resultsRemoteInput').checked; }
function refreshResultsControls() {
  $('resultsHelp').textContent = RESULTS_HELP[remoteMode() ? 'remote' : 'local'];
  $('uploadKeyField').hidden = !remoteMode();
}
$('uploadKeyInput').value = stored(UPLOAD_KEY_STORAGE_KEY, '');
$('uploadKeyInput').addEventListener('change', () => store(UPLOAD_KEY_STORAGE_KEY, $('uploadKeyInput').value.trim()));
$('resultsRemoteInput').disabled = !REMOTE_ALLOWED;
if (REMOTE_ALLOWED && stored(RESULTS_MODE_STORAGE_KEY, 'local') === 'remote') $('resultsRemoteInput').checked = true;
for (const input of [$('resultsLocalInput'), $('resultsRemoteInput')]) {
  input.addEventListener('change', () => {
    store(RESULTS_MODE_STORAGE_KEY, remoteMode() ? 'remote' : 'local');
    refreshResultsControls();
    refreshIndex();
    if (fanout) selectionChanged();
  });
}
refreshResultsControls();
function storageConfig() {
  if (!remoteMode()) return null;
  const key = $('uploadKeyInput').value.trim();
  if (!key) throw new Error('Results are set to Remote but the upload key is empty: run tools/new-upload-key.sh and paste it');
  return { url: new URL(key, UPLOAD_BASE), run: new Date().toISOString() };
}

// ---------------------------------------------------------------------------
// Catalogue + per-galaxy state
// ---------------------------------------------------------------------------
let galaxies = [];
const byId = new Map();
let storedIndex = null;            // index.json from the results server, or null
const states = new Map();          // id -> run state, shared with the fan-out view
const overrides = new Map();       // id -> { pa, inc } edited in the table
const selection = new Set();
const fitCache = new Map();        // id -> initial-fit report
const reportCache = new Map();     // id -> Map(k -> realization report)

function startGeometry(g) {
  const o = overrides.get(g.id) || {};
  return { pa: Number.isFinite(o.pa) ? o.pa : g.paEst, inc: Number.isFinite(o.inc) ? o.inc : g.incEst };
}

function storedCount(s) { return s ? (Array.isArray(s.bootstraps) ? s.bootstraps.length : Number(s.bootstraps) || 0) : 0; }
function storedEntry(id) { return storedIndex && storedIndex.galaxies ? storedIndex.galaxies[id] : null; }

// One status for the table, from this run if the galaxy is in it, else from storage.
function tableStatus(id) {
  const st = states.get(id);
  if (st && st.status !== 'queued') return st.status;
  const s = storedEntry(id);
  if (!s || !s.fit) return st ? 'queued' : 'none';
  if (!s.fit.FITAchieved) return 'fit-failed';
  const n = storedCount(s);
  return n === 0 ? 'fit' : 'stored:' + n;
}

async function loadCatalogue() {
  const meta = await fetchJson(new URL('catalogue.json', DATASET_BASE));
  galaxies = meta.galaxies.map((g) => Object.assign(g, { short: g.id.replace('WALLABY_', '') }));
  for (const g of galaxies) byId.set(g.id, g);
  await refreshIndex();
  const fields = [...new Set(galaxies.map((g) => g.field))];
  $('catalogueSummary').textContent = `${galaxies.length.toLocaleString()} galaxies in ${fields.length} fields `
    + `(${meta.releases.length} releases), built ${new Date(meta.created).toLocaleDateString()}. `
    + 'PA and inclination start from the catalogue (kin_pa and the beam-corrected ellipse); edit them per galaxy in the table.';
  setupSky(fields);
  applyFilters();
}

// ---------------------------------------------------------------------------
// Sky map: every galaxy, coloured by field; brush to filter
// ---------------------------------------------------------------------------
let skyBox = null;
const hiddenFields = new Set();
let skyScales = null;
let skyBrush = null;
let fieldColor = null;

function setupSky(fields) {
  fieldColor = d3.scaleOrdinal().domain(fields).range(['#6aa9ff', '#ffd166', '#43d6c4', '#ff8fab', '#a08bff', '#7ee081', '#ffa94d', '#5ce1ff']);
  const legend = $('skyLegend');
  legend.replaceChildren(...fields.map((f) => {
    const span = el('span', { title: 'click to show or hide' }, el('i'), `${f} (${galaxies.filter((g) => g.field === f).length})`);
    span.querySelector('i').style.background = fieldColor(f);
    span.addEventListener('click', () => {
      if (hiddenFields.has(f)) hiddenFields.delete(f); else hiddenFields.add(f);
      span.classList.toggle('off', hiddenFields.has(f));
      applyFilters();
    });
    return span;
  }));
  const svg = d3.select('#skySvg');
  skyBrush = d3.brush().extent(() => [[0, 0], [$('skyCanvas').clientWidth, $('skyCanvas').clientHeight]]).on('end', (ev) => {
    if (!ev.selection) { if (skyBox) { skyBox = null; applyFilters(); } return; }
    const [[x0, y0], [x1, y1]] = ev.selection;
    const { x, y } = skyScales;
    skyBox = { ra: [x.invert(x1), x.invert(x0)], dec: [y.invert(y1), y.invert(y0)] };
    applyFilters();
  });
  svg.append('g').attr('class', 'brush');
  svg.on('dblclick', () => { svg.select('.brush').call(skyBrush.move, null); });
  window.addEventListener('resize', () => drawSky());
}

function drawSky() {
  const canvas = $('skyCanvas');
  const w = canvas.clientWidth, h = canvas.clientHeight, dpr = window.devicePixelRatio || 1;
  canvas.width = w * dpr; canvas.height = h * dpr;
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.fillStyle = '#070a12';
  ctx.fillRect(0, 0, w, h);
  const ras = d3.extent(galaxies, (g) => g.ra), decs = d3.extent(galaxies, (g) => g.dec);
  const x = d3.scaleLinear().domain([ras[1] + 4, ras[0] - 4]).range([40, w - 12]); // east to the left
  const y = d3.scaleLinear().domain([decs[0] - 3, decs[1] + 3]).range([h - 22, 10]);
  skyScales = { x, y };
  ctx.strokeStyle = 'rgba(255,255,255,0.06)';
  ctx.fillStyle = 'rgba(160,172,200,0.75)';
  ctx.font = '10px ui-monospace, Menlo, monospace';
  for (const t of x.ticks(10)) { ctx.beginPath(); ctx.moveTo(x(t), 10); ctx.lineTo(x(t), h - 22); ctx.stroke(); ctx.fillText(`${t}°`, x(t) - 8, h - 8); }
  for (const t of y.ticks(5)) { ctx.beginPath(); ctx.moveTo(40, y(t)); ctx.lineTo(w - 12, y(t)); ctx.stroke(); ctx.fillText(`${t}°`, 4, y(t) + 3); }
  const matching = new Set(filtered.map((g) => g.id));
  for (const g of galaxies) {
    if (selection.has(g.id)) continue;
    const on = matching.has(g.id);
    ctx.globalAlpha = on ? 0.85 : 0.12;
    ctx.fillStyle = fieldColor(g.field);
    ctx.beginPath();
    ctx.arc(x(g.ra), y(g.dec), on ? 1.6 : 1.1, 0, 7);
    ctx.fill();
  }
  ctx.globalAlpha = 1;
  for (const id of selection) {
    const g = byId.get(id);
    ctx.fillStyle = '#ffffff';
    ctx.shadowColor = '#ffffff';
    ctx.shadowBlur = 8;
    ctx.beginPath();
    ctx.arc(x(g.ra), y(g.dec), 3, 0, 7);
    ctx.fill();
  }
  ctx.shadowBlur = 0;
  const svg = d3.select('#skySvg').attr('viewBox', `0 0 ${w} ${h}`);
  if (skyBrush && (svg.attr('data-w') !== String(w) || svg.attr('data-h') !== String(h))) {
    svg.attr('data-w', w).attr('data-h', h).select('.brush').call(skyBrush);
    svg.selectAll('.brush .selection').attr('fill', '#6aa9ff').attr('fill-opacity', 0.12).attr('stroke', '#6aa9ff');
  }
}

// ---------------------------------------------------------------------------
// Filters, sorting, table
// ---------------------------------------------------------------------------
let filtered = [];
let sortKey = 'ra', sortDir = 1, shown = TABLE_PAGE;

const COLUMNS = [
  { key: 'sel', label: '' },
  { key: 'name', label: 'Galaxy' },
  { key: 'field', label: 'Field' },
  { key: 'ra', label: 'RA', num: true },
  { key: 'dec', label: 'Dec', num: true },
  { key: 'vsys', label: 'Vsys km/s', num: true },
  { key: 'sizeBeams', label: 'Size (beams)', num: true },
  { key: 'paEst', label: 'PA°', num: true },
  { key: 'incEst', label: 'Inc°', num: true },
  { key: 'status', label: 'Status' },
];

function statusOrder(id) {
  const s = tableStatus(id);
  return { none: 0, queued: 1, fitting: 2, fit: 3, bootstrapping: 4, done: 6, 'fit-failed': -1, 'bs-failed': -1 }[s] ?? 5;
}

function applyFilters() {
  const q = $('searchInput').value.trim().toUpperCase().replace(/\s+/g, '_').replace(/^WALLABY_?/, '');
  const minSize = parseFloat($('minSizeInput').value) || 0;
  const flag0 = $('flagInput').checked;
  const status = $('statusFilter').value;
  filtered = galaxies.filter((g) => {
    if (hiddenFields.has(g.field)) return false;
    if (flag0 && g.flag !== 0) return false;
    if (g.sizeBeams < minSize) return false;
    if (q && !g.id.toUpperCase().includes(q)) return false;
    if (skyBox && !(g.ra >= skyBox.ra[0] && g.ra <= skyBox.ra[1] && g.dec >= skyBox.dec[0] && g.dec <= skyBox.dec[1])) return false;
    if (status) {
      const s = tableStatus(g.id);
      if (status === 'none' && !['none', 'queued'].includes(s)) return false;
      if (status === 'fit' && s !== 'fit') return false;
      if (status === 'partial' && !(s.startsWith('stored:') || s === 'bootstrapping')) return false;
      if (status === 'done' && s !== 'done' && !(s.startsWith('stored:') && Number(s.slice(7)) >= Number($('nBootstrapsInput').value))) return false;
    }
    return true;
  });
  sortRows();
  shown = TABLE_PAGE;
  renderTable();
  drawSky();
}

function sortRows() {
  const val = (g) => (sortKey === 'status' ? statusOrder(g.id) : sortKey === 'paEst' ? startGeometry(g).pa : sortKey === 'incEst' ? startGeometry(g).inc : g[sortKey]);
  filtered.sort((a, b) => {
    const va = val(a), vb = val(b);
    return (va < vb ? -1 : va > vb ? 1 : 0) * sortDir;
  });
}

function renderHead() {
  $('catHead').replaceChildren(...COLUMNS.map((c) => {
    const th = el('th', { text: c.label, class: sortKey === c.key ? 'sorted' : null, 'data-dir': sortDir > 0 ? '▲' : '▼' });
    if (c.num) th.style.textAlign = 'right';
    if (c.key !== 'sel') th.addEventListener('click', () => {
      sortDir = sortKey === c.key ? -sortDir : 1;
      sortKey = c.key;
      sortRows(); renderHead(); renderTable();
    });
    return th;
  }));
}

function statusCell(id) {
  const s = tableStatus(id);
  const st = states.get(id);
  let label, colorStatus = s;
  if (s === 'none') return el('span', { class: 'muted', text: '–' });
  if (s.startsWith('stored:')) { label = `stored, ${s.slice(7)} realizations`; colorStatus = 'fit'; }
  else label = Fanout.STATUS_LABEL[s] || s;
  if (st && st.nBoot && ['bootstrapping', 'done', 'bs-failed'].includes(s)) label += ` ${st.done.size + st.failed.size}/${st.nBoot}`;
  if (s === 'fit' && !st) label = 'fit stored';
  const dot = el('i');
  dot.style.background = Fanout.statusColor(colorStatus);
  return el('span', { class: 'gstatus' }, dot, label);
}

let focusedId = null;
function renderTable() {
  const body = $('catBody');
  const rows = filtered.slice(0, shown).map((g) => {
    const geo = startGeometry(g);
    const cb = el('input', { type: 'checkbox' });
    cb.checked = selection.has(g.id);
    cb.addEventListener('change', () => { if (cb.checked) selection.add(g.id); else selection.delete(g.id); selectionChanged(); });
    const geoInput = (key, value) => {
      const input = el('input', { type: 'number', step: '0.1' });
      input.value = value.toFixed(1);
      input.addEventListener('change', () => {
        const o = Object.assign({}, overrides.get(g.id));
        o[key] = parseFloat(input.value);
        overrides.set(g.id, o);
        if (selection.has(g.id)) selectionChanged();
      });
      return input;
    };
    const name = el('a', { href: '#', text: g.name, onclick: (ev) => { ev.preventDefault(); showDetail(g.id); } });
    const tr = el('tr', { class: [selection.has(g.id) ? 'selected' : '', focusedId === g.id ? 'focus' : ''].join(' ').trim() || null, 'data-id': g.id },
      el('td', {}, cb),
      el('td', {}, name),
      el('td', { text: g.field }),
      el('td', { class: 'num', text: g.ra.toFixed(3) }),
      el('td', { class: 'num', text: g.dec.toFixed(3) }),
      el('td', { class: 'num', text: g.vsys.toFixed(0) }),
      el('td', { class: 'num', text: g.sizeBeams.toFixed(1) }),
      el('td', { class: 'num' }, geoInput('pa', geo.pa)),
      el('td', { class: 'num' }, geoInput('inc', geo.inc)),
      el('td', {}, statusCell(g.id)),
    );
    return tr;
  });
  body.replaceChildren(...rows);
  $('showMoreBtn').hidden = filtered.length <= shown;
  $('showMoreBtn').textContent = `Show more (${(filtered.length - shown).toLocaleString()} hidden)`;
  $('matchCount').textContent = filtered.length.toLocaleString();
}

let tableRefreshQueued = false;
function refreshStatuses() {
  if (tableRefreshQueued) return;
  tableRefreshQueued = true;
  setTimeout(() => {
    tableRefreshQueued = false;
    for (const tr of $('catBody').children) {
      const id = tr.getAttribute('data-id');
      tr.lastChild.replaceChildren(statusCell(id));
    }
  }, 250);
}

['searchInput', 'minSizeInput', 'flagInput', 'statusFilter'].forEach((id) => $(id).addEventListener('input', applyFilters));
$('showMoreBtn').addEventListener('click', () => { shown += TABLE_PAGE; renderTable(); });
$('selectMatchingBtn').addEventListener('click', () => { for (const g of filtered) selection.add(g.id); selectionChanged(); });
$('selectRandomBtn').addEventListener('click', () => {
  const pool = filtered.filter((g) => !selection.has(g.id));
  for (let i = 0; i < 12 && pool.length; i++) selection.add(pool.splice(Math.floor(Math.random() * pool.length), 1)[0].id);
  selectionChanged();
});
$('clearSelectionBtn').addEventListener('click', () => { if (!running()) { selection.clear(); selectionChanged(); } });

// ---------------------------------------------------------------------------
// Selection -> estimate + fan-out preview
// ---------------------------------------------------------------------------
let fanout = null;
function selectedGalaxies() {
  return [...selection].map((id) => byId.get(id)).sort((a, b) => a.ra - b.ra);
}

function selectionChanged() {
  $('selectedCount').textContent = selection.size.toLocaleString();
  const n = selection.size, b = Math.max(1, parseInt($('nBootstrapsInput').value, 10) || 50);
  $('estimate').replaceChildren(
    n ? 'This deploys ' : 'Select galaxies to deploy.',
    ...(n ? [el('strong', { text: `1 initial-fit job (${n} slice${n > 1 ? 's' : ''})` }), ' and up to ',
      el('strong', { text: `${n} bootstrap job${n > 1 ? 's' : ''} (${(n * b).toLocaleString()} slices)` }),
      `. Each slice takes about 15 s on a laptop-class worker, so about ${Math.round((n * (b + 1) * 15) / 60).toLocaleString()} worker-minutes in all, `,
      `at most ${(n * (b + 1) * getSlicePaymentOffer()).toLocaleString()} compute credits at your offer, `,
      `and about ${(resultBytes(selectedGalaxies(), b) / 1e9).toFixed(1)} GB of results`,
      !remoteMode() && resultBytes(selectedGalaxies(), b) > LOCAL_RESULTS_LIMIT ? ' (too much for Local: use Remote).' : '.'] : []),
  );
  $('deployBtn').disabled = !n || running();
  $('bootstrapNowBtn').disabled = !n || running();
  if (!running()) {
    const list = selectedGalaxies();
    for (const id of [...states.keys()]) if (!selection.has(id)) states.delete(id);
    for (const g of list) if (!states.has(g.id)) states.set(g.id, newState());
    fanout.setGalaxies(list.map(vizGalaxy), states);
  }
  renderTable();
  drawSky();
}
$('nBootstrapsInput').addEventListener('change', selectionChanged);

function newState() { return { status: 'queued', nBoot: 0, done: new Set(), failed: new Set(), inFlight: 0, error: null, fitSummary: null }; }
function vizGalaxy(g) {
  const geo = startGeometry(g);
  return { id: g.id, name: g.name, short: g.short, field: g.field, vsys: g.vsys, sizeBeams: g.sizeBeams, pa: geo.pa, inc: geo.inc };
}

// ---------------------------------------------------------------------------
// Backends: real DCP, or a simulation with the same interface
// ---------------------------------------------------------------------------
let identityReady = null;
// dcp-client.js (a script tag in catalogue.html) is a loader: `dcp` appears some time after
// it runs, announced by a 'dcpclientready' event.
const dcpReady = new Promise((resolve, reject) => {
  if (typeof dcp !== 'undefined') { resolve(); return; }
  window.addEventListener('dcpclientready', () => resolve(), { once: true });
  setTimeout(() => (typeof dcp === 'undefined' ? reject(new Error('dcp-client did not load; see the browser console')) : resolve()), 30000);
});
dcpReady.catch(() => {});
function loadDcpClient() { return dcpReady; }
function ensureIdentity() {
  if (!identityReady) {
    identityReady = (async () => {
      await loadDcpClient();
      if (!dcp.identity.check()) await dcp.identity.set(apiKeyInput.value.trim());
      await dcp.wallet.add(await dcp.wallet.get());
    })().catch((error) => { identityReady = null; throw error; });
  }
  return identityReady;
}

let sofiaParTemplate = null;
async function getSofiaParTemplate() {
  if (!sofiaParTemplate) {
    const res = await fetch('sofia-template-par-file.par');
    if (!res.ok) throw new Error('could not load sofia-template-par-file.par');
    sofiaParTemplate = await res.text();
  }
  return sofiaParTemplate;
}

// Shared argument of the initial-fit job: the fitting options (built the same way as a
// single-galaxy payload, minus the per-galaxy fields the slices carry) and each galaxy's
// starting PA/Inc.
let builderSource = null;
async function sharedInitialFitArg(list) {
  const shared = BuildFitPayloads.buildInitialFitPayload({
    cubeBytes: new Uint8Array(1), maskBytes: new Uint8Array(1), paEstDeg: 0, incEstDeg: 0,
    options: DefaultFittingOptions.defaultFittingOptions(),
  });
  for (const k of ['observedCubeRawFitsB64', 'observedMaskRawFitsB64', 'paEstDeg', 'incEstDeg']) delete shared[k];
  shared.galaxies = {};
  for (const g of list) {
    const geo = startGeometry(g);
    shared.galaxies[g.id] = { paEstDeg: geo.pa, incEstDeg: geo.inc };
  }
  // Each fit slice also builds its galaxy's bootstrap payload (see runInitialFit), so this tab
  // fetches ~1 MB per galaxy instead of the fit and the cube.
  if (!builderSource) {
    const res = await fetch(document.querySelector('script[src*="buildFitPayloads.js"]').src);
    if (!res.ok) throw new Error('could not load buildFitPayloads.js');
    builderSource = await res.text();
  }
  shared.bootstrapBuilder = { source: builderSource, sofiaParTemplate: await getSofiaParTemplate(), options: DefaultFittingOptions.defaultFittingOptions() };
  return shared;
}

function wireJob(job, h) {
  job.on('readystatechange', (s) => h.state && h.state(s));
  job.on('accepted', () => h.accepted && h.accepted(job.id));
  job.on('status', (ev) => h.status && h.status(ev));
  job.on('result', (ev) => h.result(parseResult(ev.result), ev.sliceNumber));
  job.on('error', (ev) => h.sliceError && h.sliceError(ev));
  job.on('nofunds', () => h.log && h.log('no funds: top up the account or lower the slice offer', 'err'));
}

const RealBackend = {
  name: 'DCP',
  async deployInitialFits(list, h, storage) {
    await ensureIdentity();
    // A RemoteDataSet makes workers fetch each URL as the slice datum; a plain array of URL
    // objects from the browser client is uploaded by value instead.
    const urls = new dcp.compute.RemoteDataSet(list.map((g) => new URL(g.file, DATASET_BASE).href));
    const job = dcp.compute.for(urls, runInitialFit, [await sharedInitialFitArg(list)]);
    job.requires(window.INITIAL_FIT_ONLY_MODULES);
    job.computeGroups = getComputeGroups();
    job.public = { name: `🌌 WALLABY initial fits (${list.length})`, description: '3KIDNAS catalogue run', link: 'https://www.candiapl.ca/' };
    if (storage) job.setResultStorage(storage.url, { kind: 'fit', run: storage.run, contentType: 'application/json' });
    wireJob(job, h);
    return { cancel: () => job.cancel(), done: job.exec(getSlicePaymentOffer()) };
  },
  async deployBootstrap(g, payload, indices, h, storage) {
    await ensureIdentity();
    const job = dcp.compute.for(indices, runBootstrapRealization, [payload]);
    job.requires(window.BOOTSTRAP_ONLY_MODULES);
    job.computeGroups = getComputeGroups();
    job.public = { name: `🌌 ${g.name} (bootstrap)`, description: '3KIDNAS catalogue run', link: 'https://www.candiapl.ca/' };
    if (storage) job.setResultStorage(storage.url, { kind: 'bootstrap', galaxyId: g.id, run: storage.run, contentType: 'application/json' });
    wireJob(job, h);
    return { cancel: () => job.cancel(), done: job.exec(getSlicePaymentOffer()) };
  },
};

// Same events, synthetic values: fits converge around the catalogue estimates, and
// realizations scatter around the fit.
const SimBackend = (() => {
  const timers = new Set();
  const later = (ms, fn) => { const t = setTimeout(() => { timers.delete(t); fn(); }, ms); timers.add(t); return t; };
  const gauss = () => { let u = 0, v = 0; while (!u) u = Math.random(); while (!v) v = Math.random(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); };
  function fakeFit(g) {
    const geo = startGeometry(g);
    const nR = Math.max(3, Math.round(g.sizeBeams * 2));
    const vmax = 50 + Math.random() * 190, rt = 10 + Math.random() * 25;
    const R = Array.from({ length: nR }, (_, i) => (i + 0.5) * 15);
    return {
      galaxyId: g.id, FITAchieved: Math.random() > 0.08,
      INCLINATION: [Math.min(88, Math.max(20, geo.inc + gauss() * 6))], POSITIONANGLE: [(geo.pa + gauss() * 4 + 360) % 360],
      VSYS: [g.vsys + gauss() * 4], VDISP: [8 + Math.random() * 6], XCENTER: [20 + gauss()], YCENTER: [20 + gauss()],
      R, VROT: R.map((r) => vmax * (1 - Math.exp(-r / rt))),
    };
  }
  function fakeRealization(fit, k) {
    const s = 1 + gauss() * 0.06;
    return {
      realizationIndex: k, sofiaFailed: Math.random() < 0.04,
      INCLINATION: [fit.INCLINATION[0] + gauss() * 3.5], POSITIONANGLE: [fit.POSITIONANGLE[0] + gauss() * 2.5],
      VSYS: [fit.VSYS[0] + gauss() * 2], VDISP: [fit.VDISP[0] + gauss()], XCENTER: [fit.XCENTER[0] + gauss() * 0.3], YCENTER: [fit.YCENTER[0] + gauss() * 0.3],
      R: fit.R, VROT: fit.VROT.map((v) => v * s * (1 + gauss() * 0.03)),
    };
  }
  // ?simSpeed=0.01 runs the simulation 100 times faster, for full-catalogue dry runs.
  const speed = Math.max(0.001, Number(query.get('simSpeed')) || 1);
  function runJob(n, makeResult, h, perSlice, workers) {
    let cancelled = false, nextSlice = 0, computed = 0;
    const results = new Array(n);
    let resolveDone, rejectDone;
    const done = new Promise((res, rej) => { resolveDone = res; rejectDone = rej; });
    later((700 + Math.random() * 900) * Math.min(1, speed * 10), () => {
      if (cancelled) return;
      h.accepted && h.accepted('sim-' + Math.random().toString(36).slice(2, 10));
      const startOne = () => {
        if (cancelled || nextSlice >= n) return;
        const i = nextSlice++;
        h.status && h.status({ total: n, distributed: nextSlice, computed });
        later(perSlice(), async () => {
          if (cancelled) return;
          try { results[i] = await makeResult(i); } catch (e) { results[i] = { error: String(e.message || e) }; }
          if (cancelled) return;
          computed++;
          h.status && h.status({ total: n, distributed: nextSlice, computed });
          h.result(results[i], i + 1);
          if (computed === n) resolveDone(results);
          else startOne();
        });
      };
      for (let w = 0; w < workers; w++) startOne();
    });
    return { cancel: () => { cancelled = true; rejectDone(new Error('cancelled')); }, done };
  }
  // What a worker does with setResultStorage: POST the result as a form, and the reply
  // becomes the slice result.
  async function post(storage, params, value) {
    const body = new URLSearchParams(Object.assign({}, params, { contentType: 'application/json', content: JSON.stringify(value) }));
    const res = await fetch(storage.url, { method: 'POST', body });
    return res.json();
  }
  return {
    name: 'simulation',
    async deployInitialFits(list, h, storage) {
      return runJob(list.length, (i) => {
        const f = fakeFit(list[i]);
        if (storage) return post(storage, { kind: 'fit', run: storage.run }, f);
        fitCache.set(list[i].id, f);
        return f;
      }, h, () => (2500 + Math.random() * 9000) * speed, Math.min(list.length, 16 / Math.sqrt(speed)));
    },
    async deployBootstrap(g, payload, indices, h, storage) {
      return runJob(indices.length, (i) => {
        const r = fakeRealization(payload.simFit, indices[i]);
        return storage ? post(storage, { kind: 'bootstrap', galaxyId: g.id, run: storage.run }, r) : r;
      }, h, () => (1200 + Math.random() * 5000) * speed, 10);
    },
    stopAll() { for (const t of timers) clearTimeout(t); timers.clear(); },
  };
})();

const backend = SIMULATE ? SimBackend : RealBackend;

// ---------------------------------------------------------------------------
// The run
// ---------------------------------------------------------------------------
let run = null;
function running() { return Boolean(run && !run.finished); }

function setRunButtons() {
  $('deployBtn').disabled = running() || !selection.size;
  $('bootstrapNowBtn').disabled = running() || !selection.size;
  $('stopBtn').disabled = !running();
  $('clearSelectionBtn').disabled = running();
  $('downloadBtn').disabled = !run;
}

function touch(id) { fanout.update(); refreshStatuses(); if (focusedId === id) scheduleDetail(); }

function startRun(list) {
  run = {
    list, finished: false, jobs: new Set(), storage: storageConfig(),
    nBoot: Math.max(1, parseInt($('nBootstrapsInput').value, 10) || 50),
    auto: $('autoBootstrapInput').checked,
    limiter: new Limiter(Math.max(1, parseInt($('deployConcurrencyInput').value, 10) || 3)),
    jobSlots: new Limiter(Math.max(1, parseInt($('maxJobsInput').value, 10) || 100)),
    handledFits: new Set(), bootstrapsStarted: new Set(), pending: 0,
  };
  for (const g of list) states.set(g.id, newState());
  holdTabAwake(run);
  fanout.setGalaxies(list.map(vizGalaxy), states);
  fanout.resetClock();
  fanout.setRunning(true);
  setRunButtons();
  logLine(`${backend.name}: ${list.length} galaxies, ${run.nBoot} realizations each, results ${run.storage ? 'stored on this server' : 'returned to this tab'}`);
}

function trackJob(handle) {
  run.jobs.add(handle);
  run.pending++;
  handle.done.catch(() => {}).finally(() => {
    run.jobs.delete(handle);
    run.pending--;
    maybeFinish();
  });
}

// A long run lives in this tab. A held Web Lock keeps Chrome from freezing or discarding it in
// the background, and the wake lock keeps the screen (and so the machine) awake while visible.
function holdTabAwake(r) {
  const released = new Promise((resolve) => { r.releaseAwake = resolve; });
  try { navigator.locks && navigator.locks.request('3kidnas-run', () => released); } catch (e) { /* unsupported */ }
  const wake = () => { try { navigator.wakeLock && navigator.wakeLock.request('screen').then((l) => { r.wakeLock = l; }, () => {}); } catch (e) { /* unsupported */ } };
  wake();
  document.addEventListener('visibilitychange', () => { if (!document.hidden && running()) wake(); });
}
window.addEventListener('beforeunload', (ev) => { if (running()) { ev.preventDefault(); ev.returnValue = ''; } });

function maybeFinish() {
  if (!run || run.finished || run.pending > 0) return;
  run.finished = true;
  if (run.releaseAwake) run.releaseAwake();
  if (run.wakeLock) run.wakeLock.release().catch(() => {});
  fanout.setRunning(false);
  setRunButtons();
  const counts = {};
  for (const g of run.list) { const s = states.get(g.id).status; counts[s] = (counts[s] || 0) + 1; }
  logLine('run finished: ' + Object.entries(counts).map(([s, n]) => `${n} ${Fanout.STATUS_LABEL[s] || s}`).join(', '), 'ok');
  refreshIndex();
}

async function refreshIndex() {
  if (!remoteMode()) { storedIndex = null; refreshStatuses(); return; }
  try { storedIndex = await fetchJson(new URL('index.json', RESULTS_BASE)); } catch (e) { storedIndex = null; }
  refreshStatuses();
}

// Rough bytes of results a selection produces: a fit is ~3.5x its galaxy file (it carries the
// model cube), a realization ~1.3 kB.
function resultBytes(list, nBoot) { return list.reduce((sum, g) => sum + 3.5 * (g.bytes || 1e6) + nBoot * 1300, 0); }
const LOCAL_RESULTS_LIMIT = 1.5e9;
function checkLocalSize(list, nBoot) {
  if (remoteMode() || SIMULATE) return;
  const bytes = resultBytes(list, nBoot);
  if (bytes > LOCAL_RESULTS_LIMIT) {
    throw new Error(`about ${(bytes / 1e9).toFixed(1)} GB of results is too much to hold in this tab: switch Results to Remote`);
  }
}

async function deploy() {
  const list = selectedGalaxies();
  if (!list.length || running()) return;
  const nBoot = Math.max(1, parseInt($('nBootstrapsInput').value, 10) || 50);
  checkLocalSize(list, nBoot);
  // Deploy always starts over: each new fit replaces the stored one and deletes its old
  // realizations, so a rerun overwrites rather than adds. Bootstrap fitted galaxies is resume.
  startRun(list);
  await deployFits(list);
  touch();
}

async function deployFits(list) {
  for (const g of list) states.get(g.id).status = 'fitting';
  fanout.update();
  fanout.dispatched();
  const t0 = Date.now();
  try {
    const handle = await backend.deployInitialFits(list, {
      state: (s) => logLine('initial fits: ' + s),
      accepted: (id) => logLine(`initial-fit job ${id} accepted (${list.length} slices)`, 'ok'),
      status: () => {},
      result: (value, sliceNumber) => onFitResult(value, sliceNumber, list),
      sliceError: (ev) => {
        const g = ev && ev.sliceNumber ? list[ev.sliceNumber - 1] : null;
        logLine(`initial fit slice ${ev && ev.sliceNumber}${g ? ' (' + g.name + ')' : ''}: ${(ev && ev.message) || ev}`, 'err');
      },
      log: logLine,
    }, run.storage);
    trackJob(handle);
    handle.done.then((results) => {
      logLine(`initial-fit job complete in ${Math.round((Date.now() - t0) / 1000)} s`, 'ok');
      // Reconcile: anything a missed or undecodable result event left behind.
      const values = results && typeof results.values === 'function' ? results.values() : Array.from(results || []);
      values.forEach((v, i) => onFitResult(parseResult(v), i + 1, list));
      for (const g of list) {
        const st = states.get(g.id);
        if (st.status === 'fitting') { st.status = 'fit-failed'; st.error = st.error || 'no result returned'; }
      }
      touch();
    }, (error) => {
      logLine('initial-fit job ended: ' + (error && error.message ? error.message : error), 'err');
      for (const g of list) { const st = states.get(g.id); if (st.status === 'fitting') st.status = 'stopped'; }
      touch();
    });
  } catch (error) {
    logLine('could not deploy the initial fits: ' + (error.message || error), 'err');
    for (const g of list) states.get(g.id).status = 'stopped';
    maybeFinish();
    touch();
  }
}

function summaryOf(fit) {
  return { pa: fixed(fit.POSITIONANGLE, 1), inc: fixed(fit.INCLINATION, 1), vsys: fixed(fit.VSYS, 0) };
}

function onFitResult(value, sliceNumber, list) {
  if (!value) return;
  const id = value.galaxyId || (list[sliceNumber - 1] && list[sliceNumber - 1].id);
  if (!id || !states.has(id) || run.handledFits.has(id)) return;
  run.handledFits.add(id);
  const g = byId.get(id);
  const st = states.get(id);
  const ok = Boolean(value.FITAchieved);
  // Stored fits are fetched again when needed: holding thousands of them would exhaust the tab.
  if (value.stored) st.fitUrl = value.url;
  else fitCache.set(id, value);
  st.status = ok ? 'fit' : 'fit-failed';
  if (!ok) st.error = value.error || 'did not converge';
  if (ok) st.fitSummary = summaryOf(value);  // the stored ack carries the same three values
  fanout.fitLanded(id, ok);
  logLine(`${g.name}: ${ok ? 'fit converged' : 'fit failed: ' + String(st.error).split('\n')[0]}`, ok ? 'ok' : 'err');
  if (ok && run.auto) enqueueBootstrap(g);
  touch(id);
}

async function getFit(id) {
  if (fitCache.has(id)) return fitCache.get(id);
  const st = states.get(id);
  const url = new URL((st && st.fitUrl) || `galaxies/${id}/initial-fit.json`, RESULTS_BASE);
  return fetchJson(url);
}

// The bootstrap job's argument. Fit workers build it (see sharedInitialFitArg); fits stored
// before that, or ones whose build failed, are built here from the fit and the cube instead.
async function getBootstrapPayload(g) {
  if (SIMULATE) return { simFit: await getFit(g.id) };
  const local = fitCache.get(g.id);
  if (local && local.bootstrapPayload) return local.bootstrapPayload;
  if (!local) {
    try { return await fetchJson(new URL(`galaxies/${g.id}/bootstrap-payload.json`, RESULTS_BASE)); } catch (e) { /* build it below */ }
  }
  const fit = local || await getFit(g.id);
  const datum = await fetchJson(new URL(g.file, DATASET_BASE));
  return BuildFitPayloads.buildBootstrapPayload({
    initialFitResult: fit, cubeBytes: b64ToBytes(datum.observedCubeRawFitsB64),
    sofiaParTemplateText: await getSofiaParTemplate(), bootstrapSeed: 0,
    options: DefaultFittingOptions.defaultFittingOptions(),
  });
}

function enqueueBootstrap(g) {
  if (run.bootstrapsStarted.has(g.id)) return;
  run.bootstrapsStarted.add(g.id);
  const st = states.get(g.id);
  run.pending++;  // held until this bootstrap's job is deployed (or fails to)
  (async () => {
    await run.jobSlots.acquire();
    if (run.finished || st.status === 'stopped') { run.jobSlots.release(); run.pending--; maybeFinish(); return; }
    await run.limiter.acquire();
    let released = false, slotFreed = false;
    const release = () => { if (!released) { released = true; run.limiter.release(); } };
    const freeSlot = () => { if (!slotFreed) { slotFreed = true; run.jobSlots.release(); } };
    try {
      const payload = await getBootstrapPayload(g);
      if (!st.fitSummary) { const s = storedEntry(g.id); if (s && s.fit) st.fitSummary = summaryOf(s.fit); }
      st.nBoot = run.nBoot;
      st.done = new Set(); st.failed = new Set(); st.notes = new Map(); st.inFlight = 0;
      // Resume: realizations already stored are not run again.
      const have = run.resume && run.storage && storedCount(storedEntry(g.id)) ? await storedRealizations(g.id) : new Map();
      for (const [k, r] of have) if (k < run.nBoot) (r.sofiaFailed ? st.failed : st.done).add(k);
      const indices = Array.from({ length: run.nBoot }, (_, k) => k).filter((k) => !have.has(k));
      st.status = 'bootstrapping';
      reportCache.set(g.id, new Map());
      fanout.relayout(g.id);
      touch(g.id);
      const handle = await backend.deployBootstrap(g, payload, indices, {
        accepted: (id) => { logLine(`${g.name}: bootstrap job ${id} accepted`); release(); },
        status: (ev) => {
          if (ev && Number.isFinite(ev.distributed) && Number.isFinite(ev.computed)) st.inFlight = Math.max(0, ev.distributed - ev.computed);
          fanout.update();
        },
        result: (value, sliceNumber) => onRealization(g, value, sliceNumber),
        sliceError: (ev) => logLine(`${g.name}: realization ${ev && ev.sliceNumber}: ${(ev && ev.message) || ev}`, 'err'),
        log: logLine,
      }, run.storage);
      trackJob(handle);
      handle.done.then((results) => {
        const values = results && typeof results.values === 'function' ? results.values() : Array.from(results || []);
        values.forEach((v, i) => onRealization(g, parseResult(v), i + 1));
        st.inFlight = 0;
        st.status = st.done.size ? 'done' : 'bs-failed';
        fanout.galaxyDone(g.id);
        logLine(`${g.name}: ${st.done.size}/${st.nBoot} realizations usable` + (st.failed.size ? `, ${st.failed.size} failed` : ''), st.done.size ? 'ok' : 'err');
        touch(g.id);
      }, (error) => {
        st.inFlight = 0;
        st.status = 'stopped';
        logLine(`${g.name}: bootstrap job ended: ${(error && error.message) || error}`, 'err');
        touch(g.id);
      }).finally(() => { release(); freeSlot(); });
    } catch (error) {
      st.status = 'bs-failed';
      st.error = error.message || String(error);
      logLine(`${g.name}: could not deploy the bootstrap: ${st.error}`, 'err');
      release();
      freeSlot();
      touch(g.id);
    } finally {
      run.pending--;
      maybeFinish();
    }
  })();
}

function onRealization(g, value, sliceNumber) {
  if (!value) return;
  const k = Number.isInteger(value.realizationIndex) ? value.realizationIndex : sliceNumber - 1;
  const st = states.get(g.id);
  if (!st || !Number.isInteger(k) || k < 0 || k >= st.nBoot || st.done.has(k) || st.failed.has(k)) return;
  // A fit that hit the iteration limit still yields a model, which 3KIDNAS keeps; only a
  // realization with no model (SoFiA found nothing, or an error) is a failure.
  const ok = !value.sofiaFailed;
  (ok ? st.done : st.failed).add(k);
  const note = value.error ? 'error: ' + String(value.error).split('\n')[0].slice(0, 140)
    : !ok ? 'SoFiA found no source in this noisy cube'
    : value.converged === false ? 'fit hit the iteration limit (kept, as in 3KIDNAS)' : '';
  if (note && st.notes) st.notes.set(k, note);
  if (!value.stored) reportCache.get(g.id).set(k, value);
  fanout.realizationLanded(g.id, k, ok);
  touch(g.id);
}

async function bootstrapFitted() {
  if (remoteMode()) await refreshIndex();
  const list = selectedGalaxies().filter((g) => {
    const s = storedEntry(g.id);
    const f = fitCache.get(g.id);
    return f ? Boolean(f.FITAchieved) : Boolean(s && s.fit && s.fit.FITAchieved);
  });
  if (!list.length) { logLine('none of the selected galaxies has a converged fit to bootstrap from', 'err'); return; }
  checkLocalSize(list, Math.max(1, parseInt($('nBootstrapsInput').value, 10) || 50));
  startRun(list);
  let complete = 0;
  run.resume = true;
  for (const g of list) {
    const st = states.get(g.id);
    run.handledFits.add(g.id);
    if (run.storage && storedCount(storedEntry(g.id)) >= run.nBoot) {
      st.status = 'done'; st.nBoot = run.nBoot; complete++;
      st.done = new Set(Array.from({ length: run.nBoot }, (_, k) => k)); st.failed = new Set(); st.notes = new Map();
      continue;
    }
    st.status = 'fit';
    enqueueBootstrap(g);
  }
  if (complete) logLine(`${complete} galaxies already have ${run.nBoot} stored realizations; the rest get only the ones they are missing`);
  maybeFinish();
  touch();
}

function stopAll() {
  if (!run) return;
  logLine('stopping every job in this run…', 'err');
  for (const h of run.jobs) { try { h.cancel(); } catch (e) { /* already finished */ } }
  if (SIMULATE) SimBackend.stopAll();
  for (const g of run.list) {
    const st = states.get(g.id);
    if (['queued', 'fitting', 'bootstrapping'].includes(st.status)) st.status = 'stopped';
  }
  run.pending = 0;
  maybeFinish();
  touch();
}

$('deployBtn').addEventListener('click', () => deploy().catch((e) => logLine(e.message || String(e), 'err')));
$('bootstrapNowBtn').addEventListener('click', () => bootstrapFitted().catch((e) => logLine(e.message || String(e), 'err')));
$('stopBtn').addEventListener('click', stopAll);

// Everything this tab holds for the run's galaxies: the only copy when results are not stored.
async function downloadResults() {
  if (!run) return;
  const out = { created: new Date().toISOString(), dataset: DATASET_BASE.href, galaxies: {} };
  if (run.storage) {
    // Remote: everything is already on the results host; pulling it through this tab would be
    // tens of GB for a big run. Save where it is instead.
    out.results = RESULTS_BASE.href;
    for (const g of run.list) {
      const st = states.get(g.id);
      out.galaxies[g.id] = {
        status: st.status, error: st.error || null, realizations: st.done ? st.done.size : 0, failed: st.failed ? st.failed.size : 0,
        fit: `${RESULTS_BASE.href}galaxies/${g.id}/initial-fit.json`, bootstrap: `${RESULTS_BASE.href}galaxies/${g.id}/bootstrap.ndjson`,
      };
    }
    return saveJson(out, 'wallaby-run-manifest');
  }
  for (const g of run.list) {
    const st = states.get(g.id);
    let fit = fitCache.get(g.id) || null;
    if (!fit && run.storage && st.fitUrl) { try { fit = await getFit(g.id); } catch (e) { fit = null; } }
    out.galaxies[g.id] = {
      status: st.status, start: startGeometry(g), error: st.error || null,
      initialFit: fit, realizations: run.storage ? await realizationReports(g.id) : [...(reportCache.get(g.id) || new Map()).values()],
    };
  }
  saveJson(out, 'wallaby-fits');
}
function saveJson(out, name) {
  const blob = new Blob([JSON.stringify(out)], { type: 'application/json' });
  const a = el('a', { href: URL.createObjectURL(blob), download: `${name}-${new Date().toISOString().slice(0, 19).replace(/:/g, '')}.json` });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 10000);
}
$('downloadBtn').addEventListener('click', () => downloadResults().catch((e) => logLine('download failed: ' + (e.message || e), 'err')));

// dcp-client can reject inside its own event handling (for example a result it cannot
// decode); log it and keep the run going rather than leave an uncaught error.
window.addEventListener('unhandledrejection', (ev) => {
  const reason = ev.reason;
  logLine('unhandled: ' + ((reason && reason.message) || reason), 'err');
  ev.preventDefault();
});

// ---------------------------------------------------------------------------
// Detail panel
// ---------------------------------------------------------------------------
let detailQueued = false;
function scheduleDetail() {
  if (detailQueued) return;
  detailQueued = true;
  setTimeout(() => { detailQueued = false; if (focusedId) renderDetail(focusedId); }, 600);
}

async function realizationReports(id) {
  const cache = reportCache.get(id) || new Map();
  reportCache.set(id, cache);
  const st = states.get(id);
  const s = storedEntry(id);
  if (remoteMode() && (storedCount(s) || (st && st.done && st.done.size))) return [...(await storedRealizations(id)).values()];
  return [...cache.values()];
}

// One line per stored realization; a retried slice can repeat an index, and the last wins.
async function storedRealizations(id) {
  const out = new Map();
  try {
    const res = await fetch(new URL(`galaxies/${id}/bootstrap.ndjson`, RESULTS_BASE), { cache: 'no-store' });
    if (!res.ok) return out;
    for (const line of (await res.text()).split('\n')) {
      try { const r = JSON.parse(line); if (Number.isInteger(r.realizationIndex)) out.set(r.realizationIndex, r); } catch (e) { /* partial line */ }
    }
  } catch (e) { /* not reachable */ }
  return out;
}

function showDetail(id) {
  focusedId = id;
  fanout.select(id);
  for (const tr of $('catBody').children) tr.classList.toggle('focus', tr.getAttribute('data-id') === id);
  renderDetail(id);
  $('detail').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
}

async function renderDetail(id) {
  const g = byId.get(id);
  const st = states.get(id);
  const geo = startGeometry(g);
  let fit = fitCache.get(id) || null;
  const s = storedEntry(id);
  if (!fit && ((st && st.fitUrl) || (s && s.fit)) && !SIMULATE) { try { fit = await getFit(id); } catch (e) { fit = null; } }
  const reports = await realizationReports(id);
  if (focusedId !== id) return;
  // Bootstrap_Error_Analysis only sees realizations that produced a model.
  const ok = reports.filter((r) => r && !r.sofiaFailed && r.FITAchieved !== false);
  const row = (k, v) => el('tr', {}, el('td', { text: k }), el('td', { text: v }));
  const param = (label, key, d) => {
    const best = fit ? Number(first(fit[key])) : NaN;
    if (ok.length >= 3 && Number.isFinite(best)) {
      const values = ok.map((r) => Number(first(r[key]))).filter(Number.isFinite);
      const { err } = medianMadUncertainty(values, best, key === 'POSITIONANGLE');
      return row(label, fmtValErr(best, err));
    }
    return row(label, Number.isFinite(best) ? best.toFixed(d) : '–');
  };
  const status = tableStatus(id);
  const rows = [
    row('Status', Fanout.STATUS_LABEL[status] || status.replace('stored:', 'stored, ') ),
    row('Start PA / inc (deg)', `${geo.pa.toFixed(1)} / ${geo.inc.toFixed(1)}`),
    param('Position angle (deg)', 'POSITIONANGLE', 1),
    param('Inclination (deg)', 'INCLINATION', 1),
    param('Vsys (km/s)', 'VSYS', 1),
    param('Vdisp (km/s)', 'VDISP', 1),
    row('Realizations', reports.length ? `${ok.length} usable of ${reports.length}` : '–'),
  ];
  if (st && st.error) rows.push(row('Error', String(st.error).split('\n')[0].slice(0, 120)));
  const links = [];
  if (!SIMULATE && (fit || (s && s.fit))) links.push(el('a', { href: new URL(`galaxies/${id}/initial-fit.json`, RESULTS_BASE).href, target: '_blank', rel: 'noopener', text: 'fit JSON' }));
  links.push(el('a', { href: new URL(g.file, DATASET_BASE).href, target: '_blank', rel: 'noopener', text: 'cube + mask JSON' }));
  const linkLine = el('div', { class: 'dsub', style: 'margin-top:8px' });
  links.forEach((a, i) => { if (i) linkLine.append(' · '); a.style.color = 'inherit'; linkLine.append(a); });
  $('detail').replaceChildren(
    el('h3', { text: 'Galaxy' }),
    el('div', { class: 'dname', text: g.name }),
    el('div', { class: 'dsub', text: `${g.field} · ${g.release} · RA ${g.ra.toFixed(3)} Dec ${g.dec.toFixed(3)} · ${g.sizeBeams} beams` }),
    el('table', {}, el('tbody', {}, ...rows)),
    linkLine,
  );
  $('fanoutFoot').classList.add('with-plots');
  // Per-ring VROT_ERR, as BootsrapErrors_ForParam computes it for the error bars.
  const vrotErr = fit && Array.isArray(fit.VROT) && ok.length >= 3 ? fit.VROT.map((v, i) => {
    const values = ok.map((r) => Number((r.VROT || [])[i])).filter(Number.isFinite);
    return values.length >= 3 ? medianMadUncertainty(values, Number(v)).err : NaN;
  }) : [];
  Fanout.drawRotationCurves($('rotCurveSvg'), fit, ok, vrotErr);
  Fanout.drawGeoCloud($('geoCloudSvg'), fit, reports);
}

// ---------------------------------------------------------------------------
// Boot
// ---------------------------------------------------------------------------
fanout = Fanout.create({
  svgEl: $('fanoutSvg'), stageEl: $('fanoutStage'), tipEl: $('fanoutTip'), emptyEl: $('fanoutEmpty'),
  statsEl: $('stats'), legendEl: $('fanoutLegend'),
  onSelect: (id) => showDetail(id),
});
$('fullscreenBtn').addEventListener('click', () => {
  const card = $('fullscreenBtn').closest('.fanout-card');
  if (document.fullscreenElement) document.exitFullscreen();
  else card.requestFullscreen().catch((e) => logLine('fullscreen refused: ' + (e.message || e), 'err'));
});
$('simBanner').hidden = !SIMULATE;
renderHead();
loadCatalogue().then(() => selectionChanged()).catch((e) => {
  $('catalogueSummary').textContent = 'Could not load the catalogue: ' + (e.message || e);
});
