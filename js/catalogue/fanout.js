'use strict';
// =============================================================================
// fanout.js -- d3 views for catalogue.html.
//
//   Fanout.create(...)        the radial "fan-out of fan-outs": the initial-fit job at the
//                             centre, one galaxy per slice on a ring (drawn at its PA and
//                             inclination), and each galaxy's bootstrap realizations fanned
//                             out beyond it. Particles travel inward as results land.
//   Fanout.drawRotationCurves the selected galaxy's best-fit rotation curve over every
//                             realization's
//   Fanout.drawGeoCloud       inclination vs position angle across realizations
//
// All text that comes from results (worker-chosen values) is set with d3's .text(), never
// as markup.
// =============================================================================

const Fanout = (() => {
  const TAU = Math.PI * 2;
  const STATUS_VAR = {
    queued: '--s-queued', fitting: '--s-fitting', fit: '--s-fit', 'fit-failed': '--s-failed',
    bootstrapping: '--s-boot', done: '--s-done', 'bs-failed': '--s-failed', stopped: '--s-queued',
  };
  const STATUS_LABEL = {
    queued: 'queued', fitting: 'fitting', fit: 'fitted', 'fit-failed': 'fit failed',
    bootstrapping: 'bootstrapping', done: 'done', 'bs-failed': 'bootstrap failed', stopped: 'stopped',
  };

  function cssColor(name) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim() || '#888';
  }
  function statusColor(status) { return cssColor(STATUS_VAR[status] || '--s-queued'); }

  // Angle 0 points up, increasing clockwise.
  function polar(angle, radius) { return [radius * Math.sin(angle), -radius * Math.cos(angle)]; }

  const BIG_RUN = 300;

  function create({ svgEl, stageEl, tipEl, emptyEl, statsEl, legendEl, onSelect }) {
    const svg = d3.select(svgEl);
    const defs = svg.append('defs');
    const sky = svg.append('g').attr('class', 'sky');
    const root = svg.append('g').attr('class', 'root');
    const gSpokes = root.append('g');
    const gFanLinks = root.append('g');
    const gDots = root.append('g');
    const gArcs = root.append('g');
    const gGalaxies = root.append('g');
    const gHub = root.append('g');
    const gParticles = root.append('g');
    const tip = d3.select(tipEl);

    let W = 0, H = 0, M = 0, R1 = 0, RB0 = 0, RBmax = 0;
    let galaxies = [], byId = new Map(), states = new Map();
    let selected = null, nodeSize = 8, arcMode = false;
    let startedAt = null, running = false;
    let zoomK = 1;

    // ---------------------------------------------------------------- defs
    function glowFilter(id, std) {
      const f = defs.append('filter').attr('id', id).attr('x', '-100%').attr('y', '-100%').attr('width', '300%').attr('height', '300%');
      f.append('feGaussianBlur').attr('stdDeviation', std).attr('result', 'b');
      const m = f.append('feMerge');
      m.append('feMergeNode').attr('in', 'b');
      m.append('feMergeNode').attr('in', 'SourceGraphic');
    }
    glowFilter('fo-glow', 2.2);
    glowFilter('fo-glow-strong', 4.5);
    for (const status of Object.keys(STATUS_VAR)) {
      const g = defs.append('radialGradient').attr('id', 'fo-gal-' + status);
      g.append('stop').attr('offset', '0%').attr('stop-color', '#ffffff').attr('stop-opacity', 0.95);
      g.append('stop').attr('offset', '28%').attr('stop-color', statusColor(status)).attr('stop-opacity', 0.95);
      g.append('stop').attr('offset', '100%').attr('stop-color', statusColor(status)).attr('stop-opacity', 0);
    }
    const hubGrad = defs.append('radialGradient').attr('id', 'fo-hub');
    hubGrad.append('stop').attr('offset', '0%').attr('stop-color', '#ffffff');
    hubGrad.append('stop').attr('offset', '35%').attr('stop-color', cssColor('--s-fitting'));
    hubGrad.append('stop').attr('offset', '100%').attr('stop-color', cssColor('--s-fitting')).attr('stop-opacity', 0);

    // ---------------------------------------------------------------- zoom
    const zoom = d3.zoom().scaleExtent([0.6, 10]).on('zoom', (ev) => {
      zoomK = ev.transform.k;
      root.attr('transform', `translate(${W / 2 + ev.transform.x},${H / 2 + ev.transform.y}) scale(${ev.transform.k})`);
      // Distant stars move less: a little parallax.
      const p = Math.pow(ev.transform.k, 0.25);
      sky.attr('transform', `translate(${W / 2 + ev.transform.x * 0.15},${H / 2 + ev.transform.y * 0.15}) scale(${p})`);
    });
    svg.call(zoom).on('dblclick.zoom', null);

    // ---------------------------------------------------------------- layout
    function resize() {
      W = stageEl.clientWidth || 900;
      // Fit the viewport, leaving room for the card's header and legend; fullscreen gets the lot.
      const head = stageEl.parentElement.querySelector('.fanout-head');
      const chrome = (head ? head.offsetHeight : 0) + legendEl.offsetHeight + 8;
      const avail = window.innerHeight - chrome - (document.fullscreenElement ? 0 : 32);
      H = Math.max(360, Math.min(Math.round(W * 0.68), avail));
      svg.attr('height', H).attr('viewBox', `0 0 ${W} ${H}`);
      M = Math.min(W, H) / 2;
      R1 = M * 0.44;
      RB0 = R1 + M * 0.09;
      RBmax = M * 0.97;
      root.attr('transform', `translate(${W / 2},${H / 2})`);
      sky.attr('transform', `translate(${W / 2},${H / 2})`);
      drawSky();
      layout();
      render();
    }

    function drawSky() {
      sky.selectAll('*').remove();
      const rnd = d3.randomLcg(42);
      const n = Math.round((W * H) / 2600);
      const span = Math.max(W, H);
      for (let i = 0; i < n; i++) {
        const r = rnd() < 0.92 ? 0.35 + rnd() * 0.6 : 1 + rnd() * 0.9;
        sky.append('circle').attr('class', 'star')
          .attr('cx', (rnd() - 0.5) * span * 1.3).attr('cy', (rnd() - 0.5) * span * 1.3).attr('r', r)
          .attr('fill', rnd() < 0.15 ? '#ffe9c4' : rnd() < 0.3 ? '#c9dcff' : '#ffffff')
          .style('animation-delay', (-rnd() * 4).toFixed(2) + 's')
          .style('animation-duration', (2.5 + rnd() * 4).toFixed(2) + 's');
      }
      // A faint band of nebulosity behind the ring.
      const neb = defs.select('#fo-neb').empty() ? defs.append('radialGradient').attr('id', 'fo-neb') : defs.select('#fo-neb');
      neb.selectAll('*').remove();
      neb.append('stop').attr('offset', '0%').attr('stop-color', '#3a4fa0').attr('stop-opacity', 0.22);
      neb.append('stop').attr('offset', '55%').attr('stop-color', '#5b2f86').attr('stop-opacity', 0.08);
      neb.append('stop').attr('offset', '100%').attr('stop-color', '#000').attr('stop-opacity', 0);
      sky.insert('circle', ':first-child').attr('r', M * 1.05).attr('fill', 'url(#fo-neb)');
    }

    function layout() {
      const n = galaxies.length;
      if (!n) return;
      const sector = TAU / n;
      // Past ~90 galaxies the fans get too small to read: draw a progress ring per galaxy
      // instead, and let the galaxy ring use the space the fans would have taken.
      arcMode = n > 90;
      R1 = M * (arcMode ? 0.82 : 0.44);
      RB0 = R1 + M * 0.09;
      nodeSize = Math.max(2.6, Math.min(14, R1 * sector * 0.34));
      galaxies.forEach((g, i) => {
        g.angle = i * sector;
        [g.x, g.y] = polar(g.angle, R1);
        g.sector = sector;
      });
      for (const g of galaxies) layoutDots(g);
    }

    // Place a galaxy's realizations in rows across its sector, row by row outward.
    function layoutDots(g) {
      const st = states.get(g.id);
      const nBoot = st && st.nBoot;
      g.dots = [];
      if (!nBoot || arcMode) return;
      const width = Math.min(g.sector * 0.86, 1.25);
      const depth = RBmax - RB0;
      let sp = Math.min(17, depth / 2.5);
      let rows;
      for (; sp >= 2.2; sp -= 0.3) {
        rows = [];
        let left = nBoot, r = RB0;
        while (left > 0 && r <= RB0 + depth) {
          const perRow = Math.max(1, Math.floor((r * width) / sp));
          rows.push({ r, count: Math.min(perRow, left) });
          left -= perRow;
          r += sp;
        }
        if (left <= 0) break;
      }
      const dotR = Math.max(0.8, Math.min(3.6, sp * 0.26));
      let k = 0;
      for (const row of rows) {
        for (let j = 0; j < row.count; j++, k++) {
          const t = row.count === 1 ? 0.5 : j / (row.count - 1);
          const a = g.angle + (t - 0.5) * width * Math.min(1, (row.count * sp) / (row.r * width));
          const [x, y] = polar(a, row.r);
          g.dots.push({ g, k, a, r: row.r, x, y, dotR });
        }
      }
    }

    // ---------------------------------------------------------------- render
    function spokePath(g) {
      const [cx, cy] = polar(g.angle - 0.42, R1 * 0.55);
      return `M0,0 Q${cx},${cy} ${g.x},${g.y}`;
    }
    function fanPath(d) {
      const [cx, cy] = polar(d.g.angle, (R1 + d.r) / 2);
      return `M${d.g.x},${d.g.y} Q${cx},${cy} ${d.x},${d.y}`;
    }

    function render() {
      emptyEl.hidden = galaxies.length > 0;
      renderHub();
      const status = (g) => (states.get(g.id) || {}).status || 'queued';

      gSpokes.selectAll('path.spoke').data(galaxies, (g) => g.id).join('path')
        .attr('class', (g) => 'spoke' + (status(g) === 'fitting' ? ' flowing' : ''))
        .attr('d', spokePath).attr('fill', 'none')
        .attr('stroke', (g) => statusColor(status(g)))
        .attr('stroke-opacity', (g) => (status(g) === 'queued' ? 0.18 : 0.45))
        .attr('stroke-width', Math.max(0.5, Math.min(1.4, nodeSize / 7)));

      const allDots = arcMode ? [] : galaxies.flatMap((g) => g.dots || []);
      const dotState = (d) => {
        const st = states.get(d.g.id);
        if (!st) return 'pending';
        if (st.failed && st.failed.has(d.k)) return 'failed';
        if (st.done && st.done.has(d.k)) return 'done';
        if (st.status === 'bootstrapping' && d.k < doneCount(st) + (st.inFlight || 0)) return 'running';
        return 'pending';
      };
      gFanLinks.selectAll('path').data(allDots, (d) => d.g.id + ':' + d.k).join('path')
        .attr('d', fanPath).attr('fill', 'none')
        .attr('stroke', (d) => (dotState(d) === 'done' ? cssColor('--s-done') : dotState(d) === 'failed' ? cssColor('--s-failed') : '#9fb4ff'))
        .attr('stroke-opacity', (d) => (dotState(d) === 'pending' ? 0.06 : dotState(d) === 'running' ? 0.2 : 0.16))
        .attr('stroke-width', 0.5);
      gDots.selectAll('circle').data(allDots, (d) => d.g.id + ':' + d.k).join(
        (enter) => enter.append('circle').attr('cx', (d) => d.g.x).attr('cy', (d) => d.g.y).attr('r', 0)
          .call((e) => e.transition().duration(900).delay((d) => d.k * 12).ease(d3.easeBackOut)
            .attr('cx', (d) => d.x).attr('cy', (d) => d.y).attr('r', (d) => d.dotR)),
        (update) => update.attr('cx', (d) => d.x).attr('cy', (d) => d.y).attr('r', (d) => d.dotR),
      )
        .attr('class', (d) => 'dot' + (dotState(d) === 'running' ? ' running' : ''))
        .attr('fill', (d) => ({ done: cssColor('--s-done'), failed: cssColor('--s-failed'), running: cssColor('--s-boot') }[dotState(d)] || '#ffffff'))
        .attr('fill-opacity', (d) => (dotState(d) === 'pending' ? 0.13 : 1))
        .attr('filter', (d) => (dotState(d) === 'done' ? 'url(#fo-glow)' : null))
        .on('mouseenter mousemove', (ev, d) => showDotTip(ev, d, dotState(d)))
        .on('mouseleave', () => tip.style('opacity', 0));

      // Arc mode: a progress ring per galaxy instead of individual dots.
      const arcData = arcMode ? galaxies.filter((g) => (states.get(g.id) || {}).nBoot) : [];
      const arcGen = d3.arc();
      gArcs.selectAll('g.prog').data(arcData, (g) => g.id).join((e) => {
        const grp = e.append('g').attr('class', 'prog');
        grp.append('path').attr('class', 'track');
        grp.append('path').attr('class', 'bar');
        return grp;
      })
        .attr('transform', (g) => `translate(${g.x},${g.y})`)
        .each(function (g) {
          const st = states.get(g.id);
          const frac = st.nBoot ? doneCount(st) / st.nBoot : 0;
          const inner = nodeSize * 1.5, outer = nodeSize * 1.5 + Math.max(1.2, nodeSize * 0.45);
          d3.select(this).select('.track').attr('d', arcGen({ innerRadius: inner, outerRadius: outer, startAngle: 0, endAngle: TAU }))
            .attr('fill', '#ffffff').attr('fill-opacity', 0.08);
          d3.select(this).select('.bar').attr('d', arcGen({ innerRadius: inner, outerRadius: outer, startAngle: 0, endAngle: TAU * frac }))
            .attr('fill', frac >= 1 ? cssColor('--s-done') : cssColor('--s-boot'));
        });

      const nodes = gGalaxies.selectAll('g.gal').data(galaxies, (g) => g.id).join((enter) => {
        const grp = enter.append('g').attr('class', 'gal').style('cursor', 'pointer');
        grp.append('circle').attr('class', 'halo');
        grp.append('circle').attr('class', 'sel-ring').attr('fill', 'none');
        grp.append('ellipse').attr('class', 'disk');
        grp.append('ellipse').attr('class', 'core');
        grp.append('text').attr('class', 'lbl');
        grp.on('mouseenter', (ev, g) => showTip(ev, g)).on('mousemove', (ev, g) => showTip(ev, g))
          .on('mouseleave', () => tip.style('opacity', 0))
          .on('click', (ev, g) => { select(g.id); if (onSelect) onSelect(g.id); });
        return grp;
      });
      nodes.attr('transform', (g) => `translate(${g.x},${g.y})`);
      nodes.select('.halo')
        .attr('class', (g) => 'halo' + (['fitting', 'bootstrapping'].includes(status(g)) ? ' active' : ''))
        .attr('r', nodeSize * 1.9).attr('fill', (g) => statusColor(status(g))).attr('fill-opacity', 0.16)
        .attr('filter', 'url(#fo-glow-strong)');
      nodes.select('.sel-ring').attr('r', nodeSize * 2.4)
        .attr('stroke', '#ffffff').attr('stroke-width', 1.2).attr('stroke-dasharray', '2 3')
        .attr('opacity', (g) => (g.id === selected ? 0.9 : 0));
      // Sky PA is east of north; SVG rotates clockwise from +x, so the major axis at PA
      // points up when PA = 0.
      nodes.select('.disk')
        .attr('rx', (g) => nodeSize * Math.max(0.28, Math.cos((g.inc || 0) * Math.PI / 180)))
        .attr('ry', nodeSize)
        .attr('transform', (g) => `rotate(${g.pa || 0})`)
        .attr('fill', (g) => `url(#fo-gal-${status(g)})`);
      nodes.select('.core')
        .attr('rx', (g) => nodeSize * 0.22 * Math.max(0.4, Math.cos((g.inc || 0) * Math.PI / 180)))
        .attr('ry', nodeSize * 0.22)
        .attr('transform', (g) => `rotate(${g.pa || 0})`)
        .attr('fill', '#ffffff').attr('filter', 'url(#fo-glow)');
      const showLabels = galaxies.length <= 40;
      nodes.select('.lbl')
        .text((g) => (showLabels ? g.short : ''))
        .attr('fill', cssColor('--space-muted')).attr('font-size', 9.5)
        .attr('font-family', 'ui-monospace, Menlo, monospace')
        .attr('text-anchor', (g) => (g.angle > Math.PI ? 'start' : 'end'))
        .attr('dominant-baseline', 'middle')
        .attr('transform', (g) => {
          const [lx, ly] = polar(g.angle, -nodeSize * 2.6);
          const deg = (g.angle * 180) / Math.PI;
          const rot = g.angle > Math.PI ? deg + 90 : deg - 90;
          return `translate(${lx},${ly}) rotate(${rot})`;
        });
      renderStats();
    }

    function doneCount(st) { return (st.done ? st.done.size : 0) + (st.failed ? st.failed.size : 0); }

    function renderHub() {
      gHub.selectAll('*').remove();
      if (!galaxies.length) return;
      const r = Math.max(14, M * 0.075);
      gHub.append('circle').attr('r', r * 2.6).attr('fill', 'url(#fo-hub)').attr('opacity', running ? 0.55 : 0.25);
      gHub.append('circle').attr('class', 'hub-ring' + (running ? ' running' : '')).attr('r', r * 1.45)
        .attr('fill', 'none').attr('stroke', cssColor('--s-fitting')).attr('stroke-opacity', 0.55)
        .attr('stroke-dasharray', '2 5').attr('stroke-width', 1.2);
      gHub.append('circle').attr('r', r).attr('fill', '#0d1322').attr('stroke', cssColor('--s-fitting'))
        .attr('stroke-width', 1.5).attr('filter', 'url(#fo-glow)');
      const counts = tally();
      gHub.append('text').attr('text-anchor', 'middle').attr('dy', '-0.15em')
        .attr('fill', cssColor('--space-ink')).attr('font-size', Math.max(11, r * 0.55)).attr('font-weight', 700)
        .attr('font-family', 'ui-monospace, Menlo, monospace')
        .text(`${counts.fitDone}/${galaxies.length}`);
      gHub.append('text').attr('text-anchor', 'middle').attr('dy', '1.25em')
        .attr('fill', cssColor('--space-muted')).attr('font-size', Math.max(8, r * 0.3))
        .text('initial fits');
    }

    // ---------------------------------------------------------------- stats
    function tally() {
      const c = { total: galaxies.length, queued: 0, fitting: 0, fit: 0, bootstrapping: 0, done: 0, failed: 0, fitDone: 0, real: 0, realTotal: 0, realFailed: 0 };
      for (const g of galaxies) {
        const st = states.get(g.id) || { status: 'queued' };
        if (st.status === 'fit-failed' || st.status === 'bs-failed') c.failed++;
        else if (c[st.status] !== undefined) c[st.status]++;
        if (!['queued', 'fitting', 'stopped'].includes(st.status)) c.fitDone++;
        if (st.nBoot) { c.realTotal += st.nBoot; c.real += doneCount(st); c.realFailed += st.failed ? st.failed.size : 0; }
      }
      return c;
    }

    function renderStats() {
      const c = tally();
      const elapsed = startedAt ? Math.round((Date.now() - startedAt) / 1000) : 0;
      const items = [
        ['Galaxies', c.total, null],
        ['Fitting', c.fitting, '--s-fitting'],
        ['Fitted', c.fit, '--s-fit'],
        ['Bootstrapping', c.bootstrapping, '--s-boot'],
        ['Done', c.done, '--s-done'],
        ['Galaxies failed', c.failed, '--s-failed'],
        ['Realizations', c.realTotal ? `${c.real}/${c.realTotal}` : '0', null],
        ['Realizations failed', c.realFailed, '--s-failed'],
        ['Elapsed', startedAt ? `${Math.floor(elapsed / 60)}:${String(elapsed % 60).padStart(2, '0')}` : '–', null],
      ];
      d3.select(statsEl).selectAll('div.stat').data(items).join((e) => {
        const s = e.append('div').attr('class', 'stat');
        s.append('div').attr('class', 'v');
        const k = s.append('div').attr('class', 'k');
        k.append('i');
        k.append('span');
        return s;
      }).each(function ([k, v, color]) {
        const s = d3.select(this);
        s.select('.v').text(v);
        s.select('.k span').text(k);
        s.select('.k i').style('display', color ? null : 'none').style('background', color ? cssColor(color) : null);
      });
    }
    setInterval(() => { if (running) renderStats(); }, 1000);

    function renderLegend() {
      const items = [['queued', 'queued'], ['fitting', 'fitting'], ['fit', 'fitted'], ['bootstrapping', 'bootstrapping'], ['done', 'done'], ['fit-failed', 'failed']];
      const sel = d3.select(legendEl).selectAll('span').data(items).join((e) => {
        const s = e.append('span');
        s.append('i');
        s.append('b').style('font-weight', 400);
        return s;
      });
      sel.select('i').style('background', ([s]) => statusColor(s));
      sel.select('b').text(([, l]) => l);
    }
    renderLegend();

    // ---------------------------------------------------------------- particles
    function particleAlong(pathD, { color, r = 2.2, duration = 900, reverse = true }) {
      if (galaxies.length > BIG_RUN) return;
      const path = gParticles.append('path').attr('d', pathD).attr('fill', 'none').attr('stroke', 'none');
      const node = path.node();
      const len = node.getTotalLength();
      const dot = gParticles.append('circle').attr('r', r / Math.sqrt(zoomK)).attr('fill', color).attr('filter', 'url(#fo-glow-strong)');
      const tail = gParticles.append('circle').attr('r', (r * 0.6) / Math.sqrt(zoomK)).attr('fill', color).attr('opacity', 0.5);
      dot.transition().duration(duration).ease(d3.easeCubicIn)
        .attrTween('transform', () => (t) => {
          const p = node.getPointAtLength(len * (reverse ? 1 - t : t));
          return `translate(${p.x},${p.y})`;
        })
        .on('end', () => { dot.remove(); path.remove(); });
      tail.transition().duration(duration).delay(70).ease(d3.easeCubicIn)
        .attrTween('transform', () => (t) => {
          const p = node.getPointAtLength(len * (reverse ? 1 - t : t));
          return `translate(${p.x},${p.y})`;
        })
        .on('end', () => tail.remove());
    }
    function burst(x, y, color) {
      if (galaxies.length > BIG_RUN) return;
      const ring = gParticles.append('circle').attr('cx', x).attr('cy', y).attr('r', nodeSize)
        .attr('fill', 'none').attr('stroke', color).attr('stroke-width', 1.5).attr('opacity', 0.9);
      ring.transition().duration(900).ease(d3.easeCubicOut).attr('r', nodeSize * 4).attr('opacity', 0).remove();
    }

    // ---------------------------------------------------------------- tooltip
    function placeTip(ev) {
      const [mx, my] = d3.pointer(ev, stageEl);
      const left = Math.min(mx + 14, (stageEl.clientWidth || W) - 270);
      tip.style('left', left + 'px').style('top', (my + 14) + 'px').style('opacity', 1);
    }
    function showDotTip(ev, d, state) {
      const st = states.get(d.g.id) || {};
      const note = st.notes && st.notes.get(d.k);
      tip.selectAll('*').remove();
      tip.append('div').attr('class', 'tname').text(`${d.g.name} · realization ${d.k + 1}`);
      tip.append('div').attr('class', 'trow').text({ done: 'done', failed: 'failed', running: 'running' }[state] || 'waiting');
      if (note) tip.append('div').attr('class', 'trow').text(note);
      placeTip(ev);
    }
    function showTip(ev, g) {
      const st = states.get(g.id) || { status: 'queued' };
      const [mx, my] = d3.pointer(ev, stageEl);
      tip.selectAll('*').remove();
      tip.append('div').attr('class', 'tname').text(g.name);
      const rows = [
        `${g.field} · Vsys ${g.vsys} km/s · ${g.sizeBeams} beams`,
        `PA ${g.pa.toFixed(1)}° · inc ${g.inc.toFixed(1)}° (start)`,
        `status: ${STATUS_LABEL[st.status] || st.status}`,
      ];
      if (st.fitSummary) rows.push(`fit: PA ${st.fitSummary.pa}° · inc ${st.fitSummary.inc}° · Vsys ${st.fitSummary.vsys}`);
      if (st.nBoot) rows.push(`realizations ${doneCount(st)}/${st.nBoot}` + (st.failed && st.failed.size ? `, ${st.failed.size} failed` : ''));
      if (st.error) rows.push(`error: ${String(st.error).slice(0, 120)}`);
      for (const r of rows) tip.append('div').attr('class', 'trow').text(r);
      const left = Math.min(mx + 14, (stageEl.clientWidth || W) - 270);
      tip.style('left', left + 'px').style('top', (my + 14) + 'px').style('opacity', 1);
    }

    // ---------------------------------------------------------------- API
    let renderQueued = false;
    function scheduleRender() {
      if (renderQueued) return;
      renderQueued = true;
      // Thousands of galaxies make a full render expensive, and results can arrive hundreds
      // a second, so big runs redraw at most once a second.
      if (galaxies.length > BIG_RUN) setTimeout(() => requestAnimationFrame(() => { renderQueued = false; render(); }), 1000);
      else requestAnimationFrame(() => { renderQueued = false; render(); });
    }

    function select(id) {
      selected = id;
      gGalaxies.selectAll('g.gal').select('.sel-ring').attr('opacity', (g) => (g.id === selected ? 0.9 : 0));
    }

    window.addEventListener('resize', () => { clearTimeout(resize.t); resize.t = setTimeout(resize, 150); });
    document.addEventListener('fullscreenchange', () => { svg.call(zoom.transform, d3.zoomIdentity); resize(); });

    return {
      setGalaxies(list, stateMap) {
        galaxies = list.map((g) => Object.assign({}, g));
        byId = new Map(galaxies.map((g) => [g.id, g]));
        states = stateMap;
        svg.call(zoom.transform, d3.zoomIdentity);
        resize();
      },
      // A galaxy's realization count became known (its bootstrap job deployed).
      relayout(id) {
        const g = byId.get(id);
        if (g) layoutDots(g);
        scheduleRender();
      },
      update: scheduleRender,
      setRunning(on) {
        running = on;
        if (on && !startedAt) startedAt = Date.now();
        if (on) scheduleRender();
        else render();  // the final state, without the big-run delay
      },
      resetClock() { startedAt = null; },
      dispatched() {
        for (const g of galaxies) particleAlong(spokePath(g), { color: cssColor('--s-fitting'), reverse: false, duration: 1100, r: 1.6 });
      },
      fitLanded(id, ok) {
        const g = byId.get(id);
        if (!g) return;
        const color = ok ? cssColor('--s-fit') : cssColor('--s-failed');
        particleAlong(spokePath(g), { color, duration: 1000, r: 2.6 });
        burst(g.x, g.y, color);
        scheduleRender();
      },
      realizationLanded(id, k, ok) {
        const g = byId.get(id);
        if (!g) return;
        const d = (g.dots || [])[k];
        const color = ok ? cssColor('--s-done') : cssColor('--s-failed');
        if (d) particleAlong(fanPath(d), { color, duration: 650, r: 1.6 });
        scheduleRender();
      },
      galaxyDone(id) {
        const g = byId.get(id);
        if (g) burst(g.x, g.y, cssColor('--s-done'));
        scheduleRender();
      },
      select,
    };
  }

  // ---------------------------------------------------------------- detail plots
  function plotFrame(svgEl, h) {
    const w = svgEl.clientWidth || 420;
    const svg = d3.select(svgEl).attr('viewBox', `0 0 ${w} ${h}`);
    svg.selectAll('*').remove();
    return { svg, w, h, m: { l: 40, r: 10, t: 10, b: 30 } };
  }
  function axes(f, x, y, xLabel, yLabel) {
    const { svg, w, h, m } = f;
    const style = (g) => {
      g.select('.domain').attr('stroke', '#ffffff').attr('stroke-opacity', 0.15);
      g.selectAll('line').attr('stroke', '#ffffff').attr('stroke-opacity', 0.07);
      g.selectAll('text').attr('fill', cssColor('--space-muted')).attr('font-size', 9.5);
    };
    svg.append('g').attr('transform', `translate(0,${h - m.b})`).call(d3.axisBottom(x).ticks(5).tickSize(-(h - m.t - m.b))).call(style);
    svg.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickSize(-(w - m.l - m.r))).call(style);
    svg.append('text').attr('x', (m.l + w - m.r) / 2).attr('y', h - 3).attr('text-anchor', 'middle')
      .attr('fill', cssColor('--space-muted')).attr('font-size', 10).text(xLabel);
    svg.append('text').attr('transform', `translate(11,${(m.t + h - m.b) / 2}) rotate(-90)`).attr('text-anchor', 'middle')
      .attr('fill', cssColor('--space-muted')).attr('font-size', 10).text(yLabel);
  }
  const finiteArr = (a) => Array.isArray(a) ? a.map(Number).filter(Number.isFinite) : [];

  // Best fit with per-ring error bars over the realizations, like BootstrapBoxPlot.py.
  function drawRotationCurves(svgEl, fit, reports, vrotErr = []) {
    const f = plotFrame(svgEl, 220);
    const curves = reports.filter((r) => r && !r.sofiaFailed && Array.isArray(r.R) && Array.isArray(r.VROT))
      .map((r) => r.R.map((R, i) => [Number(R), Number(r.VROT[i])]).filter(([a, b]) => Number.isFinite(a) && Number.isFinite(b)));
    const best = fit && Array.isArray(fit.R) && Array.isArray(fit.VROT)
      ? fit.R.map((R, i) => [Number(R), Number(fit.VROT[i]), Number(vrotErr[i])]).filter(([a, b]) => Number.isFinite(a) && Number.isFinite(b)) : [];
    const pts = curves.flat().concat(best.map(([r, v, e]) => [r, Number.isFinite(e) ? v + e : v]));
    if (!pts.length) {
      f.svg.append('text').attr('x', f.w / 2).attr('y', f.h / 2).attr('text-anchor', 'middle').attr('fill', cssColor('--space-muted')).attr('font-size', 11).text('No rotation curve yet');
      return;
    }
    const x = d3.scaleLinear().domain([0, d3.max(pts, (p) => p[0]) * 1.05]).nice().range([f.m.l, f.w - f.m.r]);
    const y = d3.scaleLinear().domain([0, d3.max(pts, (p) => p[1]) * 1.12]).nice().range([f.h - f.m.b, f.m.t]);
    axes(f, x, y, 'radius (arcsec)', 'Vrot (km/s)');
    const line = d3.line().x((p) => x(p[0])).y((p) => y(p[1])).curve(d3.curveMonotoneX);
    f.svg.append('g').selectAll('path').data(curves).join('path').attr('d', line).attr('fill', 'none')
      .attr('stroke', cssColor('--s-done')).attr('stroke-opacity', Math.max(0.08, Math.min(0.5, 6 / Math.max(1, curves.length)))).attr('stroke-width', 1);
    if (best.length) {
      const bars = f.svg.append('g').selectAll('g').data(best.filter((p) => Number.isFinite(p[2]))).join('g')
        .attr('stroke', '#ffffff').attr('stroke-width', 1.4);
      bars.append('line').attr('x1', (p) => x(p[0])).attr('x2', (p) => x(p[0]))
        .attr('y1', (p) => y(Math.max(0, p[1] - p[2]))).attr('y2', (p) => y(p[1] + p[2]));
      for (const sign of [-1, 1]) {
        bars.append('line').attr('x1', (p) => x(p[0]) - 4).attr('x2', (p) => x(p[0]) + 4)
          .attr('y1', (p) => y(Math.max(0, p[1] + sign * p[2]))).attr('y2', (p) => y(Math.max(0, p[1] + sign * p[2])));
      }
      f.svg.append('path').datum(best).attr('d', line).attr('fill', 'none').attr('stroke', '#ffffff').attr('stroke-width', 2.2);
      f.svg.append('g').selectAll('circle').data(best).join('circle').attr('cx', (p) => x(p[0])).attr('cy', (p) => y(p[1])).attr('r', 3)
        .attr('fill', '#ffffff');
    }
  }

  function drawGeoCloud(svgEl, fit, reports) {
    const f = plotFrame(svgEl, 220);
    const pts = reports.filter((r) => r && !r.sofiaFailed)
      .map((r) => [Number(finiteArr(r.POSITIONANGLE)[0]), Number(finiteArr(r.INCLINATION)[0])])
      .filter(([a, b]) => Number.isFinite(a) && Number.isFinite(b));
    const best = fit ? [finiteArr(fit.POSITIONANGLE)[0], finiteArr(fit.INCLINATION)[0]] : null;
    const all = pts.concat(best && best.every(Number.isFinite) ? [best] : []);
    if (!all.length) {
      f.svg.append('text').attr('x', f.w / 2).attr('y', f.h / 2).attr('text-anchor', 'middle').attr('fill', cssColor('--space-muted')).attr('font-size', 11).text('No realizations yet');
      return;
    }
    const pad = (ext, min) => { const s = Math.max(min, (ext[1] - ext[0]) * 0.15); return [ext[0] - s, ext[1] + s]; };
    const x = d3.scaleLinear().domain(pad(d3.extent(all, (p) => p[0]), 2)).nice().range([f.m.l, f.w - f.m.r]);
    const y = d3.scaleLinear().domain(pad(d3.extent(all, (p) => p[1]), 2)).nice().range([f.h - f.m.b, f.m.t]);
    axes(f, x, y, 'position angle (deg)', 'inclination (deg)');
    if (pts.length >= 6) {
      const density = d3.contourDensity().x((p) => x(p[0])).y((p) => y(p[1])).size([f.w, f.h]).bandwidth(14).thresholds(6)(pts);
      const color = d3.scaleSequential(d3.interpolateRgb('rgba(67,214,196,0.04)', 'rgba(255,209,102,0.4)')).domain([0, d3.max(density, (d) => d.value)]);
      f.svg.append('g').selectAll('path').data(density).join('path').attr('d', d3.geoPath()).attr('fill', (d) => color(d.value))
        .attr('stroke', cssColor('--s-done')).attr('stroke-opacity', 0.18).attr('stroke-width', 0.6);
    }
    f.svg.append('g').selectAll('circle').data(pts).join('circle').attr('cx', (p) => x(p[0])).attr('cy', (p) => y(p[1])).attr('r', 2.4)
      .attr('fill', cssColor('--s-done')).attr('fill-opacity', 0.85);
    if (best && best.every(Number.isFinite)) {
      const [bx, by] = [x(best[0]), y(best[1])];
      f.svg.append('line').attr('x1', bx - 9).attr('x2', bx + 9).attr('y1', by).attr('y2', by).attr('stroke', '#fff').attr('stroke-width', 1.6);
      f.svg.append('line').attr('x1', bx).attr('x2', bx).attr('y1', by - 9).attr('y2', by + 9).attr('stroke', '#fff').attr('stroke-width', 1.6);
      f.svg.append('circle').attr('cx', bx).attr('cy', by).attr('r', 4).attr('fill', 'none').attr('stroke', '#fff').attr('stroke-width', 1.2);
    }
  }

  return { create, drawRotationCurves, drawGeoCloud, statusColor, STATUS_LABEL };
})();
