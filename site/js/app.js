/* HOHC — static web app.
 * Filtering logic is a line-by-line port of sankey.R from the original
 * Shiny app, operating on pre-aggregated mean TPM per (gene, organ, sex).
 */

const MALE_ONLY = ['Testis', 'Prostate'];
const FEMALE_ONLY = ['Fallopian Tube', 'Vagina', 'Ovary', 'Uterus', 'Cervix Uteri'];

const COMMON_ORGANS = [
  'Thyroid', 'Stomach', 'Spleen', 'Small Intestine', 'Skin', 'Salivary Gland',
  'Pituitary', 'Pancreas', 'Nerve', 'Muscle', 'Lung', 'Liver', 'Kidney', 'Heart',
  'Esophagus', 'Colon', 'Breast', 'Brain', 'Blood Vessel', 'Blood', 'Bladder',
  'Adrenal Gland', 'Adipose Tissue'
];

const STAGES = ['Hormone class', 'Hormone organ', 'HR organ', 'HR type'];

const CLASS_COLORS = {
  'Peptide': '#8b5cf6',
  'Lipid-derived': '#f59e0b',
  'Amino acid-derived': '#10b981',
  'Unknown': '#94a3b8'
};
const TYPE_COLORS = { 'Membrane HR': '#e11d48', 'Nuclear HR': '#2563eb' };

// Consistent qualitative palette for the 30 organs
const ORGAN_PALETTE = [
  '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
  '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78',
  '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d',
  '#9edae5', '#393b79', '#637939', '#8c6d31', '#843c39', '#7b4173',
  '#5254a3', '#bd9e39', '#ad494a', '#a55194', '#6b6ecf', '#b5cf6b'
];

let DATA = null;          // loaded JSON
let exprByGene = null;    // sex -> gene -> [{organ, tpm}]
let organColor = {};
let chart = null;

const state = {
  sex: 'Male',
  hormoneClass: 'All',
  hormoneNames: [],  // empty = show the full atlas by default
  hormoneTissues: [],
  receptorTissues: [],
  receptorClass: 'All',
  receptorNames: [],
  topN: 3,
  tpmLo: 0.1,
  tpmHi: 6.35
};
const DEFAULTS = JSON.parse(JSON.stringify(state));

// ---------- data ----------

async function loadData() {
  const resp = await fetch('data/hohc_data.json');
  DATA = await resp.json();

  exprByGene = { Male: new Map(), Female: new Map() };
  for (const e of DATA.expression) {
    const m = exprByGene[e.sex];
    if (!m) continue;
    if (!m.has(e.gene)) m.set(e.gene, []);
    m.get(e.gene).push({ organ: e.organ, tpm: e.tpm });
  }

  const organs = [...new Set([...COMMON_ORGANS, ...MALE_ONLY, ...FEMALE_ONLY, ...DATA.tissues])];
  organs.forEach((o, i) => { organColor[o] = ORGAN_PALETTE[i % ORGAN_PALETTE.length]; });
}

// Port of the core of sankey.R: rows after join + filters.
// Each row: {hclass, hormone, horgan, gene, rtype, rorgan, tpm}
function computeRows(s) {
  const oppose = s.sex === 'Male' ? FEMALE_ONLY : MALE_ONLY;
  const expr = exprByGene[s.sex];

  // Top-N / TPM-range passing organs per gene (thresholds computed on the
  // gene's full organ vector, as in the original group_by + filter).
  const passing = new Map();
  for (const [gene, arr] of expr) {
    const sorted = [...arr].sort((a, b) => b.tpm - a.tpm);
    if (s.topN > sorted.length) { passing.set(gene, []); continue; }
    const thr = sorted[s.topN - 1].tpm;
    passing.set(gene, arr.filter(e => {
      const l = Math.log1p(e.tpm);
      return l > s.tpmLo && l < s.tpmHi && e.tpm >= thr;
    }));
  }

  const rows = [];
  for (const r of DATA.receptors) {
    if (oppose.includes(r.hormone_organ)) continue;
    const hits = passing.get(r.gene);
    if (hits === undefined) {
      // gene absent from GTEx: left_join keeps the row with NA organ
      rows.push({ hclass: r.hormone_class, hormone: r.hormone, horgan: r.hormone_organ,
                  gene: r.gene, rtype: r.receptor_class, rorgan: null, tpm: null });
    } else {
      for (const h of hits) {
        rows.push({ hclass: r.hormone_class, hormone: r.hormone, horgan: r.hormone_organ,
                    gene: r.gene, rtype: r.receptor_class, rorgan: h.organ, tpm: h.tpm });
      }
    }
  }

  return rows.filter(row =>
    (s.hormoneClass === 'All' || row.hclass === s.hormoneClass) &&
    (s.hormoneNames.length === 0 || s.hormoneNames.includes(row.hormone)) &&
    (s.hormoneTissues.length === 0 || s.hormoneTissues.includes(row.horgan)) &&
    (s.receptorClass === 'All' || row.rtype === s.receptorClass) &&
    (s.receptorNames.length === 0 || s.receptorNames.includes(row.gene)) &&
    (s.receptorTissues.length === 0 || s.receptorTissues.includes(row.rorgan))
  );
}

// ---------- sankey ----------

function nodeColor(stage, name) {
  if (stage === 0) return CLASS_COLORS[name] || '#94a3b8';
  if (stage === 3) return TYPE_COLORS[name] || '#94a3b8';
  return organColor[name] || '#94a3b8';
}

function isLight() {
  return document.documentElement.dataset.theme === 'light';
}

function render() {
  const rows = computeRows(state);
  const msg = document.getElementById('message');

  const values = r => [r.hclass, r.horgan, r.rorgan, r.rtype];

  // links between adjacent stages (count rows where both ends are present)
  const linkCount = new Map();
  const nodeSet = [new Set(), new Set(), new Set(), new Set()];
  let chains = 0;
  for (const r of rows) {
    const v = values(r);
    if (v[0] && v[1] && v[2] && v[3]) chains++;
    for (let i = 0; i < 4; i++) if (v[i]) nodeSet[i].add(v[i]);
    for (let i = 0; i < 3; i++) {
      if (v[i] && v[i + 1]) {
        const key = i + '|' + v[i] + '|' + v[i + 1];
        linkCount.set(key, (linkCount.get(key) || 0) + 1);
      }
    }
  }

  if (chains === 0) {
    msg.hidden = false;
    msg.textContent = 'No results found. You may want to reset the filters or try different options.';
  } else {
    msg.hidden = true;
  }

  const nodes = [];
  for (let i = 0; i < 4; i++) {
    for (const name of nodeSet[i]) {
      nodes.push({
        name: i + ':' + name,
        depth: i,
        itemStyle: { color: nodeColor(i, name) },
        label: { formatter: name }
      });
    }
  }
  const links = [...linkCount.entries()].map(([key, n]) => {
    const [i, a, b] = key.split('|');
    return { source: i + ':' + a, target: (Number(i) + 1) + ':' + b, value: n };
  });

  chart.setOption({
    tooltip: {
      trigger: 'item',
      backgroundColor: isLight() ? '#ffffff' : '#131c31',
      borderColor: isLight() ? 'rgba(15,23,42,0.15)' : 'rgba(255,255,255,0.15)',
      textStyle: { color: isLight() ? '#16213a' : '#e8edf7', fontFamily: 'Inter' },
      formatter: p => {
        if (p.dataType === 'edge') {
          const a = p.data.source.slice(2), b = p.data.target.slice(2);
          return `<b>${a}</b> → <b>${b}</b><br>${p.data.value} connection${p.data.value > 1 ? 's' : ''}`;
        }
        return `<b>${p.name.slice(2)}</b>`;
      }
    },
    series: [{
      type: 'sankey',
      data: nodes,
      links: links,
      left: 10, right: chart.getWidth() < 520 ? 72 : 110, top: 14, bottom: 14,
      nodeWidth: chart.getWidth() < 520 ? 8 : 10,
      nodeGap: links.length > 60 ? (chart.getWidth() < 520 ? 5 : 8) : 18,
      nodeAlign: 'justify',
      draggable: true,
      emphasis: { focus: 'adjacency' },
      blur: { lineStyle: { opacity: 0.06 }, itemStyle: { opacity: 0.2 } },
      itemStyle: { borderRadius: 3 },
      lineStyle: { color: 'gradient', opacity: 0.38, curveness: 0.55 },
      label: {
        fontSize: chart.getWidth() < 520 ? 10 : 12.5,
        fontFamily: 'Inter', color: isLight() ? '#16213a' : '#e8edf7',
        textBorderColor: 'transparent'
      },
      layoutIterations: links.length > 150 ? 24 : 48
    }]
  }, true);

  // header stats
  const genes = new Set(rows.filter(r => r.rorgan).map(r => r.gene));
  const hormones = new Set(rows.filter(r => r.rorgan).map(r => r.hormone));
  document.getElementById('headerStats').innerHTML =
    `<span class="stat-chip"><b>${hormones.size}</b> hormones</span>` +
    `<span class="stat-chip"><b>${genes.size}</b> receptors</span>` +
    `<span class="stat-chip"><b>${links.length}</b> links</span>` +
    `<span class="stat-chip">${state.sex}</span>`;
}

// ---------- download ----------

function downloadTSV() {
  const rows = computeRows(state).filter(r => r.gene);
  const header = ['Hormone class', 'Hormone name', 'Hormone organ',
                  'HR gene name', 'HR type', 'HR organ', 'Mean TPM'];
  const lines = [header.join('\t')];
  for (const r of rows) {
    lines.push([r.hclass, r.hormone, r.horgan, r.gene, r.rtype,
                r.rorgan ?? 'NA', r.tpm ?? 'NA'].join('\t'));
  }
  const blob = new Blob([lines.join('\n')], { type: 'text/tab-separated-values' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `HOHC_data-${state.sex}.tsv`;
  a.click();
  URL.revokeObjectURL(a.href);
}

// ---------- UI: multi-select chips ----------

function makeMulti(elId, optionsFn, stateKey, onChange) {
  const root = document.getElementById(elId);
  const input = document.createElement('input');
  input.placeholder = 'Type to search…';
  const list = document.createElement('div');
  list.className = 'multi-list';
  list.hidden = true;

  function redraw() {
    root.querySelectorAll('.chip').forEach(c => c.remove());
    for (const v of state[stateKey]) {
      const chip = document.createElement('span');
      chip.className = 'chip';
      chip.textContent = v;
      const x = document.createElement('button');
      x.textContent = '×';
      x.onclick = e => {
        e.stopPropagation();
        state[stateKey] = state[stateKey].filter(s => s !== v);
        redraw(); onChange();
      };
      chip.appendChild(x);
      root.insertBefore(chip, input);
    }
  }

  function showList() {
    const q = input.value.toLowerCase();
    const opts = optionsFn().filter(o =>
      !state[stateKey].includes(o) && o.toLowerCase().includes(q));
    list.innerHTML = '';
    for (const o of opts.slice(0, 200)) {
      const d = document.createElement('div');
      d.textContent = o;
      d.onmousedown = e => {
        e.preventDefault();
        state[stateKey] = [...state[stateKey], o];
        input.value = '';
        redraw(); showList(); onChange();
      };
      list.appendChild(d);
    }
    list.hidden = opts.length === 0;
  }

  input.addEventListener('focus', showList);
  input.addEventListener('input', showList);
  input.addEventListener('blur', () => setTimeout(() => { list.hidden = true; }, 150));
  root.addEventListener('click', () => input.focus());
  root.appendChild(input);
  root.appendChild(list);
  root._redraw = redraw;
  redraw();
}

// ---------- wiring ----------

function sexTissues() {
  const extra = state.sex === 'Male' ? MALE_ONLY : FEMALE_ONLY;
  return [...COMMON_ORGANS, ...extra].sort();
}

// Coalesce rapid input events (slider drags, fast typing) into one render
// per animation frame to keep the UI responsive.
let renderQueued = false;
function scheduleRender() {
  if (renderQueued) return;
  renderQueued = true;
  requestAnimationFrame(() => {
    renderQueued = false;
    render();
  });
}

// ---------- zoom / pan ----------

let zoom = 1;
const ZOOM_MIN = 0.5, ZOOM_MAX = 4;

function applyZoom(z, cx, cy) {
  const vp = document.getElementById('chartViewport');
  const el = document.getElementById('chart');
  const prev = zoom;
  zoom = Math.min(ZOOM_MAX, Math.max(ZOOM_MIN, z));
  // explicit pixel sizes: percentage heights are unreliable inside flex.
  // The unzoomed base is captured when zooming starts, otherwise reading
  // vp.clientWidth while already zoomed would compound the scale.
  if (prev === 1 || !applyZoom.base) {
    applyZoom.base = { w: vp.clientWidth, h: vp.clientHeight };
  }
  const baseW = applyZoom.base.w, baseH = applyZoom.base.h;
  cx = cx ?? vp.clientWidth / 2;
  cy = cy ?? vp.clientHeight / 2;
  // keep the focal point (cx, cy in viewport coords) stable while resizing
  const fx = (vp.scrollLeft + cx) / prev;
  const fy = (vp.scrollTop + cy) / prev;
  el.style.width = Math.round(baseW * zoom) + 'px';
  el.style.height = Math.round(baseH * zoom) + 'px';
  chart.resize();
  vp.scrollLeft = fx * zoom - cx;
  vp.scrollTop = fy * zoom - cy;
  vp.classList.toggle('pannable', zoom > 1);
  document.getElementById('zoomLevel').textContent = Math.round(zoom * 100) + '%';
}

function initZoom() {
  const vp = document.getElementById('chartViewport');
  document.getElementById('zoomIn').onclick = () => applyZoom(zoom * 1.25);
  document.getElementById('zoomOut').onclick = () => applyZoom(zoom / 1.25);
  document.getElementById('zoomReset').onclick = () => applyZoom(1);

  // wheel / trackpad pinch zooms toward the cursor
  vp.addEventListener('wheel', e => {
    e.preventDefault();
    const rect = vp.getBoundingClientRect();
    const factor = e.ctrlKey || e.metaKey ? 1.18 : 1.12; // pinch feels stronger
    applyZoom(zoom * (e.deltaY < 0 ? factor : 1 / factor),
              e.clientX - rect.left, e.clientY - rect.top);
  }, { passive: false });

  // --- touch: one finger pans when zoomed in, two fingers pinch-zoom ---
  let touchState = null;
  const startPan = t => ({ x: t.clientX, y: t.clientY, sl: vp.scrollLeft, st: vp.scrollTop });
  vp.addEventListener('touchstart', e => {
    if (e.touches.length === 2) {
      // claim the gesture before Safari starts a page pinch
      e.preventDefault();
      const [a, b] = e.touches;
      const rect = vp.getBoundingClientRect();
      touchState = {
        pinch: Math.hypot(a.clientX - b.clientX, a.clientY - b.clientY),
        zoom0: zoom,
        cx: (a.clientX + b.clientX) / 2 - rect.left,
        cy: (a.clientY + b.clientY) / 2 - rect.top
      };
    } else if (e.touches.length === 1 && zoom > 1) {
      touchState = startPan(e.touches[0]);
    }
  }, { passive: false });
  vp.addEventListener('touchmove', e => {
    if (!touchState) return;
    if (touchState.pinch && e.touches.length >= 2) {
      e.preventDefault();
      const [a, b] = e.touches;
      const d = Math.hypot(a.clientX - b.clientX, a.clientY - b.clientY);
      applyZoom(touchState.zoom0 * (d / touchState.pinch), touchState.cx, touchState.cy);
    } else if (!touchState.pinch && e.touches.length === 1) {
      e.preventDefault();
      vp.scrollLeft = touchState.sl - (e.touches[0].clientX - touchState.x);
      vp.scrollTop = touchState.st - (e.touches[0].clientY - touchState.y);
    }
  }, { passive: false });
  vp.addEventListener('touchend', e => {
    // pinch released into a single remaining finger -> continue as a pan
    touchState = (e.touches.length === 1 && zoom > 1) ? startPan(e.touches[0]) : null;
  }, { passive: true });
  // block Safari's own gesture (page) zoom over the chart
  vp.addEventListener('gesturestart', e => e.preventDefault());

  // --- mouse: drag empty space to pan; dragging a node is left to ECharts ---
  let drag = null;
  let onNode = false;
  chart.on('mousedown', p => { if (p.dataType === 'node') onNode = true; });
  chart.getZr().on('mouseup', () => { onNode = false; });
  vp.addEventListener('mousedown', e => {
    if (zoom <= 1 || onNode) return;
    drag = { x: e.clientX, y: e.clientY, sl: vp.scrollLeft, st: vp.scrollTop };
    vp.classList.add('panning');
  });
  window.addEventListener('mousemove', e => {
    if (!drag || onNode) return;
    vp.scrollLeft = drag.sl - (e.clientX - drag.x);
    vp.scrollTop = drag.st - (e.clientY - drag.y);
  });
  window.addEventListener('mouseup', () => {
    drag = null;
    onNode = false;
    vp.classList.remove('panning');
  });
}

function init() {
  chart = echarts.init(document.getElementById('chart'));
  let resizeTimer = null;
  window.addEventListener('resize', () => {
    const el = document.getElementById('chart');
    const vp = document.getElementById('chartViewport');
    if (zoom === 1) {
      el.style.width = ''; el.style.height = '';
      applyZoom.base = null;
    } else {
      // re-measure the unzoomed base and keep the zoomed size consistent
      el.style.width = ''; el.style.height = '';
      applyZoom.base = { w: vp.clientWidth, h: vp.clientHeight };
      el.style.width = Math.round(applyZoom.base.w * zoom) + 'px';
      el.style.height = Math.round(applyZoom.base.h * zoom) + 'px';
    }
    chart.resize();
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(render, 200); // re-pick responsive sizes
  });
  initZoom();

  const hormones = [...new Set(DATA.receptors.map(r => r.hormone))].sort();
  const genes = [...new Set(DATA.receptors.map(r => r.gene))].sort();

  makeMulti('hormoneName', () => hormones, 'hormoneNames', scheduleRender);
  makeMulti('hormoneTissue', sexTissues, 'hormoneTissues', scheduleRender);
  makeMulti('receptorTissue', sexTissues, 'receptorTissues', scheduleRender);
  makeMulti('receptorName', () => genes, 'receptorNames', scheduleRender);

  document.querySelectorAll('#sexSeg .seg-btn').forEach(btn => {
    btn.onclick = () => {
      document.querySelectorAll('#sexSeg .seg-btn').forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      state.sex = btn.dataset.value;
      // drop organ selections that belong to the other sex
      const valid = sexTissues();
      state.hormoneTissues = state.hormoneTissues.filter(t => valid.includes(t));
      state.receptorTissues = state.receptorTissues.filter(t => valid.includes(t));
      document.getElementById('hormoneTissue')._redraw();
      document.getElementById('receptorTissue')._redraw();
      render();
    };
  });

  document.getElementById('hormoneClass').onchange = e => { state.hormoneClass = e.target.value; render(); };
  document.getElementById('receptorClass').onchange = e => { state.receptorClass = e.target.value; render(); };

  const topN = document.getElementById('topN');
  topN.oninput = () => {
    state.topN = Number(topN.value);
    document.getElementById('topNVal').textContent = topN.value;
    scheduleRender();
  };

  const lo = document.getElementById('tpmLo'), hi = document.getElementById('tpmHi');
  const syncTpm = () => {
    let a = Number(lo.value), b = Number(hi.value);
    if (a > b) [a, b] = [b, a];
    state.tpmLo = a; state.tpmHi = b;
    document.getElementById('tpmVal').textContent = `${a.toFixed(2)} – ${b.toFixed(2)}`;
    scheduleRender();
  };
  lo.oninput = syncTpm;
  hi.oninput = syncTpm;

  document.getElementById('downloadBtn').onclick = downloadTSV;
  document.getElementById('resetBtn').onclick = () => {
    Object.assign(state, JSON.parse(JSON.stringify(DEFAULTS)));
    document.getElementById('hormoneClass').value = 'All';
    document.getElementById('receptorClass').value = 'All';
    topN.value = 3; document.getElementById('topNVal').textContent = '3';
    lo.value = 0.1; hi.value = 6.35;
    document.getElementById('tpmVal').textContent = '0.10 – 6.35';
    document.querySelectorAll('#sexSeg .seg-btn').forEach(b =>
      b.classList.toggle('active', b.dataset.value === 'Male'));
    ['hormoneName', 'hormoneTissue', 'receptorTissue', 'receptorName']
      .forEach(id => document.getElementById(id)._redraw());
    render();
  };

  // light/dark toggle (dark by default; the choice is remembered locally)
  const themeBtn = document.getElementById('themeToggle');
  const setThemeIcon = () => { themeBtn.textContent = isLight() ? '☀' : '☾'; };
  try {
    document.documentElement.dataset.theme = localStorage.getItem('hohc-theme') || 'dark';
  } catch (e) { document.documentElement.dataset.theme = 'dark'; }
  setThemeIcon();
  themeBtn.addEventListener('click', () => {
    const next = isLight() ? 'dark' : 'light';
    document.documentElement.dataset.theme = next;
    try { localStorage.setItem('hohc-theme', next); } catch (e) {}
    setThemeIcon();
    render(); // refresh chart text/tooltip colors
  });

  // logo click: pulse the waves, reset zoom and replay the sankey animation
  const logoBtn = document.getElementById('logoBtn');
  logoBtn.addEventListener('click', () => {
    logoBtn.classList.remove('pulse');
    void logoBtn.offsetWidth; // restart the CSS animation
    logoBtn.classList.add('pulse');
    applyZoom(1);
    chart.clear();
    render();
  });

  document.getElementById('legend').innerHTML =
    STAGES.map(s => `<span>${s}</span>`).join('');

  render();
}

loadData().then(init);
