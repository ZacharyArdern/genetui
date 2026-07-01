// ── Gene hover/click info bar ─────────────────────────────────────────────────
let selectedGene = null;
window.getSelectedGene = () => selectedGene;
let _statusActive = false;
let _dmndPickerActive = false;
let _dmndCompletions = [], _dmndCompIdx = 0;

function _getOrCreateDmndDropdown() {
  let el = document.getElementById('dmnd-completions-dropdown');
  if (!el) {
    el = document.createElement('div');
    el.id = 'dmnd-completions-dropdown';
    el.style.cssText = 'position:fixed;z-index:10000;background:#1c2128;border:1px solid #30363d;'
      + 'border-radius:4px;max-height:200px;overflow-y:auto;font-size:12px;font-family:monospace;display:none;';
    document.body.appendChild(el);
  }
  return el;
}

function _renderDmndCompletions(inp, res) {
  res.innerHTML = _dmndCompletions.map((p, i) => {
    const bg = i === _dmndCompIdx ? '#1f6feb' : 'transparent';
    return `<div data-ci="${i}" style="padding:3px 10px;cursor:pointer;background:${bg};color:#c9d1d9">${esc(p)}</div>`;
  }).join('');
  res.querySelectorAll('[data-ci]').forEach(el => {
    el.addEventListener('mousedown', ev => {
      ev.preventDefault();
      inp.value = _dmndCompletions[parseInt(el.dataset.ci)];
      res.style.display = 'none';
      inp.focus();
    });
  });
}

function _showDmndCompletions(list) {
  const inp = document.getElementById('dmnd-path-input');
  if (!inp) return;
  const res = _getOrCreateDmndDropdown();
  _dmndCompletions = list;
  _dmndCompIdx = 0;
  if (!list.length) { res.style.display = 'none'; return; }
  const rect = inp.getBoundingClientRect();
  res.style.left  = rect.left + 'px';
  res.style.top   = (rect.bottom + 2) + 'px';
  res.style.width = Math.max(rect.width, 300) + 'px';
  res.style.display = 'block';
  _renderDmndCompletions(inp, res);
}

function renderDmndPicker() {
  const bar = document.getElementById('gene-info-bar');
  bar.className = 'status';
  bar.innerHTML =
    `<span class="gi-status" style="color:#f0b429;white-space:nowrap">No DIAMOND DB \u2014 enter path:</span>`
    + `<input id="dmnd-path-input" type="text" placeholder="path/to/db.dmnd" autocomplete="off"
        style="flex:1;margin:0 6px;background:#161b22;color:#c9d1d9;border:1px solid #30363d;`
    + `border-radius:3px;padding:2px 6px;font-size:12px;font-family:monospace;min-width:0">`
    + `<span class="gi-hint" style="color:#6e7681;white-space:nowrap">Tab: complete \u00b7 Enter: confirm \u00b7 Esc: cancel</span>`;

  const inp = document.getElementById('dmnd-path-input');
  const res = _getOrCreateDmndDropdown();
  inp.focus();

  inp.addEventListener('input', () => {
    if (ws && ws.readyState === WebSocket.OPEN)
      ws.send(JSON.stringify({cmd: 'complete_path', prefix: inp.value}));
  });

  inp.addEventListener('keydown', ev => {
    if (ev.key === 'Tab') {
      ev.preventDefault();
      if (_dmndCompletions.length) {
        inp.value = _dmndCompletions[_dmndCompIdx];
        res.style.display = 'none';
        _dmndCompletions = [];
      }
      if (ws && ws.readyState === WebSocket.OPEN)
        ws.send(JSON.stringify({cmd: 'complete_path', prefix: inp.value}));
    } else if (ev.key === 'Enter') {
      ev.preventDefault();
      const path = (_dmndCompletions.length ? _dmndCompletions[_dmndCompIdx] : null) || inp.value.trim();
      res.style.display = 'none';
      if (path && ws && ws.readyState === WebSocket.OPEN)
        ws.send(JSON.stringify({cmd: 'set_dmnd', path}));
    } else if (ev.key === 'Escape') {
      ev.preventDefault();
      res.style.display = 'none';
      _dmndPickerActive = false;
      if (lastState) lastState.needs_dmnd_db = false;
      renderGeneInfo(selectedGene);
    } else if (ev.key === 'ArrowDown' && _dmndCompletions.length) {
      _dmndCompIdx = Math.min(_dmndCompIdx + 1, _dmndCompletions.length - 1);
      _renderDmndCompletions(inp, res); ev.preventDefault();
    } else if (ev.key === 'ArrowUp' && _dmndCompletions.length) {
      _dmndCompIdx = Math.max(_dmndCompIdx - 1, 0);
      _renderDmndCompletions(inp, res); ev.preventDefault();
    }
  });

  inp.addEventListener('blur', () => {
    setTimeout(() => { res.style.display = 'none'; }, 150);
  });
}

function onStatusUpdate() {
  const needsDb = !!(lastState && lastState.needs_dmnd_db);
  const msg     = (lastState && lastState.status_msg) || '';
  const msaErr  = (lastState && lastState.msa_error)  || '';
  const comps   = (lastState && lastState.path_completions) || [];

  if (needsDb) {
    _statusActive = false;
    if (!_dmndPickerActive) {
      _dmndPickerActive = true;
      renderDmndPicker();
    } else if (comps.length) {
      _showDmndCompletions(comps);
    }
    return;
  }

  _dmndPickerActive = false;
  const dd = document.getElementById('dmnd-completions-dropdown');
  if (dd) dd.style.display = 'none';

  if (msg) {
    _statusActive = true;
    const bar = document.getElementById('gene-info-bar');
    bar.className = 'status';
    bar.innerHTML = `<span class="gi-status">\u23F3 ${esc(msg)}</span>`;
  } else if (msaErr) {
    _statusActive = false;
    const bar = document.getElementById('gene-info-bar');
    bar.className = 'status';
    bar.innerHTML = `<span class="gi-status" style="color:#f85149">\u26a0 MSA error: ${esc(msaErr)}</span>`;
  } else {
    _statusActive = false;
    renderGeneInfo(selectedGene);
  }
}

function geneFromIdx(idx) {
  if (!lastState) return null;
  return lastState.features.find(f => f._idx === idx)
    || (lastState.blast_features || []).find(f => f._idx === idx) || null;
}

function renderGeneInfo(f) {
  const bar = document.getElementById('gene-info-bar');
  if (!f) {
    bar.className = 'empty';
    bar.innerHTML =
      '<span class="gi-hint" style="color:#484f58"><b style="color:#8b949e">/:search</b> &nbsp; <b style="color:#8b949e">d:panels</b> &nbsp; <b style="color:#8b949e">q:quit</b></span>'
      + '<span class="gi-hint">hover a gene for info</span>';
    return;
  }
  const len = f.end - f.start;
  const lenStr = len >= 1000 ? (len/1000).toFixed(1)+'kb' : len+'bp';
  const locus  = f.locus_tag && f.locus_tag !== f.name ? f.locus_tag : '';
  const kind   = f.is_orf ? 'ORF' : (f.noncoding ? 'nc' : 'gene');
  const strandGlyph = f.strand === '+' ? '＋' : '−';
  const canFold = !f.noncoding;

  bar.className = '';
  bar.innerHTML =
    `<span class="gi-name">${esc(f.name)}</span>`
    + (locus ? `<span class="gi-locus">${esc(locus)}</span>` : '')
    + `<span class="gi-strand">${strandGlyph}</span>`
    + `<span class="gi-coords">${fmt(f.start)}–${fmt(f.end)}</span>`
    + `<span class="gi-len">${lenStr}</span>`
    + `<span class="gi-kind">${kind}</span>`
    + `<span class="gi-hint" style="margin-left:auto">`
    + (canFold ? 'press <b>f</b> to fold · <b>m</b> for MSA' : 'non-coding')
    + '</span>';
}

// SVG event delegation — hover and click
(function() {
  const svgEl = document.getElementById('svg');

  svgEl.addEventListener('mousemove', e => {
    if (_statusActive) return;
    const g = e.target.closest('[data-idx]');
    if (!g) return;
    const idx = parseInt(g.dataset.idx, 10);
    const f = geneFromIdx(idx);
    if (f && (!selectedGene || selectedGene._idx !== idx)) {
      renderGeneInfo(f);
    }
  });

  svgEl.addEventListener('mouseleave', () => {
    if (_statusActive) return;
    if (!selectedGene) renderGeneInfo(null);
    else renderGeneInfo(selectedGene);
  });

  svgEl.addEventListener('click', e => {
    const g = e.target.closest('[data-idx]');
    if (!g) {
      selectedGene = null;
      renderGeneInfo(null);
      if (typeof updateRunMsaBtn === 'function') updateRunMsaBtn();
      return;
    }
    const idx = parseInt(g.dataset.idx, 10);
    const f = geneFromIdx(idx);
    if (!f) return;
    if (selectedGene && selectedGene._idx === f._idx) {
      selectedGene = null;
      renderGeneInfo(null);
    } else {
      selectedGene = f;
      renderGeneInfo(f);
    }
    if (typeof updateRunMsaBtn === 'function') updateRunMsaBtn();
  });
})();

// 'f' → fold, 'm' → MSA
document.addEventListener('keydown', e => {
  const tag = document.activeElement && document.activeElement.tagName;
  if (tag === 'INPUT' || tag === 'SELECT' || tag === 'TEXTAREA') return;
  const target = selectedGene;
  if (!target || target.noncoding) return;
  if (!ws || ws.readyState !== WebSocket.OPEN) return;
  if (e.key === 'f' || e.key === 'F') {
    ws.send(JSON.stringify({cmd: 'fold', gene: target.name}));
  } else if (e.key === 'm' || e.key === 'M') {
    ws.send(JSON.stringify({cmd: 'msa',  gene: target.name}));
  }
});

// Clear selection when new state arrives (genome changed)
(function() {
  const _orig = window.renderLocal;
  // Hook into state changes: reset selection if genome name changes
  let lastGenomeName = '';
  const orig_renderLocal = window.renderLocal;
  window.renderLocal = function() {
    if (lastState && lastState.genome_name !== lastGenomeName) {
      lastGenomeName = lastState.genome_name;
      selectedGene = null;
    }
    if (orig_renderLocal) orig_renderLocal.apply(this, arguments);
  };
})();
