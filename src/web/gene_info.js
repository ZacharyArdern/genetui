// ── Gene hover/click info bar ─────────────────────────────────────────────────
let selectedGene = null;
let _statusActive = false;

function onStatusUpdate() {
  const msg = (lastState && lastState.status_msg) || '';
  _statusActive = !!msg;
  if (msg) {
    const bar = document.getElementById('gene-info-bar');
    bar.className = 'status';
    bar.innerHTML = `<span class="gi-status">\u23F3 ${esc(msg)}</span>`;
  } else {
    renderGeneInfo(selectedGene);
  }
}

function geneFromIdx(idx) {
  if (!lastState) return null;
  return lastState.features.find(f => f._idx === idx) || null;
}

function renderGeneInfo(f) {
  const bar = document.getElementById('gene-info-bar');
  if (!f) {
    bar.className = 'empty';
    bar.innerHTML = '<span class="gi-name">—</span>'
      + '<span class="gi-hint">hover or click a gene in the track above</span>';
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
