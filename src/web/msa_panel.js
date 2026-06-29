// ── MSA panel — canvas-based alignment viewer ─────────────────────────────────
// Amino-acid colors: biotite Taylor color scheme, pastelled for dark UI
// (biotite.sequence.graphics.ProteinColoringScheme / Taylor 2002)
const MSA_TAYLOR = {
  A:'#CCFF00', R:'#0000FF', N:'#CC00FF', D:'#FF0000', C:'#FFFF00',
  Q:'#FF00CC', E:'#FF0066', G:'#FF9900', H:'#0066FF', I:'#66FF00',
  L:'#33FF00', K:'#6600FF', M:'#00FF00', F:'#00FF66', P:'#FFAA00',
  S:'#FF3300', T:'#FF6600', W:'#00CCFF', Y:'#00FFCC', V:'#99FF00',
};

// Mix each Taylor color 45% with white → pastel, readable on dark background
function pastelize(hex) {
  const r = parseInt(hex.slice(1,3),16), g = parseInt(hex.slice(3,5),16), b = parseInt(hex.slice(5,7),16);
  return `rgb(${Math.round(r*.55+255*.45)},${Math.round(g*.55+255*.45)},${Math.round(b*.55+255*.45)})`;
}
const MSA_COLORS = {};
for (const [aa, col] of Object.entries(MSA_TAYLOR)) MSA_COLORS[aa] = pastelize(col);

const MSA_CHAR_W  = 9;
const MSA_CHAR_H  = 14;
const MSA_NAME_W  = 160;
const MSA_GAP_COL = '#404855';

let msaRowOff = 0, msaColOff = 0;
let msaLastGene = '';

// ── Panel hover tracking ──────────────────────────────────────────────────────
(function(){
  const msaWrap = document.getElementById('msa-canvas-wrap');
  msaWrap.addEventListener('mouseenter', () => { window.hoveredPanel = 'msa'; });
  msaWrap.addEventListener('mouseleave', () => { window.hoveredPanel = 'genome'; });
})();

// ── Horizontal drag resize between msa-panel and fig-panel ───────────────────
(function(){
  const h = document.getElementById('msa-h-handle');
  const mp = document.getElementById('msa-panel');
  let drag = false, sx = 0, sw = 0;
  h.addEventListener('mousedown', e => {
    drag = true; sx = e.clientX; sw = mp.offsetWidth; e.preventDefault();
  });
  document.addEventListener('mousemove', e => {
    if (!drag) return;
    const row = document.getElementById('bottom-row');
    const maxW = row.clientWidth - 100;
    mp.style.width = Math.max(180, Math.min(sw + e.clientX - sx, maxW)) + 'px';
    renderMsa();
  });
  document.addEventListener('mouseup', () => { drag = false; });
})();

function msaVisible() {
  return document.getElementById('msa-panel').classList.contains('visible');
}

function renderMsa() {
  if (!msaVisible() || !lastState || !lastState.msa_sequences || lastState.msa_sequences.length === 0) return;

  const seqs  = lastState.msa_sequences; // [[id, seq], …]
  const nSeqs = seqs.length;
  const alnLen = seqs[0] ? seqs[0][1].length : 0;

  const wrap   = document.getElementById('msa-canvas-wrap');
  const canvas = document.getElementById('msa-canvas');
  const W = wrap.clientWidth  || 600;
  const H = wrap.clientHeight || 300;

  canvas.width  = W;
  canvas.height = H;
  canvas.style.width  = W + 'px';
  canvas.style.height = H + 'px';

  const ctx = canvas.getContext('2d');
  ctx.fillStyle = '#0d1117';
  ctx.fillRect(0, 0, W, H);

  const seqW    = W - MSA_NAME_W;
  const visCols = Math.floor(seqW / MSA_CHAR_W);
  const visRows = Math.floor((H - MSA_CHAR_H) / MSA_CHAR_H); // reserve bottom row for ruler

  msaColOff = Math.max(0, Math.min(msaColOff, Math.max(0, alnLen - visCols)));
  msaRowOff = Math.max(0, Math.min(msaRowOff, Math.max(0, nSeqs  - visRows)));

  ctx.font = `${MSA_CHAR_H - 3}px monospace`;
  ctx.textBaseline = 'top';

  for (let ri = 0; ri < visRows && ri + msaRowOff < nSeqs; ri++) {
    const [id, seq] = seqs[ri + msaRowOff];
    const y = ri * MSA_CHAR_H;

    // Name column
    ctx.fillStyle = ri % 2 ? '#0a0e18' : '#090d16';
    ctx.fillRect(0, y, MSA_NAME_W, MSA_CHAR_H);
    const nameStr = id.length > 22 ? id.slice(0, 21) + '…' : id;
    ctx.fillStyle = '#6888b8';
    ctx.fillText(nameStr, 4, y + 1);

    // Separator
    ctx.fillStyle = '#21262d';
    ctx.fillRect(MSA_NAME_W - 1, y, 1, MSA_CHAR_H);

    // Sequence characters
    for (let ci = 0; ci < visCols && ci + msaColOff < seq.length; ci++) {
      const c   = seq[ci + msaColOff].toUpperCase();
      const col = MSA_COLORS[c];
      const x   = MSA_NAME_W + ci * MSA_CHAR_W;
      if (col) {
        ctx.fillStyle = col;
        ctx.fillRect(x, y, MSA_CHAR_W, MSA_CHAR_H);
        ctx.fillStyle = '#111';
      } else {
        ctx.fillStyle = MSA_GAP_COL;
      }
      ctx.fillText(c, x + 1, y + 1);
    }
  }

  // Column ruler
  const rulerY = H - MSA_CHAR_H;
  ctx.fillStyle = '#0a0e18';
  ctx.fillRect(0, rulerY, W, MSA_CHAR_H);
  ctx.fillStyle = '#484f58';
  ctx.font = '9px monospace';
  for (let ci = 0; ci < visCols; ci++) {
    const absCol = ci + msaColOff + 1;
    if (absCol % 10 === 0) {
      const x = MSA_NAME_W + ci * MSA_CHAR_W;
      ctx.fillStyle = '#21262d';
      ctx.fillRect(x, rulerY, 1, MSA_CHAR_H);
      ctx.fillStyle = '#484f58';
      ctx.fillText(absCol, x + 1, rulerY + 2);
    }
  }

  document.getElementById('msa-status').textContent =
    `${nSeqs} seqs · ${alnLen} cols · row ${msaRowOff + 1} · col ${msaColOff + 1}`;
  document.getElementById('msa-gene-name').textContent = lastState.msa_gene || '';
}

// ── Scroll handling — wheel ───────────────────────────────────────────────────
document.getElementById('msa-canvas-wrap').addEventListener('wheel', e => {
  e.preventDefault();
  e.stopPropagation();
  if (e.shiftKey || Math.abs(e.deltaX) > Math.abs(e.deltaY)) {
    msaColOff += Math.sign(e.deltaX || e.deltaY) * Math.max(1, Math.round(Math.abs(e.deltaX || e.deltaY) / MSA_CHAR_W));
  } else {
    msaRowOff += Math.sign(e.deltaY) * Math.max(1, Math.round(Math.abs(e.deltaY) / MSA_CHAR_H));
  }
  renderMsa();
}, {passive: false});

// ── Scroll handling — keyboard (only when mouse is over MSA panel) ────────────
document.addEventListener('keydown', e => {
  if (window.hoveredPanel !== 'msa') return;
  if (!msaVisible()) return;
  const tag = document.activeElement && document.activeElement.tagName;
  if (tag === 'INPUT' || tag === 'SELECT' || tag === 'TEXTAREA') return;
  let handled = true;
  switch(e.key) {
    case 'ArrowDown':  msaRowOff += 1;  break;
    case 'ArrowUp':    msaRowOff -= 1;  break;
    case 'ArrowRight': msaColOff += 5;  break;
    case 'ArrowLeft':  msaColOff -= 5;  break;
    case 'PageDown':   msaRowOff += 10; break;
    case 'PageUp':     msaRowOff -= 10; break;
    case 'Home':       msaColOff  = 0;  break;
    case 'End':        msaColOff  = 1e9; break;
    default: handled = false;
  }
  if (handled) { e.preventDefault(); e.stopPropagation(); renderMsa(); }
});

// ── State update hook (called from websocket.js) ──────────────────────────────
function onMsaStateUpdate() {
  if (!lastState) return;
  const gene = lastState.msa_gene || '';
  const seqs = lastState.msa_sequences;
  if (gene && seqs && seqs.length > 0) {
    if (gene !== msaLastGene) { msaRowOff = 0; msaColOff = 0; msaLastGene = gene; }
    const panel = document.getElementById('msa-panel');
    const handle = document.getElementById('msa-h-handle');
    panel.classList.add('visible');
    handle.classList.add('visible');
    if (window._panelCbs && window._panelCbs['msa-panel']) {
      window._panelCbs['msa-panel'].disabled = false;
      window._panelCbs['msa-panel'].checked  = true;
    }
    renderMsa();
  }
}

// ── Resize observer ───────────────────────────────────────────────────────────
(function(){
  if (!window.ResizeObserver) return;
  new ResizeObserver(() => renderMsa())
    .observe(document.getElementById('msa-canvas-wrap'));
})();
