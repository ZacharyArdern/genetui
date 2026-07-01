// ── Fancy arrows: dna_features_viewer SVG tile system ─────────────────────────
// Each tile covers ~2-3× the viewport span.  Rendered by Pyodide at the actual
// display pixel width (capped at 1600 px for speed), so fonts are always 1:1.
// Cache is keyed by genome + tile range + gencode + scale-band + mode.

const FANCY_TILE_MULT = 2;    // tile covers ≥ MULT × viewport
const FANCY_RENDER_MAX = 900; // max pixel width sent to matplotlib (keep SVG small)
const FANCY_CACHE_MAX  = 50;
// Matches drawSixFrame totalBottom (without nucleotide tracks):
// 3*(FH+GH)*2 strands + SEP + 2*FH(nc) = 3*(19+13)*2 + 4 + 38 = 234
const FANCY_TILE_H = 3*(19+13)*2 + 4 + 19*2; // = 234

const _fancyCache = new Map();
let _fancyTicket  = 0;
let _fancyLastKey = ''; // composite invalidation key

// ── Helpers ───────────────────────────────────────────────────────────────────

function getFancyMode() { return false; }
window.getFancyMode = getFancyMode;

// Quantise to nearest 12% step in log2 space (~2 steps per octave)
function _scaleBand(scale) {
  return Math.round(Math.log2(Math.max(scale, 1e-9)) * 8);
}

function _snapTile(vs, ve) {
  const span = Math.max(ve - vs, 1);
  const raw  = span * FANCY_TILE_MULT;
  const mag  = Math.pow(10, Math.floor(Math.log10(raw)));
  let ts = mag * 10;
  for (const m of [1, 2, 5, 10]) { if (mag * m >= raw) { ts = mag * m; break; } }
  return { tileStart: Math.floor(vs / ts) * ts, tileEnd: Math.floor(vs / ts) * ts + ts, tileSpan: ts };
}

function _tileKey(g, ts, te, gc, sb, sf) {
  return `${g}|${ts}|${te}|${gc}|${sb}|${sf}`;
}

// ── Cache ─────────────────────────────────────────────────────────────────────

function _evict() {
  if (_fancyCache.size <= FANCY_CACHE_MAX) return;
  [..._fancyCache.entries()]
    .sort((a, b) => a[1].usedAt - b[1].usedAt)
    .slice(0, _fancyCache.size - FANCY_CACHE_MAX)
    .forEach(([k, v]) => { v.el && v.el.parentNode && v.el.parentNode.removeChild(v.el); _fancyCache.delete(k); });
}

function _clearFancyCache() {
  document.getElementById('fancy-layer') && (document.getElementById('fancy-layer').innerHTML = '');
  _fancyCache.clear();
  _fancyTicket++;
}

// ── Positioning ───────────────────────────────────────────────────────────────

function _positionTile(entry, scale, vs, panelW) {
  const el = entry.el; if (!el) return;
  const left = (entry.tileStart - vs) * scale;
  const wid  = (entry.tileEnd - entry.tileStart) * scale;
  const top  = typeof window._fancyYGeneTop === 'number' ? window._fancyYGeneTop : RULER_H;
  el.style.left    = left + 'px';
  el.style.width   = wid  + 'px';
  el.style.top     = top  + 'px';
  el.style.height  = FANCY_TILE_H + 'px';
  el.style.display = (left > panelW || left + wid < 0) ? 'none' : '';
  // Re-sync SVG width whenever tile is repositioned (display width changes with zoom)
  const svgEl = el.querySelector('svg');
  if (svgEl) svgEl.style.width = wid + 'px';
}

function _repositionAll() {
  const panel = document.getElementById('live-panel'); if (!panel) return;
  const W = panel.clientWidth || window.innerWidth;
  const scale = W / Math.max(localVE - localVS, 1);
  for (const e of _fancyCache.values())
    if (e.el && e.el.parentNode) _positionTile(e, scale, localVS, W);
}

// ── Python builder ────────────────────────────────────────────────────────────

function _buildPy(tileStart, tileEnd, renderW, features, sixFrame, sequence, seqStart, gencode) {
  const seqLen = tileEnd - tileStart;
  const figW   = (renderW / 96).toFixed(3);
  const gc     = gencode || '1';

  const visC = features.filter(f => f.end > tileStart && f.start < tileEnd && !f.noncoding);
  const visN = features.filter(f => f.end > tileStart && f.start < tileEnd &&  f.noncoding);

  const fd  = JSON.stringify(visC.map(f => ({
    s:  Math.max(0,      f.start - tileStart),
    e:  Math.min(seqLen, f.end   - tileStart),
    d:  f.strand === '+' ? 1 : -1,
    c:  f.color || '#4a7fbf',
    l:  f.name,
    fr: f.strand === '+' ? (f.start - 1) % 3 : (f.end - 1) % 3,
  })));
  const nd  = JSON.stringify(visN.map(f => ({
    s: Math.max(0,      f.start - tileStart),
    e: Math.min(seqLen, f.end   - tileStart),
    d: f.strand === '+' ? 1 : -1,
    c: '#6e7681', l: f.name,
  })));

  // Sequence slice within tile bounds
  const seqEnd   = seqStart + (sequence ? sequence.length : 0);
  const oStart   = Math.max(tileStart, seqStart);
  const oEnd     = Math.min(tileEnd,   seqEnd);
  const seqSlice = (oEnd > oStart && sequence) ? sequence.slice(oStart - seqStart, oEnd - seqStart) : '';
  const seqOff   = oStart - tileStart;

  const stopArr  = JSON.stringify([...getStopSet()]);

  // Each JS frame row = FH=19 px = 19/96 ≈ 0.198 in; match that height so tile ≈ JS track size.
  // Labels are skipped (too slow; gene names visible on hover). Empty rows skip plot() entirely.
  const FH_IN = (19 / 96).toFixed(4); // inches per frame row

  if (sixFrame) {
    return `
import json as _j,io as _io,matplotlib as _ml;_ml.use('Agg')
import matplotlib.pyplot as _plt,matplotlib.gridspec as _gs
from dna_features_viewer import GraphicFeature,GraphicRecord as _GR
_fd=${JSON.stringify(fd)};_nd=${JSON.stringify(nd)}
_seq=${JSON.stringify(seqSlice)};_soff=${seqOff};_slen=${seqLen};_W=${figW}
_STOP=set(_j.loads(${JSON.stringify(stopArr)}))
_CMP={'A':'T','T':'A','G':'C','C':'G','a':'t','t':'a','g':'c','c':'g'}
def _rc(s): return ''.join(_CMP.get(c,'N') for c in reversed(s))
_fs=[[],[],[]];_rs=[[],[],[]]
if _seq:
    _u=_seq.upper();_rcu=_rc(_seq)
    for _i in range(len(_u)-2):
        if _u[_i:_i+3] in _STOP: _gp=_soff+_i;_fs[_gp%3].append(_gp)
    for _i in range(len(_rcu)-2):
        if _rcu[_i:_i+3] in _STOP: _ge=_soff+len(_seq)-1-_i;_rs[_ge%3].append(_ge-2)
_feats=_j.loads(_fd);_nc=_j.loads(_nd)
_ROWS=[('+1',1,0),('+2',1,1),('+3',1,2),('-1',-1,0),('-2',-1,1),('-3',-1,2)]
_FH=${FH_IN}
_nc_h=_FH*0.85
_all_h=[_nc_h,_FH,_FH,_FH,_FH*0.15,_FH,_FH,_FH,_nc_h]
_fig=_plt.figure(figsize=(_W,sum(_all_h)));_fig.patch.set_alpha(0)
_gsp=_gs.GridSpec(9,1,figure=_fig,hspace=0,height_ratios=_all_h,top=1.0,bottom=0.0,left=0.07,right=1.0)
for _ri,(_lbl,_si,_fr) in enumerate(_ROWS):
    _idx=_ri+1 if _ri<3 else _ri+2
    _ax=_fig.add_subplot(_gsp[_idx])
    _gfs=[GraphicFeature(start=f['s'],end=f['e'],strand=f['d'],color=f['c'])
          for f in _feats if f['d']==_si and f['fr']==_fr]
    if _gfs:
        _GR(sequence_length=_slen,features=_gfs).plot(ax=_ax,with_ruler=False,draw_line=True)
    else:
        _ax.set_xlim(0,_slen);_ax.axhline(0,color='#30363d',lw=0.5)
    _ax.set_ylim(-0.55,0.55);_ax.patch.set_alpha(0);_ax.axis('off')
    _ax.text(0,0.5,_lbl,transform=_ax.transAxes,fontsize=5,color='#6e7681',va='top')
    for _x in (_fs[_fr] if _si==1 else _rs[_fr]): _ax.axvline(_x,color='#f85149',alpha=0.7,lw=0.6)
_sep=_fig.add_subplot(_gsp[4]);_sep.set_visible(False)
for _gi,_si,_lbl in [(0,1,'nc+'),(8,-1,'nc-')]:
    _ax=_fig.add_subplot(_gsp[_gi])
    _gfs=[GraphicFeature(start=f['s'],end=f['e'],strand=f['d'],color=f['c']) for f in _nc if f['d']==_si]
    if _gfs:
        _GR(sequence_length=_slen,features=_gfs).plot(ax=_ax,with_ruler=False,draw_line=True)
    else:
        _ax.set_xlim(0,_slen);_ax.axhline(0,color='#30363d',lw=0.5)
    _ax.set_ylim(-0.55,0.55);_ax.patch.set_alpha(0);_ax.axis('off')
    _ax.text(0,0.5,_lbl,transform=_ax.transAxes,fontsize=5,color='#6e7681',va='top')
_b=_io.BytesIO()
_fig.savefig(_b,format='svg',bbox_inches='tight',pad_inches=0,transparent=True)
_plt.close('all');_b.getvalue().decode()
`.trim();
  }

  // Simple mode — single pass, no labels
  return `
import json as _j,io as _io,matplotlib as _ml;_ml.use('Agg')
import matplotlib.pyplot as _plt
from dna_features_viewer import GraphicFeature,GraphicRecord as _GR
_fd=${JSON.stringify(fd)};_nd=${JSON.stringify(nd)};_slen=${seqLen};_W=${figW}
_fc=[GraphicFeature(start=f['s'],end=f['e'],strand=f['d'],color=f['c']) for f in _j.loads(_fd)]
_fn=[GraphicFeature(start=f['s'],end=f['e'],strand=f['d'],color=f['c']) for f in _j.loads(_nd)]
_all=_fc+_fn
_fig_h=${FH_IN}*3 if _all else ${FH_IN}
_ax,_=_GR(sequence_length=_slen,features=_all).plot(figure_width=_W,figure_height=_fig_h)
_ax.patch.set_alpha(0);_ax.figure.patch.set_alpha(0);_ax.axis('off')
_b=_io.BytesIO()
_ax.figure.savefig(_b,format='svg',bbox_inches='tight',pad_inches=0,transparent=True)
_plt.close('all');_b.getvalue().decode()
`.trim();
}

// ── Fetch a tile ──────────────────────────────────────────────────────────────

async function _fetchTile(ts, te, gName, gc, scale, sb, sixFrame, feats, seq, seqStart) {
  if (te <= ts) return;
  const key = _tileKey(gName, ts, te, gc, sb, sixFrame ? 1 : 0);
  if (_fancyCache.has(key)) { _fancyCache.get(key).usedAt = Date.now(); return; }

  // Render at actual display width, capped for speed
  const displayW = Math.round((te - ts) * scale);
  const renderW  = Math.min(FANCY_RENDER_MAX, Math.max(200, displayW));

  const ticket  = _fancyTicket;
  let svgStr;
  try {
    svgStr = await pyodide.runPythonAsync(
      _buildPy(ts, te, renderW, feats, sixFrame, seq, seqStart, gc));
  } catch (e) { console.warn('fancy tile:', e); return; }
  if (_fancyTicket !== ticket || !svgStr) return;

  const div = document.createElement('div');
  div.className = 'ftile';
  div.innerHTML = svgStr;
  const svgEl = div.querySelector('svg');
  if (svgEl) {
    // Fix height explicitly so SVG doesn't scale vertically when displayed wider than rendered.
    // preserveAspectRatio=none allows horizontal stretch without distorting height.
    svgEl.setAttribute('preserveAspectRatio', 'none');
    svgEl.style.cssText = `width:${displayW}px;height:${FANCY_TILE_H}px;display:block;`;
  }

  const entry = { el: div, tileStart: ts, tileEnd: te, usedAt: Date.now() };
  _fancyCache.set(key, entry);
  _evict();

  const layer = document.getElementById('fancy-layer');
  if (!layer) return;
  layer.appendChild(div);

  const panel = document.getElementById('live-panel');
  const W = panel ? panel.clientWidth || window.innerWidth : window.innerWidth;
  _positionTile(entry, scale, localVS, W);
}

// ── Public API ────────────────────────────────────────────────────────────────

window.updateFancyLayer = async function() {
  const layer = document.getElementById('fancy-layer');
  if (!layer) return;
  if (!getFancyMode()) { layer.style.display = 'none'; return; }
  layer.style.display = '';

  _repositionAll();
  if (!lastState || !pyReady) return;

  const panel    = document.getElementById('live-panel');
  const W        = panel ? panel.clientWidth || window.innerWidth : window.innerWidth;
  const H        = panel ? panel.clientHeight || 300 : 300;
  const vs       = localVS, ve = localVE;
  const scale    = W / Math.max(ve - vs, 1);
  const sb       = _scaleBand(scale);
  const gName    = lastState.genome_name || '';
  const gc       = (document.getElementById('opt-gencode') || {}).value || '1';
  const sixFrame = !!(document.getElementById('opt-sixframe') || {}).checked;
  const gs       = lastState.genome_size || 0;
  const feats    = lastState.features    || [];
  const seq      = lastState.sequence    || '';
  const seqStart = lastState.seq_start   || 0;

  // Only clear cache when genome or frame mode changes (not on zoom — old tiles stay
  // visible while new ones render, preventing the "disappears on zoom" flicker).
  const ck = `${gName}|${sixFrame}`;
  if (ck !== _fancyLastKey) { _clearFancyCache(); _fancyLastKey = ck; }

  _repositionAll();

  const { tileStart, tileEnd, tileSpan } = _snapTile(vs, ve);
  const cs = Math.max(0, tileStart);
  const ce = gs > 0 ? Math.min(gs, tileEnd) : tileEnd;

  // Fetch current tile (await — ensures it's visible before prefetch fires)
  await _fetchTile(cs, ce, gName, gc, scale, sb, sixFrame, feats, seq, seqStart);

  // Prefetch neighbours
  const ps = Math.max(0, cs - tileSpan);
  if (cs > ps) _fetchTile(ps, cs, gName, gc, scale, sb, sixFrame, feats, seq, seqStart);
  const ne = gs > 0 ? Math.min(gs, ce + tileSpan) : ce + tileSpan;
  if (ne > ce) _fetchTile(ce, ne, gName, gc, scale, sb, sixFrame, feats, seq, seqStart);
};

window.invalidateFancyCache = _clearFancyCache;

// ── Main browser options dropdown ─────────────────────────────────────────────
(function() {
  const btn = document.getElementById('btn-main-browser-opts');
  const box = document.getElementById('main-browser-opts-box');
  if (!btn || !box) return;
  let open = false;
  function pos() {
    const r = btn.getBoundingClientRect();
    box.style.left = Math.max(4, r.right - box.offsetWidth) + 'px';
    box.style.top  = (r.bottom + 4) + 'px';
  }
  btn.addEventListener('click', e => {
    e.stopPropagation(); open = !open;
    box.style.display = open ? 'block' : 'none';
    btn.classList.toggle('active', open);
    if (open) pos();
  });
  document.addEventListener('click', e => {
    if (open && !box.contains(e.target) && e.target !== btn) {
      open = false; box.style.display = 'none'; btn.classList.remove('active');
    }
  });
})();

// ── Arrow style / six-frame change → clear cache ──────────────────────────────
(function() {
  ['opt-arrow-style', 'opt-sixframe'].forEach(id => {
    const el = document.getElementById(id);
    if (!el) return;
    el.addEventListener('change', () => {
      _clearFancyCache();
      if (typeof renderLocal === 'function') renderLocal();
      if (window.updateFancyLayer) window.updateFancyLayer();
    });
  });
})();
