// ── Anchored draggable dropdown factory ──────────────────────────────────────
function _makeDropdown(btnId, boxId) {
  const btn=document.getElementById(btnId);
  const box=document.getElementById(boxId);
  if (!btn||!box) return;
  let open=false, dragging=false, ox=0, oy=0, mx=0, my=0;
  function positionBox() {
    const r=btn.getBoundingClientRect();
    let left=r.right-box.offsetWidth;
    if (left<4) left=4;
    box.style.left=left+'px'; box.style.top=(r.bottom+4)+'px';
    box.style.right='auto'; box.style.bottom='auto';
  }
  btn.addEventListener('click', e=>{
    e.stopPropagation();
    open=!open; box.style.display=open?'block':'none';
    btn.classList.toggle('active',open);
    if (open) positionBox();
  });
  document.addEventListener('click', e=>{
    if (open && !box.contains(e.target) && e.target!==btn) {
      open=false; box.style.display='none'; btn.classList.remove('active');
    }
  });
  const h=box.querySelector('h3');
  if (h) h.addEventListener('mousedown', e=>{
    dragging=true; mx=e.clientX; my=e.clientY;
    const r=box.getBoundingClientRect(); ox=r.left; oy=r.top;
    e.preventDefault();
  });
  document.addEventListener('mousemove', e=>{
    if (!dragging) return;
    box.style.left=(ox+e.clientX-mx)+'px'; box.style.top=(oy+e.clientY-my)+'px';
  });
  document.addEventListener('mouseup', ()=>{dragging=false;});
}

_makeDropdown('btn-opts',  'opts-box');
_makeDropdown('btn-files', 'files-box');

// ── BAM / codon colour state ──────────────────────────────────────────────────
window._bamColorOverrides = {};  // idx -> hex (user set)
window._codonColors = {};        // codon -> hex (from codon-lines field)
window._codonAutoColors = {};    // codon -> hex (auto-assigned for BAM-name codons not in codon-lines)

const _DEFAULT_PAL = ['#3fc9d0','#56d364','#a371f7','#f85149','#e3b341','#d66adc'];

function _nextPalColor() {
  const used = new Set([...Object.values(_codonColors), ...Object.values(_codonAutoColors)]);
  return _DEFAULT_PAL.find(c => !used.has(c)) || _DEFAULT_PAL[Object.keys(_codonAutoColors).length % _DEFAULT_PAL.length];
}

function _codonFromLabel(label) {
  const m = (label||'').match(/(?:^|[^A-Za-z])([ACGTacgt]{3})(?:[^A-Za-z]|$)/);
  return m ? m[1].toUpperCase() : null;
}

window.getBamEffectiveColor = function(idx, label) {
  if (_bamColorOverrides[idx] != null) return _bamColorOverrides[idx];
  const codon = _codonFromLabel(label);
  if (codon) {
    if (_codonColors[codon]) return _codonColors[codon];
    // Auto-assign a stable distinct colour for this codon
    if (!_codonAutoColors[codon]) _codonAutoColors[codon] = _nextPalColor();
    return _codonAutoColors[codon];
  }
  return _DEFAULT_PAL[idx % _DEFAULT_PAL.length];
};

function _syncCodonColsFromBams(labels) {
  let changed = false;
  (labels||[]).forEach(lbl => {
    const codon = _codonFromLabel(lbl);
    if (codon && !_codonColors[codon]) {
      const col = _codonAutoColors[codon] || _nextPalColor();
      _codonAutoColors[codon] = col;
      _codonColors[codon] = col;
      changed = true;
    }
  });
  if (changed) {
    document.getElementById('fig-codon-lines').value =
      Object.entries(_codonColors).map(([c,v]) => `${c}:${v}`).join(', ');
    _renderCodonColPickers();
    if (typeof scheduleAutoRender === 'function') scheduleAutoRender();
  }
}

function _parseFigCodonCols() {
  const txt = document.getElementById('fig-codon-lines').value;
  const entries = [];
  txt.split(/[,;]+/).forEach(e => {
    const m = e.trim().match(/^([A-Za-z]{3})(?:[:\s]+([#\w]+))?$/);
    if (m) entries.push([m[1].toUpperCase(), m[2] || null]);
  });
  // Assign distinct palette colours to codons without explicit colour
  window._codonColors = {};
  let palIdx = 0;
  entries.forEach(([cdn, col]) => {
    if (col) {
      window._codonColors[cdn] = col;
    } else {
      // Pick next palette color not already used by other codon lines entries
      const usedSoFar = new Set(Object.values(window._codonColors));
      const auto = _DEFAULT_PAL.find(c => !usedSoFar.has(c)) || _DEFAULT_PAL[palIdx % _DEFAULT_PAL.length];
      window._codonColors[cdn] = auto;
      palIdx++;
    }
    // If this codon was auto-assigned before, promote it to explicit
    delete _codonAutoColors[cdn];
  });
  _renderCodonColPickers();
}

function _renderCodonColPickers() {
  const el = document.getElementById('fig-codon-col-pickers');
  if (!el) return;
  const entries = Object.entries(window._codonColors);
  if (!entries.length) { el.innerHTML = ''; return; }
  el.innerHTML = '<div style="display:flex;flex-wrap:wrap;gap:6px 14px;margin-top:4px">'
    + entries.map(([cdn, col]) =>
        `<div style="display:flex;align-items:center;gap:4px">`
        + `<input type="color" value="${col}" data-cdn="${cdn}"`
        + ` style="width:20px;height:16px;padding:0;border:none;background:none;cursor:pointer">`
        + `<span style="font-size:11px;font-family:monospace;color:#c9d1d9">${cdn}</span>`
        + `</div>`
      ).join('')
    + '</div>';
  el.querySelectorAll('input[type=color]').forEach(inp => {
    inp.addEventListener('input', function() {
      window._codonColors[this.dataset.cdn] = this.value;
      document.getElementById('fig-codon-lines').value =
        Object.entries(window._codonColors).map(([c,v]) => `${c}:${v}`).join(', ');
      scheduleAutoRender();
    });
  });
}

// ── Option change listeners ───────────────────────────────────────────────────
function optionChanged() { renderLocal(); scheduleAutoRender(); }
['opt-sixframe','opt-nucleotides','opt-stopcodons','fig-ruler','fig-white','fig-show-labels','fig-hide-long-labels'].forEach(id=>{
  document.getElementById(id).addEventListener('change', optionChanged);
});
document.getElementById('opt-gencode').addEventListener('change', optionChanged);
document.getElementById('opt-max-span').addEventListener('change', () => {
  clampLocal(); renderLocal(); scheduleAutoRender();
});
[['fig-width','lbl-width',1],['fig-fs','lbl-fs',0],['fig-overflow-thresh','lbl-overflow-thresh',3]].forEach(([id,lb,d])=>{
  document.getElementById(id).addEventListener('input', function(){
    document.getElementById(lb).textContent=parseFloat(this.value).toFixed(d);
    scheduleAutoRender();
  });
});
document.getElementById('fig-hide-long-labels').addEventListener('change', function(){
  document.getElementById('overflow-thresh-row').style.display = this.checked ? '' : 'none';
});
document.getElementById('cov-style-fig').addEventListener('change', function(){
  scheduleAutoRender();
});
document.getElementById('fig-cov-layout').addEventListener('change', function(){
  scheduleAutoRender();
});
document.getElementById('fig-cov-log').addEventListener('change', function(){
  scheduleAutoRender();
});
document.getElementById('fig-cov-ymax').addEventListener('input', function(){
  scheduleAutoRender();
});
document.getElementById('cov-style-svg').addEventListener('change', function(){
  const hr = document.getElementById('cov-height-row');
  if (hr) hr.style.display = this.value !== 'none' ? '' : 'none';
  renderLocal();
});
// Apply initial state for cov-height-row since default is now "reads"
(function(){ const hr=document.getElementById('cov-height-row'); if(hr) hr.style.display=''; })();
document.getElementById('cov-height').addEventListener('input', function(){
  const lbl = document.getElementById('lbl-cov-height');
  if (lbl) lbl.textContent = this.value;
  renderLocal();
});
document.getElementById('opt-codons').addEventListener('input', function(){
  renderLocal();
});
document.getElementById('fig-codon-lines').addEventListener('input', function(){
  _parseFigCodonCols();
  scheduleAutoRender();
});
_parseFigCodonCols();

// ── Files panel ──────────────────────────────────────────────────────────────
function renderFilesPanel() {
  if (lastState) _syncCodonColsFromBams(lastState.bam_labels);
  const el = document.getElementById('loaded-files-list');
  if (!el || !lastState) return;
  let html = '';
  if (lastState.genome_label)
    html += `<div style="color:#8b949e">🧬 ${esc(lastState.genome_label)}</div>`;
  if (lastState.gff_label)
    html += `<div style="color:#8b949e">➡️ ${esc(lastState.gff_label)}</div>`;
  window._extraGffColors = window._extraGffColors || {};
  const egLabels = lastState.extra_gff_labels || [];
  egLabels.forEach((lbl, i) => {
    if (!window._extraGffColors[i]) window._extraGffColors[i] = '#3fb950';
    const col = window._extraGffColors[i];
    html += `<div style="display:flex;align-items:center;gap:4px;margin:2px 0">`
      + `<input type="color" value="${col}" title="Pick colour" style="width:22px;height:18px;padding:0;border:none;background:none;cursor:pointer"`
      + ` onchange="window._extraGffColors[${i}]=this.value;scheduleAutoRender&&scheduleAutoRender()">`
      + `<input type="text" value="${col}" title="Hex or matplotlib colour name" placeholder="#hex or name"`
      + ` style="width:80px;font-size:10px;background:#161b22;color:#c9d1d9;border:1px solid #30363d;border-radius:3px;padding:1px 4px"`
      + ` onchange="window._extraGffColors[${i}]=this.value;scheduleAutoRender&&scheduleAutoRender()" oninput="this.previousElementSibling&&(this.previousElementSibling.value=(/^#[0-9a-fA-F]{6}$/.test(this.value)?this.value:this.previousElementSibling.value))">`
      + `<span style="color:#c9d1d9;font-size:11px">➡️ ${esc(lbl)}</span></div>`;
  });
  if (egLabels.length > 0) {
    const ov = lastState.gff_overlay !== false;
    html += `<div style="margin-left:0;margin-top:2px">`
      + `<select style="font-size:10px;background:#161b22;color:#c9d1d9;border:1px solid #30363d;border-radius:3px;padding:1px 4px"`
      + ` onchange="ws&&ws.readyState===1&&ws.send(JSON.stringify({cmd:'set_gff_overlay',overlay:this.value==='overlay'}))">`
      + `<option value="overlay"${ov?' selected':''}>Overlay</option>`
      + `<option value="track"${!ov?' selected':''}>Separate track</option>`
      + `</select></div>`;
  }
  const bamLabels = lastState.bam_labels || [];
  const mask = lastState.active_bam_mask || bamLabels.map(() => true);
  bamLabels.forEach((lbl, i) => {
    if (i > 0) html += `<hr style="border:none;border-top:1px solid #30363d;margin:3px 0">`;
    const chk = mask[i] !== false ? 'checked' : '';
    const col = window.getBamEffectiveColor(i, lbl);
    html += `<div style="display:flex;align-items:center;gap:4px">`
      + `<input type="color" value="${col}"`
      + ` style="width:20px;height:16px;padding:0;border:none;background:none;cursor:pointer"`
      + ` onchange="window._bamColorOverrides[${i}]=this.value;renderLocal&&renderLocal();scheduleAutoRender&&scheduleAutoRender()">`
      + `<label style="color:#c9d1d9;cursor:pointer;flex:1">`
      + `<input type="checkbox" ${chk} style="margin-right:4px"`
      + ` onchange="ws&&ws.readyState===1&&ws.send(JSON.stringify({cmd:'toggle_bam',idx:${i}}))">`
      + `💥 ${esc(lbl)}</label></div>`;
  });
  el.innerHTML = html;
}
window.renderFilesPanel = renderFilesPanel;

function openAddFileOverlay() {
  let overlay = document.getElementById('add-file-overlay');
  if (!overlay) {
    overlay = document.createElement('div');
    overlay.id = 'add-file-overlay';
    overlay.style.cssText = 'display:none;position:fixed;inset:0;z-index:600;background:rgba(0,0,0,0.55);align-items:flex-start;justify-content:center';
    overlay.innerHTML =
      `<div style="margin-top:120px;width:520px;background:#161b22;border:1px solid #58a6ff;`
      + `border-radius:8px;overflow:visible;box-shadow:0 8px 32px rgba(0,0,0,0.7);padding:16px">`
      + `<div style="color:#c9d1d9;font-size:13px;font-family:monospace;margin-bottom:10px">Add file (BAM / GFF)</div>`
      + `<input id="add-file-input" type="text" placeholder="path/to/file.bam or .gff" autocomplete="off"`
      + ` style="width:100%;box-sizing:border-box;background:#0d1117;color:#e6edf3;border:1px solid #58a6ff;`
      + `border-radius:4px;padding:6px 10px;font-size:13px;font-family:monospace">`
      + `<div id="add-file-comps" style="position:absolute;left:16px;right:16px;background:#1c2128;`
      + `border:1px solid #58a6ff;border-radius:4px;max-height:200px;overflow-y:auto;`
      + `font-size:12px;font-family:monospace;display:none;z-index:10"></div>`
      + `<div style="color:#6e7681;font-size:11px;font-family:monospace;margin-top:8px">`
      + `Tab: complete · Enter: load · Esc: cancel · supports .bam .gff .gff3</div>`
      + `</div>`;
    document.body.appendChild(overlay);
    let _comps = [], _ci = 0;
    const inp = overlay.querySelector('#add-file-input');
    const compsEl = overlay.querySelector('#add-file-comps');
    function renderC() {
      compsEl.innerHTML = _comps.map((p,i)=>`<div data-ci="${i}" style="padding:3px 10px;cursor:pointer;background:${i===_ci?'#1f6feb':'transparent'};color:#e6edf3">${esc(p)}</div>`).join('');
      compsEl.querySelectorAll('[data-ci]').forEach(el2=>el2.addEventListener('mousedown',ev=>{ev.preventDefault();inp.value=_comps[parseInt(el2.dataset.ci)];compsEl.style.display='none';_comps=[];inp.focus();}));
      compsEl.style.display=_comps.length?'block':'none';
    }
    inp.addEventListener('input',()=>{ if(ws&&ws.readyState===1) ws.send(JSON.stringify({cmd:'complete_path',prefix:inp.value})); });
    inp.addEventListener('keydown',ev=>{
      if(ev.key==='Tab'){ev.preventDefault();if(_comps.length){inp.value=_comps[_ci];_comps=[];compsEl.style.display='none';}if(ws&&ws.readyState===1)ws.send(JSON.stringify({cmd:'complete_path',prefix:inp.value}));}
      else if(ev.key==='Enter'){ev.preventDefault();const p=(_comps.length?_comps[_ci]:null)||inp.value.trim();if(p&&ws&&ws.readyState===1)ws.send(JSON.stringify({cmd:'upload_file',path:p}));closeAddFileOverlay();}
      else if(ev.key==='Escape'){ev.preventDefault();closeAddFileOverlay();}
      else if(ev.key==='ArrowDown'&&_comps.length){_ci=Math.min(_ci+1,_comps.length-1);renderC();ev.preventDefault();}
      else if(ev.key==='ArrowUp'&&_comps.length){_ci=Math.max(_ci-1,0);renderC();ev.preventDefault();}
    });
    inp.addEventListener('blur',()=>setTimeout(()=>{compsEl.style.display='none';},150));
    overlay.addEventListener('mousedown',ev=>{ if(ev.target===overlay) closeAddFileOverlay(); });
    window._updateAddFileComps = list => { _comps=list; _ci=0; renderC(); };
  }
  overlay.style.display='flex';
  const inp=overlay.querySelector('#add-file-input');
  inp.value=''; inp.focus();
}
function closeAddFileOverlay(){
  const o=document.getElementById('add-file-overlay');
  if(o) o.style.display='none';
}
window.openAddFileOverlay = openAddFileOverlay;

// ── Panel selector (d key) ────────────────────────────────────────────────────
(function(){
  const PANELS=[
    {id:'live-panel',   label:'Synced browser'},
    {id:'fig-panel',    label:'Custom gene plot'},
    {id:'circ-panel',   label:'Circular map'},
    {id:'struct-panel', label:'Structure viewer'},
    {id:'msa-panel',    label:'MSA viewer'},
  ];
  const panelCbs={};
  const ov=document.createElement('div');
  ov.id='panel-sel-overlay';
  ov.style.cssText='display:none;position:fixed;inset:0;z-index:300;';
  const box=document.createElement('div');
  box.style.cssText='position:fixed;top:50%;left:50%;transform:translate(-50%,-50%);'
    +'background:#161b22;border:1px solid #30363d;border-radius:8px;padding:16px 20px;'
    +'min-width:240px;box-shadow:0 8px 32px rgba(0,0,0,.6);z-index:301;';
  box.innerHTML='<div style="font-size:11px;color:#8b949e;text-transform:uppercase;'
    +'letter-spacing:.06em;margin-bottom:10px">Panels  <span style="font-size:9px;color:#484f58">(d)</span></div>';
  PANELS.forEach(p=>{
    const row=document.createElement('div');
    row.style.cssText='display:flex;align-items:center;gap:8px;margin:6px 0;font-size:12px;color:#c9d1d9;';
    const cb=document.createElement('input'); cb.type='checkbox';
    cb.style.accentColor='#58a6ff';
    const isStruct=(p.id==='struct-panel');
    const isMsa=(p.id==='msa-panel');
    const isCirc=(p.id==='circ-panel');
    cb.checked=!isStruct&&!isMsa; cb.disabled=isStruct;
    cb.addEventListener('change', ()=>{
      const el=document.getElementById(p.id);
      if (!el) return;
      if (p.id==='live-panel') {
        el.style.display=cb.checked?'':'none';
        document.getElementById('drag-handle').style.display=cb.checked?'':'none';
      } else if (p.id==='struct-panel') {
        el.classList.toggle('visible', cb.checked);
        document.getElementById('struct-h-handle').classList.toggle('visible', cb.checked);
      } else if (p.id==='msa-panel') {
        el.classList.toggle('visible', cb.checked);
        document.getElementById('msa-h-handle').classList.toggle('visible', cb.checked);
      } else if (p.id==='circ-panel') {
        el.classList.toggle('visible', cb.checked);
        document.getElementById('circ-h-handle').classList.toggle('visible', cb.checked);
        if (cb.checked && typeof onCircStateUpdate === 'function') onCircStateUpdate();
      } else {
        el.style.display=cb.checked?'':'none';
      }
    });
    panelCbs[p.id]=cb;
    const lbl=document.createElement('label'); lbl.textContent=p.label;
    row.append(cb,lbl); box.appendChild(row);
  });


  ov.appendChild(box); document.body.appendChild(ov);
  window._panelCbs=panelCbs;
  function toggleSel() {
    const open=ov.style.display==='none';
    ov.style.display=open?'block':'none';
  }
  ov.addEventListener('click', e=>{if(e.target===ov) ov.style.display='none';});
  document.addEventListener('keydown', e=>{
    if (e.key==='d'||e.key==='D') {
      const tag=document.activeElement&&document.activeElement.tagName;
      if (tag==='INPUT'||tag==='SELECT'||tag==='TEXTAREA') return;
      toggleSel();
    }
  });
})();
