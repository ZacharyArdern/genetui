// ── Panel focus tracking (set by msa_panel.js, gene_info.js, etc.) ───────────
window.hoveredPanel = 'genome';

// ── Local viewport helpers ────────────────────────────────────────────────────
function renderLocal() {
  if (!lastState) return;
  drawSVG(lastState, localVS, localVE);
  if (typeof window.updateFancyLayer === 'function') window.updateFancyLayer();
}

function getMaxSpan() {
  const el = document.getElementById('opt-max-span');
  const v = el ? parseInt(el.value) : 100000;
  return (v > 0) ? v : Infinity;
}

function clampLocal() {
  const gs = lastState ? (lastState.genome_size||0) : 0;
  // Enforce max zoom-out before clamping to genome bounds
  const maxSpan = getMaxSpan();
  if (localVE - localVS > maxSpan) {
    const center = (localVS + localVE) / 2;
    localVS = Math.round(center - maxSpan / 2);
    localVE = localVS + maxSpan;
  }
  const span = localVE-localVS;
  if (gs>0) {
    if (localVS<0) { localVS=0; localVE=span; }
    if (localVE>gs) { localVE=gs; localVS=Math.max(0,gs-span); }
  }
  localVS = Math.max(0, localVS);
}

let vpSendTimer = null;
function sendViewport() {
  if (ws && ws.readyState===WebSocket.OPEN)
    ws.send(JSON.stringify({start: localVS, end: localVE}));
}
function enterLocalMode() {
  localMode = true;
  clearTimeout(localModeTimer);
  localModeTimer = setTimeout(() => { localMode=false; }, 8000);
  clearTimeout(vpSendTimer);
  vpSendTimer = setTimeout(sendViewport, 200);
}

// ── Wheel: pinch/spread → zoom, two-finger scroll → pan ──────────────────────
// On macOS: pinch fires wheel with ctrlKey=true; two-finger swipe fires with ctrlKey=false.
// On Linux/Windows: Ctrl+scroll also zooms; plain scroll pans.
document.addEventListener('wheel', e => {
  const panel = document.getElementById('live-panel');
  if (!e.composedPath().some(el => el===panel)) return;
  e.preventDefault();
  if (!lastState) return;
  enterLocalMode();
  const span = localVE - localVS;
  const W = panel.clientWidth || window.innerWidth;

  if (e.ctrlKey) {
    // Pinch / Ctrl+scroll → zoom around cursor position
    const rect = panel.getBoundingClientRect();
    const frac = (e.clientX - rect.left) / W;
    const factor = e.deltaY > 0 ? 1.18 : 1 / 1.18;
    const newSpan = Math.max(10, span * factor);
    const center = localVS + frac * span;
    localVS = Math.round(center - frac * newSpan);
    localVE = Math.round(localVS + newSpan);
  } else {
    // Two-finger scroll → pan (use dominant axis)
    const raw = Math.abs(e.deltaX) >= Math.abs(e.deltaY) ? e.deltaX : e.deltaY;
    const shift = Math.round((raw / W) * span);
    localVS += shift;
    localVE += shift;
  }

  clampLocal(); renderLocal(); scheduleAutoRender();
  if (typeof drawViewportMarker === 'function') drawViewportMarker();
}, {passive:false});

// ── Drag pan ──────────────────────────────────────────────────────────────────
document.getElementById('live-panel').addEventListener('mousedown', e => {
  if (e.button!==0) return;
  dragStart=e.clientX; dragVS0=localVS; dragVE0=localVE;
  document.getElementById('live-panel').classList.add('dragging');
});
document.addEventListener('mousemove', e => {
  if (dragStart===null) return;
  enterLocalMode();
  const p=document.getElementById('live-panel');
  const dx=e.clientX-dragStart, span=dragVE0-dragVS0;
  localVS=dragVS0+Math.round(-(dx/(p.clientWidth||1))*span);
  localVE=dragVE0+Math.round(-(dx/(p.clientWidth||1))*span);
  clampLocal(); renderLocal();
  if (typeof drawViewportMarker === 'function') drawViewportMarker();
});
document.addEventListener('mouseup', () => {
  if (dragStart!==null) {
    dragStart=null;
    document.getElementById('live-panel').classList.remove('dragging');
    scheduleAutoRender();
  }
});

// ── Keyboard navigation ───────────────────────────────────────────────────────
document.addEventListener('keydown', e => {
  if (!lastState) return;
  const tag=document.activeElement&&document.activeElement.tagName;
  if (tag==='INPUT'||tag==='SELECT'||tag==='TEXTAREA') return;
  if (window.hoveredPanel && window.hoveredPanel !== 'genome') return;
  if (e.key === '/') {
    e.preventDefault();
    const inp = document.getElementById('browser-search');
    if (inp) { inp.focus(); inp.select(); }
    return;
  }
  const span=localVE-localVS, step=Math.max(1,Math.round(span*.12));
  let changed=true;
  switch(e.key) {
    case 'ArrowLeft':  localVS-=step; localVE-=step; break;
    case 'ArrowRight': localVS+=step; localVE+=step; break;
    case 'ArrowUp':   {const d=Math.round(span*.15); localVS+=d; localVE-=d; break;}
    case 'ArrowDown': {const d=Math.round(span*.15); localVS-=d; localVE+=d; break;}
    case '+': case '=': {const d=Math.round(span*.15); localVS+=d; localVE-=d; break;}
    case '-':           {const d=Math.round(span*.15); localVS-=d; localVE+=d; break;}
    default: changed=false;
  }
  if (!changed) return;
  e.preventDefault();
  enterLocalMode(); clampLocal(); renderLocal(); scheduleAutoRender();
  if (typeof drawViewportMarker === 'function') drawViewportMarker();
});

window.addEventListener('resize', renderLocal);

// ── Browser gene / coordinate search ─────────────────────────────────────────
(function() {
  const inp = document.getElementById('browser-search');
  const btn = document.getElementById('btn-browser-search');
  const res = document.getElementById('browser-search-results');
  let searchHits = [], searchIdx = 0;

  function allSearchFeats() {
    if (!lastState) return [];
    return [...(lastState.features||[]), ...(lastState.blast_features||[])];
  }

  function parseCoord(q) {
    const s = q.replace(/[,_\s]/g,'');
    const m = s.match(/^(\d+)[.\-]+(\d+)$/);
    if (m) return [parseInt(m[1]), parseInt(m[2])];
    const n = parseInt(s);
    if (!isNaN(n) && s.match(/^\d+$/)) return [Math.max(0,n-5000), n+5000];
    return null;
  }

  function navigateTo(f) {
    const span = Math.max(f.end - f.start, 1);
    const pad  = Math.max(span * 2, 5000);
    localVS = Math.max(0, f.start - pad);
    localVE = f.end + pad;
    clampLocal(); renderLocal(); sendViewport();
    enterLocalMode();
  }

  function showDropdown(hits) {
    searchHits = hits;
    searchIdx  = 0;
    if (!hits.length) { res.style.display='none'; return; }
    const rect = inp.getBoundingClientRect();
    res.style.left  = rect.left + 'px';
    res.style.top   = (rect.bottom + 2) + 'px';
    res.style.display = 'block';
    renderDropdown();
  }

  function renderDropdown() {
    res.innerHTML = searchHits.map((f, i) => {
      const bg = i === searchIdx ? '#1f6feb' : 'transparent';
      const locus = f.locus_tag && f.locus_tag !== f.name ? ` <span style="color:#6e7681">${esc(f.locus_tag)}</span>` : '';
      return `<div data-i="${i}" style="padding:3px 10px;cursor:pointer;background:${bg};color:#c9d1d9">
        <span style="color:#58a6ff">${esc(f.name)}</span>${locus}
        <span style="color:#484f58;float:right">${fmt(f.start)}</span></div>`;
    }).join('');
    res.querySelectorAll('[data-i]').forEach(el => {
      el.addEventListener('mousedown', e => {
        e.preventDefault();
        const i = parseInt(el.dataset.i);
        navigateTo(searchHits[i]);
        res.style.display='none';
        inp.value='';
      });
    });
  }

  function doSearch(commit) {
    const q = inp.value.trim();
    if (!q) { res.style.display='none'; return; }
    const coord = parseCoord(q);
    if (coord) {
      if (commit) {
        localVS=coord[0]; localVE=coord[1];
        clampLocal(); renderLocal(); sendViewport(); enterLocalMode();
        res.style.display='none'; inp.value=''; inp.blur();
      }
      return;
    }
    const ql = q.toLowerCase();
    const feats = allSearchFeats();
    const exact  = feats.filter(f => f.name.toLowerCase()===ql || (f.locus_tag||'').toLowerCase()===ql);
    const prefix = feats.filter(f => !exact.includes(f) && (f.name.toLowerCase().startsWith(ql) || (f.locus_tag||'').toLowerCase().startsWith(ql)));
    const hits = [...exact, ...prefix].slice(0, 20);
    if (commit && hits.length > 0) {
      navigateTo(hits[0]);
      res.style.display='none'; inp.value=''; inp.blur();
    } else {
      showDropdown(hits);
    }
  }

  inp.addEventListener('input', () => doSearch(false));
  inp.addEventListener('keydown', e => {
    if (e.key === 'Enter') { e.preventDefault(); doSearch(true); }
    else if (e.key === 'Escape') { res.style.display='none'; inp.blur(); }
    else if (e.key === 'ArrowDown' && searchHits.length) {
      searchIdx = Math.min(searchIdx+1, searchHits.length-1); renderDropdown(); e.preventDefault();
    } else if (e.key === 'ArrowUp' && searchHits.length) {
      searchIdx = Math.max(searchIdx-1, 0); renderDropdown(); e.preventDefault();
    }
  });
  inp.addEventListener('focus', () => { window.hoveredPanel='search'; });
  inp.addEventListener('blur',  () => { setTimeout(()=>{ res.style.display='none'; window.hoveredPanel='genome'; },150); });
  btn.addEventListener('click', () => doSearch(true));
  document.getElementById('browser-search-form').addEventListener('submit', () => doSearch(true));
})();

// ── Vertical drag-resize (live-panel height) ──────────────────────────────────
(function(){
  const h=document.getElementById('drag-handle');
  const p=document.getElementById('live-panel');
  let drag=false, sy=0, sh=0;
  h.addEventListener('mousedown', e=>{drag=true;sy=e.clientY;sh=p.offsetHeight;e.preventDefault();});
  document.addEventListener('mousemove', e=>{
    if (!drag) return;
    p.style.height=Math.max(60,Math.min(sh+e.clientY-sy,window.innerHeight-160))+'px';
    renderLocal();
  });
  document.addEventListener('mouseup', ()=>{drag=false;});
})();
