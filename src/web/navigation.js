// ── Panel focus tracking ──────────────────────────────────────────────────────
window.hoveredPanel = 'genome';
window.focusedPanel = 'genome';

(function() {
  const panelIds = ['live-panel','circ-panel','fig-panel','struct-panel','msa-panel','circ-body'];
  function setFocus(id) {
    window.focusedPanel = id;
    ['live-panel','circ-panel','fig-panel','struct-panel','msa-panel'].forEach(pid => {
      const el = document.getElementById(pid);
      if (el) el.classList.toggle('panel-focused', pid === id);
    });
  }
  panelIds.forEach(id => {
    const el = document.getElementById(id);
    if (!el) return;
    el.addEventListener('mousedown', () => {
      const canonical = id === 'circ-body' ? 'circ-panel' : id;
      setFocus(canonical);
    });
  });
  // Default focus on genome track
  setFocus('live-panel');
})();

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
    openSearchOverlay();
    return;
  }
  if (e.key === 'q' || e.key === 'Q') {
    e.preventDefault();
    window.close();
    return;
  }
  if (e.key === 'u' || e.key === 'U') {
    e.preventDefault();
    openUploadBamOverlay();
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

// ── Upload BAM overlay (u key) ────────────────────────────────────────────────
let _uploadBamComps = [], _uploadBamCompIdx = 0;

function openUploadBamOverlay() {
  let overlay = document.getElementById('upload-bam-overlay');
  if (!overlay) {
    overlay = document.createElement('div');
    overlay.id = 'upload-bam-overlay';
    overlay.style.cssText = 'display:none;position:fixed;inset:0;z-index:600;background:rgba(0,0,0,0.55);'
      + 'align-items:flex-start;justify-content:center';
    overlay.innerHTML =
      `<div style="margin-top:120px;width:520px;background:#161b22;border:1px solid #3fb950;`
      + `border-radius:8px;overflow:visible;box-shadow:0 8px 32px rgba(0,0,0,0.7);padding:16px">`
      + `<div style="color:#3fb950;font-size:13px;font-family:monospace;margin-bottom:10px">Upload BAM / SAM / CRAM</div>`
      + `<input id="upload-bam-input" type="text" placeholder="path/to/file.bam" autocomplete="off"`
      + ` style="width:100%;box-sizing:border-box;background:#0d1117;color:#e6edf3;border:1px solid #3fb950;`
      + `border-radius:4px;padding:6px 10px;font-size:13px;font-family:monospace">`
      + `<div id="upload-bam-comps" style="position:absolute;left:16px;right:16px;background:#1c2128;`
      + `border:1px solid #3fb950;border-radius:4px;max-height:200px;overflow-y:auto;`
      + `font-size:12px;font-family:monospace;display:none;z-index:10"></div>`
      + `<div style="color:#6e7681;font-size:11px;font-family:monospace;margin-top:8px">`
      + `Tab: complete · Enter: load · Esc: cancel</div>`
      + `</div>`;
    document.body.appendChild(overlay);

    const inp = overlay.querySelector('#upload-bam-input');
    const compsEl = overlay.querySelector('#upload-bam-comps');

    function renderComps() {
      compsEl.innerHTML = _uploadBamComps.map((p, i) => {
        const bg = i === _uploadBamCompIdx ? '#1f6feb' : 'transparent';
        return `<div data-ci="${i}" style="padding:3px 10px;cursor:pointer;background:${bg};color:#e6edf3">${esc(p)}</div>`;
      }).join('');
      compsEl.querySelectorAll('[data-ci]').forEach(el => {
        el.addEventListener('mousedown', ev => {
          ev.preventDefault();
          inp.value = _uploadBamComps[parseInt(el.dataset.ci)];
          compsEl.style.display = 'none';
          _uploadBamComps = [];
          inp.focus();
        });
      });
      compsEl.style.display = _uploadBamComps.length ? 'block' : 'none';
    }

    inp.addEventListener('input', () => {
      if (ws && ws.readyState === WebSocket.OPEN)
        ws.send(JSON.stringify({cmd: 'complete_path', prefix: inp.value}));
    });

    inp.addEventListener('keydown', ev => {
      if (ev.key === 'Tab') {
        ev.preventDefault();
        if (_uploadBamComps.length) {
          inp.value = _uploadBamComps[_uploadBamCompIdx];
          _uploadBamComps = []; compsEl.style.display = 'none';
        }
        if (ws && ws.readyState === WebSocket.OPEN)
          ws.send(JSON.stringify({cmd: 'complete_path', prefix: inp.value}));
      } else if (ev.key === 'Enter') {
        ev.preventDefault();
        const path = (_uploadBamComps.length ? _uploadBamComps[_uploadBamCompIdx] : null) || inp.value.trim();
        if (path && ws && ws.readyState === WebSocket.OPEN)
          ws.send(JSON.stringify({cmd: 'upload_bam', path}));
        closeUploadBamOverlay();
      } else if (ev.key === 'Escape') {
        ev.preventDefault();
        closeUploadBamOverlay();
      } else if (ev.key === 'ArrowDown' && _uploadBamComps.length) {
        _uploadBamCompIdx = Math.min(_uploadBamCompIdx + 1, _uploadBamComps.length - 1);
        renderComps(); ev.preventDefault();
      } else if (ev.key === 'ArrowUp' && _uploadBamComps.length) {
        _uploadBamCompIdx = Math.max(_uploadBamCompIdx - 1, 0);
        renderComps(); ev.preventDefault();
      }
    });

    inp.addEventListener('blur', () => {
      setTimeout(() => { compsEl.style.display = 'none'; }, 150);
    });

    overlay.addEventListener('mousedown', ev => {
      if (ev.target === overlay) closeUploadBamOverlay();
    });

    window._updateUploadBamComps = function(list) {
      _uploadBamComps = list; _uploadBamCompIdx = 0; renderComps();
    };
  }
  _uploadBamComps = []; _uploadBamCompIdx = 0;
  overlay.style.display = 'flex';
  const inp = overlay.querySelector('#upload-bam-input');
  inp.value = '';
  inp.focus();
}

function closeUploadBamOverlay() {
  const overlay = document.getElementById('upload-bam-overlay');
  if (overlay) overlay.style.display = 'none';
  _uploadBamComps = [];
}

// ── Search overlay (gene/position + DIAMOND blast) ────────────────────────────
// Tab switching
window.switchSearchTab = function(tab) {
  const isGene = tab === 'gene';
  const tG = document.getElementById('stab-gene');
  const tD = document.getElementById('stab-dmnd');
  tG.style.background = isGene ? '#1f6feb' : 'transparent';
  tG.style.color      = isGene ? '#fff'    : '#8b949e';
  tD.style.background = isGene ? 'transparent' : '#1f6feb';
  tD.style.color      = isGene ? '#8b949e'     : '#fff';
  document.getElementById('search-pane-gene').style.display = isGene ? '' : 'none';
  document.getElementById('search-pane-dmnd').style.display = isGene ? 'none' : '';
  if (isGene) {
    const i = document.getElementById('browser-search'); i.focus(); i.select();
  } else {
    document.getElementById('dmnd-query-input').focus();
  }
};

(function() {
  const overlay = document.getElementById('search-overlay');
  const inp     = document.getElementById('browser-search');
  const res     = document.getElementById('browser-search-results');
  let searchHits = [], searchIdx = 0;

  function closeSearch() {
    overlay.style.display='none';
    res.style.display='none';
    inp.value='';
    searchHits=[];
    const qR = document.getElementById('dmnd-query-completions');
    if (qR) qR.style.display='none';
    window.hoveredPanel='genome';
  }

  window.openSearchOverlay = function() {
    overlay.style.display='flex';
    switchSearchTab('gene');
    window.hoveredPanel='search';
  };

  // Click outside the inner box closes overlay
  overlay.addEventListener('mousedown', e => {
    if (e.target===overlay) closeSearch();
  });

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

  function renderDropdown() {
    res.innerHTML = searchHits.map((f, i) => {
      const bg = i === searchIdx ? '#1f6feb' : 'transparent';
      const locus = f.locus_tag && f.locus_tag !== f.name ? ` <span style="color:#6e7681">${esc(f.locus_tag)}</span>` : '';
      return `<div data-i="${i}" style="padding:4px 10px;cursor:pointer;background:${bg};color:#c9d1d9">
        <span style="color:#58a6ff">${esc(f.name)}</span>${locus}
        <span style="color:#484f58;float:right">${fmt(f.start)}</span></div>`;
    }).join('');
    res.querySelectorAll('[data-i]').forEach(el => {
      el.addEventListener('mousedown', e => {
        e.preventDefault();
        navigateTo(searchHits[parseInt(el.dataset.i)]);
        closeSearch();
      });
    });
  }

  function doSearch(commit) {
    const q = inp.value.trim();
    if (!q) { res.style.display='none'; searchHits=[]; return; }
    const coord = parseCoord(q);
    if (coord) {
      if (commit) {
        localVS=coord[0]; localVE=coord[1];
        clampLocal(); renderLocal(); sendViewport(); enterLocalMode();
        closeSearch();
      }
      return;
    }
    const ql = q.toLowerCase();
    const feats = allSearchFeats();
    const exact  = feats.filter(f => f.name.toLowerCase()===ql || (f.locus_tag||'').toLowerCase()===ql);
    const prefix = feats.filter(f => !exact.includes(f) && (f.name.toLowerCase().startsWith(ql) || (f.locus_tag||'').toLowerCase().startsWith(ql)));
    searchHits = [...exact, ...prefix].slice(0, 20);
    searchIdx  = 0;
    if (commit && searchHits.length > 0) {
      navigateTo(searchHits[0]);
      closeSearch();
    } else if (searchHits.length) {
      res.style.display='block';
      renderDropdown();
    } else {
      res.style.display='none';
    }
  }

  inp.addEventListener('input', () => doSearch(false));
  inp.addEventListener('keydown', e => {
    if (e.key === 'Enter') { e.preventDefault(); doSearch(true); }
    else if (e.key === 'Escape') { closeSearch(); }
    else if (e.key === 'ArrowDown' && searchHits.length) {
      searchIdx = Math.min(searchIdx+1, searchHits.length-1); renderDropdown(); e.preventDefault();
    } else if (e.key === 'ArrowUp' && searchHits.length) {
      searchIdx = Math.max(searchIdx-1, 0); renderDropdown(); e.preventDefault();
    }
  });

  // ── DIAMOND blast pane ──────────────────────────────────────────────────────
  const qInp = document.getElementById('dmnd-query-input');
  const qRes = document.getElementById('dmnd-query-completions');
  let qComps = [], qCompIdx = 0;

  function renderQComps() {
    qRes.innerHTML = qComps.map((p,i) => {
      const bg = i===qCompIdx ? '#1f6feb' : 'transparent';
      return `<div data-qi="${i}" style="padding:3px 10px;cursor:pointer;background:${bg};color:#c9d1d9;font-family:monospace;font-size:11px">${esc(p)}</div>`;
    }).join('');
    qRes.querySelectorAll('[data-qi]').forEach(el => {
      el.addEventListener('mousedown', ev => {
        ev.preventDefault();
        qInp.value = qComps[parseInt(el.dataset.qi)];
        qRes.style.display='none'; qComps=[];
      });
    });
  }

  // Called by onStatusUpdate in gene_info.js when path_completions arrive
  window._updateDmndQueryComps = function(list) {
    if (document.getElementById('search-pane-dmnd').style.display === 'none') return;
    qComps = list; qCompIdx = 0;
    if (!list.length) { qRes.style.display='none'; return; }
    qRes.style.display='block';
    renderQComps();
  };

  function runDmnd() {
    const query = (qComps.length ? qComps[qCompIdx] : null) || qInp.value.trim();
    if (!query) return;
    const use6ft = document.getElementById('dmnd-t0').checked;
    if (ws && ws.readyState === WebSocket.OPEN)
      ws.send(JSON.stringify({cmd: 'run_diamond', query, use_6ft: use6ft}));
    closeSearch();
  }

  qInp.addEventListener('input', () => {
    if (ws && ws.readyState === WebSocket.OPEN)
      ws.send(JSON.stringify({cmd: 'complete_path', prefix: qInp.value}));
  });
  qInp.addEventListener('keydown', e => {
    if (e.key === 'Tab') {
      e.preventDefault();
      if (qComps.length) { qInp.value = qComps[qCompIdx]; qRes.style.display='none'; qComps=[]; }
      if (ws && ws.readyState === WebSocket.OPEN)
        ws.send(JSON.stringify({cmd: 'complete_path', prefix: qInp.value}));
    } else if (e.key === 'Enter') { e.preventDefault(); runDmnd(); }
    else if (e.key === 'Escape') { closeSearch(); }
    else if (e.key === 'ArrowDown' && qComps.length) {
      qCompIdx = Math.min(qCompIdx+1, qComps.length-1); renderQComps(); e.preventDefault();
    } else if (e.key === 'ArrowUp' && qComps.length) {
      qCompIdx = Math.max(qCompIdx-1, 0); renderQComps(); e.preventDefault();
    }
  });
  document.getElementById('btn-run-dmnd').addEventListener('click', runDmnd);
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
