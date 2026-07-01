// ── Horizontal drag resize ────────────────────────────────────────────────────
(function(){
  const h=document.getElementById('circ-h-handle');
  const cp=document.getElementById('circ-panel');
  let drag=false, sx=0, sw=0;
  h.addEventListener('mousedown',e=>{drag=true;sx=e.clientX;sw=cp.offsetWidth;e.preventDefault();});
  document.addEventListener('mousemove',e=>{
    if(!drag)return;
    const row=document.getElementById('bottom-row');
    cp.style.width=Math.max(180,Math.min(sw+e.clientX-sx,row.clientWidth-100))+'px';
  });
  document.addEventListener('mouseup',()=>{drag=false;});
})();

// ── State ─────────────────────────────────────────────────────────────────────
let circLastFetchKey  = null;
let circDebounceTimer = null;
let circDragState     = null;  // active drag of "you are here" arc
let circDragMoved     = false; // suppress click after drag
let circZoom          = 0.65;  // fraction of panel width; wheel zooms this

function applyCircZoom() {
  document.querySelectorAll('.circ-map-item svg').forEach(svg => {
    svg.style.width  = Math.round(circZoom * 100) + '%';
    svg.style.height = 'auto';
  });
}

// ── Wheel zoom on circle panel ────────────────────────────────────────────────
document.getElementById('circ-body').addEventListener('wheel', e => {
  e.preventDefault();
  const factor = e.deltaY > 0 ? 0.92 : 1 / 0.92;
  circZoom = Math.max(0.25, Math.min(1.5, circZoom * factor));
  applyCircZoom();
}, { passive: false });

function setCircStatus(msg){const el=document.getElementById('circ-status');if(el)el.textContent=msg;}
function circOpt(id){return document.getElementById(id);}

function circFetchKey(state) {
  if (!state) return null;
  return JSON.stringify([state.genome_size, state.main_size, (state.plasmid_sizes||[]), (state.blast_features||[]).length]);
}

// ── Fetch and render all genome maps ─────────────────────────────────────────
async function fetchAllCircMaps() {
  if (!lastState || lastState.genome_size === 0) return;
  const dark   = !circOpt('circ-white').checked ? 1 : 0;
  const nc     = circOpt('circ-show-nc').checked ? 1 : 0;
  const legend = circOpt('circ-show-legend').checked ? 1 : 0;
  const titleInp = document.getElementById('circ-title-input');
  const titleVal = titleInp && titleInp.value.trim();
  const titleParam = titleVal ? `&title=${encodeURIComponent(titleVal)}` : '';

  // Build list of genome indices to render
  const nPlasmids = (lastState.plasmid_names && circOpt('circ-show-plasmids').checked)
    ? lastState.plasmid_names.length : 0;
  const genomeCount = 1 + nPlasmids;

  setCircStatus('rendering\u2026');

  // Fetch all SVGs in parallel
  let svgTexts;
  try {
    svgTexts = await Promise.all(
      Array.from({length: genomeCount}, (_, i) =>
        fetch(`/circular-map.svg?genome=${i}&dark=${dark}&nc=${nc}&legend=${legend}${titleParam}`)
          .then(r => r.ok ? r.text() : Promise.reject(r.status))
      )
    );
  } catch(err) {
    setCircStatus('error: ' + err);
    return;
  }

  const body = document.getElementById('circ-body');
  body.innerHTML = '';

  svgTexts.forEach((svgText, genomeIdx) => {
    const wrapper = document.createElement('div');
    wrapper.className = 'circ-map-item';
    wrapper.dataset.genome = genomeIdx;
    if (genomeIdx > 0) {
      // Plasmid maps rendered at 60% width, centred
      wrapper.style.cssText = 'width:60%;margin:0 auto;';
    }
    wrapper.innerHTML = svgText;
    const svgEl = wrapper.querySelector('svg');
    if (svgEl) {
      svgEl.removeAttribute('width');
      svgEl.removeAttribute('height');
      svgEl.style.cssText = `width:${Math.round(circZoom*100)}%;height:auto;display:block;margin:0 auto;`;

      // Title double-click → focus options
      const titleEl = svgEl.getElementById('circ-title-svg');
      if (titleEl) {
        titleEl.addEventListener('dblclick', () => {
          const box = document.getElementById('circ-opts-box');
          const inp = document.getElementById('circ-title-input');
          if (!box||!inp) return;
          box.style.display='block';
          const btn=document.getElementById('btn-circ-opts');
          if(btn){const r=btn.getBoundingClientRect();box.style.right=(document.body.clientWidth-r.right)+'px';box.style.top=(r.bottom+4)+'px';}
          setTimeout(()=>inp.focus(),30);
        });
      }

      // Click to navigate
      svgEl.style.cursor = 'crosshair';
      svgEl.addEventListener('click', e => {
        if (circDragMoved) { circDragMoved = false; return; }
        if (!lastState || lastState.genome_size === 0) return;
        const cx = parseFloat(svgEl.getAttribute('data-cx'));
        const cy = parseFloat(svgEl.getAttribute('data-cy'));
        const rOuter = parseFloat(svgEl.getAttribute('data-r-outer'));
        const rInner = parseFloat(svgEl.getAttribute('data-r-inner'));
        if (isNaN(cx)||isNaN(cy)||isNaN(rOuter)||isNaN(rInner)) return;

        const pt = svgEl.createSVGPoint();
        pt.x = e.clientX; pt.y = e.clientY;
        const svgPt = pt.matrixTransform(svgEl.getScreenCTM().inverse());
        const dx = svgPt.x - cx, dy = svgPt.y - cy;
        const dist = Math.sqrt(dx*dx + dy*dy);
        if (dist < rInner * 0.85 || dist > rOuter * 1.08) return;

        let angle = Math.atan2(dy, dx) + Math.PI/2;
        if (angle < 0) angle += 2*Math.PI;

        // Genome size for this map
        const gsize = genomeIdx === 0
          ? (lastState.main_size || lastState.genome_size)
          : (lastState.plasmid_sizes && lastState.plasmid_sizes[genomeIdx-1]) || 0;
        if (gsize === 0) return;

        const clickPos = Math.round(angle / (2*Math.PI) * gsize);
        const viewSpan = lastState.view_end - lastState.view_start;
        const half = Math.max(5000, Math.round(viewSpan / 2));
        const vs = Math.max(0, clickPos - half);
        const ve = Math.min(gsize, vs + 2*half);

        if (genomeIdx !== (lastState.active_genome || 0)) {
          // Switch genome and navigate
          if (ws && ws.readyState === WebSocket.OPEN) {
            ws.send(JSON.stringify({cmd:'switch_genome', genome:genomeIdx, start:vs, end:ve}));
          }
        } else {
          // Same genome — just navigate
          localVS = vs; localVE = ve;
          sendViewport();
        }
      });
    }
    body.appendChild(wrapper);
  });

  circLastFetchKey = circFetchKey(lastState);
  const activeName = lastState.genome_name;
  const activeMb   = (lastState.genome_size/1e6).toFixed(2);
  setCircStatus(`${activeName} — ${activeMb} Mb${nPlasmids > 0 ? ` + ${nPlasmids} plasmid(s)` : ''}`);
  drawViewportMarker();
}

// ── "you are here" marker — drawn into each SVG's viewport layer ─────────────
function drawViewportMarker() {
  if (!lastState || !circOpt('circ-show-viewport').checked) {
    document.querySelectorAll('.circ-map-item svg #circ-viewport-layer').forEach(el => el.remove());
    return;
  }

  // Use localVS/localVE for immediate responsiveness (globals from navigation.js)
  const vs = (typeof localVS !== 'undefined') ? localVS : lastState.view_start;
  const ve = (typeof localVE !== 'undefined') ? localVE : lastState.view_end;

  const ns = 'http://www.w3.org/2000/svg';
  document.querySelectorAll('.circ-map-item').forEach(wrapper => {
    const genomeIdx = parseInt(wrapper.dataset.genome, 10);
    const svgEl = wrapper.querySelector('svg');
    if (!svgEl) return;

    const old = svgEl.getElementById('circ-viewport-layer');
    if (old) old.remove();

    const activeIdx = lastState.active_genome || 0;
    if (genomeIdx !== activeIdx) return;

    const cx = parseFloat(svgEl.getAttribute('data-cx'));
    const cy = parseFloat(svgEl.getAttribute('data-cy'));
    const rOuter = parseFloat(svgEl.getAttribute('data-r-outer'));
    const rInner = parseFloat(svgEl.getAttribute('data-r-inner'));
    if (isNaN(cx)||isNaN(cy)||isNaN(rOuter)||isNaN(rInner)) return;

    const gsize = lastState.genome_size;
    if (!gsize) return;

    const layer = document.createElementNS(ns, 'g');
    layer.id = 'circ-viewport-layer';
    svgEl.appendChild(layer);

    function posToAngle(p) { return -Math.PI/2 + 2*Math.PI * p / gsize; }
    const a1   = posToAngle(vs);
    const a2   = posToAngle(ve);
    const span = 2*Math.PI*(ve - vs)/gsize;
    const laf  = span > Math.PI ? 1 : 0;
    const fmt  = v => v.toFixed(2);

    const x1o=cx+rOuter*Math.cos(a1), y1o=cy+rOuter*Math.sin(a1);
    const x2o=cx+rOuter*Math.cos(a2), y2o=cy+rOuter*Math.sin(a2);
    const x1i=cx+rInner*Math.cos(a1), y1i=cy+rInner*Math.sin(a1);
    const x2i=cx+rInner*Math.cos(a2), y2i=cy+rInner*Math.sin(a2);
    const d=`M ${fmt(x1o)} ${fmt(y1o)} A ${fmt(rOuter)} ${fmt(rOuter)} 0 ${laf} 1 ${fmt(x2o)} ${fmt(y2o)} L ${fmt(x2i)} ${fmt(y2i)} A ${fmt(rInner)} ${fmt(rInner)} 0 ${laf} 0 ${fmt(x1i)} ${fmt(y1i)} Z`;

    const fill=document.createElementNS(ns,'path');
    fill.setAttribute('d',d);
    fill.setAttribute('fill','rgba(88,166,255,0.22)');
    fill.setAttribute('stroke','#58a6ff');
    fill.setAttribute('stroke-width','1.5');
    fill.setAttribute('style','cursor:grab');
    fill.setAttribute('pointer-events','all');
    layer.appendChild(fill);

    // Drag on arc: mousedown starts drag
    fill.addEventListener('mousedown', e => {
      if (e.button !== 0) return;
      e.stopPropagation();
      const pt = svgEl.createSVGPoint();
      pt.x = e.clientX; pt.y = e.clientY;
      const sp = pt.matrixTransform(svgEl.getScreenCTM().inverse());
      const startAngle = Math.atan2(sp.y - cy, sp.x - cx);
      circDragState = { svgEl, cx, cy, gsize, rInner, startAngle, startVS: vs, startVE: ve };
      circDragMoved = false;
      e.preventDefault();
    });

    // "you are here" label
    const midA=(a1+a2)/2;
    const rLabel = rOuter + (cx * 0.055);
    const lx=cx+rLabel*Math.cos(midA), ly=cy+rLabel*Math.sin(midA);
    const rot=((midA+Math.PI/2)*180/Math.PI)%360;
    const lbl=document.createElementNS(ns,'text');
    lbl.setAttribute('x',lx.toFixed(1));lbl.setAttribute('y',ly.toFixed(1));
    lbl.setAttribute('text-anchor','middle');lbl.setAttribute('dominant-baseline','middle');
    lbl.setAttribute('font-size','8.5');lbl.setAttribute('font-family','sans-serif');
    lbl.setAttribute('fill','#58a6ff');lbl.setAttribute('font-weight','bold');
    lbl.setAttribute('transform',`rotate(${rot.toFixed(1)},${lx.toFixed(1)},${ly.toFixed(1)})`);
    lbl.setAttribute('pointer-events','none');
    lbl.textContent='you are here';
    layer.appendChild(lbl);
  });
}

// ── Drag "you are here" arc ───────────────────────────────────────────────────
document.addEventListener('mousemove', e => {
  if (!circDragState) return;
  const { svgEl, cx, cy, gsize, rInner, startAngle, startVS, startVE } = circDragState;
  const pt = svgEl.createSVGPoint();
  pt.x = e.clientX; pt.y = e.clientY;
  const sp = pt.matrixTransform(svgEl.getScreenCTM().inverse());
  // Project cursor onto ring if inside inner radius to prevent degenerate arcs
  const dx = sp.x - cx, dy = sp.y - cy;
  const dist = Math.sqrt(dx*dx + dy*dy);
  const safeR = Math.max(dist, rInner * 1.1);
  const px = cx + dx * safeR / dist;
  const py = cy + dy * safeR / dist;
  let delta = Math.atan2(py - cy, px - cx) - startAngle;
  // Normalise to -π..π
  if (delta > Math.PI) delta -= 2*Math.PI;
  if (delta < -Math.PI) delta += 2*Math.PI;
  const posDelta = Math.round(delta / (2*Math.PI) * gsize);
  localVS = Math.max(0, startVS + posDelta);
  localVE = Math.min(gsize, startVE + posDelta);
  circDragMoved = true;
  enterLocalMode();   // keeps localMode=true so server updates don't override
  renderLocal();      // update genome browser track immediately
  drawViewportMarker();
  if (typeof scheduleAutoRender === 'function') scheduleAutoRender(); // trigger gene plot
});

document.addEventListener('mouseup', e => {
  if (!circDragState) return;
  circDragState = null;
  if (circDragMoved) { sendViewport(); if (typeof scheduleAutoRender === 'function') scheduleAutoRender(); }
});

// ── PNG download — downloads all maps as one tall PNG ────────────────────────
async function downloadCircPng() {
  const svgEls = [...document.querySelectorAll('.circ-map-item svg')];
  if (!svgEls.length) return;
  setCircStatus('preparing PNG\u2026');
  try {
    const sz = 2000;
    const c = document.createElement('canvas');
    c.width = sz; c.height = sz * svgEls.length;
    const ctx = c.getContext('2d');
    for (let i = 0; i < svgEls.length; i++) {
      await new Promise((res, rej) => {
        const svgStr = new XMLSerializer().serializeToString(svgEls[i]);
        const blob = new Blob([svgStr],{type:'image/svg+xml'});
        const url = URL.createObjectURL(blob);
        const img = new Image();
        img.onload = () => { ctx.drawImage(img,0,i*sz,sz,sz); URL.revokeObjectURL(url); res(); };
        img.onerror = rej;
        img.src = url;
      });
    }
    const a = document.createElement('a');
    a.href = c.toDataURL('image/png');
    a.download = `${(lastState&&lastState.genome_name||'genome').replace(/\s+/g,'_')}_circular.png`;
    a.click();
    setCircStatus('PNG downloaded');
  } catch(e) { setCircStatus('PNG error'); }
}

// ── Options dropdown ──────────────────────────────────────────────────────────
(function(){
  const btn=document.getElementById('btn-circ-opts');
  const box=document.getElementById('circ-opts-box');
  let open=false;
  btn.addEventListener('click',e=>{
    e.stopPropagation();open=!open;
    box.style.display=open?'block':'none';
    btn.classList.toggle('active',open);
    if(open){const r=btn.getBoundingClientRect();box.style.right=(document.body.clientWidth-r.right)+'px';box.style.top=(r.bottom+4)+'px';}
  });
  document.addEventListener('click',e=>{
    if(open&&!box.contains(e.target)&&e.target!==btn){open=false;box.style.display='none';btn.classList.remove('active');}
  });
})();

['circ-white','circ-show-nc','circ-show-legend','circ-show-plasmids'].forEach(id=>{
  document.getElementById(id).addEventListener('change',()=>fetchAllCircMaps());
});
document.getElementById('circ-show-viewport').addEventListener('change',()=>drawViewportMarker());
(function(){
  const inp=document.getElementById('circ-title-input');
  if(!inp)return;
  inp.addEventListener('keydown',e=>{if(e.key==='Enter')fetchAllCircMaps();});
  inp.addEventListener('blur',()=>fetchAllCircMaps());
})();
document.getElementById('btn-circ-render').addEventListener('click',()=>fetchAllCircMaps());
document.getElementById('btn-circ-png').addEventListener('click',()=>downloadCircPng());
document.getElementById('btn-circ-zoom-in').addEventListener('click',()=>{
  circZoom=Math.min(1.5,circZoom*1.15); applyCircZoom();
});
document.getElementById('btn-circ-zoom-out').addEventListener('click',()=>{
  circZoom=Math.max(0.25,circZoom/1.15); applyCircZoom();
});

// ── Hooks called from websocket.js / panels.js ────────────────────────────────
function onCircStateUpdate() {
  const panel=document.getElementById('circ-panel');
  if(!panel.classList.contains('visible'))return;
  if(!lastState)return;
  const key = circFetchKey(lastState);
  if (key !== circLastFetchKey) {
    clearTimeout(circDebounceTimer);
    circDebounceTimer = setTimeout(fetchAllCircMaps, 300);
  } else {
    drawViewportMarker();
  }
}
