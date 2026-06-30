// ── Panel focus tracking (set by msa_panel.js, gene_info.js, etc.) ───────────
window.hoveredPanel = 'genome';

// ── Local viewport helpers ────────────────────────────────────────────────────
function renderLocal() { if (lastState) drawSVG(lastState, localVS, localVE); }

function clampLocal() {
  const gs = lastState ? (lastState.genome_size||0) : 0;
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
