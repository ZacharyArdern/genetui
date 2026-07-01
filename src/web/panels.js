// ── Options dropdown (anchored to btn-opts, draggable) ───────────────────────
(function(){
  const btn=document.getElementById('btn-opts');
  const box=document.getElementById('opts-box');
  let open=false, dragging=false, ox=0, oy=0, mx=0, my=0;

  function positionBox() {
    const r=btn.getBoundingClientRect();
    let left=r.right-box.offsetWidth;
    if (left<4) left=4;
    const top=r.bottom+4;
    box.style.left=left+'px'; box.style.top=top+'px';
    box.style.right='auto'; box.style.bottom='auto';
  }

  btn.addEventListener('click', function(e) {
    e.stopPropagation();
    open=!open;
    box.style.display=open?'block':'none';
    btn.classList.toggle('active',open);
    if (open) positionBox();
  });

  document.addEventListener('click', e=>{
    if (open && !box.contains(e.target) && e.target!==btn) {
      open=false; box.style.display='none'; btn.classList.remove('active');
    }
  });

  box.querySelector('h3').style.cursor='move';
  box.querySelector('h3').addEventListener('mousedown', e=>{
    dragging=true; mx=e.clientX; my=e.clientY;
    const r=box.getBoundingClientRect(); ox=r.left; oy=r.top;
    e.preventDefault();
  });
  document.addEventListener('mousemove', e=>{
    if (!dragging) return;
    box.style.left=(ox+e.clientX-mx)+'px';
    box.style.top=(oy+e.clientY-my)+'px';
  });
  document.addEventListener('mouseup', ()=>{dragging=false;});
})();

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
document.getElementById('cov-style-svg').addEventListener('change', function(){
  const hr = document.getElementById('cov-height-row');
  if (hr) hr.style.display = this.value !== 'none' ? '' : 'none';
  renderLocal();
});
document.getElementById('cov-height').addEventListener('input', function(){
  const lbl = document.getElementById('lbl-cov-height');
  if (lbl) lbl.textContent = this.value;
  renderLocal();
});
document.getElementById('opt-codons').addEventListener('input', function(){
  renderLocal();
});

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
