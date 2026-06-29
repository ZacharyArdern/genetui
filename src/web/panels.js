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
['opt-sixframe','opt-stopcodons','fig-ruler','fig-white','fig-show-labels','fig-hide-long-labels'].forEach(id=>{
  document.getElementById(id).addEventListener('change', optionChanged);
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
    cb.checked=!isStruct&&!isMsa; cb.disabled=isStruct||isMsa;
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

  // Coverage style row (live SVG track only)
  const covSep = document.createElement('div');
  covSep.style.cssText='border-top:1px solid #21262d;margin:10px 0 8px;';
  box.appendChild(covSep);
  const covRow = document.createElement('div');
  covRow.style.cssText='display:flex;align-items:center;gap:8px;font-size:12px;color:#c9d1d9;';
  const covLbl = document.createElement('span'); covLbl.textContent='Coverage (live)';
  covLbl.style.cssText='min-width:80px;color:#8b949e;';
  const covSel = document.createElement('select');
  covSel.id = 'cov-style-svg';
  covSel.style.cssText='background:#0d1117;color:#c9d1d9;border:1px solid #30363d;border-radius:4px;padding:2px 4px;font-size:11px;';
  [['none','None'],['histogram','Histogram'],['kernel','Kernel'],['reads','Raw reads']].forEach(([v,l])=>{
    const o=document.createElement('option'); o.value=v; o.textContent=l; covSel.appendChild(o);
  });
  const covHRow = document.createElement('div');
  covHRow.id='cov-height-row';
  covHRow.style.cssText='display:none;align-items:center;gap:6px;font-size:12px;color:#c9d1d9;margin-top:6px;';
  const covHLbl = document.createElement('span'); covHLbl.textContent='Height';
  covHLbl.style.cssText='min-width:80px;color:#8b949e;';
  const covHIn = document.createElement('input');
  covHIn.type='range'; covHIn.id='cov-height'; covHIn.min=30; covHIn.max=150; covHIn.step=5; covHIn.value=70;
  covHIn.style.cssText='width:90px;';
  const covHVal = document.createElement('span'); covHVal.id='lbl-cov-height'; covHVal.textContent='70';
  covHVal.style.cssText='min-width:24px;font-size:11px;color:#8b949e;';
  covHRow.append(covHLbl, covHIn, covHVal, Object.assign(document.createElement('span'),{textContent:'px',style:'font-size:10px;color:#8b949e;'}));
  covSel.addEventListener('change', ()=>{
    covHRow.style.display = covSel.value !== 'none' ? 'flex' : 'none';
    renderLocal();
  });
  covHIn.addEventListener('input', ()=>{
    covHVal.textContent = covHIn.value;
    renderLocal();
  });
  covRow.append(covLbl, covSel); box.appendChild(covRow); box.appendChild(covHRow);

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
