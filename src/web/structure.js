// ── 3Dmol.js structure viewer ─────────────────────────────────────────────────
let mol3dViewer=null, mol3dLoaded=false, mol3dStyle='plddt', lastPdbName='';
let pendingPdb=null, pendingPdbName='';

// ── Horizontal drag resize between struct-panel and fig-panel ─────────────────
(function(){
  const h=document.getElementById('struct-h-handle');
  const sp=document.getElementById('struct-panel');
  let drag=false, sx=0, sw=0;
  h.addEventListener('mousedown', e=>{drag=true;sx=e.clientX;sw=sp.offsetWidth;e.preventDefault();});
  document.addEventListener('mousemove', e=>{
    if (!drag) return;
    const row=document.getElementById('bottom-row');
    const maxW=row.clientWidth-100;
    sp.style.width=Math.max(150,Math.min(sw+e.clientX-sx,maxW))+'px';
    if (mol3dViewer) mol3dViewer.resize();
  });
  document.addEventListener('mouseup', ()=>{drag=false;});
})();

function load3DmolLib(cb) {
  if (window.$3Dmol) { cb(); return; }
  const s=document.createElement('script');
  s.src='https://3Dmol.org/build/3Dmol-min.js';
  s.onload=cb;
  s.onerror=()=>{ document.getElementById('struct-ph').textContent='3Dmol.js failed to load (offline?)'; };
  document.head.appendChild(s);
}

function getBgColor() {
  return document.getElementById('chk-white-bg').checked ? 0xffffff : 0x060a10;
}

function applyAfPlddt(opacity) {
  mol3dViewer.setStyle({},{cartoon:{
    colorscheme:{prop:'b', gradient:'roygb', min:0, max:100},
    opacity
  }});
}

function applyStyle() {
  if (!mol3dViewer) return;
  mol3dViewer.removeAllSurfaces();
  mol3dViewer.setStyle({});
  const plddt={colorscheme:{prop:'b',gradient:'roygb',min:0,max:100}};
  switch(mol3dStyle) {
    case 'plddt':
      applyAfPlddt(0.72); break;
    case 'spectrum':
      mol3dViewer.setStyle({},{cartoon:{color:'spectrum', opacity:0.88}}); break;
    case 'chain':
      mol3dViewer.setStyle({},{cartoon:{colorscheme:'chainHetatm', opacity:0.88}}); break;
    case 'secondary':
      mol3dViewer.setStyle({},{cartoon:{colorscheme:'ssJmol', opacity:0.88}}); break;
    case 'surface_plddt':
      mol3dViewer.setStyle({},{cartoon:{color:'white', opacity:0.18}});
      mol3dViewer.addSurface(window.$3Dmol.SurfaceType.MS,{...plddt, opacity:0.78}); break;
    case 'surface_white':
      mol3dViewer.setStyle({},{cartoon:{color:'white', opacity:0.18}});
      mol3dViewer.addSurface(window.$3Dmol.SurfaceType.MS,{color:'white', opacity:0.72}); break;
    case 'stick':
      mol3dViewer.setStyle({},{stick:{colorscheme:'Jmol', radius:0.12}}); break;
    case 'sphere':
      mol3dViewer.setStyle({},{sphere:{colorscheme:'Jmol', radius:0.35}}); break;
  }
  const showBar = mol3dStyle==='plddt' || mol3dStyle==='surface_plddt';
  document.getElementById('plddt-bar').classList.toggle('hidden',
    !showBar || !document.getElementById('chk-plddt-bar').checked);
  mol3dViewer.render();
}

async function downloadStructurePng(transparent) {
  if (!mol3dViewer) return;
  if (transparent) {
    mol3dViewer.setBackgroundColor(0x000000, 0);
    mol3dViewer.render();
  }
  const uri=mol3dViewer.pngURI();
  const blob=await (await fetch(uri)).blob();
  dlBlob(blob, (lastPdbName||'structure')+(transparent?'_transparent':'')+'.png');
  if (transparent) {
    mol3dViewer.setBackgroundColor(getBgColor());
    mol3dViewer.render();
  }
}

function initViewer(pdb, name) {
  const el=document.getElementById('struct-view');
  el.innerHTML='';
  mol3dViewer=window.$3Dmol.createViewer(el,{backgroundColor:getBgColor(),antialias:true});
  mol3dViewer.addModel(pdb,'pdb');
  applyStyle();
  mol3dViewer.zoomTo();
  mol3dViewer.render();
  lastPdbName=name;
  document.getElementById('struct-name').textContent=name;
  mol3dLoaded=true;
}

function updateStructure(pdb, name) {
  const panel=document.getElementById('struct-panel');
  panel.classList.add('visible');
  document.getElementById('struct-h-handle').classList.add('visible');
  if (window._panelCbs && window._panelCbs['struct-panel']) {
    window._panelCbs['struct-panel'].disabled=false;
    window._panelCbs['struct-panel'].checked=true;
  }
  if (!window.$3Dmol) {
    pendingPdb=pdb; pendingPdbName=name;
    load3DmolLib(()=>{ if (pendingPdb) { initViewer(pendingPdb,pendingPdbName); pendingPdb=null; } });
    return;
  }
  if (mol3dLoaded && name===lastPdbName) return;
  initViewer(pdb, name);
}

// ── Structure button handlers ─────────────────────────────────────────────────
document.getElementById('struct-style-sel').addEventListener('change', function() {
  mol3dStyle=this.value;
  if (mol3dLoaded) applyStyle();
});
document.getElementById('chk-plddt-bar').addEventListener('change', ()=>{ if (mol3dLoaded) applyStyle(); });
document.getElementById('chk-white-bg').addEventListener('change', ()=>{
  if (!mol3dViewer) return;
  mol3dViewer.setBackgroundColor(getBgColor());
  document.getElementById('struct-view').style.background = document.getElementById('chk-white-bg').checked ? '#fff' : '#060a10';
  mol3dViewer.render();
});
document.getElementById('btn-struct-png').addEventListener('click', ()=>downloadStructurePng(false));
document.getElementById('btn-struct-png-t').addEventListener('click', ()=>downloadStructurePng(true));
document.getElementById('btn-struct-dl').addEventListener('click', ()=>{
  if (!lastState||!lastState.protein_pdb) return;
  dlBlob(new Blob([lastState.protein_pdb],{type:'chemical/x-pdb'}),lastPdbName+'.pdb');
});
