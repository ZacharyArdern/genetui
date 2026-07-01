// ── WebSocket connection ──────────────────────────────────────────────────────
function connect() {
  ws = new WebSocket(`ws://${location.host}/ws`);
  ws.onopen = () => {
    document.getElementById('status').className='badge connected';
    document.getElementById('status').textContent='live';
    clearTimeout(reconnTimer);
  };
  ws.onmessage = ev => {
    try {
      const state=JSON.parse(ev.data);
      state.features.forEach((f,i)=>{f._idx=i;});
      if (state.blast_features) {
        state.blast_features.forEach((f,i)=>{f._idx=state.features.length+i;});
      }
      lastState=state;
      if (!localMode) {localVS=state.view_start; localVE=state.view_end;}
      renderLocal();
      scheduleAutoRender();
      if (state.protein_pdb) updateStructure(state.protein_pdb, state.protein_name||'protein');
      if (typeof onMsaStateUpdate  === 'function') onMsaStateUpdate();
      if (typeof onStatusUpdate    === 'function') onStatusUpdate();
      if (typeof onCircStateUpdate === 'function') onCircStateUpdate();
    } catch(_) {}
  };
  ws.onclose = () => {
    document.getElementById('status').className='badge disconnected';
    document.getElementById('status').textContent='reconnecting\u2026';
    reconnTimer=setTimeout(connect,2000);
  };
  ws.onerror = () => ws.close();
}
connect();
