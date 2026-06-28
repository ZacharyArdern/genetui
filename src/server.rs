use std::sync::{Arc, Mutex};
use tokio::sync::broadcast;
use axum::{
    Router,
    extract::{
        ws::{Message, WebSocket, WebSocketUpgrade},
        State,
    },
    response::{Html, IntoResponse},
    routing::get,
};
use futures_util::{sink::SinkExt, stream::StreamExt};
use serde::{Deserialize, Serialize};

#[derive(Clone, Serialize, Deserialize, Default, Debug)]
pub struct BrowserState {
    pub view_start: u64,
    pub view_end:   u64,
    pub genome_size: u64,
    pub genome_name: String,
    pub features: Vec<BrowserFeature>,
    pub sequence: String,   // DNA bases for seq_start..(seq_start+sequence.len()) (may be empty)
    pub seq_start: u64,    // genomic offset of sequence[0]
    pub protein_pdb: String,   // raw PDB text of selected+folded protein (may be empty)
    pub protein_name: String,  // gene name for the PDB
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct BrowserFeature {
    pub name:      String,
    pub start:     u64,
    pub end:       u64,
    pub strand:    char,
    pub color:     String,
    pub noncoding: bool,
}

pub struct WebServer {
    state:        Arc<Mutex<BrowserState>>,
    tx:           broadcast::Sender<String>,
    viewport_cmd: Arc<Mutex<Option<(u64, u64)>>>,
}

impl WebServer {
    pub fn new() -> Arc<Self> {
        let (tx, _) = broadcast::channel(32);
        Arc::new(Self {
            state:        Arc::new(Mutex::new(BrowserState::default())),
            tx,
            viewport_cmd: Arc::new(Mutex::new(None)),
        })
    }

    pub fn take_viewport_cmd(&self) -> Option<(u64, u64)> {
        self.viewport_cmd.lock().ok()?.take()
    }

    /// Push new state to all connected browser clients (sync-safe, callable from event loop).
    pub fn push(&self, state: BrowserState) {
        if let Ok(json) = serde_json::to_string(&state) {
            if let Ok(mut s) = self.state.lock() { *s = state; }
            let _ = self.tx.send(json);
        }
    }

    pub async fn run(self: Arc<Self>, port: u16) {
        let app = Router::new()
            .route("/",   get(serve_html))
            .route("/ws", get(ws_handler))
            .with_state(self);

        match tokio::net::TcpListener::bind(format!("127.0.0.1:{}", port)).await {
            Ok(listener) => { let _ = axum::serve(listener, app).await; }
            Err(e) => eprintln!("genetui web server: could not bind port {port}: {e}"),
        }
    }
}

async fn serve_html() -> Html<&'static str> { Html(HTML) }

async fn ws_handler(
    ws: WebSocketUpgrade,
    State(server): State<Arc<WebServer>>,
) -> impl IntoResponse {
    ws.on_upgrade(move |socket| handle_ws(socket, server))
}

async fn handle_ws(socket: WebSocket, server: Arc<WebServer>) {
    let (mut sink, mut stream) = socket.split();

    // Send current state immediately on connect (drop lock before await).
    let initial = server.state.lock().ok()
        .and_then(|s| serde_json::to_string(&*s).ok());
    if let Some(json) = initial {
        let _ = sink.send(Message::Text(json.into())).await;
    }

    let mut rx = server.tx.subscribe();
    loop {
        tokio::select! {
            Ok(json) = rx.recv() => {
                if sink.send(Message::Text(json.into())).await.is_err() { break; }
            }
            Some(msg) = stream.next() => {
                match msg {
                    Ok(Message::Close(_)) | Err(_) => break,
                    Ok(Message::Text(txt)) => {
                        if let Ok(v) = serde_json::from_str::<serde_json::Value>(&txt) {
                            if let (Some(s), Some(e)) = (v["start"].as_u64(), v["end"].as_u64()) {
                                if let Ok(mut cmd) = server.viewport_cmd.lock() {
                                    *cmd = Some((s, e));
                                }
                            }
                        }
                    }
                    _ => {}
                }
            }
            else => break,
        }
    }
}

// ── Embedded web app ─────────────────────────────────────────────────────────

const HTML: &str = r#####"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>genetui</title>
<style>
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
  background: #0d1117; color: #c9d1d9;
  font-family: 'Segoe UI', system-ui, sans-serif;
  height: 100vh; display: flex; flex-direction: column; overflow: hidden;
}

/* ── Header ──────────────────────────────────────────────────────────────── */
#header {
  height: 34px; padding: 0 12px; background: #161b22;
  border-bottom: 1px solid #30363d;
  display: flex; align-items: center; gap: 10px; flex-shrink: 0;
}
#genome-name { font-size: 13px; font-weight: 600; color: #58a6ff; }
#position    { font-size: 11px; color: #8b949e; font-family: monospace; }
.badge { font-size: 10px; padding: 2px 7px; border-radius: 10px; white-space: nowrap; }
#status { margin-left: auto; }
.connected    { background: #1a7f37; color: #fff; }
.disconnected { background: #6e7681; color: #fff; }
#py-badge { background: #6e3694; color: #fff; display: none; }
#py-badge.ready { background: #1a7f37; }
.hbtn {
  padding: 2px 9px; border: 1px solid #30363d; background: #21262d;
  color: #c9d1d9; border-radius: 5px; font-size: 11px; cursor: pointer;
}
.hbtn:hover { background: #30363d; }
.hbtn.active { background: #1f6feb; border-color: #1f6feb; color: #fff; }

/* ── Live SVG track ──────────────────────────────────────────────────────── */
#live-panel { flex-shrink: 0; overflow: hidden; position: relative; cursor: grab; }
#live-panel.dragging { cursor: grabbing; }
#svg { width: 100%; height: 100%; display: block; }
.gene:hover polygon { opacity: 0.75; }

/* ── Drag handle ─────────────────────────────────────────────────────────── */
#drag-handle {
  height: 5px; background: #21262d; cursor: ns-resize; flex-shrink: 0;
  border-top: 1px solid #30363d; border-bottom: 1px solid #30363d;
}
#drag-handle:hover { background: #58a6ff; }

/* ── Bottom row (structure left + gene plot right) ───────────────────────── */
#bottom-row { flex: 1; display: flex; flex-direction: row; min-height: 0; overflow: hidden; }

/* ── Structure panel (left, hidden until protein folded) ─────────────────── */
#struct-panel { width: 38%; min-width: 200px; display: none; flex-direction: column; overflow: hidden; flex-shrink: 0; }
#struct-panel.visible { display: flex; }
#struct-h-handle {
  width: 5px; background: #21262d; cursor: ew-resize; flex-shrink: 0;
  border-left: 1px solid #30363d; border-right: 1px solid #30363d; display: none;
}
#struct-h-handle.visible { display: block; }
#struct-h-handle:hover { background: #58a6ff; }
#struct-bar {
  height: 30px; padding: 0 10px; background: #161b22; border-bottom: 1px solid #30363d;
  display: flex; align-items: center; gap: 8px; flex-shrink: 0;
}
#struct-bar-label { font-size: 11px; color: #8b949e; }
#struct-name { font-size: 11px; color: #58a6ff; font-family: monospace; }
#struct-view-wrap { flex: 1; display: flex; flex-direction: row; min-height: 0; }
#struct-view { flex: 1; position: relative; background: #060a10; }
#struct-ph { color: #484f58; font-size: 12px; text-align: center; margin-top: 40px; }
#plddt-bar {
  width: 22px; flex-shrink: 0; display: flex; flex-direction: column; align-items: center;
  padding: 8px 2px; gap: 2px; background: transparent; pointer-events: none;
}
#plddt-bar.hidden { display: none; }
#plddt-bar-gradient {
  width: 10px; flex: 1; border-radius: 5px;
  background: linear-gradient(to bottom, #0000ff 0%, #00ff00 40%, #ffff00 60%, #ff7700 80%, #ff0000 100%);
}
#plddt-bar-labels {
  position: absolute; right: 20px; top: 8px; bottom: 8px; display: flex; flex-direction: column;
  justify-content: space-between; pointer-events: none;
}
#plddt-bar-labels span { font-size: 8px; color: #8b949e; line-height: 1; }

/* ── Custom gene plot panel (right) ──────────────────────────────────────── */
#fig-panel { flex: 1; min-width: 0; min-height: 0; display: flex; flex-direction: column; overflow: hidden; }
#fig-bar {
  height: 30px; padding: 0 10px; background: #161b22; border-bottom: 1px solid #30363d;
  display: flex; align-items: center; gap: 8px; flex-shrink: 0;
}
#fig-bar-label { font-size: 11px; color: #8b949e; }
#fig-status    { font-size: 10px; color: #6e7681; font-family: monospace; }
.dbtn {
  padding: 2px 8px; border: 1px solid #30363d; background: #21262d;
  color: #8b949e; border-radius: 4px; font-size: 10px; cursor: pointer;
}
.dbtn:hover:not(:disabled) { color: #c9d1d9; }
.dbtn:disabled { opacity: 0.35; cursor: default; }
#fig-output {
  flex: 1; overflow: auto; padding: 10px;
  display: flex; justify-content: center; align-items: flex-start;
}
#fig-output img { max-width: 100%; border-radius: 4px; }
.fig-ph { color: #484f58; font-size: 12px; text-align: center; margin-top: 20px; line-height: 1.8; }
.fig-err { color: #f85149; font-size: 11px; font-family: monospace; white-space: pre-wrap; word-break: break-all; padding: 8px; }

/* ── Options dropdown ────────────────────────────────────────────────────── */
#opts-box {
  position: fixed; z-index: 300;
  background: #161b22; border: 1px solid #30363d; border-radius: 8px;
  padding: 14px 16px; width: 380px; box-shadow: 0 8px 32px rgba(0,0,0,0.6);
  cursor: default; user-select: none;
}
#opts-box h3 { font-size: 11px; color: #8b949e; text-transform: uppercase; letter-spacing: .06em; margin-bottom: 10px; }
.og { display: grid; grid-template-columns: 1fr 1fr; gap: 8px 18px; }
.cr { display: flex; align-items: center; gap: 6px; font-size: 11px; color: #8b949e; }
.cr label { white-space: nowrap; }
.cr input[type=range]  { flex: 1; min-width: 70px; accent-color: #58a6ff; }
.cr select { background: #0d1117; color: #c9d1d9; border: 1px solid #30363d; border-radius: 4px; padding: 2px 5px; font-size: 11px; }
.cr input[type=number] { width: 60px; background: #0d1117; color: #c9d1d9; border: 1px solid #30363d; border-radius: 4px; padding: 2px 5px; font-size: 11px; }
.cr input[type=checkbox] { accent-color: #58a6ff; }
.cv { color: #c9d1d9; font-family: monospace; min-width: 26px; }
.osep { border: none; border-top: 1px solid #21262d; margin: 8px 0; }
</style>
</head>
<body>

<!-- Header -->
<div id="header">
  <span id="genome-name">genetui</span>
  <span id="position">—</span>
  <span id="py-badge" class="badge">pyodide</span>
  <span id="status" class="badge disconnected">disconnected</span>
</div>

<!-- Live SVG track -->
<div id="live-panel" style="height:55vh">
  <svg id="svg" xmlns="http://www.w3.org/2000/svg"></svg>
</div>

<!-- Resize handle -->
<div id="drag-handle"></div>

<!-- Bottom row: structure (left) + custom gene plot (right) -->
<div id="bottom-row">

  <!-- Structure panel (appears when protein is folded in TUI) -->
  <div id="struct-panel">
    <div id="struct-bar">
      <span id="struct-bar-label">Structure</span>
      <span id="struct-name"></span>
      <div style="margin-left:auto;display:flex;gap:5px;align-items:center">
        <select id="struct-style-sel" class="dbtn" style="padding:2px 4px;cursor:pointer">
          <option value="plddt">pLDDT colouring</option>
          <option value="spectrum">Spectrum (N→C)</option>
          <option value="chain">Chain colours</option>
          <option value="secondary">Secondary structure</option>
          <option value="surface_plddt">Surface — pLDDT</option>
          <option value="surface_white">Surface — white</option>
          <option value="stick">Stick</option>
          <option value="sphere">Sphere</option>
        </select>
        <label style="font-size:10px;color:#8b949e;display:flex;align-items:center;gap:3px;cursor:pointer">
          <input type="checkbox" id="chk-plddt-bar" checked style="cursor:pointer"> pLDDT bar
        </label>
        <label style="font-size:10px;color:#8b949e;display:flex;align-items:center;gap:3px;cursor:pointer">
          <input type="checkbox" id="chk-white-bg" style="cursor:pointer"> White bg
        </label>
        <button class="dbtn" id="btn-struct-png">PNG</button>
        <button class="dbtn" id="btn-struct-png-t">PNG transp.</button>
        <button class="dbtn" id="btn-struct-dl">PDB</button>
      </div>
    </div>
    <div id="struct-view-wrap">
      <div id="struct-view"><div id="struct-ph">Fold a protein in the TUI (select gene &#8594; f) to view its structure here.</div></div>
      <div id="plddt-bar">
        <span style="font-size:8px;color:#8b949e">100</span>
        <div id="plddt-bar-gradient"></div>
        <span style="font-size:8px;color:#8b949e"> 0 </span>
      </div>
    </div>
  </div>

  <!-- Horizontal resize handle between structure and gene plot -->
  <div id="struct-h-handle"></div>

  <!-- Custom gene plot panel -->
  <div id="fig-panel">
    <div id="fig-bar">
      <span id="fig-bar-label">Custom gene plot</span>
      <span id="fig-status">loading Pyodide…</span>
      <div style="margin-left:auto;display:flex;gap:5px">
        <button class="dbtn" id="btn-opts">&#9881; Options</button>
        <button class="dbtn" id="btn-render" disabled>Render</button>
        <button class="dbtn" id="btn-svg"   disabled>SVG</button>
        <button class="dbtn" id="btn-png"   disabled>PNG 600dpi</button>
      </div>
    </div>
    <div id="fig-output">
      <div class="fig-ph">Custom gene plot appears here automatically when scrolling stops.<br>Pyodide + dna_features_viewer is loading in the background.</div>
    </div>
  </div>

</div>

<!-- Options dropdown (anchored to btn-opts, not an overlay) -->
<div id="opts-box" style="display:none">
  <h3>Custom gene plot options</h3>
  <div class="og">
    <div class="cr">
      <label>Width</label>
      <input type="range" id="fig-width" min="4" max="24" step="0.5" value="10">
      <span class="cv" id="lbl-width">10</span><span style="font-size:10px">in</span>
    </div>
    <div class="cr">
      <label>Font</label>
      <input type="range" id="fig-fs" min="5" max="16" step="1" value="11">
      <span class="cv" id="lbl-fs">11</span><span style="font-size:10px">pt</span>
    </div>
    <div class="cr">
      <input type="checkbox" id="fig-ruler" checked>
      <label for="fig-ruler">Ruler</label>
    </div>
    <div class="cr">
      <input type="checkbox" id="fig-white" checked>
      <label for="fig-white">White background</label>
    </div>
    <div class="cr">
      <input type="checkbox" id="fig-show-labels" checked>
      <label for="fig-show-labels">Labels</label>
    </div>
    <div class="cr">
      <input type="checkbox" id="fig-hide-long-labels">
      <label for="fig-hide-long-labels">Hide overflow labels</label>
    </div>
    <div class="cr" id="overflow-thresh-row" style="display:none">
      <label>Overflow sensitivity</label>
      <input type="range" id="fig-overflow-thresh" min="0.07" max="0.18" step="0.005" value="0.105">
      <span class="cv" id="lbl-overflow-thresh">0.105</span>
    </div>
  </div>
  <hr class="osep">
  <div class="og">
    <div class="cr">
      <input type="checkbox" id="opt-sixframe" checked>
      <label for="opt-sixframe">Six-frame layout (live + figure)</label>
    </div>
    <div class="cr">
      <input type="checkbox" id="opt-stopcodons" checked>
      <label for="opt-stopcodons">Stop codons (live + figure)</label>
    </div>
  </div>
</div>

<script>
const NC_COLOR = '#6e7681';
const RULER_H  = 28;
const FH       = 19;   // frame-row height (six-frame mode)
const LEVEL_GAP = 24;  // gene-level gap (simple mode)
const ARROW_H  = 14;
const STOP_COL = '#f85149';
const STOP_SET = new Set(['TAA','TAG','TGA']);
const COMP     = {A:'T',T:'A',G:'C',C:'G',a:'t',t:'a',g:'c',c:'g'};

let lastState = null, localVS = 0, localVE = 0;
let localMode = false, localModeTimer = null;
let pyodide = null, pyReady = false, pyRendering = false;
let lastSvgData = null, debounceTimer = null;
let ws, reconnTimer, dragStart = null, dragVS0 = 0, dragVE0 = 0;

// ── Helpers ───────────────────────────────────────────────────────────────────
function fmt(n) {
  if (n >= 1e6) return (n/1e6).toFixed(2)+' Mb';
  if (n >= 1e3) return (n/1e3).toFixed(1)+' kb';
  return n+' bp';
}
function trunc(s,n)  { return s.length<=n ? s : s.slice(0,n-1)+'\u2026'; }
function esc(s)      { return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
function setFigStatus(m) { document.getElementById('fig-status').textContent = m; }
function niceTick(span, W) {
  const raw = span / Math.max(Math.floor(W/90),1);
  const mag = Math.pow(10, Math.floor(Math.log10(Math.max(raw,1))));
  for (const m of [1,2,5,10]) { if (mag*m >= raw) return mag*m; }
  return mag*10;
}

// ── Level assignment ──────────────────────────────────────────────────────────
function assignLevels(feats, scale, vs) {
  const ends = [];
  return feats.map(f => {
    const x1 = (f.start-vs)*scale, x2 = (f.end-vs)*scale;
    let lv = 0;
    while (ends[lv] !== undefined && ends[lv] > x1-6) lv++;
    ends[lv] = x2+6;
    return {...f, x1, x2, lv};
  });
}

// ── Arrow polygon ─────────────────────────────────────────────────────────────
function arrowPts(x1,x2,cy,h,strand) {
  const w=x2-x1, ah=Math.min(h*.9,w*.25,10), t=cy-h/2, b=cy+h/2;
  if (strand==='+') {
    if (w<3) return `${x1},${t} ${x2},${cy} ${x1},${b}`;
    return `${x1},${t} ${x2-ah},${t} ${x2},${cy} ${x2-ah},${b} ${x1},${b}`;
  } else {
    if (w<3) return `${x2},${t} ${x1},${cy} ${x2},${b}`;
    return `${x2},${t} ${x1+ah},${t} ${x1},${cy} ${x1+ah},${b} ${x2},${b}`;
  }
}

// ── Six-frame helpers ─────────────────────────────────────────────────────────
function rcSeq(s) {
  let r=''; for (let i=s.length-1;i>=0;i--) r += COMP[s[i]]||'N'; return r;
}
function stopPositions(seq, frame) {
  const out=[], u=seq.toUpperCase();
  for (let i=frame; i+2<u.length; i+=3) { if (STOP_SET.has(u.slice(i,i+3))) out.push(i); }
  return out;
}

// ── SVG render ────────────────────────────────────────────────────────────────
function drawSVG(state, vs, ve) {
  const svgEl = document.getElementById('svg');
  const panel = document.getElementById('live-panel');
  const W = panel.clientWidth  || window.innerWidth;
  const H = panel.clientHeight || 300;
  const {genome_size, genome_name, features, sequence='', seq_start: seqOff=0} = state;
  const span  = Math.max(ve-vs,1), scale = W/span;
  const sixFrame = document.getElementById('opt-sixframe').checked;

  document.getElementById('genome-name').textContent = genome_name || 'genome';
  document.getElementById('position').textContent =
    `${fmt(vs)} \u2013 ${fmt(ve)}  (${fmt(span)})`;

  const out = [];

  // Ruler
  const tick=niceTick(span,W), first=Math.ceil(vs/tick)*tick;
  out.push(`<line x1="0" y1="${RULER_H}" x2="${W}" y2="${RULER_H}" stroke="#30363d" stroke-width="1"/>`);
  for (let pos=first; pos<=ve; pos+=tick) {
    const x=(pos-vs)*scale;
    out.push(`<line x1="${x|0}" y1="${RULER_H-7}" x2="${x|0}" y2="${RULER_H}" stroke="#484f58" stroke-width="1"/>`);
    out.push(`<text x="${(x+3)|0}" y="${RULER_H-10}" font-size="10" fill="#8b949e" font-family="monospace">${fmt(pos)}</text>`);
  }

  if (sixFrame) {
    drawSixFrame(features, sequence, seqOff, vs, ve, W, H, scale, out);
  } else {
    drawSimple(features, vs, ve, W, H, scale, out);
  }

  // Minimap
  if (genome_size > 0) {
    const mw=Math.min(W-20,200), mh=4, mx=W-mw-8, my=H-12;
    out.push(`<rect x="${mx}" y="${my}" width="${mw}" height="${mh}" fill="#21262d" rx="2"/>`);
    const vx1=mx+(vs/genome_size)*mw, vx2=mx+(ve/genome_size)*mw;
    out.push(`<rect x="${vx1|0}" y="${my}" width="${Math.max(2,(vx2-vx1))|0}" height="${mh}" fill="#58a6ff" rx="1"/>`);
  }

  svgEl.innerHTML = out.join('');
}

// ── Simple +/- track ──────────────────────────────────────────────────────────
function drawSimple(features, vs, ve, W, H, scale, out) {
  const SGAP = 8;
  const plus  = assignLevels(features.filter(f=>f.strand==='+'), scale, vs);
  const minus = assignLevels(features.filter(f=>f.strand==='-'), scale, vs);
  const nPlus = plus.length  ? Math.max(0,...plus.map(f=>f.lv))+1 : 0;
  const divY  = RULER_H + nPlus*LEVEL_GAP + SGAP;
  out.push(`<line x1="0" y1="${divY}" x2="${W}" y2="${divY}" stroke="#21262d" stroke-width="1"/>`);
  const drawG = (f, cy, strand) => {
    const col = f.noncoding ? NC_COLOR : f.color;
    const x1c=Math.max(0,f.x1), x2c=Math.min(W,f.x2);
    if (x2c<=x1c) return;
    const pts=arrowPts(x1c,x2c,cy,ARROW_H,strand), lw=x2c-x1c;
    out.push(`<g class="gene"><title>${esc(f.name)}</title><polygon points="${pts}" fill="${col}" fill-opacity="0.88"/>`);
    if (lw>24) {
      const lbl=trunc(f.name,Math.floor(lw/6.5));
      const ly = strand==='+'? (cy-ARROW_H/2-3):(cy+ARROW_H/2+11);
      out.push(`<text x="${((x1c+x2c)/2)|0}" y="${ly|0}" text-anchor="middle" font-size="10" fill="#c9d1d9" font-family="monospace">${esc(lbl)}</text>`);
    }
    out.push('</g>');
  };
  plus.forEach(f  => drawG(f, divY-SGAP-ARROW_H/2-f.lv*LEVEL_GAP, '+'));
  minus.forEach(f => drawG(f, divY+SGAP+ARROW_H/2+f.lv*LEVEL_GAP, '-'));
}

// ── Six-frame translation track ───────────────────────────────────────────────
function drawSixFrame(features, seq, seqOff, vs, ve, W, H, scale, out) {
  const rc      = rcSeq(seq);
  const seqLen  = seq.length;
  const showSC  = document.getElementById('opt-stopcodons').checked;
  // Bin stop codons by ABSOLUTE GENOMIC frame (0/1/2) so they align with gene frames.
  // Gene frame (1-based coords): (f.start-1)%3 (+) or (f.end-1)%3 (-).
  // seqOff is 0-based, so (seqOff+p)%3 == (1-based-1)%3 — already consistent.
  const fwdStops = [[],[],[]];
  const revStops = [[],[],[]];
  if (showSC && seq) {
    for (let rcFr = 0; rcFr < 3; rcFr++) {
      for (const p of stopPositions(seq, rcFr)) {
        const gp = seqOff + p;          // genomic start of codon
        fwdStops[gp % 3].push(gp);
      }
      for (const p of stopPositions(rc, rcFr)) {
        const gEnd = seqOff + seqLen - p - 1;   // genomic end (3′) of minus codon
        revStops[gEnd % 3].push(seqOff + seqLen - p - 2);  // display at codon midpoint
      }
    }
  }

  const divY = RULER_H + FH*3;
  const ncYp = divY + FH*3 + 3;   // non-coding + row
  const ncYm = ncYp + FH + 2;     // non-coding - row

  // Frame row backgrounds
  for (let fr=0;fr<3;fr++) {
    const bg = fr%2 ? '#0f1520' : '#0d1117';
    out.push(`<rect x="0" y="${RULER_H+FH*fr}" width="${W}" height="${FH}" fill="${bg}"/>`);
    out.push(`<rect x="0" y="${divY+FH*fr}" width="${W}" height="${FH}" fill="${bg}"/>`);
  }

  // Strand divider
  out.push(`<line x1="0" y1="${divY}" x2="${W}" y2="${divY}" stroke="#30363d" stroke-width="1"/>`);

  // Stop codon ticks — forward
  for (let fr=0;fr<3;fr++) {
    const y1=RULER_H+FH*fr+2, y2=RULER_H+FH*(fr+1)-2;
    for (const gp of fwdStops[fr]) {
      const x=(gp-vs)*scale;
      if (x<-1||x>W+1) continue;
      out.push(`<line x1="${x|0}" y1="${y1}" x2="${x|0}" y2="${y2}" stroke="${STOP_COL}" stroke-width="1" opacity="0.6"/>`);
    }
    out.push(`<text x="3" y="${(RULER_H+FH*fr+FH/2+3)|0}" font-size="8" fill="#484f58" font-family="monospace">+${fr+1}</text>`);
  }
  // Stop codon ticks — reverse
  for (let fr=0;fr<3;fr++) {
    const y1=divY+FH*fr+2, y2=divY+FH*(fr+1)-2;
    for (const gp of revStops[fr]) {
      const x=(gp-vs)*scale;
      if (x<-1||x>W+1) continue;
      out.push(`<line x1="${x|0}" y1="${y1}" x2="${x|0}" y2="${y2}" stroke="${STOP_COL}" stroke-width="1" opacity="0.6"/>`);
    }
    out.push(`<text x="3" y="${(divY+FH*fr+FH/2+3)|0}" font-size="8" fill="#484f58" font-family="monospace">-${fr+1}</text>`);
  }

  // Coding genes in their frame row
  const gh = FH-6;
  for (const f of features) {
    if (f.noncoding) continue;
    const x1=(f.start-vs)*scale, x2=(f.end-vs)*scale;
    const x1c=Math.max(0,x1), x2c=Math.min(W,x2);
    if (x2c<=x1c) continue;
    const col = f.color;
    const fr  = f.strand==='+' ? (f.start-1)%3 : (f.end-1)%3;
    const cy  = f.strand==='+' ? RULER_H+FH*fr+FH/2 : divY+FH*fr+FH/2;
    const pts = arrowPts(x1c, x2c, cy, gh, f.strand);
    const lw  = x2c-x1c;
    out.push(`<g class="gene"><title>${esc(f.name)}</title><polygon points="${pts}" fill="${col}" fill-opacity="0.9"/>`);
    if (lw>18) {
      const lbl = trunc(f.name, Math.floor(lw/5.5));
      out.push(`<text x="${((x1c+x2c)/2)|0}" y="${(cy+3)|0}" text-anchor="middle" font-size="9" fill="#fff" fill-opacity="0.85" font-family="monospace">${esc(lbl)}</text>`);
    }
    out.push('</g>');
  }

  // Non-coding track
  out.push(`<line x1="0" y1="${ncYp-2}" x2="${W}" y2="${ncYp-2}" stroke="#21262d" stroke-width="1"/>`);
  out.push(`<text x="3" y="${(ncYp+FH/2+3)|0}" font-size="8" fill="#484f58" font-family="monospace">nc+</text>`);
  out.push(`<text x="3" y="${(ncYm+FH/2+3)|0}" font-size="8" fill="#484f58" font-family="monospace">nc-</text>`);
  for (const f of features) {
    if (!f.noncoding) continue;
    const x1=(f.start-vs)*scale, x2=(f.end-vs)*scale;
    const x1c=Math.max(0,x1), x2c=Math.min(W,x2);
    if (x2c<=x1c) continue;
    const cy = f.strand==='+' ? ncYp+FH/2 : ncYm+FH/2;
    const pts = arrowPts(x1c, x2c, cy, FH-6, f.strand);
    out.push(`<g class="gene"><title>${esc(f.name)}</title><polygon points="${pts}" fill="${NC_COLOR}" fill-opacity="0.7"/></g>`);
  }
}

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

// ── Navigation — wheel (on document, path-checked) ────────────────────────────
document.addEventListener('wheel', e => {
  const panel = document.getElementById('live-panel');
  if (!e.composedPath().some(el => el===panel)) return;
  e.preventDefault();
  if (!lastState) return;
  enterLocalMode();
  const rect=panel.getBoundingClientRect(), frac=(e.clientX-rect.left)/panel.clientWidth;
  const span=localVE-localVS, factor=e.deltaY>0?1.18:1/1.18;
  const newSpan=Math.max(10,span*factor), center=localVS+frac*span;
  localVS=Math.round(center-frac*newSpan);
  localVE=Math.round(localVS+newSpan);
  clampLocal(); renderLocal(); scheduleAutoRender();
}, {passive:false});

// Drag pan
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
});
document.addEventListener('mouseup', () => {
  if (dragStart!==null) {
    dragStart=null;
    document.getElementById('live-panel').classList.remove('dragging');
    scheduleAutoRender();
  }
});

// Arrow keys
document.addEventListener('keydown', e => {
  if (!lastState) return;
  const tag=document.activeElement&&document.activeElement.tagName;
  if (tag==='INPUT'||tag==='SELECT'||tag==='TEXTAREA') return;
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
});

window.addEventListener('resize', renderLocal);

// Drag-resize handle
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

// ── WebSocket ─────────────────────────────────────────────────────────────────
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
      lastState=state;
      if (!localMode) {localVS=state.view_start; localVE=state.view_end;}
      renderLocal();
      scheduleAutoRender();
      if (state.protein_pdb) updateStructure(state.protein_pdb, state.protein_name||'protein');
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

// ── Options dropdown (anchored to btn-opts, draggable) ───────────────────────
(function(){
  const btn=document.getElementById('btn-opts');
  const box=document.getElementById('opts-box');
  let open=false, dragging=false, ox=0, oy=0, mx=0, my=0;

  function positionBox() {
    const r=btn.getBoundingClientRect();
    // Prefer anchoring below-right, flip if near right edge
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

  // Draggable by header
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

// Auto-render on any figure option change
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

// ── Panel selector (d key) ────────────────────────────────────────────────────
(function(){
  // struct-panel is only included when protein data is present — always listed but
  // its checkbox stays disabled until structure data arrives.
  const PANELS=[
    {id:'live-panel',   label:'Synced browser'},
    {id:'fig-panel',    label:'Custom gene plot'},
    {id:'struct-panel', label:'Structure viewer'},
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
    // struct-panel starts hidden & unchecked until protein data arrives
    const isStruct=(p.id==='struct-panel');
    cb.checked=!isStruct; cb.disabled=isStruct;
    cb.addEventListener('change', ()=>{
      const el=document.getElementById(p.id);
      if (!el) return;
      if (p.id==='live-panel') {
        el.style.display=cb.checked?'':'none';
        document.getElementById('drag-handle').style.display=cb.checked?'':'none';
      } else if (p.id==='struct-panel') {
        el.classList.toggle('visible', cb.checked);
        document.getElementById('struct-h-handle').classList.toggle('visible', cb.checked);
      } else {
        el.style.display=cb.checked?'':'none';
      }
    });
    panelCbs[p.id]=cb;
    const lbl=document.createElement('label'); lbl.textContent=p.label;
    row.append(cb,lbl); box.appendChild(row);
  });
  ov.appendChild(box); document.body.appendChild(ov);
  window._panelCbs=panelCbs;  // expose so struct viewer can enable checkbox
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

// ── Pyodide ───────────────────────────────────────────────────────────────────
async function initPyodide() {
  const badge=document.getElementById('py-badge');
  badge.style.display=''; badge.textContent='py\u2026';
  setFigStatus('loading Pyodide\u2026');
  try {
    pyodide = await loadPyodide({indexURL:'https://cdn.jsdelivr.net/pyodide/v0.26.4/full/'});
    setFigStatus('installing packages\u2026');
    await pyodide.loadPackage(['micropip','matplotlib','numpy']);
    await pyodide.runPythonAsync('import micropip\nawait micropip.install("dna_features_viewer", keep_going=True)');
    badge.textContent='py ready'; badge.className='badge ready';
    setFigStatus('ready \u2014 renders when scrolling stops');
    pyReady=true;
    document.getElementById('btn-render').disabled=false;
    if (lastState) scheduleAutoRender();
  } catch(err) {
    badge.textContent='py error'; badge.style.background='#da3633';
    setFigStatus('Pyodide error');
    document.getElementById('fig-output').innerHTML=`<pre class="fig-err">${esc(String(err))}</pre>`;
  }
}
(function(){
  const s=document.createElement('script');
  s.src='https://cdn.jsdelivr.net/pyodide/v0.26.4/full/pyodide.js';
  s.onload=()=>initPyodide();
  s.onerror=()=>{
    document.getElementById('py-badge').style.display='';
    document.getElementById('py-badge').textContent='offline';
    setFigStatus('CDN unavailable');
  };
  document.head.appendChild(s);
})();

// ── Auto-render debounce (always on) ─────────────────────────────────────────
function scheduleAutoRender() {
  if (!pyReady) return;
  clearTimeout(debounceTimer);
  debounceTimer=setTimeout(()=>triggerRender(false), 700);
}

// ── Python code builder ───────────────────────────────────────────────────────
function buildPyCode(vs,ve,forPng) {
  if (!lastState) return null;
  const width     = parseFloat(document.getElementById('fig-width').value);
  const fs        = parseInt(document.getElementById('fig-fs').value);
  const showLbls  = document.getElementById('fig-show-labels').checked;
  const hideLong  = document.getElementById('fig-hide-long-labels').checked;
  const ruler     = document.getElementById('fig-ruler').checked;
  const white     = document.getElementById('fig-white').checked;
  const sixframe  = document.getElementById('opt-sixframe').checked;
  const stopCod   = document.getElementById('opt-stopcodons').checked;
  const seqLen    = Math.max(1,ve-vs);
  const feats     = lastState.features.filter(f=>f.start<=ve&&f.end>=vs);
  const bg        = white?"'white'":"'#0d1117'";
  const fg        = white?"'#111'":"'#c9d1d9'";
  const dpi       = forPng?600:96;
  const fmt_      = forPng?'png':'svg';

  // Approximate character width at current font size (inches per char, tuned to dna_features_viewer defaults)
  const charWidthIn = (fs / 11) * parseFloat(document.getElementById('fig-overflow-thresh').value);
  function labelFitsInFeature(f) {
    const visibleBp = Math.min(f.end, ve) - Math.max(f.start, vs);
    const featWidthIn = visibleBp / seqLen * width;
    return featWidthIn >= f.name.length * charWidthIn;
  }

  // Build feature tuples: (start_rel, end_rel, strand_int, label, color, noncoding, frame)
  const pyFeats = feats.map(f=>{
    const lbl = showLbls && !(hideLong && !labelFitsInFeature(f)) ? f.name.replace(/\\/g,'\\\\').replace(/'/g,"\\'") : '';
    const col = f.noncoding ? NC_COLOR : f.color;
    const s  = Math.max(0, f.start-vs);
    const ee = Math.min(seqLen, f.end-vs);
    const si = f.strand==='+'?1:-1;
    const fr = f.strand==='+'? (f.start-1)%3 : (f.end-1)%3;
    return `    (${s},${ee},${si},${lbl?`'${lbl}'`:'None'},'${col}',${f.noncoding?'True':'False'},${fr})`;
  }).join(',\n');

  if (sixframe) {
    return `
import io,base64,traceback,matplotlib,matplotlib.gridspec as gridspec
matplotlib.use('Agg')
import matplotlib.pyplot as plt
_result=None
try:
    from dna_features_viewer import GraphicRecord,GraphicFeature
    _bg=${bg}; _fg=${fg}
    _seqLen=${seqLen}; _width=${width}; _ruler=${ruler?'True':'False'}
    _stopCod=${stopCod?'True':'False'}
    _feats_raw=[
${pyFeats}
    ]
    _STOP={'TAA','TAG','TGA'}
    _COMP={'A':'T','T':'A','G':'C','C':'G'}
    def _rc(s):
        return ''.join(_COMP.get(c,'N') for c in reversed(s.upper()))
    _seq_sub=globals().get('_seq_sub',''); _seq_sub_off=globals().get('_seq_sub_off',0)
    # Genomic start of _seq_sub (used to compute absolute genomic frame of each codon)
    _geno_off=${vs}+_seq_sub_off
    # Pre-bin stop codons by absolute genomic frame so they align with gene frames.
    # Gene frame: start%3 (+) or end%3 (-) — both absolute genomic coords.
    _fwd_stops=[[],[],[]]
    _rev_stops=[[],[],[]]
    if _stopCod and _seq_sub:
        _u=_seq_sub.upper()
        _rc_u=_rc(_seq_sub)
        for _i in range(len(_u)-2):
            if _u[_i:_i+3] in _STOP:
                _gfr=(_geno_off+_i)%3
                _fwd_stops[_gfr].append(_seq_sub_off+_i)  # viewport-relative x
        for _i in range(len(_rc_u)-2):
            if _rc_u[_i:_i+3] in _STOP:
                _gend=_geno_off+len(_seq_sub)-_i-1
                _gfr=_gend%3
                _rev_stops[_gfr].append(_seq_sub_off+len(_seq_sub)-_i-2)
    # 6 rows: +1 +2 +3  -1 -2 -3  (strand_int, absolute_genomic_frame 0/1/2)
    _ROWS=[('+1',1,0),('+2',1,1),('+3',1,2),('-1',-1,0),('-2',-1,1),('-3',-1,2)]
    # 7 gridspec rows: 3 fwd + thin separator + 3 rev
    fig=plt.figure(figsize=(_width,6*0.38+0.6))
    fig.patch.set_facecolor(_bg)
    gs=gridspec.GridSpec(7,1,figure=fig,hspace=0.0,
        height_ratios=[1,1,1,0.25,1,1,1],
        top=0.96,bottom=0.06,left=0.07,right=0.98)
    _axes=[]
    for _ri,(_lbl,_si,_fr) in enumerate(_ROWS):
        _gsi = _ri if _ri<3 else _ri+1   # skip gridspec row 3 (separator)
        ax=fig.add_subplot(gs[_gsi])
        _gfs=[]
        for (s,e,st,lb,col,noncoding,frame) in _feats_raw:
            if st!=_si: continue
            _effective_fr=0 if noncoding else frame
            if _effective_fr!=_fr: continue
            _gfs.append(GraphicFeature(start=s,end=e,strand=_si,label=lb,color=col))
        rec=GraphicRecord(sequence_length=_seqLen,features=_gfs)
        rec.plot(ax=ax,with_ruler=(_ri==5 and _ruler),draw_line=True)
        ax.set_ylim(-0.55,0.55)   # keep baseline centred; half-height tracks
        ax.set_facecolor(_bg)
        ax.set_ylabel(_lbl,fontsize=7,rotation=0,labelpad=22,va='center',color=_fg)
        for sp in ax.spines.values(): sp.set_color(_fg)
        ax.tick_params(colors=_fg,labelsize=${fs-2})
        if _ri<5: ax.set_xticklabels([])
        if _stopCod and _seq_sub:
            _sc_list=_fwd_stops[_fr] if _si==1 else _rev_stops[_fr]
            for _x in _sc_list:
                ax.axvline(_x,color='black',alpha=0.7,linewidth=0.8)
        _axes.append(ax)
    # Draw strand-separator axis in the blank gridspec row between +3 and -1
    _ax_sep=fig.add_subplot(gs[3])
    _ax_sep.set_visible(False)
    _pos3 =_axes[2].get_position()
    _pos4 =_axes[3].get_position()
    _ymid =(_pos3.y0+_pos4.y1)/2
    from matplotlib.lines import Line2D
    fig.add_artist(Line2D([0.07,0.98],[_ymid,_ymid],transform=fig.transFigure,
        color=_fg,linewidth=0.8,linestyle='-'))
    buf=io.BytesIO()
    fig.savefig(buf,format='${fmt_}',dpi=${dpi},bbox_inches='tight',facecolor=_bg)
    plt.close(fig); buf.seek(0)
    _result='OK:'+base64.b64encode(buf.read()).decode()
except Exception:
    _result='ERR:'+traceback.format_exc()
`;
  }

  // Simple mode
  const pyFeatsSingle = feats.map(f=>{
    const lbl=showLbls&&!(hideLong&&!labelFitsInFeature(f))?f.name.replace(/\\/g,'\\\\').replace(/'/g,"\\'"):'';
    const col=f.noncoding?NC_COLOR:f.color;
    const s=Math.max(0,f.start-vs), ee=Math.min(seqLen,f.end-vs);
    const si=f.strand==='+'?1:-1;
    return `    GraphicFeature(start=${s},end=${ee},strand=${si},label=${lbl?`'${lbl}'`:'None'},color='${col}')`;
  }).join(',\n');
  return `
import io,base64,traceback,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
_result=None
try:
    from dna_features_viewer import GraphicRecord,GraphicFeature
    _f=[
${pyFeatsSingle}
    ]
    rec=GraphicRecord(sequence_length=${seqLen},features=_f)
    ax,_=rec.plot(figure_width=${width},with_ruler=${ruler?'True':'False'},draw_line=True)
    fig=ax.get_figure()
    fig.patch.set_facecolor(${bg})
    ax.set_facecolor(${bg})
    for sp in ax.spines.values(): sp.set_color(${fg})
    ax.tick_params(colors=${fg})
    for t in ax.texts: t.set_fontsize(${fs})
    buf=io.BytesIO()
    fig.savefig(buf,format='${fmt_}',dpi=${dpi},bbox_inches='tight',facecolor=${bg})
    plt.close(fig); buf.seek(0)
    _result='OK:'+base64.b64encode(buf.read()).decode()
except Exception:
    _result='ERR:'+traceback.format_exc()
`;
}

// ── Render ────────────────────────────────────────────────────────────────────
async function triggerRender(forPng) {
  if (!pyReady||!lastState||pyRendering) return;
  pyRendering=true;
  document.getElementById('btn-render').disabled=true;
  setFigStatus(forPng?'rendering PNG\u2026':'rendering\u2026');
  // Pass sequence window to Python for stop codon computation
  const seq      = lastState.sequence || '';
  const seqStart = lastState.seq_start || 0;
  const relVS    = Math.max(0, localVS - seqStart);
  const relVE    = Math.min(seq.length, localVE - seqStart);
  const subseq   = relVS < relVE ? seq.slice(relVS, relVE) : '';
  const subOff   = Math.max(0, seqStart - localVS);
  pyodide.globals.set('_seq_sub', subseq);
  pyodide.globals.set('_seq_sub_off', subOff);
  try {
    await pyodide.runPythonAsync(buildPyCode(localVS,localVE,forPng));
    const result=pyodide.globals.get('_result');
    if (!result||typeof result!=='string') {
      setFig('<div class="fig-ph">Python returned null — check console.</div>');
      setFigStatus('error');
    } else if (result.startsWith('ERR:')) {
      setFig(`<pre class="fig-err">${esc(result.slice(4))}</pre>`);
      setFigStatus('Python error \u2014 see figure panel');
    } else if (result.startsWith('OK:')) {
      const b64=result.slice(3);
      if (forPng) {
        const blob=await (await fetch('data:image/png;base64,'+b64)).blob();
        dlBlob(blob,'genetui_figure.png'); setFigStatus('PNG downloaded');
      } else {
        lastSvgData=atob(b64);
        setFig(`<img src="data:image/svg+xml;base64,${b64}" style="max-width:100%">`);
        document.getElementById('btn-svg').disabled=false;
        document.getElementById('btn-png').disabled=false;
        setFigStatus('ready');
      }
    }
  } catch(err) {
    setFig(`<pre class="fig-err">${esc(String(err))}</pre>`);
    setFigStatus('JS error');
  }
  document.getElementById('btn-render').disabled=false;
  pyRendering=false;
}

function setFig(html) { document.getElementById('fig-output').innerHTML=html; }
function dlBlob(blob,name) {
  const url=URL.createObjectURL(blob);
  Object.assign(document.createElement('a'),{href:url,download:name}).click();
  URL.revokeObjectURL(url);
}

document.getElementById('btn-render').addEventListener('click',()=>triggerRender(false));
document.getElementById('btn-svg').addEventListener('click',()=>{
  if (lastSvgData) dlBlob(new Blob([lastSvgData],{type:'image/svg+xml'}),'genetui_figure.svg');
});
document.getElementById('btn-png').addEventListener('click',()=>triggerRender(true));

// ── 3Dmol.js structure viewer ─────────────────────────────────────────────────
let mol3dViewer=null, mol3dLoaded=false, mol3dStyle='plddt', lastPdbName='';
let pendingPdb=null, pendingPdbName='';

// Horizontal drag resize between struct-panel and fig-panel
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
// AlphaFold-like pLDDT gradient: orange(low)→yellow→cyan→blue(high)
// 3Dmol built-in 'roygb' goes red→yellow→green→blue; we reverse it (min>max) to
// get blue(high)→green→yellow→red(low), then shift hue by using rwb as a fallback.
// Best approximation available without custom gradient API: use 'roygb' with min/max
// swapped so high b-factor = blue end.
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
</script>
</body>
</html>"#####;
