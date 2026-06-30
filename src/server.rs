use std::sync::{Arc, Mutex};
use std::collections::HashMap;
use std::f64::consts::PI;
use tokio::sync::broadcast;
use axum::{
    Router,
    extract::{
        ws::{Message, WebSocket, WebSocketUpgrade},
        State, Query,
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
    pub sequence: String,        // DNA bases for seq_start..(seq_start+sequence.len()) (may be empty)
    pub seq_start: u64,         // genomic offset of sequence[0]
    pub protein_pdb: String,    // raw PDB text of selected+folded protein (may be empty)
    pub protein_name: String,   // gene name for the PDB
    pub coverage_plus:      Vec<u32>,   // per-bin read depth, plus strand
    pub coverage_minus:     Vec<u32>,   // per-bin read depth, minus strand
    pub coverage_bin_size:  u64,        // bp per bin
    pub coverage_bin_start: u64,        // genomic position of first bin
    pub input_files: Vec<String>,       // source file paths (gff, fasta, bam)
    pub msa_gene: String,               // gene name for the active MSA (empty = none)
    pub msa_sequences: Vec<[String; 2]>,// [(id, aligned_sequence), …]
    pub status_msg: String,             // progress message for slow tasks (empty = idle)
    pub gc_skew: Vec<f32>,              // downsampled GC skew values (-1..1) for circular map (legacy, kept for compatibility)
    pub active_genome: usize,
    pub main_name: String,
    pub main_size: u64,
    pub main_features: Vec<BrowserFeature>,
    pub main_gc_skew: Vec<f32>,
    pub plasmid_names: Vec<String>,
    pub plasmid_sizes: Vec<u64>,
    pub plasmid_features: Vec<Vec<BrowserFeature>>,
    pub plasmid_gc_skew: Vec<Vec<f32>>,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct BrowserFeature {
    pub name:      String,
    pub locus_tag: String,
    pub start:     u64,
    pub end:       u64,
    pub strand:    char,
    pub color:     String,
    pub noncoding: bool,
    pub is_orf:    bool,
}

pub struct WebServer {
    state:              Arc<Mutex<BrowserState>>,
    tx:                 broadcast::Sender<String>,
    viewport_cmd:       Arc<Mutex<Option<(u64, u64)>>>,
    fold_gene_cmd:      Arc<Mutex<Option<String>>>,
    msa_gene_cmd:       Arc<Mutex<Option<String>>>,
    switch_genome_cmd:  Arc<Mutex<Option<(usize, u64, u64)>>>,
}

impl WebServer {
    pub fn new() -> Arc<Self> {
        let (tx, _) = broadcast::channel(32);
        Arc::new(Self {
            state:             Arc::new(Mutex::new(BrowserState::default())),
            tx,
            viewport_cmd:      Arc::new(Mutex::new(None)),
            fold_gene_cmd:     Arc::new(Mutex::new(None)),
            msa_gene_cmd:      Arc::new(Mutex::new(None)),
            switch_genome_cmd: Arc::new(Mutex::new(None)),
        })
    }

    pub fn take_viewport_cmd(&self) -> Option<(u64, u64)> {
        self.viewport_cmd.lock().ok()?.take()
    }

    pub fn take_fold_cmd(&self) -> Option<String> {
        self.fold_gene_cmd.lock().ok()?.take()
    }

    pub fn take_msa_cmd(&self) -> Option<String> {
        self.msa_gene_cmd.lock().ok()?.take()
    }

    pub fn take_switch_cmd(&self) -> Option<(usize, u64, u64)> {
        self.switch_genome_cmd.lock().ok()?.take()
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
            .route("/static/genome_track.js", get(|| async { ([("content-type", "application/javascript")], JS_GENOME_TRACK) }))
            .route("/static/coverage.js",     get(|| async { ([("content-type", "application/javascript")], JS_COVERAGE) }))
            .route("/static/navigation.js",   get(|| async { ([("content-type", "application/javascript")], JS_NAVIGATION) }))
            .route("/static/panels.js",       get(|| async { ([("content-type", "application/javascript")], JS_PANELS) }))
            .route("/static/gene_plot.js",    get(|| async { ([("content-type", "application/javascript")], JS_GENE_PLOT) }))
            .route("/static/structure.js",    get(|| async { ([("content-type", "application/javascript")], JS_STRUCTURE) }))
            .route("/static/websocket.js",    get(|| async { ([("content-type", "application/javascript")], JS_WEBSOCKET) }))
            .route("/static/gene_info.js",    get(|| async { ([("content-type", "application/javascript")], JS_GENE_INFO) }))
            .route("/static/msa_panel.js",      get(|| async { ([("content-type", "application/javascript")], JS_MSA_PANEL) }))
            .route("/static/circular_plot.js", get(|| async { ([("content-type", "application/javascript")], JS_CIRCULAR_PLOT) }))
            .route("/circular-map.svg",        get(circular_map_svg))
            .with_state(self);

        match tokio::net::TcpListener::bind(format!("127.0.0.1:{}", port)).await {
            Ok(listener) => { let _ = axum::serve(listener, app).await; }
            Err(e) => eprintln!("genetui web server: could not bind port {port}: {e}"),
        }
    }
}

async fn serve_html() -> Html<&'static str> { Html(HTML) }

async fn circular_map_svg(
    State(srv): State<Arc<WebServer>>,
    Query(params): Query<HashMap<String, String>>,
) -> impl IntoResponse {
    let dark        = params.get("dark").map(|v| v != "0").unwrap_or(true);
    let show_nc     = params.get("nc").map(|v| v != "0").unwrap_or(false);
    let show_legend = params.get("legend").map(|v| v != "0").unwrap_or(true);
    let title       = params.get("title").filter(|s| !s.is_empty()).map(|s| s.as_str());
    let genome_idx  = params.get("genome").and_then(|s| s.parse::<usize>().ok()).unwrap_or(0);
    let state = match srv.state.lock() {
        Ok(s) => s,
        Err(_) => return ([("content-type", "image/svg+xml")], "<svg/>".to_string()),
    };

    let (name, size, features, gc_skew) = if genome_idx == 0 {
        let n = if state.main_name.is_empty() { &state.genome_name } else { &state.main_name };
        let s = if state.main_size > 0 { state.main_size } else { state.genome_size };
        let f = if !state.main_features.is_empty() { &state.main_features } else { &state.features };
        let g = if !state.main_gc_skew.is_empty() { &state.main_gc_skew } else { &state.gc_skew };
        (n.clone(), s, f.clone(), g.clone())
    } else {
        let pi = genome_idx - 1;
        if pi >= state.plasmid_names.len() {
            return ([("content-type", "image/svg+xml")], "<svg/>".to_string());
        }
        let pfeats = state.plasmid_features.get(pi).cloned().unwrap_or_default();
        let mut sz = state.plasmid_sizes.get(pi).copied().unwrap_or(0);
        if sz == 0 {
            sz = pfeats.iter().map(|f| f.end).max().unwrap_or(0);
        }
        (
            state.plasmid_names[pi].clone(),
            sz,
            pfeats,
            state.plasmid_gc_skew.get(pi).cloned().unwrap_or_default(),
        )
    };

    let svg = build_circular_svg(&name, size, &features, &gc_skew, dark, show_nc, show_legend, title, genome_idx);
    ([("content-type", "image/svg+xml")], svg)
}

fn svg_esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

fn lerp_color(t: f64, r0: u8, g0: u8, b0: u8, r1: u8, g1: u8, b1: u8) -> String {
    let t = t.clamp(0.0, 1.0);
    let r = (r0 as f64 + t * (r1 as f64 - r0 as f64)).round() as u8;
    let g = (g0 as f64 + t * (g1 as f64 - g0 as f64)).round() as u8;
    let b = (b0 as f64 + t * (b1 as f64 - b0 as f64)).round() as u8;
    format!("#{:02x}{:02x}{:02x}", r, g, b)
}

fn build_circular_svg(
    name: &str,
    size: u64,
    features: &[BrowserFeature],
    gc_skew: &[f32],
    dark: bool, show_nc: bool, show_legend: bool, custom_title: Option<&str>,
    id_sfx: usize,
) -> String {
    if size == 0 { return "<svg/>".to_string(); }

    let svg_w  = 1000.0f64;
    let svg_h  = 800.0f64;
    let cx     = 400.0f64;   // circle center X
    let cy     = 400.0f64;   // circle center Y
    let gsize  = size as f64;

    let bg         = if dark { "#0d1117" } else { "#ffffff" };
    let fg         = if dark { "#c9d1d9" } else { "#333333" };
    let ring_color = if dark { "#30363d" } else { "#cccccc" };

    // Track radii as fractions of cx (working inward)
    let r_lbl_f      = 0.97f64;
    let r_tick_o_f   = 0.94f64;
    let r_plus_o_f   = 0.93f64;  // +strand density outer
    let r_plus_i_f   = 0.81f64;  // +strand density inner  (12%)
    let r_minus_o_f  = 0.80f64;  // -strand density outer
    let r_minus_i_f  = 0.68f64;  // -strand density inner  (12%)
    let r_nc_o_f     = 0.67f64;
    let r_nc_i_f     = 0.59f64;
    let (r_skew_o_f, r_skew_i_f) = if show_nc { (0.58f64, 0.49f64) } else { (0.66f64, 0.57f64) };
    let r_inner_f    = r_skew_i_f;

    // Absolute pixel values in viewBox units for JS
    let r_outer_abs = cx * r_plus_o_f;
    let r_inner_abs = cx * r_inner_f;

    let rf = |f: f64| cx * f;

    let pos_to_angle = |pos: u64| -> f64 { -PI / 2.0 + 2.0 * PI * pos as f64 / gsize };

    let arc_path = |start: u64, end: u64, r_in: f64, r_out: f64| -> String {
        let a1   = pos_to_angle(start);
        let a2   = pos_to_angle(end);
        let span = 2.0 * PI * (end - start) as f64 / gsize;
        let laf  = if span > PI { 1 } else { 0 };
        let (x1o,y1o) = (cx + r_out*a1.cos(), cy + r_out*a1.sin());
        let (x2o,y2o) = (cx + r_out*a2.cos(), cy + r_out*a2.sin());
        let (x1i,y1i) = (cx + r_in *a1.cos(), cy + r_in *a1.sin());
        let (x2i,y2i) = (cx + r_in *a2.cos(), cy + r_in *a2.sin());
        format!("M {x1o:.2} {y1o:.2} A {r_out:.2} {r_out:.2} 0 {laf} 1 {x2o:.2} {y2o:.2} \
                 L {x2i:.2} {y2i:.2} A {r_in:.2} {r_in:.2} 0 {laf} 0 {x1i:.2} {y1i:.2} Z")
    };

    // ── Sliding-window density (per strand) ───────────────────────────────────
    let win  = 100_000u64.min(size / 4 + 1);
    let step = (win / 10).max(1);
    let n_w  = ((gsize / step as f64).ceil() as usize).max(1);
    let mut plus_c  = vec![0u64; n_w];
    let mut minus_c = vec![0u64; n_w];
    for f in features {
        if f.noncoding || f.end <= f.start { continue; }
        let first_w = if f.start >= win { ((f.start - win) / step + 1) as usize } else { 0 };
        let last_w  = (f.end / step) as usize;
        let counts  = if f.strand == '+' { &mut plus_c } else { &mut minus_c };
        for wi in first_w..=last_w.min(n_w.saturating_sub(1)) {
            let ws = wi as u64 * step;
            let os = f.start.max(ws);
            let oe = f.end.min(ws + win);
            if oe > os { counts[wi] += oe - os; }
        }
    }
    // Per-strand min/max for full-range normalisation
    let plus_min  = *plus_c.iter().min().unwrap_or(&0) as f64;
    let plus_max  = (*plus_c.iter().max().unwrap_or(&1) as f64).max(plus_min + 1.0);
    let minus_min = *minus_c.iter().min().unwrap_or(&0) as f64;
    let minus_max = (*minus_c.iter().max().unwrap_or(&1) as f64).max(minus_min + 1.0);

    // Strand density colour ramps
    let (p_r0,p_g0,p_b0, p_r1,p_g1,p_b1) = if dark {
        (0x0du8,0x11u8,0x17u8, 0x39u8,0xd3u8,0x53u8)  // dark bg → green (+ strand)
    } else {
        (0xf0u8,0xf5u8,0xf0u8, 0x00u8,0x7au8,0x33u8)
    };
    let (m_r0,m_g0,m_b0, m_r1,m_g1,m_b1) = if dark {
        (0x0du8,0x11u8,0x17u8, 0xf9u8,0x73u8,0x16u8)  // dark bg → orange (- strand)
    } else {
        (0xf5u8,0xf0u8,0xeeu8, 0xc0u8,0x50u8,0x00u8)
    };

    // ── GC skew min/max ───────────────────────────────────────────────────────
    let (skew_min, skew_max) = if gc_skew.is_empty() {
        (-1.0f64, 1.0f64)
    } else {
        let mn = gc_skew.iter().cloned().fold(f32::INFINITY,  f32::min) as f64;
        let mx = gc_skew.iter().cloned().fold(f32::NEG_INFINITY, f32::max) as f64;
        let pad = ((mx - mn) * 0.05).max(1e-6);
        (mn - pad, mx + pad)
    };
    let skew_col_lo  = if dark { "#61afef" } else { "#2780b9" };
    let skew_col_mid = if dark { "#21262d" } else { "#f0f0f0" };
    let skew_col_hi  = if dark { "#e06c75" } else { "#c0392b" };
    // Per-SVG unique IDs so multiple inline SVGs don't clash
    let gp = format!("lg-plus-{id_sfx}");
    let gm = format!("lg-minus-{id_sfx}");
    let gs = format!("lg-skew-{id_sfx}");
    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {svg_w} {svg_h}\" \
         data-cx=\"{cx}\" data-cy=\"{cy}\" data-r-outer=\"{r_outer_abs:.0}\" data-r-inner=\"{r_inner_abs:.0}\">\n\
<defs>\n\
  <linearGradient id=\"{gp}\" x1=\"0\" y1=\"0\" x2=\"1\" y2=\"0\">\
    <stop offset=\"0%\" stop-color=\"#{p_r0:02x}{p_g0:02x}{p_b0:02x}\"/>\
    <stop offset=\"100%\" stop-color=\"#{p_r1:02x}{p_g1:02x}{p_b1:02x}\"/>\
  </linearGradient>\n\
  <linearGradient id=\"{gm}\" x1=\"0\" y1=\"0\" x2=\"1\" y2=\"0\">\
    <stop offset=\"0%\" stop-color=\"#{m_r0:02x}{m_g0:02x}{m_b0:02x}\"/>\
    <stop offset=\"100%\" stop-color=\"#{m_r1:02x}{m_g1:02x}{m_b1:02x}\"/>\
  </linearGradient>\n\
  <linearGradient id=\"{gs}\" x1=\"0\" y1=\"0\" x2=\"1\" y2=\"0\">\
    <stop offset=\"0%\" stop-color=\"{skew_col_lo}\"/>\
    <stop offset=\"50%\" stop-color=\"{skew_col_mid}\"/>\
    <stop offset=\"100%\" stop-color=\"{skew_col_hi}\"/>\
  </linearGradient>\n\
</defs>\n\
<rect width=\"{svg_w}\" height=\"{svg_h}\" fill=\"{bg}\"/>\n",
    );

    // ── Track ring outlines ────────────────────────────────────────────────────
    for &f in &[r_plus_o_f, r_plus_i_f, r_minus_o_f, r_minus_i_f] {
        let r = rf(f);
        svg.push_str(&format!(
            r#"<circle cx="{cx:.1}" cy="{cy:.1}" r="{r:.1}" fill="none" stroke="{ring_color}" stroke-width="1.0"/>
"#));
    }
    if show_nc {
        for &f in &[r_nc_o_f, r_nc_i_f] {
            let r = rf(f);
            svg.push_str(&format!(
                r#"<circle cx="{cx:.1}" cy="{cy:.1}" r="{r:.1}" fill="none" stroke="{ring_color}" stroke-width="1.0"/>
"#));
        }
    }
    for &f in &[r_skew_o_f, r_skew_i_f] {
        let r = rf(f);
        svg.push_str(&format!(
            r#"<circle cx="{cx:.1}" cy="{cy:.1}" r="{r:.1}" fill="none" stroke="{ring_color}" stroke-width="1.0"/>
"#));
    }

    // ── +strand density heatmap ────────────────────────────────────────────────
    for wi in 0..n_w {
        let ws = wi as u64 * step;
        if ws >= size { break; }
        let we    = (ws + step).min(size);
        let t     = (plus_c[wi] as f64 - plus_min) / (plus_max - plus_min);
        let color = lerp_color(t, p_r0,p_g0,p_b0, p_r1,p_g1,p_b1);
        let path  = arc_path(ws, we, rf(r_plus_i_f), rf(r_plus_o_f));
        svg.push_str(&format!(r#"<path d="{path}" fill="{color}" stroke="none"/>
"#));
    }

    // ── -strand density heatmap ────────────────────────────────────────────────
    for wi in 0..n_w {
        let ws = wi as u64 * step;
        if ws >= size { break; }
        let we    = (ws + step).min(size);
        let t     = (minus_c[wi] as f64 - minus_min) / (minus_max - minus_min);
        let color = lerp_color(t, m_r0,m_g0,m_b0, m_r1,m_g1,m_b1);
        let path  = arc_path(ws, we, rf(r_minus_i_f), rf(r_minus_o_f));
        svg.push_str(&format!(r#"<path d="{path}" fill="{color}" stroke="none"/>
"#));
    }

    // ── NC feature arcs (optional) ────────────────────────────────────────────
    if show_nc {
        for f in features {
            if !f.noncoding || f.end <= f.start { continue; }
            let path = arc_path(f.start, f.end, rf(r_nc_i_f), rf(r_nc_o_f));
            svg.push_str(&format!(r#"<path d="{path}" fill="{}" stroke="none"/>
"#, f.color));
        }
    }

    // ── GC skew heatmap (innermost band) ──────────────────────────────────────
    if !gc_skew.is_empty() {
        let n = gc_skew.len();
        for i in 0..n {
            let ws = (i as u64 * size) / n as u64;
            let we = ((i + 1) as u64 * size) / n as u64;
            if we <= ws { continue; }
            let v = gc_skew[i] as f64;
            // Normalise within actual data range; mid-point = 0
            let color = if v >= 0.0 {
                let t = (v / skew_max.max(1e-9)).clamp(0.0, 1.0);
                if dark { lerp_color(t, 0x21,0x26,0x2d, 0xe0,0x6c,0x75) }
                else    { lerp_color(t, 0xf0,0xf0,0xf0, 0xc0,0x39,0x2b) }
            } else {
                let t = (v / skew_min.min(-1e-9)).clamp(0.0, 1.0);
                if dark { lerp_color(t, 0x21,0x26,0x2d, 0x61,0xaf,0xef) }
                else    { lerp_color(t, 0xf0,0xf0,0xf0, 0x27,0x80,0xb9) }
            };
            let path = arc_path(ws, we, rf(r_skew_i_f), rf(r_skew_o_f));
            svg.push_str(&format!(r#"<path d="{path}" fill="{color}" stroke="none"/>
"#));
        }
    }

    // ── Tick marks + labels ────────────────────────────────────────────────────
    let raw      = size / 8;
    let mag      = 10u64.pow(raw.max(1).ilog10());
    let interval = (raw / mag) * mag;
    let mut pos  = 0u64;
    while pos < size {
        let a = pos_to_angle(pos);
        let (ca, sa) = (a.cos(), a.sin());
        svg.push_str(&format!(
            r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="{ring_color}" stroke-width="0.8"/>
"#,
            cx + rf(r_inner_f)*ca, cy + rf(r_inner_f)*sa,
            cx + rf(r_tick_o_f)*ca, cy + rf(r_tick_o_f)*sa,
        ));
        let label = if pos == 0 { "0".into() }
            else if pos >= 1_000_000 { format!("{:.1}M", pos as f64 / 1e6) }
            else { format!("{}k", pos / 1_000) };
        svg.push_str(&format!(
            r#"<text x="{:.2}" y="{:.2}" text-anchor="middle" dominant-baseline="middle" font-size="11" font-family="sans-serif" fill="{fg}">{label}</text>
"#,
            cx + rf(r_lbl_f)*ca, cy + rf(r_lbl_f)*sa,
        ));
        pos += interval;
    }

    // ── Title in centre (word-wrapped) ────────────────────────────────────────
    let title_str = custom_title.unwrap_or(name);
    let words: Vec<&str> = title_str.split_whitespace().collect();
    let mut lines: Vec<String> = Vec::new();
    let mut cur = String::new();
    for w in &words {
        if cur.is_empty() { cur.push_str(w); }
        else if cur.len() + 1 + w.len() <= 18 { cur.push(' '); cur.push_str(w); }
        else { lines.push(cur.clone()); cur = w.to_string(); }
    }
    if !cur.is_empty() { lines.push(cur); }
    if lines.is_empty() { lines.push(title_str.to_string()); }
    lines.truncate(4);

    let line_h  = 15.0f64;
    let n_lines = lines.len() as f64;
    let y_start = cy - (n_lines - 1.0) * line_h / 2.0;
    svg.push_str(&format!(
        r#"<text id="circ-title-svg" font-family="sans-serif" fill="{fg}" font-weight="bold" font-size="12" text-anchor="middle" style="cursor:pointer">"#
    ));
    for (i, line) in lines.iter().enumerate() {
        svg.push_str(&format!(
            r#"<tspan x="{cx:.1}" y="{:.1}">{}</tspan>"#,
            y_start + i as f64 * line_h, svg_esc(line)
        ));
    }
    svg.push_str("</text>\n");

    let mb = size as f64 / 1e6;
    svg.push_str(&format!(
        r#"<text x="{cx:.1}" y="{:.1}" text-anchor="middle" dominant-baseline="middle" font-size="10" font-family="sans-serif" fill="{ring_color}">{mb:.2} Mb</text>
"#,
        y_start + n_lines * line_h,
    ));

    // ── Colour-bar legend (right strip, x=808..1000) ──────────────────────────
    if show_legend {
        let bw   = 70.0f64;      // bar width
        let bh   = 8.0f64;       // bar height
        let gap  = 55.0f64;      // row pitch
        let lx   = 820.0f64;     // left edge of bars in right strip
        let _lbl_x = lx - 4.0;  // label right-edge (unused; label goes above bar)

        // Separator background for right strip
        svg.push_str(&format!(
            r#"<rect x="808" y="0" width="192" height="{svg_h}" fill="{bg}"/>
"#));

        // Rows: +strand, -strand, GC skew, [nc if shown]
        struct LegRow { label: &'static str, grad: String, lo: String, hi: String }
        let rows = {
            let mut v = vec![
                LegRow {
                    label: "+ strand",
                    grad:  gp.clone(),
                    lo:    format!("{:.0}%", plus_min / win as f64 * 100.0),
                    hi:    format!("{:.0}%", plus_max / win as f64 * 100.0),
                },
                LegRow {
                    label: "- strand",
                    grad:  gm.clone(),
                    lo:    format!("{:.0}%", minus_min / win as f64 * 100.0),
                    hi:    format!("{:.0}%", minus_max / win as f64 * 100.0),
                },
                LegRow {
                    label: "GC skew",
                    grad:  gs.clone(),
                    lo:    format!("{:.3}", skew_min),
                    hi:    format!("{:.3}", skew_max),
                },
            ];
            if show_nc { v.push(LegRow { label: "NC feat.", grad: gp.clone(), lo: String::new(), hi: String::new() }); }
            v
        };

        // Legend title
        svg.push_str(&format!(
            r#"<text x="{:.1}" y="35" text-anchor="middle" font-size="9" font-family="sans-serif" fill="{fg}" font-weight="bold">Legend</text>
"#,
            lx + bw / 2.0
        ));

        for (ri, row) in rows.iter().enumerate() {
            let ry = 60.0 + ri as f64 * gap;
            // row label above bar
            svg.push_str(&format!(
                r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" dominant-baseline="auto" font-size="9" font-family="sans-serif" fill="{fg}">{}</text>
"#,
                lx + bw / 2.0, ry - 3.0, row.label
            ));
            if row.label == "NC feat." {
                // Solid swatch
                let nc_col = if dark { "#8b949e" } else { "#6e7681" };
                svg.push_str(&format!(
                    r#"<rect x="{lx:.1}" y="{ry:.1}" width="{bw:.1}" height="{bh:.1}" rx="1" fill="{nc_col}"/>
"#));
            } else {
                // Gradient bar
                svg.push_str(&format!(
                    r#"<rect x="{lx:.1}" y="{ry:.1}" width="{bw:.1}" height="{bh:.1}" rx="1" fill="url(#{})"/>
"#, row.grad));
                // border
                svg.push_str(&format!(
                    r#"<rect x="{lx:.1}" y="{ry:.1}" width="{bw:.1}" height="{bh:.1}" rx="1" fill="none" stroke="{ring_color}" stroke-width="0.4"/>
"#));
                // lo/hi values
                if !row.lo.is_empty() {
                    svg.push_str(&format!(
                        r#"<text x="{lx:.1}" y="{:.1}" text-anchor="start" dominant-baseline="auto" font-size="7" font-family="monospace" fill="{fg}">{}</text>
"#,
                        ry + bh + 3.0, svg_esc(&row.lo)
                    ));
                    svg.push_str(&format!(
                        r#"<text x="{:.1}" y="{:.1}" text-anchor="end" dominant-baseline="auto" font-size="7" font-family="monospace" fill="{fg}">{}</text>
"#,
                        lx + bw, ry + bh + 3.0, svg_esc(&row.hi)
                    ));
                }
            }
        }
    }

    svg.push_str("</svg>");
    svg
}

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
                            if v["cmd"].as_str() == Some("fold") {
                                if let Some(gene) = v["gene"].as_str() {
                                    if let Ok(mut cmd) = server.fold_gene_cmd.lock() {
                                        *cmd = Some(gene.to_string());
                                    }
                                }
                            }
                            if v["cmd"].as_str() == Some("msa") {
                                if let Some(gene) = v["gene"].as_str() {
                                    if let Ok(mut cmd) = server.msa_gene_cmd.lock() {
                                        *cmd = Some(gene.to_string());
                                    }
                                }
                            }
                            if v["cmd"].as_str() == Some("switch_genome") {
                                if let (Some(g), Some(s), Some(e)) = (v["genome"].as_u64(), v["start"].as_u64(), v["end"].as_u64()) {
                                    if let Ok(mut cmd) = server.switch_genome_cmd.lock() {
                                        *cmd = Some((g as usize, s, e));
                                    }
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

// ── Static JS assets (embedded at compile time) ───────────────────────────────
const JS_GENOME_TRACK: &str = include_str!("web/genome_track.js");
const JS_COVERAGE:     &str = include_str!("web/coverage.js");
const JS_NAVIGATION:   &str = include_str!("web/navigation.js");
const JS_PANELS:       &str = include_str!("web/panels.js");
const JS_GENE_PLOT:    &str = include_str!("web/gene_plot.js");
const JS_STRUCTURE:    &str = include_str!("web/structure.js");
const JS_WEBSOCKET:    &str = include_str!("web/websocket.js");
const JS_GENE_INFO:    &str = include_str!("web/gene_info.js");
const JS_MSA_PANEL:    &str = include_str!("web/msa_panel.js");
const JS_CIRCULAR_PLOT: &str = include_str!("web/circular_plot.js");

// ── Embedded HTML shell ───────────────────────────────────────────────────────
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

/* ── Gene info bar ───────────────────────────────────────────────────────── */
#gene-info-bar {
  height: 28px; padding: 0 12px; background: #0d1117; flex-shrink: 0;
  border-top: 1px solid #21262d; border-bottom: 1px solid #21262d;
  display: flex; align-items: center; gap: 14px; font-size: 11px;
  font-family: monospace; overflow: hidden; transition: opacity .15s;
}
#gene-info-bar.empty { opacity: 0.3; }
#gene-info-bar .gi-name  { color: #58a6ff; font-weight: 600; }
#gene-info-bar .gi-locus { color: #8b949e; }
#gene-info-bar .gi-coords{ color: #c9d1d9; }
#gene-info-bar .gi-strand{ color: #e3b341; }
#gene-info-bar .gi-len   { color: #7ee787; }
#gene-info-bar .gi-kind  { color: #f0883e; font-size: 9px; padding: 1px 5px;
  border: 1px solid #30363d; border-radius: 3px; }
#gene-info-bar .gi-hint   { color: #484f58; margin-left: auto; font-size: 10px; }
#gene-info-bar.status    { opacity: 1; }
#gene-info-bar .gi-status{ color: #e3b341; font-style: italic; }

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
  display: flex; justify-content: center; align-items: center;
}
#fig-output img { max-width: 100%; border-radius: 4px; }
.fig-ph { color: #484f58; font-size: 12px; text-align: center; margin-top: 20px; line-height: 1.8; }
.fig-err { color: #f85149; font-size: 11px; font-family: monospace; white-space: pre-wrap; word-break: break-all; padding: 8px; }

/* ── MSA panel ───────────────────────────────────────────────────────────── */
#msa-panel { width: 42%; min-width: 200px; display: none; flex-direction: column; overflow: hidden; flex-shrink: 0; }
#msa-panel.visible { display: flex; }
#msa-h-handle {
  width: 5px; background: #21262d; cursor: ew-resize; flex-shrink: 0;
  border-left: 1px solid #30363d; border-right: 1px solid #30363d; display: none;
}
#msa-h-handle.visible { display: block; }
#msa-h-handle:hover { background: #58a6ff; }
#msa-bar {
  height: 30px; padding: 0 10px; background: #161b22; border-bottom: 1px solid #30363d;
  display: flex; align-items: center; gap: 8px; flex-shrink: 0;
}
#msa-bar-label { font-size: 11px; color: #8b949e; }
#msa-gene-name { font-size: 11px; color: #58a6ff; font-family: monospace; }
#msa-status { font-size: 10px; color: #6e7681; font-family: monospace; margin-left: auto; }
#msa-canvas-wrap { flex: 1; display: flex; flex-direction: column; min-height: 0; overflow: hidden; position: relative; }
#msa-canvas { display: block; cursor: default; image-rendering: pixelated; }
#msa-names-canvas { position: absolute; left: 0; top: 0; pointer-events: none; }

/* ── Circular genome map panel ───────────────────────────────────────────── */
#circ-panel { width: 45%; min-width: 220px; display: none; flex-direction: column; overflow: hidden; flex-shrink: 0; }
#circ-panel.visible { display: flex; }
#circ-h-handle {
  width: 5px; background: #21262d; cursor: ew-resize; flex-shrink: 0;
  border-left: 1px solid #30363d; border-right: 1px solid #30363d; display: none;
}
#circ-h-handle.visible { display: block; }
#circ-h-handle:hover { background: #58a6ff; }
#circ-bar {
  height: 30px; padding: 0 10px; background: #161b22; border-bottom: 1px solid #30363d;
  display: flex; align-items: center; gap: 8px; flex-shrink: 0;
}
#circ-bar-label { font-size: 11px; color: #8b949e; }
#circ-status { font-size: 10px; color: #6e7681; font-family: monospace; }
#circ-body {
  flex: 1; overflow-y: auto; overflow-x: hidden; background: #0d1117;
}
.circ-map-item { width: 100%; border-bottom: 1px solid #21262d; }
.circ-map-item svg { width: 100%; height: auto; display: block; }
#circ-opts-box {
  position: fixed; z-index: 300;
  background: #161b22; border: 1px solid #30363d; border-radius: 8px;
  padding: 12px 14px; width: 230px; box-shadow: 0 8px 32px rgba(0,0,0,0.6);
  display: none;
}

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

<!-- Gene info bar -->
<div id="gene-info-bar" class="empty">
  <span class="gi-name">—</span>
  <span class="gi-hint">hover or click a gene in the track above</span>
</div>

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

  <!-- Horizontal resize handle between structure and MSA -->
  <div id="struct-h-handle"></div>

  <!-- MSA panel -->
  <div id="msa-panel">
    <div id="msa-bar">
      <span id="msa-bar-label">MSA</span>
      <span id="msa-gene-name"></span>
      <span id="msa-status">click a gene then Run MSA, or press m in TUI</span>
      <div style="margin-left:auto">
        <button class="dbtn" id="btn-run-msa" disabled>Run MSA</button>
      </div>
    </div>
    <div id="msa-canvas-wrap">
      <canvas id="msa-canvas"></canvas>
    </div>
  </div>

  <!-- Horizontal resize handle between MSA and gene plot -->
  <div id="msa-h-handle"></div>

  <!-- Circular genome map panel -->
  <div id="circ-panel" class="visible">
    <div id="circ-bar">
      <span id="circ-bar-label">Circular map</span>
      <span id="circ-status">enable panel (d) to render</span>
      <div style="margin-left:auto;display:flex;gap:5px">
        <button class="dbtn" id="btn-circ-opts">&#9881; Options</button>
        <button class="dbtn" id="btn-circ-render">Render</button>
        <button class="dbtn" id="btn-circ-png">PNG</button>
      </div>
    </div>
    <div id="circ-body">
      <div class="fig-ph">Circular genome map — click Render to generate.</div>
    </div>
  </div>

  <!-- Circular map options dropdown -->
  <div id="circ-opts-box">
    <div style="font-size:11px;color:#8b949e;text-transform:uppercase;letter-spacing:.06em;margin-bottom:10px">Circular map options</div>
    <div class="cr" style="margin-top:4px">
      <label for="circ-title-input" style="color:#8b949e;min-width:36px;font-size:11px">Title</label>
      <input type="text" id="circ-title-input" placeholder="(genome name)" style="flex:1;background:#0d1117;color:#c9d1d9;border:1px solid #30363d;border-radius:4px;padding:2px 5px;font-size:11px">
    </div>
    <div class="cr" style="margin-top:6px"><input type="checkbox" id="circ-white"> <label for="circ-white">White background</label></div>
    <div class="cr" style="margin-top:6px"><input type="checkbox" id="circ-show-nc"> <label for="circ-show-nc">Show NC features</label></div>
    <div class="cr" style="margin-top:6px"><input type="checkbox" id="circ-show-legend" checked> <label for="circ-show-legend">Show colour legend</label></div>
    <div class="cr" style="margin-top:6px"><input type="checkbox" id="circ-show-viewport" checked> <label for="circ-show-viewport">Show viewport marker</label></div>
    <div class="cr" style="margin-top:6px"><input type="checkbox" id="circ-show-plasmids" checked> <label for="circ-show-plasmids">Show plasmids</label></div>
  </div>

  <!-- Horizontal resize handle between circular map and gene plot -->
  <div id="circ-h-handle" class="visible"></div>

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
        <button class="dbtn" id="btn-py"    disabled>Download code</button>
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
    <div class="cr" style="grid-column:1/-1">
      <label for="opt-gencode" style="min-width:90px">Genetic code</label>
      <select id="opt-gencode">
        <option value="1">Standard</option>
        <option value="11">Bacterial / Archaeal / Plastid</option>
        <option value="2">Vertebrate Mitochondrial</option>
        <option value="4">Mycoplasma / Spiroplasma</option>
      </select>
    </div>
  </div>
  <hr class="osep">
  <div class="og">
    <div class="cr">
      <label for="cov-style-fig">Coverage (figure)</label>
      <select id="cov-style-fig">
        <option value="none">None</option>
        <option value="histogram">Histogram</option>
        <option value="kernel">Kernel (smoothed)</option>
        <option value="reads">Raw reads</option>
      </select>
    </div>
  </div>
</div>

<script src="/static/genome_track.js"></script>
<script src="/static/coverage.js"></script>
<script src="/static/navigation.js"></script>
<script src="/static/panels.js"></script>
<script src="/static/gene_plot.js"></script>
<script src="/static/structure.js"></script>
<script src="/static/websocket.js"></script>
<script src="/static/gene_info.js"></script>
<script src="/static/msa_panel.js"></script>
<script src="/static/circular_plot.js"></script>
</body>
</html>"#####;
