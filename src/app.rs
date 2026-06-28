use std::collections::HashMap;
use ratatui::layout::Rect;
use crate::core::Feature;

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum ActivePanel { Genome, Protein, Msa }

pub struct MsaPanel {
    pub gene_name:    String,
    pub sequences:    Vec<(String, String)>,  // (id, aligned_sequence)
    pub loading:      bool,
    pub error:        Option<String>,
    pub viewport_row: usize,
    pub viewport_col: usize,
}

pub struct ProteinPanel {
    pub gene_name: String,
    pub nt_seq: String,
    pub aa_seq: String,
    pub atoms: Vec<crate::protein::ProteinAtom>,
    pub camera: crate::pv::camera::Camera,
    pub drag_last: Option<(u16, u16)>,
    pub folding: bool,
    pub error: Option<String>,
    pub img_cache: Option<image::DynamicImage>,
    pub pdb_raw: String,
}

/// Which UI panels are currently visible.
pub struct DisplayOpts {
    pub show_chr_map:      bool,
    pub show_plasmid_maps: bool,
    pub show_legend:       bool,
    pub show_gene_tracks:  bool,
    pub show_coverage:     bool,
}

impl Default for DisplayOpts {
    fn default() -> Self {
        DisplayOpts {
            show_chr_map:      true,
            show_plasmid_maps: true,
            show_legend:       true,
            show_gene_tracks:  true,
            show_coverage:     true,
        }
    }
}

/// (row_idx, col_start, col_end, feature_idx) — tracks where each gene was rendered.
pub type HitMap = Vec<(usize, usize, usize, usize)>;

/// Data for a single plasmid, stored with local (plasmid-relative) coordinates.
pub struct PlasmidData {
    pub name: String,
    pub features: Vec<Feature>,
    pub genome_size: u64,
    pub gc_skew: Vec<f64>,
    pub view_start: u64,
    pub view_end: u64,
    pub sequence: String,
}

pub struct App {
    pub features: Vec<Feature>,
    pub genome_size: u64,
    pub genome_seq: Option<String>,
    pub genome_name: String,
    pub stop_codons: HashMap<(char, u8), Vec<usize>>,
    pub view_start: u64,
    pub view_end: u64,
    pub selected: Option<usize>,
    pub hovered: Option<usize>,
    pub status_msg: String,
    pub status_msg_ticks: u32,  // countdown; 0 = message expired
    pub on_click_cmd: Option<String>,
    /// (row_idx, col_start, col_end, feature_idx)
    pub hit_map: HitMap,
    pub gc_skew: Vec<f64>,
    pub plasmids: Vec<PlasmidData>,
    pub active_genome: usize,
    pub minimap_rect: Rect,
    pub legend_rect: Rect,
    pub plasmid_rects: Vec<Rect>,
    /// Which circle map is currently hovered: 0 = chromosome, 1+ = plasmid index+1, None = none
    pub hovered_map: Option<usize>,
    pub hovered_legend: bool,
    // Gene search state
    pub search_mode: bool,
    pub search_query: String,
    pub search_results: Vec<usize>,  // indices into active_features()
    pub search_idx: usize,
    pub coverage: Option<crate::core::StrandCoverage>,
    pub display_opts: DisplayOpts,
    pub display_menu_open: bool,
    pub display_menu_idx: usize,
    pub display_menu_rect: Rect,
    // Search results popup (shown when multiple hits exist)
    pub search_popup_open: bool,
    pub search_popup_idx: usize,
    pub search_popup_rect: Rect,
    // Protein structure panel
    pub protein: Option<ProteinPanel>,
    pub protein_panel_rect: Rect,
    pub active_panel: ActivePanel,
    pub gene_track_rect: Rect,
    /// out_dir from --minifold_mlx, used to locate output PDB after folding.
    pub fold_out_dir: Option<String>,
    // MSA panel
    pub msa: Option<MsaPanel>,
    pub msa_panel_rect: Rect,
    pub dmnd_db: Option<String>,
    pub famsa_bin: Option<String>,
    /// Incremented each tick; used by the UI for spinner animation.
    pub anim_tick: u64,
    /// Set to true to request a full terminal clear on the next frame (clears external writes).
    pub needs_clear: bool,
    pub source_files: Vec<String>,
}

impl App {
    pub fn new(
        features: Vec<Feature>,
        genome_size: u64,
        genome_seq: Option<String>,
        genome_name: String,
        stop_codons: HashMap<(char, u8), Vec<usize>>,
        on_click_cmd: Option<String>,
        gc_skew: Vec<f64>,
        plasmids: Vec<PlasmidData>,
        coverage: Option<crate::core::StrandCoverage>,
        fold_out_dir: Option<String>,
        dmnd_db: Option<String>,
        famsa_bin: Option<String>,
    ) -> Self {
        let view_end = genome_size.min(20_000).max(1);
        App {
            features,
            genome_size,
            genome_seq,
            genome_name,
            stop_codons,
            view_start: 1,
            view_end,
            selected: None,
            hovered: None,
            status_msg: String::new(),
            status_msg_ticks: 0,
            on_click_cmd,
            hit_map: Vec::new(),
            gc_skew,
            plasmids,
            active_genome: 0,
            minimap_rect: Rect::default(),
            legend_rect: Rect::default(),
            plasmid_rects: Vec::new(),
            hovered_map: None,
            hovered_legend: false,
            search_mode: false,
            search_query: String::new(),
            search_results: Vec::new(),
            search_idx: 0,
            coverage,
            display_opts: DisplayOpts::default(),
            display_menu_open: false,
            display_menu_idx: 0,
            display_menu_rect: Rect::default(),
            search_popup_open: false,
            search_popup_idx: 0,
            search_popup_rect: Rect::default(),
            protein: None,
            protein_panel_rect: Rect::default(),
            active_panel: ActivePanel::Genome,
            gene_track_rect: Rect::default(),
            fold_out_dir,
            msa: None,
            msa_panel_rect: Rect::default(),
            dmnd_db,
            famsa_bin,
            anim_tick: 0,
            needs_clear: false,
            source_files: Vec::new(),
        }
    }

    pub fn open_protein_panel(&mut self, name: String, nt_seq: String, aa_seq: String) {
        self.protein = Some(ProteinPanel {
            gene_name: name,
            nt_seq,
            aa_seq,
            atoms: Vec::new(),
            camera: crate::pv::camera::Camera::default(),
            drag_last: None,
            folding: false,
            error: None,
            img_cache: None,
            pdb_raw: String::new(),
        });
    }

    /// Set a transient status message that auto-expires after ~5 s (100 ticks at ~50 ms each).
    pub fn set_status(&mut self, msg: impl Into<String>) {
        self.status_msg = msg.into().replace('\n', " ").replace('\r', "");
        self.status_msg_ticks = 100;
    }

    /// Call once per tick; clears the message when its TTL expires.
    pub fn tick_status(&mut self) {
        self.anim_tick = self.anim_tick.wrapping_add(1);
        if self.status_msg_ticks > 0 {
            self.status_msg_ticks -= 1;
            if self.status_msg_ticks == 0 {
                self.status_msg.clear();
            }
        }
    }

    /// Toggle the display option currently highlighted in the menu.
    #[allow(unused_assignments)]
    pub fn toggle_display_opt(&mut self) {
        let has_coverage  = self.coverage.is_some();
        let has_plasmids  = !self.plasmids.is_empty();
        // Build same item list as the UI menu to resolve the index.
        let mut idx = 0;
        macro_rules! item {
            ($field:expr, $avail:expr) => {
                if $avail {
                    if idx == self.display_menu_idx { $field = !$field; return; }
                    idx += 1;
                }
            };
        }
        item!(self.display_opts.show_chr_map,      true);
        item!(self.display_opts.show_plasmid_maps, has_plasmids);
        item!(self.display_opts.show_legend,       true);
        item!(self.display_opts.show_gene_tracks,  true);
        item!(self.display_opts.show_coverage,     has_coverage);
        // "Structure panel" item: closing it = Esc on the protein panel.
        if self.protein.is_some() {
            if idx == self.display_menu_idx {
                self.protein = None;
                self.active_panel = ActivePanel::Genome;
                return;
            }
            idx += 1;
        }
        // "MSA panel" item: closing it clears the msa.
        if self.msa.is_some() {
            if idx == self.display_menu_idx {
                self.msa = None;
                self.active_panel = ActivePanel::Genome;
                return;
            }
        }
    }

    pub fn display_menu_item_count(&self) -> usize {
        let mut n = 3; // chr_map, legend, gene_tracks always present
        if !self.plasmids.is_empty() { n += 1; }
        if self.coverage.is_some()   { n += 1; }
        if self.protein.is_some()    { n += 1; }
        if self.msa.is_some()        { n += 1; }
        n
    }

    pub fn active_features(&self) -> &[Feature] {
        if self.active_genome == 0 {
            &self.features
        } else {
            self.plasmids.get(self.active_genome - 1)
                .map(|p| p.features.as_slice())
                .unwrap_or(&self.features)
        }
    }

    pub fn active_genome_size(&self) -> u64 {
        if self.active_genome == 0 { self.genome_size }
        else { self.plasmids.get(self.active_genome - 1).map(|p| p.genome_size).unwrap_or(self.genome_size) }
    }

    pub fn active_view_start(&self) -> u64 {
        if self.active_genome == 0 { self.view_start }
        else { self.plasmids.get(self.active_genome - 1).map(|p| p.view_start).unwrap_or(self.view_start) }
    }

    pub fn active_view_end(&self) -> u64 {
        if self.active_genome == 0 { self.view_end }
        else { self.plasmids.get(self.active_genome - 1).map(|p| p.view_end).unwrap_or(self.view_end) }
    }

    pub fn switch_genome(&mut self, idx: usize) {
        if idx == 0 || idx <= self.plasmids.len() {
            self.active_genome = idx;
            self.hovered = None;
            self.selected = None;
            self.search_results.clear();
            self.search_popup_open = false;
            self.search_query.clear();
            self.search_mode = false;
        }
    }

    /// Pan the viewport by `fraction` of its current span.
    /// Uses smooth circular wrap: view_start wraps modulo genome_size;
    /// view_end may exceed genome_size to show a wrap-around view.
    pub fn pan(&mut self, fraction: f64) {
        let cur_start = self.active_view_start();
        let cur_end   = self.active_view_end();
        let g_size    = self.active_genome_size();
        // Span is based on the raw stored end (which may already exceed g_size from a wrap)
        let span = cur_end.saturating_sub(cur_start).min(g_size.saturating_sub(1));
        let delta = (span as f64 * fraction).round() as i64;
        // Circular wrap of view_start: modulo genome_size (1-based)
        let raw = cur_start as i64 + delta;
        let new_vs = ((raw - 1).rem_euclid(g_size as i64) + 1) as u64;
        let new_ve = new_vs + span; // may exceed g_size — display handles the wrap
        if self.active_genome == 0 {
            self.view_start = new_vs;
            self.view_end   = new_ve;
        } else if let Some(p) = self.plasmids.get_mut(self.active_genome - 1) {
            p.view_start = new_vs;
            p.view_end   = new_ve;
        }
    }

    /// Zoom the viewport around its center by `factor`. Minimum 100bp window.
    pub fn zoom(&mut self, factor: f64) {
        let cur_start = self.active_view_start();
        let cur_end   = self.active_view_end();
        let g_size    = self.active_genome_size();
        let centre = (cur_start + cur_end) / 2;
        let span = cur_end.saturating_sub(cur_start);
        let new_span = ((span as f64 * factor).round() as u64)
            .max(100)
            .min(g_size);
        let ns = centre.saturating_sub(new_span / 2).max(1);
        let ne = (ns + new_span).min(g_size);
        if self.active_genome == 0 {
            self.view_start = ns;
            self.view_end = ne;
        } else if let Some(p) = self.plasmids.get_mut(self.active_genome - 1) {
            p.view_start = ns;
            p.view_end = ne;
        }
    }

    /// Jump to a specific range.
    pub fn goto(&mut self, start: u64, end: u64) {
        self.view_start = start.max(1);
        self.view_end = end.min(self.genome_size);
    }

    /// Jump the active genome's viewport to a specific range.
    pub fn goto_active(&mut self, start: u64, end: u64) {
        let g = self.active_genome_size();
        if self.active_genome == 0 {
            self.view_start = start.max(1);
            self.view_end   = end.min(g);
        } else if let Some(p) = self.plasmids.get_mut(self.active_genome - 1) {
            p.view_start = start.max(1);
            p.view_end   = end.min(g);
        }
    }

    /// Jump viewport to show a feature (by its index in active_features()).
    pub fn jump_to_feature(&mut self, feat_idx: usize) {
        let (fs, fe) = {
            let feats = self.active_features();
            match feats.get(feat_idx) {
                Some(f) => (f.start, f.end),
                None => return,
            }
        };
        let span = fe.saturating_sub(fs) + 1;
        let pad  = (span * 2).max(5000);
        let start = fs.saturating_sub(pad).max(1);
        let end   = fe.saturating_add(pad);
        self.goto_active(start, end);
    }

    /// Parse a query as a genomic coordinate: "1234567" or "1234-5678" / "1234..5678".
    /// Returns Some((start, end)) if it parses as a coordinate, None otherwise.
    fn parse_coord(query: &str) -> Option<(u64, u64)> {
        let q = query.replace(',', "").replace('_', "");
        // Range formats: "start-end" or "start..end"
        let sep = if q.contains("..") { Some("..") }
                  else if let Some(i) = q[1..].find('-').map(|j| j + 1) {
                      if q[i..].parse::<u64>().is_ok() { Some(&q[i..i+1]) } else { None }
                  } else { None };
        if let Some(sep) = sep {
            let parts: Vec<&str> = q.splitn(2, sep).collect();
            if parts.len() == 2 {
                if let (Ok(s), Ok(e)) = (parts[0].parse::<u64>(), parts[1].parse::<u64>()) {
                    return Some((s, e));
                }
            }
        }
        // Single position
        if let Ok(pos) = q.parse::<u64>() {
            return Some((pos.saturating_sub(1000).max(1), pos + 1000));
        }
        None
    }

    /// Run a case-insensitive prefix search over name and locus_tag in active features.
    /// Populates search_results and jumps to the first match.
    /// If query looks like a coordinate, jumps there instead.
    pub fn do_search(&mut self) {
        let query = self.search_query.trim().to_string();
        if query.is_empty() { return; }

        // Coordinate search: if query parses as a position, jump to it directly
        if let Some((start, end)) = Self::parse_coord(&query) {
            self.search_results.clear();
            self.search_popup_open = false;
            self.goto_active(start, end);
            self.status_msg = format!("→ {}–{}", start, end);
            return;
        }

        let ql = query.to_ascii_lowercase();
        // Collect matching indices (exact first, then prefix)
        let mut exact:  Vec<usize> = Vec::new();
        let mut prefix: Vec<usize> = Vec::new();
        for (i, f) in self.active_features().iter().enumerate() {
            let n = f.name.to_ascii_lowercase();
            let l = f.locus_tag.to_ascii_lowercase();
            if n == ql || l == ql {
                exact.push(i);
            } else if n.starts_with(&ql) || l.starts_with(&ql) {
                prefix.push(i);
            }
        }
        // Merge, sort by start position, deduplicate
        exact.extend(prefix);
        exact.dedup();
        {
            let feats = self.active_features();
            exact.sort_unstable_by_key(|&i| feats[i].start);
        }

        if exact.is_empty() {
            self.status_msg = format!("\"{}\" not found", query);
            self.search_results.clear();
            self.search_popup_open = false;
        } else if exact.len() == 1 {
            self.search_results = exact;
            self.search_idx = 0;
            self.search_popup_open = false;
            self.jump_search_result();
        } else {
            self.search_results = exact;
            self.search_idx = 0;
            self.search_popup_open = true;
            self.search_popup_idx = 0;
            self.jump_search_result();
        }
    }

    /// Jump to the search result currently highlighted in the popup.
    pub fn search_popup_select(&mut self) {
        if self.search_popup_idx < self.search_results.len() {
            self.search_idx = self.search_popup_idx;
            self.search_popup_open = false;
            self.jump_search_result();
        }
    }

    /// Jump to the current search result and update status_msg.
    fn jump_search_result(&mut self) {
        let idx = self.search_idx;
        let n   = self.search_results.len();
        if n == 0 { return; }
        let feat_idx = self.search_results[idx];
        self.jump_to_feature(feat_idx);
        let (name, locus, strand, start, end) = {
            let f = &self.active_features()[feat_idx];
            (f.name.clone(), f.locus_tag.clone(), f.strand, f.start, f.end)
        };
        let counter = if n > 1 { format!("  [{}/{}]", idx + 1, n) } else { String::new() };
        self.status_msg = format!(
            "{} {} {}  {}–{}{}",
            name, locus, strand, start, end, counter
        );
    }

    /// Advance to the next search result.
    pub fn search_next(&mut self) {
        if self.search_results.is_empty() { return; }
        self.search_idx = (self.search_idx + 1) % self.search_results.len();
        self.jump_search_result();
    }

    /// Visible features in the current viewport.
    pub fn visible_features(&self) -> Vec<(usize, &Feature)> {
        self.features
            .iter()
            .enumerate()
            .filter(|(_, f)| f.end >= self.view_start && f.start <= self.view_end)
            .collect()
    }
}
