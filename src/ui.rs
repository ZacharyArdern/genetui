use std::f64::consts::PI;

use ratatui::{
    Frame,
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Style},
    text::{Line, Span},
    widgets::{
        Block, Borders, Clear, Paragraph,
        canvas::{Canvas, Points},
    },
};

use crate::app::App;
use crate::core::{format_bp_tick, nice_tick_spacing};

// Plus strand gene colours — bright teal/cyan family
const PLUS_COLORS: &[(u8, u8, u8)] = &[
    (72,  210, 180),
    (80,  170, 230),
    (100, 220, 150),
    (55,  225, 200),
    (95,  185, 240),
    (120, 215, 170),
];

// Minus strand gene colours — bright orange/amber family
const MINUS_COLORS: &[(u8, u8, u8)] = &[
    (240, 148,  50),
    (225, 118,  70),
    (245, 175,  45),
    (230, 132,  80),
    (238, 160,  55),
    (220, 105,  75),
];

/// 9 logical track rows: +nc, +fr0, +fr1, +fr2, ruler, -fr0, -fr1, -fr2, -nc
const TRACK_ROWS: usize = 9;
const RULER_ROW: usize = 4;
/// Physical display lines: 2 nc rows × 1 + 6 coding rows × 2 + 1 ruler = 15
const TRACK_PHYS: usize = 15;

/// Frame backgrounds — dark but distinct, easy on the eye.
const FRAME_BG: &[(u8, u8, u8)] = &[
    (18,  25,  48),  // frame 0 — dark navy
    (38,  30,  16),  // frame 1 — dark brown
    (30,  16,  40),  // frame 2 — dark purple
];

const MINIMAP_W: u16 = 44;
const MINIMAP_H: u16 = 22; // W ≈ 2*H keeps braille dots square (cells are ~2:1 tall:wide)

pub fn draw(f: &mut Frame, app: &mut App) {
    let size = f.area();

    const COV_H: u16 = 4;
    // Coverage is chromosome-only; don't show it when a plasmid is active.
    let cov_shown   = app.coverage.is_some() && app.display_opts.show_coverage && app.active_genome == 0;
    let cov_extra   = if cov_shown { COV_H * 2 } else { 0 };
    let fixed_rows  = TRACK_PHYS as u16 + 1 + 2 + 1 + cov_extra;
    let minimap_h   = MINIMAP_H.min(size.height.saturating_sub(fixed_rows));
    let track_min   = TRACK_PHYS as u16 + 1 + cov_extra;

    let msa_open = app.msa.is_some();
    let (track_constraint, msa_constraint) = if msa_open {
        (Constraint::Length(track_min), Constraint::Min(10))
    } else {
        (Constraint::Min(track_min), Constraint::Length(0))
    };

    let vertical = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            track_constraint,
            msa_constraint,
            Constraint::Length(minimap_h),
            Constraint::Length(2),
            Constraint::Length(1),
        ])
        .split(size);

    let track_area  = vertical[0];
    app.gene_track_rect = track_area;
    let msa_area    = vertical[1];
    let bottom_area = vertical[2];
    let info_area   = vertical[3];
    let status_area = vertical[4];

    // Protein panel splits the bottom row horizontally when open.
    const LEGEND_W: u16 = 16;
    const PLASMID_W: u16 = 22; // one plasmid map column width
    // maps_min_w must fit legend + minimap + at least one plasmid map.
    let maps_min_w = LEGEND_W + MINIMAP_W + PLASMID_W;
    let (protein_area_opt, maps_area) = if app.protein.is_some()
        && bottom_area.width > maps_min_w + 30
    {
        // Give maps area exactly maps_min_w; protein panel gets the rest.
        let maps_w = maps_min_w;
        let protein_w = bottom_area.width.saturating_sub(maps_w);
        let pa = Rect { x: bottom_area.x, y: bottom_area.y,
            width: protein_w, height: bottom_area.height };
        let ma = Rect { x: bottom_area.x + protein_w, y: bottom_area.y,
            width: maps_w, height: bottom_area.height };
        app.protein_panel_rect = pa;
        (Some(pa), ma)
    } else {
        app.protein_panel_rect = Rect::default();
        (None, bottom_area)
    };

    // Bottom row: legend on left, maps (minimap + plasmids) centred in remaining space.
    let legend_w  = if app.display_opts.show_legend  { LEGEND_W } else { 0 };
    let minimap_w = if app.display_opts.show_chr_map { MINIMAP_W.min(maps_area.width) } else { 0 };

    let legend_area  = Rect { x: maps_area.x, y: maps_area.y,
        width: legend_w, height: maps_area.height };

    // Centre the minimap + plasmid block within the space after the legend.
    let avail_x      = maps_area.x + legend_w;
    let avail_w      = (maps_area.x + maps_area.width).saturating_sub(avail_x);
    let n_plasmids   = if app.display_opts.show_plasmid_maps { app.plasmids.len() } else { 0 };
    let plasmid_total_w = (n_plasmids as u16).saturating_mul(PLASMID_W).min(avail_w.saturating_sub(minimap_w));
    let maps_total_w = minimap_w + plasmid_total_w;
    let maps_offset  = avail_w.saturating_sub(maps_total_w) / 2;
    let maps_start_x = avail_x + maps_offset;

    let minimap_area = Rect { x: maps_start_x, y: maps_area.y,
        width: minimap_w, height: maps_area.height };
    let plasmid_x    = maps_start_x + minimap_w;
    let plasmid_area = Rect {
        x: plasmid_x, y: maps_area.y,
        width: plasmid_total_w,
        height: maps_area.height,
    };

    app.minimap_rect = if app.display_opts.show_chr_map { minimap_area } else { Rect::default() };
    app.legend_rect  = if app.display_opts.show_legend  { legend_area  } else { Rect::default() };

    // Split track area: cov+ | gene tracks | cov- | spacer
    if cov_shown && app.display_opts.show_gene_tracks {
        let track_split = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(COV_H),
                Constraint::Length(TRACK_PHYS as u16),
                Constraint::Length(COV_H),
                Constraint::Min(0),
            ])
            .split(track_area);
        draw_coverage_track(f, app, track_split[0], '+');
        draw_gene_tracks(f, app, track_split[1]);
        draw_coverage_track(f, app, track_split[2], '-');
    } else if cov_shown {
        let track_split = Layout::default()
            .direction(Direction::Vertical)
            .constraints([Constraint::Length(COV_H), Constraint::Length(COV_H), Constraint::Min(0)])
            .split(track_area);
        draw_coverage_track(f, app, track_split[0], '+');
        draw_coverage_track(f, app, track_split[1], '-');
    } else if app.display_opts.show_gene_tracks {
        draw_gene_tracks(f, app, track_area);
    }

    // Fill entire bottom row with black
    f.render_widget(Block::default().style(Style::default().bg(Color::Black)), bottom_area);

    // White border around gene track area when it's the active panel
    if app.protein.is_some() && app.active_panel == crate::app::ActivePanel::Genome {
        let border_block = Block::default()
            .borders(Borders::ALL)
            .border_style(Style::default().fg(Color::White));
        f.render_widget(border_block, track_area);
    }

    if let Some(pa) = protein_area_opt { draw_protein_panel(f, app, pa); }
    if msa_open                           { draw_msa_panel(f, app, msa_area); }
    if app.display_opts.show_chr_map      { draw_minimap(f, app, minimap_area); }
    if app.display_opts.show_legend       { draw_legend(f, app, legend_area); }
    if app.display_opts.show_plasmid_maps { draw_plasmid_maps(f, app, plasmid_area); }
    draw_info_panel(f, app, info_area);
    draw_status(f, app, status_area);
    if app.search_popup_open  { draw_search_popup(f, app); }
    if app.display_menu_open  { draw_display_menu(f, app); }
}

/// Collected feature data for rendering, avoids borrow conflicts with hit_map write.
struct VFeat {
    idx: usize,
    cs: usize,
    ce: usize,
    fr: usize,
    strand: char,
    name: String,
    color_idx: usize,
    #[allow(dead_code)]
    is_orf: bool,
    noncoding: bool,
}

fn draw_gene_tracks(f: &mut Frame, app: &mut App, area: Rect) {
    let width = area.width as usize;
    let height = area.height as usize;

    if width == 0 || height == 0 {
        return;
    }

    // Clear area first so stale cells from previous frames don't ghost through.
    f.render_widget(Clear, area);

    const LABEL_W: usize = 4;
    let feat_w = if width > LABEL_W { width - LABEL_W } else { 1 };

    let view_start  = app.active_view_start();
    let view_end    = app.active_view_end();
    let genome_size = app.active_genome_size();
    let span        = (view_end - view_start) as f64;
    // wrap_end > 0 means the view wraps past genome end; features in [1, wrap_end] are shown
    let wrap_end: u64 = if view_end > genome_size { view_end - genome_size } else { 0 };

    // bp_to_col: handles both the normal zone [view_start, genome_size]
    // and the wrapped zone [1, wrap_end] (added genome_size to reach past the seam).
    let bp_to_col = |bp: u64| -> usize {
        if span <= 0.0 { return 0; }
        let eff = if bp >= view_start {
            bp as f64
        } else {
            bp as f64 + genome_size as f64  // wrapped position
        };
        let col = ((eff - view_start as f64) / span * (feat_w - 1) as f64).round() as i64;
        col.max(0).min(feat_w as i64 - 1) as usize
    };

    // Collect visible features: each feature may appear once (normal) or twice (both zones).
    let vfeats: Vec<VFeat> = {
        let feats = app.active_features();
        let mut out = Vec::new();
        for (idx, f) in feats.iter().enumerate() {
            let fr = (if f.strand == '+' {
                f.start.saturating_sub(1) % 3
            } else {
                f.end.saturating_sub(1) % 3
            }) as usize;
            // Normal zone: feature overlaps [view_start, genome_size]
            if f.end >= view_start && f.start <= genome_size {
                let cs = bp_to_col(f.start.max(view_start));
                let ce = bp_to_col(f.end.min(genome_size));
                out.push(VFeat { idx, cs, ce: ce.max(cs), fr, strand: f.strand,
                    name: f.name.clone(), color_idx: f.color_idx, is_orf: f.is_orf,
                    noncoding: f.noncoding });
            }
            // Wrapped zone: feature overlaps [1, wrap_end]
            if wrap_end > 0 && f.start <= wrap_end {
                let cs = bp_to_col(f.start.max(1));
                let ce = bp_to_col(f.end.min(wrap_end));
                out.push(VFeat { idx, cs, ce: ce.max(cs), fr, strand: f.strand,
                    name: f.name.clone(), color_idx: f.color_idx, is_orf: f.is_orf,
                    noncoding: f.noncoding });
            }
        }
        out
    };

    // Build buckets: (strand, frame) -> Vec<index into vfeats> (coding only)
    // nc_plus/nc_minus: non-coding features per strand
    let mut buckets: std::collections::HashMap<(char, u8), Vec<usize>> =
        std::collections::HashMap::new();
    let mut nc_plus:  Vec<usize> = Vec::new();
    let mut nc_minus: Vec<usize> = Vec::new();
    for &s in &['+', '-'] {
        for fr in 0u8..3 {
            buckets.insert((s, fr), Vec::new());
        }
    }
    for (vi, vf) in vfeats.iter().enumerate() {
        if vf.noncoding {
            if vf.strand == '+' { nc_plus.push(vi); } else { nc_minus.push(vi); }
        } else {
            buckets.entry((vf.strand, vf.fr as u8)).or_default().push(vi);
        }
    }

    // Build stop codon columns per (strand, frame)
    let stop_cols: std::collections::HashMap<(char, u8), Vec<usize>> = if app.active_genome == 0 {
        let mut map = std::collections::HashMap::new();
        for &s in &['+', '-'] {
            for fr in 0u8..3 {
                let positions = app.stop_codons.get(&(s, fr)).cloned().unwrap_or_default();
                let cols: Vec<usize> = positions
                    .iter()
                    .filter(|&&p| {
                        // Normal zone or wrapped zone
                        (p as u64 >= view_start && p as u64 <= genome_size) ||
                        (wrap_end > 0 && p as u64 <= wrap_end)
                    })
                    .map(|&p| bp_to_col(p as u64))
                    .filter(|&c| c < feat_w)
                    .collect();
                map.insert((s, fr), cols);
            }
        }
        map
    } else {
        let mut m = std::collections::HashMap::new();
        for &s in &['+', '-'] { for fr in 0u8..3 { m.insert((s, fr), Vec::<usize>::new()); } }
        m
    };

    let top_pad = if height > TRACK_PHYS { (height - TRACK_PHYS) / 2 } else { 0 };

    let hovered_feat  = app.hovered;
    let selected_feat = app.selected;

    let mut new_hit_map: crate::app::HitMap = Vec::new();
    let mut lines: Vec<Line> = Vec::new();

    for _ in 0..top_pad {
        lines.push(Line::from(vec![Span::raw(" ".repeat(width))]));
    }

    let mut phys_row = top_pad;

    // Helper: render a non-coding row (rRNA, tRNA etc.) — one line, name inside white bar.
    let render_nc = |lines: &mut Vec<Line>,
                     new_hit_map: &mut crate::app::HitMap,
                     nc_list: &[usize],
                     strand: char,
                     phys_row: usize| {
        let label_color = if strand == '+' { Color::Rgb(137, 180, 250) } else { Color::Rgb(243, 139, 168) };
        let nc_bg = Color::Rgb(15, 15, 25);
        // Single cell array: bg=White inside genes, name char where it fits.
        let mut cells: Vec<(char, Color, Option<Color>)> = vec![(' ', Color::Reset, None); feat_w];
        for &vi in nc_list {
            let vf = &vfeats[vi];
            let (cs, ce) = (vf.cs, vf.ce);
            let fw = (ce + 1).saturating_sub(cs).max(1);
            // Fill extent with white background
            for j in 0..fw {
                if cs + j < feat_w {
                    cells[cs + j] = (' ', Color::Rgb(10, 10, 10), Some(Color::White));
                }
            }
            // Name centered inside the bar
            let name_chars: Vec<char> = vf.name.chars().collect();
            if name_chars.len() <= fw {
                let off = (fw - name_chars.len()) / 2;
                for (j, &ch) in name_chars.iter().enumerate() {
                    let col = cs + off + j;
                    if col < feat_w { cells[col] = (ch, Color::Rgb(10, 10, 10), Some(Color::White)); }
                }
            }
            let abs_row = area.y as usize + phys_row;
            let abs_cs  = area.x as usize + LABEL_W + cs;
            let abs_ce  = area.x as usize + LABEL_W + ce;
            new_hit_map.push((abs_row, abs_cs, abs_ce, vf.idx));
        }
        lines.push(cells_to_line(cells, " nc ", label_color, true, nc_bg));
    };

    for row_idx in 0..TRACK_ROWS {
        // Non-coding rows (single line each)
        if row_idx == 0 {
            let nc = nc_plus.clone();
            render_nc(&mut lines, &mut new_hit_map, &nc, '+', phys_row);
            phys_row += 1;
            continue;
        }
        if row_idx == TRACK_ROWS - 1 {
            let nc = nc_minus.clone();
            render_nc(&mut lines, &mut new_hit_map, &nc, '-', phys_row);
            phys_row += 1;
            continue;
        }
        if row_idx == RULER_ROW {
            lines.push(make_ruler(view_start, view_end, genome_size, feat_w, LABEL_W));
            phys_row += 1;
            continue;
        }

        let (strand, fr) = if row_idx < RULER_ROW {
            ('+', (row_idx - 1) as u8)   // rows 1,2,3 → fr 0,1,2
        } else {
            ('-', (row_idx - RULER_ROW - 1) as u8)  // rows 5,6,7 → fr 0,1,2
        };

        let label_color = if strand == '+' {
            Color::Rgb(137, 180, 250)
        } else {
            Color::Rgb(243, 139, 168)
        };

        let vi_list = buckets.get(&(strand, fr)).cloned().unwrap_or_default();
        let stops = stop_cols.get(&(strand, fr)).cloned().unwrap_or_default();

        let (fbr, fbg, fbb) = FRAME_BG[fr as usize % FRAME_BG.len()];
        let frame_bg = Color::Rgb(fbr, fbg, fbb);

        let mut cells_arrow: Vec<(char, Color, Option<Color>)> = vec![(' ', Color::Reset, None); feat_w];
        let mut cells_gap:   Vec<(char, Color, Option<Color>)> = vec![(' ', Color::Reset, None); feat_w];

        for &col in &stops {
            if col < feat_w { cells_gap[col] = ('│', Color::Rgb(204, 68, 68), None); }
        }

        for &vi in &vi_list {
            let vf = &vfeats[vi];
            let cs = vf.cs;
            let ce = vf.ce;
            let fw = (ce + 1).saturating_sub(cs).max(1);
            let (r, g, b) = if vf.strand == '+' {
                PLUS_COLORS[vf.color_idx % PLUS_COLORS.len()]
            } else {
                MINUS_COLORS[vf.color_idx % MINUS_COLORS.len()]
            };
            let is_selected = selected_feat == Some(vf.idx);
            let is_hovered  = hovered_feat  == Some(vf.idx);
            let gene_col = if is_selected {
                // White fill with dark text for selected gene
                Color::White
            } else if is_hovered {
                // Brighten toward white
                let bri = |x: u8| -> u8 { 255u8.min(x.saturating_add(70)) };
                Color::Rgb(bri(r), bri(g), bri(b))
            } else {
                Color::Rgb(r, g, b)
            };
            let tip = if vf.strand == '+' { '▶' } else { '◀' };
            for j in 0..fw {
                if cs + j < feat_w {
                    let ch = if vf.strand == '+' && j == fw - 1 { tip }
                             else if vf.strand == '-' && j == 0  { tip }
                             else { ' ' };
                    let fg = if ch == tip { Color::White } else { Color::Rgb(15, 15, 25) };
                    cells_arrow[cs + j] = (ch, fg, Some(gene_col));
                }
            }

            let name_chars: Vec<char> = vf.name.chars().collect();
            if name_chars.len() <= fw {
                let offset = (fw - name_chars.len()) / 2;
                let name_start = cs + offset;
                let name_end = name_start + name_chars.len();
                let blocked = (name_start..name_end).any(|col| col < feat_w && cells_gap[col].0 == '│');
                if !blocked {
                    for (j, &ch) in name_chars.iter().enumerate() {
                        let col = cs + offset + j;
                        if col < feat_w { cells_gap[col] = (ch, gene_col, None); }
                    }
                }
            }

            let abs_row = area.y as usize + phys_row;
            let abs_cs  = area.x as usize + LABEL_W + cs;
            let abs_ce  = area.x as usize + LABEL_W + ce;
            new_hit_map.push((abs_row,     abs_cs, abs_ce, vf.idx));
            new_hit_map.push((abs_row + 1, abs_cs, abs_ce, vf.idx));
        }

        let label = format!("{}{} ", strand, fr);
        lines.push(cells_to_line(cells_arrow, &label, label_color, true,  frame_bg));
        lines.push(cells_to_line(cells_gap,   "   ", Color::DarkGray, false, frame_bg));
        phys_row += 2;
    }

    app.hit_map = new_hit_map;
    f.render_widget(Paragraph::new(lines), area);
}

/// Convert a cell row to a ratatui Line, with label and separator prepended.
fn cells_to_line(
    cells: Vec<(char, Color, Option<Color>)>,
    label: &str,
    label_color: Color,
    bold: bool,
    row_bg: Color,
) -> Line<'static> {
    let label_style = {
        let s = Style::default().fg(label_color).bg(row_bg);
        if bold { s.add_modifier(ratatui::style::Modifier::BOLD) } else { s }
    };
    let mut spans: Vec<Span<'static>> = vec![
        Span::styled(label.to_owned(), label_style),
        Span::styled("│", Style::default().fg(Color::DarkGray).bg(row_bg)),
    ];

    let mut i = 0;
    while i < cells.len() {
        let (ch, fg, bg) = cells[i];
        let mut j = i + 1;
        while j < cells.len()
            && cells[j].1 == fg
            && cells[j].2 == bg
            && ch != '│'
            && cells[j].0 != '│'
        {
            j += 1;
        }
        let text: String = cells[i..j].iter().map(|(c, _, _)| *c).collect();
        let style = match bg {
            Some(bg_col) => Style::default().fg(fg).bg(bg_col),
            // No explicit bg: use frame bg (covers empty space, stop codons, gene name text)
            None => if fg == Color::Reset {
                Style::default().bg(row_bg)
            } else {
                Style::default().fg(fg).bg(row_bg)
            },
        };
        spans.push(Span::styled(text, style));
        i = j;
    }
    Line::from(spans)
}


fn make_ruler(view_start: u64, view_end: u64, genome_size: u64, feat_w: usize, label_w: usize) -> Line<'static> {
    let span = view_end - view_start;
    let tick_bp = nice_tick_spacing(span, feat_w);
    let wrap_end: u64 = if view_end > genome_size { view_end - genome_size } else { 0 };

    let bp_to_col = |bp: u64| -> usize {
        if span == 0 || feat_w == 0 { return 0; }
        let eff = if bp >= view_start { bp as f64 } else { bp as f64 + genome_size as f64 };
        let col = ((eff - view_start as f64) / span as f64 * (feat_w - 1) as f64).round() as i64;
        col.max(0).min(feat_w as i64 - 1) as usize
    };

    let mut chars: Vec<char> = vec!['─'; feat_w];

    // Place a tick and label at the given position
    let place_tick = |tick: u64, chars: &mut Vec<char>| {
        let col = bp_to_col(tick);
        if col < feat_w {
            chars[col] = '┼';
            let lbl = format_bp_tick(tick, tick_bp);
            for (j, ch) in lbl.chars().enumerate() {
                let c = col + j + 1;
                if c < feat_w { chars[c] = ch; }
            }
        }
    };

    // Normal zone: [view_start, min(view_end, genome_size)]
    if tick_bp > 0 {
        let first_tick = (view_start / tick_bp + 1) * tick_bp;
        let mut tick = first_tick;
        let normal_end = view_end.min(genome_size);
        while tick <= normal_end {
            place_tick(tick, &mut chars);
            tick += tick_bp;
        }
        // Wrapped zone: [1, wrap_end]
        if wrap_end > 0 {
            // Genome boundary marker
            let boundary_col = bp_to_col(genome_size);
            if boundary_col > 0 && boundary_col < feat_w { chars[boundary_col] = '│'; }
            // Ticks from 0
            let mut tick_w = tick_bp;
            while tick_w <= wrap_end {
                place_tick(tick_w, &mut chars);
                tick_w += tick_bp;
            }
        }
    }

    let ruler_str: String = chars.into_iter().collect();
    let label_dashes = "─".repeat(label_w);

    Line::from(vec![
        Span::styled(label_dashes, Style::default().fg(Color::DarkGray)),
        Span::styled(ruler_str, Style::default().fg(Color::Rgb(170, 170, 204))),
    ])
}

/// Map a normalised GC skew value [-1,1] to RGB.
/// t=+1 → bright green (G surplus); t=-1 → bright purple (C surplus); t=0 → mid-grey.
fn skew_to_color(t: f64) -> Color {
    let t = t.clamp(-1.0, 1.0);
    if t >= 0.0 {
        // grey → green
        let r = (80.0 - t * 50.0).round() as u8;
        let g = (80.0 + t * 150.0).round() as u8;
        let b = (80.0 - t * 40.0).round() as u8;
        Color::Rgb(r, g, b)
    } else {
        let s = -t;
        // grey → purple
        let r = (80.0 + s * 140.0).round() as u8;
        let g = (80.0 - s * 55.0).round() as u8;
        let b = (80.0 + s * 175.0).round() as u8;
        Color::Rgb(r, g, b)
    }
}

fn draw_minimap(f: &mut Frame, app: &App, area: Rect) {
    let genome_size = app.genome_size as f64;
    let w = area.width  as f64;
    let h = area.height as f64;
    if w == 0.0 || h == 0.0 { return; }

    const CELL_RATIO: f64 = 2.15;
    let x_range = 1.1;
    let y_range = x_range * (CELL_RATIO * h) / w;

    // Radii (data units)
    const R_OUT:      f64 = 0.82;                    // outer edge of density bands
    const R_IN:       f64 = R_OUT * 0.70;            // ≈ 0.574 — inner edge of density bands
    const R_MID:      f64 = (R_OUT + R_IN) * 0.5;   // plus/minus band boundary ≈ 0.697
    const R_GC_INNER: f64 = R_OUT + 0.02;            // GC skew ring inner edge (just outside bands)
    const R_GC_OUTER: f64 = R_OUT + 0.07;            // GC skew ring outer edge (0.05 wide)
    const GENE_OUT:   f64 = R_OUT + 0.14;            // viewport marker radius (outside GC ring)

    let dot_size = x_range / (2.0 * w);
    let ang_step = (dot_size / R_OUT).max(0.004);
    let rad_step = dot_size.max(0.008);

    let pos_angle = |pos: f64| -> f64 { (pos / genome_size) * 2.0 * PI };

    // ── Gene density per angular bin ──────────────────────────────────────────
    let n_dens = 360usize;
    let mut plus_dens  = vec![0u32; n_dens];
    let mut minus_dens = vec![0u32; n_dens];
    for feat in &app.features {
        let mid = (feat.start as f64 + feat.end as f64) / 2.0;
        let frac = (mid / genome_size).clamp(0.0, 1.0 - 1e-9);
        let bin = (frac * n_dens as f64) as usize;
        if feat.strand == '+' { plus_dens[bin]  += 1; }
        else                  { minus_dens[bin] += 1; }
    }
    let max_plus  = (*plus_dens .iter().max().unwrap_or(&1)).max(1) as f64;
    let max_minus = (*minus_dens.iter().max().unwrap_or(&1)).max(1) as f64;

    // ── Layer A: GC skew ring ─────────────────────────────────────────────────
    // Always rendered: grey when no FASTA loaded, coloured otherwise.
    // Sweep by angle (not by window) so density is uniform regardless of n_win.
    let mut skew_map: std::collections::HashMap<(u8,u8,u8), Vec<(f64,f64)>> = Default::default();
    {
        let n_win = app.gc_skew.len();
        // Normalise using actual observed min/max so the full colour range is used
        let (gc_min, gc_max) = if n_win > 0 {
            let mn = app.gc_skew.iter().cloned().fold(f64::INFINITY,     f64::min);
            let mx = app.gc_skew.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            (mn, mx)
        } else { (-1.0, 1.0) };
        let gc_range = (gc_max - gc_min).max(1e-9);

        let n_ang_gc = (2.0 * PI / ang_step).ceil() as usize;
        for i in 0..n_ang_gc {
            let a    = i as f64 * ang_step;
            let frac = a / (2.0 * PI);
            // t in [-1, +1]: min → -1 (C-rich / purple), max → +1 (G-rich / green)
            let t = if n_win > 0 {
                let win = ((frac * n_win as f64) as usize).min(n_win - 1);
                ((app.gc_skew[win] - gc_min) / gc_range * 2.0 - 1.0).clamp(-1.0, 1.0)
            } else {
                0.0  // grey ring when no data
            };
            let (rv, gv, bv) = match skew_to_color(t) {
                Color::Rgb(r, g, b) => (r, g, b),
                _ => (80, 80, 80),
            };
            let (sa, ca) = (a.sin(), a.cos());
            let entry = skew_map.entry((rv, gv, bv)).or_default();
            let mut r = R_GC_INNER;
            while r <= R_GC_OUTER {
                entry.push((r * sa, r * ca));
                r += rad_step;
            }
        }
    }
    let mut skew_batches: Vec<(Color, Vec<(f64,f64)>)> = skew_map.into_iter()
        .map(|((r,g,b), pts)| (Color::Rgb(r,g,b), pts))
        .collect();
    skew_batches.sort_unstable_by_key(|(c, _)| match c { Color::Rgb(r,g,b) => (*r,*g,*b), _ => (0,0,0) });

    // ── Layer B: Gene density bands (plus outer, minus inner) ─────────────────
    let mut density_map: std::collections::HashMap<(u8,u8,u8), Vec<(f64,f64)>> = Default::default();
    for i in 0..n_dens {
        let a_start = i       as f64 / n_dens as f64 * 2.0 * PI;
        let a_end   = (i + 1) as f64 / n_dens as f64 * 2.0 * PI;
        let tp = plus_dens[i]  as f64 / max_plus;
        let tm = minus_dens[i] as f64 / max_minus;
        // Plus: dark navy → bright teal
        let cp = (
            (15.0 + tp * 57.0 ).round() as u8,
            (20.0 + tp * 190.0).round() as u8,
            (20.0 + tp * 160.0).round() as u8,
        );
        // Minus: dark brown → bright amber
        let cm = (
            (15.0 + tm * 225.0).round() as u8,
            (10.0 + tm * 138.0).round() as u8,
            ( 5.0 + tm *  45.0).round() as u8,
        );

        let mut plus_pts_bin:  Vec<(f64,f64)> = Vec::new();
        let mut minus_pts_bin: Vec<(f64,f64)> = Vec::new();
        let mut a = a_start;
        while a < a_end {
            let (sa, ca) = (a.sin(), a.cos());
            let mut r = R_MID;
            while r <= R_OUT { plus_pts_bin.push((r * sa, r * ca)); r += rad_step; }
            let mut r = R_IN;
            while r <  R_MID { minus_pts_bin.push((r * sa, r * ca)); r += rad_step; }
            a += ang_step;
        }
        density_map.entry(cp).or_default().extend(plus_pts_bin);
        density_map.entry(cm).or_default().extend(minus_pts_bin);
    }
    let mut density_batches: Vec<(Color, Vec<(f64,f64)>)> = density_map.into_iter()
        .map(|((r,g,b), pts)| (Color::Rgb(r,g,b), pts))
        .collect();
    density_batches.sort_unstable_by_key(|(c, _)| match c { Color::Rgb(r,g,b) => (*r,*g,*b), _ => (0,0,0) });

    // ── Layer C: 1 Mb position ticks and labels ───────────────────────────────
    struct MbLabel { x: f64, y: f64, text: String }
    let mut tick_pts: Vec<(f64,f64)> = Vec::new();
    let mut mb_labels: Vec<MbLabel> = Vec::new();
    let mb_size = 1_000_000u64;
    let mut pos_mb = 0u64;
    while pos_mb < app.genome_size {
        let a = pos_angle(pos_mb as f64);
        let (sa, ca) = (a.sin(), a.cos());
        for k in 0..=4usize {
            let r = GENE_OUT + 0.03 + 0.07 * k as f64 / 4.0;
            tick_pts.push((r * sa, r * ca));
        }
        mb_labels.push(MbLabel {
            x: (GENE_OUT + 0.16) * sa,
            y: (GENE_OUT + 0.16) * ca,
            text: format!("{}Mb", pos_mb / mb_size),
        });
        pos_mb += mb_size;
    }

    // ── Layer D: Track circle + viewport position marker ─────────────────────
    // Thin circle at GENE_OUT that the marker rides along
    let mut track_pts: Vec<(f64,f64)> = Vec::new();
    {
        let n_track = (2.0 * PI / ang_step).ceil() as usize;
        for i in 0..n_track {
            let a = i as f64 * ang_step;
            track_pts.push((GENE_OUT * a.sin(), GENE_OUT * a.cos()));
        }
    }

    // ── Layer D: track circle + viewport position marker ─────────────────────
    // Marker is only shown when the chromosome is the active genome.
    let chromosome_active = app.active_genome == 0;

    let view_centre = (app.view_start as f64 + app.view_end as f64) / 2.0;
    let ma = pos_angle(view_centre.min(genome_size));
    let (msa, mca) = (ma.sin(), ma.cos());
    let mx = GENE_OUT * msa;
    let my = GENE_OUT * mca;
    let ddx = x_range / w;

    // ── Text helpers ──────────────────────────────────────────────────────────
    let char_w = 2.0 * x_range / w;
    let line_h = 2.0 * y_range / h;

    // Word-wrap genome name to fit inside inner circle (no word breaking, max 2 lines).
    // Skip leading accession-style words (contain a digit: NZ_CP015855.1, NC_000001, etc.)
    let inner_usable = (2.0 * R_IN * 0.85 / char_w).floor() as usize;
    let max_chars = inner_usable.max(6);
    let name_lines: Vec<String> = {
        let all_words: Vec<&str> = app.genome_name.split_whitespace().collect();
        let words: Vec<&str> = {
            let skip = all_words.iter().take_while(|w| w.chars().any(|c| c.is_ascii_digit())).count();
            if skip < all_words.len() { all_words[skip..].to_vec() } else { all_words }
        };
        let mut lines: Vec<String> = Vec::new();
        let mut cur = String::new();
        for w in &words {
            if lines.len() >= 2 { break; }
            if cur.is_empty() {
                cur = w.to_string();
            } else if cur.len() + 1 + w.len() <= max_chars {
                cur.push(' ');
                cur.push_str(w);
            } else {
                lines.push(cur.clone());
                cur = w.to_string();
            }
        }
        if !cur.is_empty() && lines.len() < 2 { lines.push(cur); }
        if lines.is_empty() { lines.push("genome".to_string()); }
        lines
    };
    let you_text = "you are";
    let you_w = you_text.len() as f64 * char_w;
    let you_x_right = mx + 3.0 * ddx;
    let you_x = if you_x_right + you_w <= x_range { you_x_right } else { mx - you_w - 3.0 * ddx };
    let you_y1 = my + line_h * 0.5;
    let you_y2 = my - line_h * 0.6;

    let canvas = Canvas::default()
        .x_bounds([-x_range, x_range])
        .y_bounds([-y_range, y_range])
        .background_color(Color::Rgb(0, 0, 0))
        .paint(|ctx| {
            for (color, pts) in &skew_batches {
                ctx.draw(&Points { coords: pts, color: *color });
            }
            for (color, pts) in &density_batches {
                ctx.draw(&Points { coords: pts, color: *color });
            }
            ctx.draw(&Points { coords: &tick_pts, color: Color::DarkGray });
            for lbl in &mb_labels {
                ctx.print(lbl.x, lbl.y, ratatui::text::Line::styled(
                    lbl.text.clone(), Style::default().fg(Color::Rgb(170, 170, 170)),
                ));
            }
            ctx.draw(&Points { coords: &track_pts, color: Color::Rgb(60, 60, 80) });
            if chromosome_active {
                ctx.print(mx - char_w * 0.5, my, ratatui::text::Line::styled(
                    "◆", Style::default().fg(Color::Yellow),
                ));
                ctx.print(you_x, you_y1, ratatui::text::Line::styled(
                    you_text, Style::default().fg(Color::Yellow),
                ));
                ctx.print(you_x, you_y2, ratatui::text::Line::styled(
                    "here", Style::default().fg(Color::Yellow),
                ));
            }
        });

    f.render_widget(canvas, area);

    // Overlay the genome name as a Paragraph centred on the inner circle.
    // ctx.print has coordinate rounding issues; Paragraph is reliable.
    {
        // Inner-circle diameter in terminal columns: R_IN / x_range * width
        let inner_cols = ((R_IN / 1.1) * area.width as f64 * 0.80) as u16;
        let inner_cols = inner_cols.max(8).min(area.width);
        let n_lines    = name_lines.len() as u16;
        let text_h     = n_lines;
        let cx = area.x + area.width  / 2;
        let cy = area.y + area.height / 2;
        let tx = cx.saturating_sub(inner_cols / 2);
        let ty = cy.saturating_sub(text_h / 2);
        let text_rect = Rect {
            x: tx, y: ty,
            width: inner_cols,
            height: text_h.max(1),
        };
        // Clip to canvas area
        let text_rect = Rect {
            x: text_rect.x.max(area.x),
            y: text_rect.y.max(area.y),
            width: text_rect.width.min(area.x + area.width - text_rect.x.max(area.x)),
            height: text_rect.height.min(area.y + area.height - text_rect.y.max(area.y)),
        };
        let para_lines: Vec<Line> = name_lines.iter().enumerate().map(|(i, s)| {
            let color = if i == 0 { Color::Rgb(200, 200, 230) } else { Color::Rgb(150, 150, 190) };
            Line::from(Span::styled(s.clone(), Style::default().fg(color).bg(Color::Black)))
        }).collect();
        f.render_widget(
            Paragraph::new(para_lines).alignment(ratatui::layout::Alignment::Center),
            text_rect,
        );
    }
}

/// Lay out and render one small circle map per plasmid in the spacer area.
fn draw_plasmid_maps(f: &mut Frame, app: &mut App, area: Rect) {
    if app.plasmids.is_empty() || area.width < 10 || area.height < 6 {
        app.plasmid_rects = Vec::new();
        return;
    }
    const MAP_W: u16 = 22;
    let n = ((area.width / MAP_W) as usize).min(app.plasmids.len());
    if n == 0 {
        app.plasmid_rects = Vec::new();
        return;
    }

    let constraints: Vec<Constraint> = (0..n)
        .map(|_| Constraint::Length(MAP_W))
        .collect();

    let rects = Layout::default()
        .direction(Direction::Horizontal)
        .constraints(constraints)
        .split(area);

    // Collect rects for click detection (only the plasmid rects, not the trailing Min(0))
    let plasmid_rects_vec: Vec<Rect> = (0..n).map(|i| rects[i]).collect();
    app.plasmid_rects = plasmid_rects_vec.clone();

    for i in 0..n {
        let is_active = app.active_genome == i + 1;
        // We need to pass plasmid data without holding a mutable borrow to app
        // So we read what we need and pass individual fields
        let rect = rects[i];
        draw_plasmid_minimap(f, app, i, rect, is_active);
    }
}

/// Small circular map for a single plasmid contig.
fn draw_plasmid_minimap(f: &mut Frame, app: &App, plasmid_idx: usize, area: Rect, is_active: bool) {
    let plasmid = &app.plasmids[plasmid_idx];
    let genome_size = plasmid.genome_size as f64;
    let w = area.width  as f64;
    let h = area.height as f64;
    if w == 0.0 || h == 0.0 || genome_size == 0.0 { return; }

    const CELL_RATIO: f64 = 2.15;
    let x_range = 1.1;
    let y_range = x_range * (CELL_RATIO * h) / w;

    const R_OUT:      f64 = 0.82;
    const R_IN:       f64 = R_OUT * 0.70;
    const R_MID:      f64 = (R_OUT + R_IN) * 0.5;
    const R_GC_INNER: f64 = R_OUT + 0.02;
    const R_GC_OUTER: f64 = R_OUT + 0.07;

    let dot_size = x_range / (2.0 * w);
    let ang_step = (dot_size / R_OUT).max(0.006);
    let rad_step = dot_size.max(0.010);

    let n_dens = 180usize;
    let mut plus_dens  = vec![0u32; n_dens];
    let mut minus_dens = vec![0u32; n_dens];
    for feat in &plasmid.features {
        let mid  = (feat.start as f64 + feat.end as f64) / 2.0;
        let frac = (mid / genome_size).clamp(0.0, 1.0 - 1e-9);
        let bin  = (frac * n_dens as f64) as usize;
        if feat.strand == '+' { plus_dens[bin]  += 1; }
        else                  { minus_dens[bin] += 1; }
    }
    let max_plus  = (*plus_dens .iter().max().unwrap_or(&1)).max(1) as f64;
    let max_minus = (*minus_dens.iter().max().unwrap_or(&1)).max(1) as f64;

    // GC skew ring
    let n_win = plasmid.gc_skew.len();
    let (gc_min, gc_max) = if n_win > 0 {
        let mn = plasmid.gc_skew.iter().cloned().fold(f64::INFINITY,     f64::min);
        let mx = plasmid.gc_skew.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (mn, mx)
    } else { (-1.0, 1.0) };
    let gc_range = (gc_max - gc_min).max(1e-9);

    let mut skew_pts: Vec<(Color, Vec<(f64,f64)>)> = {
        let n_ang = (2.0 * PI / ang_step).ceil() as usize;
        let mut map: std::collections::HashMap<(u8,u8,u8), Vec<(f64,f64)>> = Default::default();
        for i in 0..n_ang {
            let a    = i as f64 * ang_step;
            let frac = a / (2.0 * PI);
            let t = if n_win > 0 {
                let win = ((frac * n_win as f64) as usize).min(n_win - 1);
                ((plasmid.gc_skew[win] - gc_min) / gc_range * 2.0 - 1.0).clamp(-1.0, 1.0)
            } else { 0.0 };
            let (rv, gv, bv) = match skew_to_color(t) {
                Color::Rgb(r, g, b) => (r, g, b),
                _ => (80, 80, 80),
            };
            let (sa, ca) = (a.sin(), a.cos());
            let entry = map.entry((rv, gv, bv)).or_default();
            let mut r = R_GC_INNER;
            while r <= R_GC_OUTER { entry.push((r * sa, r * ca)); r += rad_step; }
        }
        let mut v: Vec<(Color, Vec<(f64,f64)>)> = map.into_iter()
            .map(|((r,g,b), pts)| (Color::Rgb(r,g,b), pts))
            .collect();
        v.sort_unstable_by_key(|(c,_)| match c { Color::Rgb(r,g,b) => (*r,*g,*b), _ => (0,0,0) });
        v
    };
    // Ensure sort is applied (redundant but explicit)
    skew_pts.sort_unstable_by_key(|(c,_)| match c { Color::Rgb(r,g,b) => (*r,*g,*b), _ => (0,0,0) });

    // Density bands
    let mut density_pts: Vec<(Color, Vec<(f64,f64)>)> = {
        let mut map: std::collections::HashMap<(u8,u8,u8), Vec<(f64,f64)>> = Default::default();
        for i in 0..n_dens {
            let a_start = i       as f64 / n_dens as f64 * 2.0 * PI;
            let a_end   = (i + 1) as f64 / n_dens as f64 * 2.0 * PI;
            let tp = plus_dens[i]  as f64 / max_plus;
            let tm = minus_dens[i] as f64 / max_minus;
            let cp = (
                (15.0 + tp * 57.0 ).round() as u8,
                (20.0 + tp * 190.0).round() as u8,
                (20.0 + tp * 160.0).round() as u8,
            );
            let cm = (
                (15.0 + tm * 225.0).round() as u8,
                (10.0 + tm * 138.0).round() as u8,
                ( 5.0 + tm *  45.0).round() as u8,
            );
            let mut a = a_start;
            while a < a_end {
                let (sa, ca) = (a.sin(), a.cos());
                let mut r = R_MID; while r <= R_OUT { map.entry(cp).or_default().push((r*sa, r*ca)); r += rad_step; }
                let mut r = R_IN;  while r <  R_MID { map.entry(cm).or_default().push((r*sa, r*ca)); r += rad_step; }
                a += ang_step;
            }
        }
        let mut v: Vec<(Color, Vec<(f64,f64)>)> = map.into_iter()
            .map(|((r,g,b), pts)| (Color::Rgb(r,g,b), pts))
            .collect();
        v.sort_unstable_by_key(|(c,_)| match c { Color::Rgb(r,g,b) => (*r,*g,*b), _ => (0,0,0) });
        v
    };
    density_pts.sort_unstable_by_key(|(c,_)| match c { Color::Rgb(r,g,b) => (*r,*g,*b), _ => (0,0,0) });

    // ── "You are here" marker for active plasmid ─────────────────────────────
    const GENE_OUT_P: f64 = 0.82 + 0.14;
    let mut p_track_pts: Vec<(f64,f64)> = Vec::new();
    let mut p_mx = 0.0_f64;
    let mut p_my = 0.0_f64;
    let mut p_you_x = 0.0_f64;
    let mut p_you_y1 = 0.0_f64;
    let mut p_you_y2 = 0.0_f64;
    let p_char_w = 2.0 * x_range / w;
    let p_ddx = x_range / w;
    if is_active {
        let n_track = (2.0 * std::f64::consts::PI / ang_step).ceil() as usize;
        for i in 0..n_track {
            let a = i as f64 * ang_step;
            p_track_pts.push((GENE_OUT_P * a.sin(), GENE_OUT_P * a.cos()));
        }
        let view_centre = (plasmid.view_start as f64 + plasmid.view_end as f64) / 2.0;
        let ma = (view_centre.min(plasmid.genome_size as f64) / plasmid.genome_size as f64) * 2.0 * std::f64::consts::PI;
        let (msa, mca) = (ma.sin(), ma.cos());
        p_mx = GENE_OUT_P * msa;
        p_my = GENE_OUT_P * mca;
        let you_w = "you are".len() as f64 * p_char_w;
        p_you_x = if p_mx + 3.0 * p_ddx + you_w <= x_range { p_mx + 3.0 * p_ddx } else { p_mx - you_w - 3.0 * p_ddx };
        let line_h_p = 2.0 * y_range / h;
        p_you_y1 = p_my + line_h_p * 0.5;
        p_you_y2 = p_my - line_h_p * 0.6;
    }

    // Short plasmid name for centre label
    let label = {
        let words: Vec<&str> = plasmid.name.split_whitespace().collect();
        // Try to find a word that looks like a plasmid name (after "plasmid")
        let mut short = String::new();
        for (i, &w) in words.iter().enumerate() {
            if w.to_ascii_lowercase() == "plasmid" {
                if let Some(&next) = words.get(i + 1) {
                    short = next.trim_end_matches(',').to_string();
                }
                break;
            }
        }
        if short.is_empty() {
            let s = words.first().copied().unwrap_or("plasmid");
            if s.len() > 10 { s[..10].to_string() } else { s.to_string() }
        } else {
            short
        }
    };
    let char_w = 2.0 * x_range / w;
    let label_x = -(label.len() as f64 * char_w) / 2.0;

    let label_color = if is_active { Color::Yellow } else { Color::Rgb(160, 160, 200) };

    // ── Position tick marks and labels around the circle ─────────────────────
    // Use ~6 ticks regardless of plasmid size; adapt unit (bp / kb / Mb).
    const POS_TICK_IN:  f64 = GENE_OUT_P + 0.01;
    const POS_TICK_OUT: f64 = GENE_OUT_P + 0.08;
    const POS_LABEL_R:  f64 = GENE_OUT_P + 0.14;

    let gsize = plasmid.genome_size;
    let pos_tick_size = nice_tick_spacing(gsize, 60);   // ~6 ticks around full circle

    // Compact label: unit chosen by tick size
    let fmt_pos = |pos: u64| -> String {
        if pos_tick_size >= 1_000_000 {
            format!("{}M", pos / 1_000_000)
        } else if pos_tick_size >= 1_000 {
            format!("{}k", pos / 1_000)
        } else {
            format!("{}", pos)
        }
    };

    let mut pos_tick_pts: Vec<(f64,f64)> = Vec::new();
    struct PosLbl { x: f64, y: f64, text: String, right_align: bool }
    let mut pos_lbls: Vec<PosLbl> = Vec::new();

    if pos_tick_size > 0 {
        let mut pos = 0u64;
        while pos < gsize {
            let a   = (pos as f64 / genome_size) * 2.0 * PI;
            let (sa, ca) = (a.sin(), a.cos());
            let mut r = POS_TICK_IN;
            while r <= POS_TICK_OUT { pos_tick_pts.push((r * sa, r * ca)); r += rad_step; }
            let text  = fmt_pos(pos);
            let lx    = POS_LABEL_R * sa;
            let ly    = POS_LABEL_R * ca;
            pos_lbls.push(PosLbl { x: lx, y: ly, text, right_align: sa < -0.15 });
            pos += pos_tick_size;
        }
    }

    let canvas = Canvas::default()
        .x_bounds([-x_range, x_range])
        .y_bounds([-y_range, y_range])
        .background_color(Color::Rgb(0, 0, 0))
        .paint(move |ctx| {
            for (color, pts) in &skew_pts {
                ctx.draw(&Points { coords: pts, color: *color });
            }
            for (color, pts) in &density_pts {
                ctx.draw(&Points { coords: pts, color: *color });
            }
            // Position ticks
            if !pos_tick_pts.is_empty() {
                ctx.draw(&Points { coords: &pos_tick_pts, color: Color::Rgb(130, 130, 160) });
            }
            // Position labels (right-aligned on left half of circle)
            for lbl in &pos_lbls {
                let x = if lbl.right_align {
                    lbl.x - lbl.text.len() as f64 * p_char_w
                } else {
                    lbl.x
                };
                ctx.print(x, lbl.y, ratatui::text::Line::styled(
                    lbl.text.clone(), Style::default().fg(Color::Rgb(130, 130, 160)),
                ));
            }
            if is_active {
                ctx.draw(&Points { coords: &p_track_pts, color: Color::Rgb(60, 60, 80) });
                ctx.print(p_mx - p_char_w * 0.5, p_my, ratatui::text::Line::styled(
                    "◆", Style::default().fg(Color::Yellow),
                ));
                ctx.print(p_you_x, p_you_y1, ratatui::text::Line::styled(
                    "you are", Style::default().fg(Color::Yellow),
                ));
                ctx.print(p_you_x, p_you_y2, ratatui::text::Line::styled(
                    "here", Style::default().fg(Color::Yellow),
                ));
            }
            ctx.print(label_x, 0.0, ratatui::text::Line::styled(
                label.clone(), Style::default().fg(label_color),
            ));
        });

    f.render_widget(canvas, area);
}

/// Draw a per-strand coverage bar chart.
/// `strand` is '+' (bars grow up, teal) or '-' (bars grow down, amber).
/// Coverage is drawn for the active genome's viewport.
fn draw_coverage_track(f: &mut Frame, app: &App, area: Rect, strand: char) {
    let cov = match &app.coverage {
        Some(c) => c,
        None => return,
    };
    let width  = area.width  as usize;
    let height = area.height as usize;
    if width == 0 || height == 0 { return; }
    f.render_widget(Clear, area);

    const LABEL_W: usize = 4;
    let feat_w = if width > LABEL_W { width - LABEL_W } else { 1 };

    let view_start  = app.active_view_start();
    let view_end    = app.active_view_end();
    let genome_size = app.active_genome_size();
    let span        = (view_end - view_start) as f64;
    if span <= 0.0 { return; }

    // Build per-column max coverage in the current viewport
    let data = if strand == '+' { &cov.plus } else { &cov.minus };
    let _bin_size = cov.bin_size as f64;

    let mut col_max = vec![0u32; feat_w];
    let wrap_end: u64 = if view_end > genome_size { view_end - genome_size } else { 0 };

    let col_for_global = |gpos: u64| -> Option<usize> {
        let eff = if gpos >= view_start {
            gpos as f64
        } else {
            gpos as f64 + genome_size as f64
        };
        let col = ((eff - view_start as f64) / span * (feat_w - 1) as f64).round() as i64;
        if col >= 0 && (col as usize) < feat_w { Some(col as usize) } else { None }
    };

    for (bin_idx, &count) in data.iter().enumerate() {
        if count == 0 { continue; }
        let gpos = bin_idx as u64 * cov.bin_size;
        // Normal zone
        if gpos >= view_start && gpos <= genome_size {
            if let Some(col) = col_for_global(gpos) {
                col_max[col] = col_max[col].max(count);
            }
        }
        // Wrapped zone
        if wrap_end > 0 && gpos <= wrap_end {
            if let Some(col) = col_for_global(gpos) {
                col_max[col] = col_max[col].max(count);
            }
        }
    }

    let max_val = col_max.iter().copied().max().unwrap_or(1).max(1) as f64;
    let log_max = (max_val + 1.0).ln();
    let levels = height;

    let (base_r, base_g, base_b, hi_r, hi_g, hi_b) = if strand == '+' {
        (15u8, 40u8, 50u8, 72u8, 210u8, 180u8) // teal
    } else {
        (40u8, 20u8, 10u8, 240u8, 148u8, 50u8)  // amber
    };

    // Format a coverage count compactly into ≤4 chars
    let fmt_cov = |n: f64| -> String {
        if n >= 1_000_000.0      { format!("{:.0}M", n / 1_000_000.0) }
        else if n >= 10_000.0    { format!("{:.0}k", n / 1_000.0) }
        else if n >= 1_000.0     { format!("{:.1}k", n / 1_000.0) }
        else                     { format!("{:.0}", n) }
    };

    // Row adjacent to gene track: show track label instead of y-axis value
    // cov+ track: bottom row (row = levels-1) is adjacent to gene track
    // cov- track: top row (row = 0) is adjacent to gene track
    let label_row = if strand == '+' { levels - 1 } else { 0 };
    let track_label = if strand == '+' { "cov+" } else { "cov-" };

    let mut lines: Vec<Line<'static>> = Vec::with_capacity(height);
    for row in 0..height {
        let mut spans: Vec<Span<'static>> = Vec::with_capacity(feat_w + 2);

        let label = if row == label_row {
            track_label.to_string()
        } else {
            let threshold = if strand == '+' {
                ((levels - row) as f64 / levels as f64 * log_max).exp() - 1.0
            } else {
                ((row + 1) as f64 / levels as f64 * log_max).exp() - 1.0
            };
            format!("{:>4}", fmt_cov(threshold.max(0.0)))
        };
        let label_color = if row == label_row {
            if strand == '+' { Color::Rgb(72, 210, 180) } else { Color::Rgb(240, 148, 50) }
        } else {
            Color::DarkGray
        };
        spans.push(Span::styled(label, Style::default().fg(label_color)));
        spans.push(Span::styled("│", Style::default().fg(Color::DarkGray)));

        for col in 0..feat_w {
            let fill_frac = (col_max[col] as f64 + 1.0).ln() / log_max;
            let filled_rows = (fill_frac * levels as f64).round() as usize;
            let is_filled = if strand == '+' {
                // bars grow up from bottom: bottom `filled_rows` rows are filled
                row >= levels.saturating_sub(filled_rows)
            } else {
                // bars grow down from top: top `filled_rows` rows are filled
                row < filled_rows
            };
            let t = col_max[col] as f64 / max_val;
            let (r, g, b) = if is_filled {
                (
                    (base_r as f64 + t * (hi_r as f64 - base_r as f64)).round() as u8,
                    (base_g as f64 + t * (hi_g as f64 - base_g as f64)).round() as u8,
                    (base_b as f64 + t * (hi_b as f64 - base_b as f64)).round() as u8,
                )
            } else {
                (base_r / 2, base_g / 2, base_b / 2)
            };
            spans.push(Span::styled(" ", Style::default().bg(Color::Rgb(r, g, b))));
        }
        lines.push(Line::from(spans));
    }

    f.render_widget(Paragraph::new(lines), area);
}

// Portrait legend: three vertical gradient bars with labels.
// Layout (14 cols): each column slot is 4 chars [3ch bar + 1ch sep], 2ch trailing.
// Top = high / dense / G-rich,  Bottom = low / sparse / C-rich.
fn draw_legend(f: &mut Frame, _app: &App, area: Rect) {
    let bg  = Color::Rgb(14, 14, 24);
    let dim = Color::Rgb(100, 100, 135);
    let total_h = area.height as usize;
    if total_h == 0 { return; }

    let legend_h = (total_h / 2).max(5);
    let gap  = || Span::styled(" ", Style::default().bg(bg));
    let tail = || Span::styled(" ", Style::default().bg(bg));
    let blank = || Line::from(Span::styled(
        " ".repeat(area.width as usize),
        Style::default().bg(bg),
    ));

    let col_p  = Color::Rgb(72,  210, 180); // + strand teal
    let col_m  = Color::Rgb(240, 148,  50); // - strand orange
    let col_gc = Color::Rgb(140, 140, 200); // GC label
    let col_g  = Color::Rgb(100, 210, 120); // G-rich green
    let col_c  = Color::Rgb(160,  80, 200); // C-rich purple

    let mut lines: Vec<Line<'static>> = Vec::with_capacity(total_h);

    let blank_rows = total_h.saturating_sub(legend_h);
    for _ in 0..blank_rows { lines.push(blank()); }

    // Layout: [4ch col][1ch gap] × 3 + [1ch tail] = 16
    // ── Header: two-row label ─────────────────────────────────────────────────
    lines.push(Line::from(vec![
        Span::styled("gene", Style::default().fg(dim   ).bg(bg)), gap(),
        Span::styled("gene", Style::default().fg(dim   ).bg(bg)), gap(),
        Span::styled("    ", Style::default().fg(dim   ).bg(bg)), tail(),
    ]));
    lines.push(Line::from(vec![
        Span::styled("cov+", Style::default().fg(col_p ).bg(bg)), gap(),
        Span::styled("cov-", Style::default().fg(col_m ).bg(bg)), gap(),
        Span::styled("skew", Style::default().fg(col_gc).bg(bg)), tail(),
    ]));
    // ── Hi labels ────────────────────────────────────────────────────────────
    lines.push(Line::from(vec![
        Span::styled(" hi ", Style::default().fg(col_p).bg(bg)), gap(),
        Span::styled(" hi ", Style::default().fg(col_m).bg(bg)), gap(),
        Span::styled(" +G ", Style::default().fg(col_g).bg(bg)), tail(),
    ]));

    // ── Gradient rows (4ch bar + 1ch gap per column) ─────────────────────────
    let bar_rows = legend_h.saturating_sub(5).max(1);
    for row in 0..bar_rows {
        let t = if bar_rows > 1 { 1.0 - row as f64 / (bar_rows - 1) as f64 } else { 1.0 };

        let (rp, gp, bp_c) = (
            (15.0 + t *  57.0).round() as u8,
            (20.0 + t * 190.0).round() as u8,
            (20.0 + t * 160.0).round() as u8,
        );
        let (rm, gm, bm) = (
            (15.0 + t * 225.0).round() as u8,
            (10.0 + t * 138.0).round() as u8,
            ( 5.0 + t *  45.0).round() as u8,
        );
        let (rs, gs, bs) = match skew_to_color(t * 2.0 - 1.0) {
            Color::Rgb(r, g, b) => (r, g, b),
            _ => (80, 80, 80),
        };

        lines.push(Line::from(vec![
            Span::styled("    ", Style::default().bg(Color::Rgb(rp, gp, bp_c))), gap(),
            Span::styled("    ", Style::default().bg(Color::Rgb(rm, gm, bm))),   gap(),
            Span::styled("    ", Style::default().bg(Color::Rgb(rs, gs, bs))),   tail(),
        ]));
    }

    // ── Lo labels ────────────────────────────────────────────────────────────
    lines.push(Line::from(vec![
        Span::styled(" lo ", Style::default().fg(dim  ).bg(bg)), gap(),
        Span::styled(" lo ", Style::default().fg(dim  ).bg(bg)), gap(),
        Span::styled(" +C ", Style::default().fg(col_c).bg(bg)), tail(),
    ]));

    while lines.len() < total_h { lines.push(blank()); }

    f.render_widget(Paragraph::new(lines).style(Style::default().bg(bg)), area);
}

fn draw_info_panel(f: &mut Frame, app: &App, area: Rect) {
    use ratatui::style::Modifier;

    let bg  = Color::Rgb(22, 22, 38);
    let dim = Color::Rgb(80, 80, 110);
    let shortcuts = Span::styled(
        "/:search  d:menu",
        Style::default().fg(Color::Rgb(60, 65, 95)).bg(bg),
    );

    let lines = if let Some(idx) = app.hovered {
        let feat = &app.active_features()[idx];
        let len_bp = feat.end.saturating_sub(feat.start) + 1;
        let kind = if feat.is_orf { "ORF" } else { "gene" };
        let locus = if feat.locus_tag.is_empty() { "—".to_string() } else { feat.locus_tag.clone() };

        let line1 = Line::from(vec![
            Span::styled(" name: ", Style::default().fg(dim).bg(bg)),
            Span::styled(feat.name.clone(), Style::default().fg(Color::Rgb(205, 214, 244)).bg(bg).add_modifier(Modifier::BOLD)),
            Span::styled("   locus: ", Style::default().fg(dim).bg(bg)),
            Span::styled(locus, Style::default().fg(Color::Rgb(180, 190, 230)).bg(bg)),
            Span::styled(format!("   {}", kind), Style::default().fg(dim).bg(bg)),
        ]);
        let line2 = Line::from(vec![
            Span::styled(" coords: ", Style::default().fg(dim).bg(bg)),
            Span::styled(
                format!("{}–{}", feat.start, feat.end),
                Style::default().fg(Color::Rgb(160, 200, 240)).bg(bg),
            ),
            Span::styled("   strand: ", Style::default().fg(dim).bg(bg)),
            Span::styled(
                feat.strand.to_string(),
                Style::default().fg(Color::Rgb(160, 200, 240)).bg(bg),
            ),
            Span::styled("   length: ", Style::default().fg(dim).bg(bg)),
            Span::styled(
                crate::core::format_bp(len_bp),
                Style::default().fg(Color::Rgb(160, 200, 240)).bg(bg),
            ),
        ]);
        vec![line1, line2]
    } else if let Some(map_idx) = app.hovered_map {
        // Hovering over a circle map — show the contig/genome name
        let name = if map_idx == 0 {
            app.genome_name.clone()
        } else {
            app.plasmids.get(map_idx - 1)
                .map(|p| p.name.clone())
                .unwrap_or_default()
        };
        let size = if map_idx == 0 {
            app.genome_size
        } else {
            app.plasmids.get(map_idx - 1).map(|p| p.genome_size).unwrap_or(0)
        };
        vec![
            Line::from(vec![
                Span::styled(" ", Style::default().bg(bg)),
                Span::styled(name, Style::default().fg(Color::Rgb(205, 214, 244)).bg(bg).add_modifier(Modifier::BOLD)),
            ]),
            Line::from(vec![
                Span::styled(
                    format!(" {}  ({})", if map_idx == 0 { "chromosome" } else { "plasmid" },
                        crate::core::format_bp(size)),
                    Style::default().fg(dim).bg(bg),
                ),
                Span::styled("  —  click to navigate", Style::default().fg(Color::Rgb(55, 60, 88)).bg(bg)),
            ]),
        ]
    } else if app.hovered_legend {
        vec![
            Line::from(Span::styled(
                " colour bars: gene density per 10 kb (+ strand, − strand) and GC skew",
                Style::default().fg(Color::Rgb(180, 190, 220)).bg(bg),
            )),
            Line::from(shortcuts),
        ]
    } else {
        let hint = if app.active_genome > 0 {
            if let Some(p) = app.plasmids.get(app.active_genome - 1) {
                let name = &p.name;
                let short = if name.len() > 40 { format!("{}…", &name[..39]) } else { name.clone() };
                format!(" viewing: {}", short)
            } else {
                " hover over a gene or feature for details".to_string()
            }
        } else {
            " hover over a gene or feature for details".to_string()
        };
        vec![
            Line::from(Span::styled(hint, Style::default().fg(dim).bg(bg))),
            Line::from(shortcuts),
        ]
    };

    f.render_widget(
        Paragraph::new(lines).style(Style::default().bg(bg)),
        area,
    );
}

fn draw_search_popup(f: &mut Frame, app: &mut App) {
    let results = &app.search_results;
    if results.is_empty() { return; }

    let n = results.len();
    let max_rows: usize = 12;
    let visible = n.min(max_rows);
    let w: u16 = 46;
    let h: u16 = visible as u16 + 2;

    let total = f.area();
    let x = total.width.saturating_sub(w + 1);
    let y = total.height.saturating_sub(h + 4); // above status + info
    let popup = Rect { x, y, width: w, height: h };
    app.search_popup_rect = popup;

    let feats = app.active_features();

    f.render_widget(Clear, popup);

    let bg  = Color::Rgb(20, 20, 38);
    let sel = Color::Rgb(60, 60, 110);
    let dim = Color::Rgb(90, 90, 130);
    let fg  = Color::Rgb(210, 215, 250);
    let hi  = Color::Rgb(137, 180, 250);

    // Scroll so selected item is visible
    let offset = if app.search_popup_idx >= max_rows {
        app.search_popup_idx + 1 - max_rows
    } else { 0 };

    let lines: Vec<Line> = (offset..offset + visible)
        .map(|i| {
            let feat_idx = results[i];
            let f = &feats[feat_idx];
            let name = if f.name.len() > 16 { format!("{:.16}", f.name) } else { f.name.clone() };
            let locus = if f.locus_tag.len() > 12 { format!("{:.12}", f.locus_tag) } else { f.locus_tag.clone() };
            let row_bg = if i == app.search_popup_idx { sel } else { bg };
            Line::from(vec![
                Span::styled(format!(" {:<16} ", name),  Style::default().fg(hi).bg(row_bg)),
                Span::styled(format!("{:<12} ", locus),  Style::default().fg(dim).bg(row_bg)),
                Span::styled(format!("{}{:>8}", f.strand, f.start), Style::default().fg(fg).bg(row_bg)),
            ])
        })
        .collect();

    let title = format!(" {} results — j/k Enter ", n);
    let block = Block::default()
        .borders(Borders::ALL)
        .title(title)
        .border_style(Style::default().fg(dim))
        .style(Style::default().bg(bg));

    f.render_widget(Paragraph::new(lines).block(block), popup);
}

fn draw_display_menu(f: &mut Frame, app: &mut App) {
    let has_coverage = app.coverage.is_some();
    let has_plasmids = !app.plasmids.is_empty();

    struct MenuItem { label: &'static str, checked: bool }
    let items: Vec<MenuItem> = {
        let o = &app.display_opts;
        let mut v = vec![
            MenuItem { label: "Chromosome map",  checked: o.show_chr_map },
            MenuItem { label: "Plasmid map(s)",  checked: o.show_plasmid_maps },
            MenuItem { label: "Colour bar",      checked: o.show_legend },
            MenuItem { label: "Gene tracks",     checked: o.show_gene_tracks },
            MenuItem { label: "Coverage tracks", checked: o.show_coverage },
        ];
        if !has_plasmids { v.remove(1); }
        if !has_coverage { v.pop(); }
        if app.protein.is_some() {
            v.push(MenuItem { label: "Structure panel", checked: true });
        }
        if app.msa.is_some() {
            v.push(MenuItem { label: "MSA panel",       checked: true });
        }
        v
    };

    let n   = items.len() as u16;
    let w: u16 = 24;
    let h: u16 = n + 2;

    let total = f.area();
    let x = total.width.saturating_sub(w + 1);
    let y = total.height.saturating_sub(h + 1);
    let popup = Rect { x, y, width: w, height: h };

    app.display_menu_rect = popup;
    f.render_widget(Clear, popup);

    let bg   = Color::Rgb(30, 30, 50);
    let sel  = Color::Rgb(70, 70, 120);
    let dim  = Color::Rgb(100, 100, 150);
    let tick = Color::Rgb(130, 210, 130);
    let fg   = Color::Rgb(210, 215, 250);

    let lines: Vec<Line> = items.iter().enumerate().map(|(i, it)| {
        let check = if it.checked { "✓" } else { " " };
        let row_bg = if i == app.display_menu_idx { sel } else { bg };
        Line::from(vec![
            Span::styled(" [", Style::default().fg(dim).bg(row_bg)),
            Span::styled(check, Style::default().fg(tick).bg(row_bg)),
            Span::styled("] ", Style::default().fg(dim).bg(row_bg)),
            Span::styled(it.label, Style::default().fg(fg).bg(row_bg)),
        ])
    }).collect();

    let block = Block::default()
        .borders(Borders::ALL)
        .title(" Display (d) ")
        .border_style(Style::default().fg(dim))
        .style(Style::default().bg(bg));

    f.render_widget(
        Paragraph::new(lines).block(block),
        popup,
    );
}

fn draw_protein_panel(f: &mut Frame, app: &App, area: Rect) {
    let panel = match &app.protein { Some(p) => p, None => return };

    f.render_widget(Clear, area);

    let border_color = if app.active_panel == crate::app::ActivePanel::Protein {
        Color::White
    } else if panel.folding {
        Color::Rgb(80, 130, 200)
    } else if panel.error.is_some() {
        Color::Rgb(180, 70, 70)
    } else if !panel.atoms.is_empty() {
        Color::Rgb(70, 110, 180)
    } else {
        Color::Rgb(45, 50, 75)
    };

    let title = format!(" {} — f:fold  a:save AA  n:save NT  Esc:close ", panel.gene_name);
    let block = Block::default()
        .borders(Borders::ALL)
        .title(title)
        .border_style(Style::default().fg(border_color))
        .style(Style::default().bg(Color::Black));

    let inner = block.inner(area);
    f.render_widget(block, area);
    if inner.width < 4 || inner.height < 2 { return; }

    if panel.folding {
        const SPIN: &[char] = &['⠋','⠙','⠹','⠸','⠼','⠴','⠦','⠧','⠇','⠏'];
        let spin = SPIN[(app.anim_tick as usize / 2) % SPIN.len()];
        let aa_len = panel.aa_seq.len();
        let line1 = Line::from(Span::styled(
            format!("  {}  Folding {}…", spin, panel.gene_name),
            Style::default().fg(Color::Rgb(100, 160, 230)).bg(Color::Black),
        ));
        let line2 = Line::from(Span::styled(
            format!("     {} amino acids", aa_len),
            Style::default().fg(Color::Rgb(60, 80, 130)).bg(Color::Black),
        ));
        let y_off = inner.height.saturating_sub(2) / 2;
        f.render_widget(
            Paragraph::new(vec![line1, line2]).style(Style::default().bg(Color::Black)),
            Rect { y: inner.y + y_off, height: inner.height.saturating_sub(y_off), ..inner },
        );
        return;
    }

    if let Some(ref err) = panel.error {
        let lines = vec![
            Line::from(Span::styled("  Error:", Style::default().fg(Color::Rgb(200, 80, 80)).bg(Color::Black))),
            Line::from(Span::styled(format!("  {}", err), Style::default().fg(Color::Rgb(220, 110, 110)).bg(Color::Black))),
        ];
        f.render_widget(Paragraph::new(lines).style(Style::default().bg(Color::Black)), inner);
        return;
    }

    if panel.atoms.is_empty() {
        let aa_len = panel.aa_seq.len();
        let lines = vec![
            Line::from(Span::styled(
                format!("  {}  ({} aa)", panel.gene_name, aa_len),
                Style::default().fg(Color::Rgb(180, 190, 230)).bg(Color::Black),
            )),
            Line::from(Span::styled(
                "  Press f to fold with minifold",
                Style::default().fg(Color::Rgb(70, 80, 120)).bg(Color::Black),
            )),
        ];
        let y_off = inner.height.saturating_sub(2) / 2;
        f.render_widget(
            Paragraph::new(lines).style(Style::default().bg(Color::Black)),
            Rect { y: inner.y + y_off, height: inner.height.min(2), ..inner },
        );
        return;
    }

    // Dark background
    let bg_rgb = Color::Rgb(18, 18, 30);
    f.render_widget(
        Paragraph::new("").style(Style::default().bg(bg_rgb)),
        inner,
    );
    // Render cached Kitty image if available
    if let Some(ref img) = panel.img_cache {
        let render_area = if inner.height > 1 {
            Rect { height: inner.height - 1, ..inner }
        } else {
            inner
        };
        if let Some(widget) = crate::pv::kitty_png::KittyPngImage::new(img, render_area) {
            f.render_widget(widget, render_area);
        }
    }
    // Status line at bottom of inner area
    let n_res = panel.atoms.len();
    let avg_plddt = if n_res > 0 {
        panel.atoms.iter().map(|a| a.plddt as f64).sum::<f64>() / n_res as f64
    } else {
        0.0
    };
    if inner.height > 1 {
        let status = format!("  {} residues   avg pLDDT {:.0}   drag to rotate", n_res, avg_plddt);
        let st_area = Rect { y: inner.y + inner.height - 1, height: 1, ..inner };
        f.render_widget(
            Paragraph::new(status).style(Style::default().fg(Color::Rgb(70, 80, 120)).bg(bg_rgb)),
            st_area,
        );
    }
}

fn draw_status(f: &mut Frame, app: &App, area: Rect) {
    // In search mode: show the query being typed
    let (msg, fg) = if app.search_mode {
        (format!(" / {}▌  (Enter to search, Esc to cancel)", app.search_query),
         Color::Rgb(250, 220, 100))
    } else {
        let view_start = app.active_view_start();
        let view_end   = app.active_view_end();
        let genome_size = app.active_genome_size();
        let span = view_end.saturating_sub(view_start);
        let span_kb = span as f64 / 1000.0;
        let n_genes = app.active_features()
            .iter()
            .filter(|f| {
                let wrap_end: u64 = if view_end > genome_size { view_end - genome_size } else { 0 };
                (f.end >= view_start && f.start <= genome_size) ||
                (wrap_end > 0 && f.start <= wrap_end)
            })
            .count();
        let genome_hint = if app.active_genome > 0 {
            format!(" [plasmid {}]", app.active_genome)
        } else {
            String::new()
        };
        let search_hint = if !app.search_results.is_empty() {
            format!("  [{}/{}] press n for next",
                app.search_idx + 1, app.search_results.len())
        } else {
            String::new()
        };
        let m = format!(
            " view: {}–{} | span: {:.1}kb | {} genes{}{}{}{}",
            view_start, view_end, span_kb, n_genes, genome_hint,
            if app.status_msg.is_empty() { "" } else { " | " },
            app.status_msg, search_hint
        );
        (m, Color::Rgb(205, 214, 244))
    };

    let status = Paragraph::new(msg)
        .style(Style::default().fg(fg).bg(Color::Rgb(49, 50, 68)));
    f.render_widget(status, area);
}

// ── MSA panel ────────────────────────────────────────────────────────────────

fn aa_color(c: char) -> Color {
    match c.to_ascii_uppercase() {
        'K' | 'R'                         => Color::Red,
        'A' | 'F' | 'I' | 'L' | 'M' | 'V' | 'W' => Color::Blue,
        'N' | 'Q' | 'S' | 'T'            => Color::Green,
        'H' | 'Y'                         => Color::Cyan,
        'C'                               => Color::LightRed,
        'D' | 'E'                         => Color::Magenta,
        'P'                               => Color::Yellow,
        'G'                               => Color::Rgb(255, 200, 150),
        '-' | '.'                         => Color::DarkGray,
        _                                 => Color::White,
    }
}

pub fn draw_msa_panel(f: &mut Frame, app: &mut App, area: Rect) {
    app.msa_panel_rect = area;
    f.render_widget(Clear, area);
    let panel = match &app.msa { Some(p) => p, None => return };
    let is_active = app.active_panel == crate::app::ActivePanel::Msa;

    let border_style = if is_active {
        Style::default().fg(Color::White)
    } else {
        Style::default().fg(Color::Rgb(50, 60, 80))
    };
    let block = Block::default()
        .title(format!(" MSA: {}  (Esc close  ←→↑↓ scroll  Tab switch) ", panel.gene_name))
        .borders(Borders::ALL)
        .border_style(border_style)
        .style(Style::default().bg(Color::Black));
    let inner = block.inner(area);
    f.render_widget(block, area);

    if panel.loading {
        const SPIN: &[char] = &['⠋','⠙','⠹','⠸','⠼','⠴','⠦','⠧','⠇','⠏'];
        let spin = SPIN[(app.anim_tick as usize / 2) % SPIN.len()];
        let line1 = Line::from(Span::styled(
            format!("  {}  Searching homologs of {}…", spin, panel.gene_name),
            Style::default().fg(Color::Rgb(80, 200, 180)).bg(Color::Black),
        ));
        let line2 = Line::from(Span::styled(
            "     DIAMOND → FAMSA alignment in progress",
            Style::default().fg(Color::Rgb(40, 100, 90)).bg(Color::Black),
        ));
        let y_off = inner.height.saturating_sub(2) / 2;
        f.render_widget(
            Paragraph::new(vec![line1, line2]).style(Style::default().bg(Color::Black)),
            Rect { y: inner.y + y_off, height: inner.height.saturating_sub(y_off), ..inner },
        );
        return;
    }
    if let Some(ref err) = panel.error {
        let lines = vec![
            Line::from(Span::styled("  Error:", Style::default().fg(Color::Red).bg(Color::Black))),
            Line::from(Span::styled(format!("  {}", err), Style::default().fg(Color::LightRed).bg(Color::Black))),
        ];
        f.render_widget(Paragraph::new(lines).style(Style::default().bg(Color::Black)), inner);
        return;
    }
    if panel.sequences.is_empty() {
        f.render_widget(
            Paragraph::new("  No homologs found").style(Style::default().fg(Color::DarkGray).bg(Color::Black)),
            inner,
        );
        return;
    }

    const NAME_W: u16 = 20;
    if inner.width <= NAME_W + 2 || inner.height < 2 { return; }
    let seq_rows  = inner.height.saturating_sub(1) as usize;
    let status_y  = inner.y + inner.height.saturating_sub(1);
    let name_area = Rect { x: inner.x, y: inner.y, width: NAME_W, height: inner.height.saturating_sub(1) };
    let seq_area  = Rect { x: inner.x + NAME_W, y: inner.y,
                           width: inner.width.saturating_sub(NAME_W), height: inner.height.saturating_sub(1) };

    // Names column
    let name_lines: Vec<Line> = panel.sequences.iter()
        .skip(panel.viewport_row)
        .take(seq_rows)
        .map(|(id, _)| {
            let s = if id.len() >= NAME_W as usize { &id[..NAME_W as usize - 1] } else { id.as_str() };
            Line::from(Span::styled(
                format!("{:<width$}", s, width = NAME_W as usize),
                Style::default().fg(Color::Rgb(100, 160, 230)).bg(Color::Black),
            ))
        })
        .collect();
    f.render_widget(Paragraph::new(name_lines).style(Style::default().bg(Color::Black)), name_area);

    // Sequence columns — Seaview-style AA colours
    let vis_cols = seq_area.width as usize;
    let seq_lines: Vec<Line> = panel.sequences.iter()
        .skip(panel.viewport_row)
        .take(seq_rows)
        .map(|(_, seq)| {
            let spans: Vec<Span> = seq.chars()
                .skip(panel.viewport_col)
                .take(vis_cols)
                .map(|c| Span::styled(
                    c.to_string(),
                    Style::default().fg(aa_color(c)).bg(Color::Black),
                ))
                .collect();
            Line::from(spans)
        })
        .collect();
    f.render_widget(Paragraph::new(seq_lines).style(Style::default().bg(Color::Black)), seq_area);

    // Status line
    let aln_len = panel.sequences.first().map(|(_, s)| s.len()).unwrap_or(0);
    let status = format!(
        "  col {}/{} | seq {}/{} | {} total",
        panel.viewport_col.saturating_add(1), aln_len,
        panel.viewport_row.saturating_add(1), panel.sequences.len(),
        panel.sequences.len(),
    );
    f.render_widget(
        Paragraph::new(Line::from(Span::styled(status, Style::default().fg(Color::Rgb(80,90,110)).bg(Color::Black)))),
        Rect { x: inner.x, y: status_y, width: inner.width, height: 1 },
    );
}
