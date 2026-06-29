use ratatui::{
    Frame,
    layout::Rect,
    style::{Color, Style},
    text::{Line, Span},
    widgets::{Clear, Paragraph},
};

use crate::app::App;
use crate::core::{format_bp_tick, nice_tick_spacing};

pub(super) const PLUS_COLORS: &[(u8, u8, u8)] = &[
    (72,  210, 180),
    (80,  170, 230),
    (100, 220, 150),
    (55,  225, 200),
    (95,  185, 240),
    (120, 215, 170),
];

pub(super) const MINUS_COLORS: &[(u8, u8, u8)] = &[
    (240, 148,  50),
    (225, 118,  70),
    (245, 175,  45),
    (230, 132,  80),
    (238, 160,  55),
    (220, 105,  75),
];

pub(super) const TRACK_ROWS: usize = 9;
pub(super) const RULER_ROW: usize = 4;
pub(super) const TRACK_PHYS: usize = 15;

pub(super) const FRAME_BG: &[(u8, u8, u8)] = &[
    (18,  25,  48),
    (38,  30,  16),
    (30,  16,  40),
];

pub(super) struct VFeat {
    pub idx: usize,
    pub cs: usize,
    pub ce: usize,
    pub fr: usize,
    pub strand: char,
    pub name: String,
    pub color_idx: usize,
    #[allow(dead_code)]
    pub is_orf: bool,
    pub noncoding: bool,
}

pub(super) fn draw_gene_tracks(f: &mut Frame, app: &mut App, area: Rect) {
    let width = area.width as usize;
    let height = area.height as usize;

    if width == 0 || height == 0 {
        return;
    }

    f.render_widget(Clear, area);

    const LABEL_W: usize = 4;
    let feat_w = if width > LABEL_W { width - LABEL_W } else { 1 };

    let view_start  = app.active_view_start();
    let view_end    = app.active_view_end();
    let genome_size = app.active_genome_size();
    let span        = (view_end - view_start) as f64;
    let wrap_end: u64 = if view_end > genome_size { view_end - genome_size } else { 0 };

    let bp_to_col = |bp: u64| -> usize {
        if span <= 0.0 { return 0; }
        let eff = if bp >= view_start {
            bp as f64
        } else {
            bp as f64 + genome_size as f64
        };
        let col = ((eff - view_start as f64) / span * (feat_w - 1) as f64).round() as i64;
        col.max(0).min(feat_w as i64 - 1) as usize
    };

    let vfeats: Vec<VFeat> = {
        let feats = app.active_features();
        let mut out = Vec::new();
        for (idx, f) in feats.iter().enumerate() {
            let fr = (if f.strand == '+' {
                f.start.saturating_sub(1) % 3
            } else {
                f.end.saturating_sub(1) % 3
            }) as usize;
            if f.end >= view_start && f.start <= genome_size {
                let cs = bp_to_col(f.start.max(view_start));
                let ce = bp_to_col(f.end.min(genome_size));
                out.push(VFeat { idx, cs, ce: ce.max(cs), fr, strand: f.strand,
                    name: f.name.clone(), color_idx: f.color_idx, is_orf: f.is_orf,
                    noncoding: f.noncoding });
            }
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

    let stop_cols: std::collections::HashMap<(char, u8), Vec<usize>> = if app.active_genome == 0 {
        let mut map = std::collections::HashMap::new();
        for &s in &['+', '-'] {
            for fr in 0u8..3 {
                let positions = app.stop_codons.get(&(s, fr)).cloned().unwrap_or_default();
                let cols: Vec<usize> = positions
                    .iter()
                    .filter(|&&p| {
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

    let basic = app.basic_mode;

    let render_nc = |lines: &mut Vec<Line>,
                     new_hit_map: &mut crate::app::HitMap,
                     nc_list: &[usize],
                     strand: char,
                     phys_row: usize| {
        let label_color = if basic {
            if strand == '+' { Color::Cyan } else { Color::Yellow }
        } else if strand == '+' { Color::Rgb(137, 180, 250) } else { Color::Rgb(243, 139, 168) };
        let nc_bg = if basic { Color::Black } else { Color::Rgb(15, 15, 25) };
        let nc_col = Color::White;
        let nc_fg  = Color::Rgb(10, 10, 10);
        let mut cells: Vec<(char, Color, Option<Color>)> = vec![(' ', Color::Reset, None); feat_w];
        for &vi in nc_list {
            let vf = &vfeats[vi];
            let (cs, ce) = (vf.cs, vf.ce);
            let fw = (ce + 1).saturating_sub(cs).max(1);
            let tip = if basic {
                if strand == '+' { '>' } else { '<' }
            } else if strand == '+' { '▶' } else { '◀' };
            for j in 0..fw {
                if cs + j < feat_w {
                    let ch = if strand == '+' && j == fw - 1 { tip }
                             else if strand == '-' && j == 0  { tip }
                             else { ' ' };
                    cells[cs + j] = (ch, nc_fg, Some(nc_col));
                }
            }
            let name_chars: Vec<char> = vf.name.chars().collect();
            if name_chars.len() <= fw {
                let off = (fw - name_chars.len()) / 2;
                for (j, &ch) in name_chars.iter().enumerate() {
                    let col = cs + off + j;
                    if col < feat_w { cells[col] = (ch, nc_fg, Some(nc_col)); }
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
            lines.push(make_ruler(view_start, view_end, genome_size, feat_w, LABEL_W, basic));
            phys_row += 1;
            continue;
        }

        let (strand, fr) = if row_idx < RULER_ROW {
            ('+', (row_idx - 1) as u8)
        } else {
            ('-', (row_idx - RULER_ROW - 1) as u8)
        };

        let label_color = if basic {
            if strand == '+' { Color::Cyan } else { Color::Yellow }
        } else if strand == '+' {
            Color::Rgb(137, 180, 250)
        } else {
            Color::Rgb(243, 139, 168)
        };

        let vi_list = buckets.get(&(strand, fr)).cloned().unwrap_or_default();
        let stops = stop_cols.get(&(strand, fr)).cloned().unwrap_or_default();

        let frame_bg = if basic {
            Color::Black
        } else {
            let (fbr, fbg, fbb) = FRAME_BG[fr as usize % FRAME_BG.len()];
            Color::Rgb(fbr, fbg, fbb)
        };

        let mut cells_arrow: Vec<(char, Color, Option<Color>)> = vec![(' ', Color::Reset, None); feat_w];
        let mut cells_gap:   Vec<(char, Color, Option<Color>)> = vec![(' ', Color::Reset, None); feat_w];

        for &col in &stops {
            let stop_color = if basic { Color::Red } else { Color::Rgb(204, 68, 68) };
            if col < feat_w { cells_gap[col] = ('|', stop_color, None); }
        }

        for &vi in &vi_list {
            let vf = &vfeats[vi];
            let cs = vf.cs;
            let ce = vf.ce;
            let fw = (ce + 1).saturating_sub(cs).max(1);
            let is_selected = selected_feat == Some(vf.idx);
            let is_hovered  = hovered_feat  == Some(vf.idx);
            let gene_col = if basic {
                if is_selected { Color::White }
                else if is_hovered {
                    if vf.strand == '+' { Color::LightCyan } else { Color::LightYellow }
                } else {
                    if vf.strand == '+' { Color::Cyan } else { Color::Yellow }
                }
            } else {
                let (r, g, b) = if vf.strand == '+' {
                    PLUS_COLORS[vf.color_idx % PLUS_COLORS.len()]
                } else {
                    MINUS_COLORS[vf.color_idx % MINUS_COLORS.len()]
                };
                if is_selected { Color::White }
                else if is_hovered {
                    let bri = |x: u8| -> u8 { 255u8.min(x.saturating_add(70)) };
                    Color::Rgb(bri(r), bri(g), bri(b))
                } else {
                    Color::Rgb(r, g, b)
                }
            };
            let tip = if basic {
                if vf.strand == '+' { '>' } else { '<' }
            } else {
                if vf.strand == '+' { '▶' } else { '◀' }
            };
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

pub(super) fn cells_to_line(
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

pub(super) fn make_ruler(view_start: u64, view_end: u64, genome_size: u64, feat_w: usize, label_w: usize, basic: bool) -> Line<'static> {
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

    if tick_bp > 0 {
        let first_tick = (view_start / tick_bp + 1) * tick_bp;
        let mut tick = first_tick;
        let normal_end = view_end.min(genome_size);
        while tick <= normal_end {
            place_tick(tick, &mut chars);
            tick += tick_bp;
        }
        if wrap_end > 0 {
            let boundary_col = bp_to_col(genome_size);
            if boundary_col > 0 && boundary_col < feat_w { chars[boundary_col] = '│'; }
            let mut tick_w = tick_bp;
            while tick_w <= wrap_end {
                place_tick(tick_w, &mut chars);
                tick_w += tick_bp;
            }
        }
    }

    let ruler_str: String = chars.into_iter().collect();
    let label_dashes = "─".repeat(label_w);

    let ruler_color = if basic { Color::White } else { Color::Rgb(170, 170, 204) };
    Line::from(vec![
        Span::styled(label_dashes, Style::default().fg(Color::DarkGray)),
        Span::styled(ruler_str, Style::default().fg(ruler_color)),
    ])
}
