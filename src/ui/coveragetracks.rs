use ratatui::{
    Frame,
    layout::Rect,
    style::{Color, Style},
    text::{Line, Span},
    widgets::{Clear, Paragraph},
};

use crate::app::App;

// (base_r, base_g, base_b, hi_r, hi_g, hi_b) per track index
const PLUS_PALETTE: &[(u8,u8,u8,u8,u8,u8)] = &[
    (15, 40, 50,  72, 210, 180),  // teal
    (10, 40, 15,  50, 200,  80),  // green
    (30, 10, 50, 140,  60, 210),  // purple
    (40, 10, 10, 220,  70,  70),  // red
];
const MINUS_PALETTE: &[(u8,u8,u8,u8,u8,u8)] = &[
    (40, 20, 10, 240, 148,  50),  // orange
    (40, 10, 30, 220,  60, 150),  // pink
    (40, 35,  5, 200, 180,  30),  // yellow
    (10, 30, 40,  40, 160, 200),  // blue
];
const BASIC_PLUS:  &[Color] = &[Color::Cyan,    Color::Green,   Color::Blue,  Color::White];
const BASIC_MINUS: &[Color] = &[Color::Yellow,  Color::Magenta, Color::Red,   Color::LightYellow];

pub(super) fn draw_coverage_track(f: &mut Frame, app: &App, area: Rect, strand: char, track_idx: usize) {
    let cov = match app.coverages.get(track_idx) {
        Some((_, c)) => c,
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

    let data = if strand == '+' { &cov.plus } else { &cov.minus };

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
        if gpos >= view_start && gpos <= genome_size {
            if let Some(col) = col_for_global(gpos) {
                col_max[col] = col_max[col].max(count);
            }
        }
        if wrap_end > 0 && gpos <= wrap_end {
            if let Some(col) = col_for_global(gpos) {
                col_max[col] = col_max[col].max(count);
            }
        }
    }

    let max_val = col_max.iter().copied().max().unwrap_or(1).max(1) as f64;
    let log_max = (max_val + 1.0).ln();
    let levels = height;

    let basic = app.basic_mode;
    let ti = track_idx % PLUS_PALETTE.len();
    let (base_r, base_g, base_b, hi_r, hi_g, hi_b) = if strand == '+' {
        PLUS_PALETTE[ti]
    } else {
        MINUS_PALETTE[ti]
    };
    let basic_hi = if strand == '+' { BASIC_PLUS[ti] } else { BASIC_MINUS[ti] };

    let fmt_cov = |n: f64| -> String {
        if n >= 1_000_000.0      { format!("{:.0}M", n / 1_000_000.0) }
        else if n >= 10_000.0    { format!("{:.0}k", n / 1_000.0) }
        else if n >= 1_000.0     { format!("{:.1}k", n / 1_000.0) }
        else                     { format!("{:.0}", n) }
    };

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
        let hi_color = if basic {
            basic_hi
        } else {
            Color::Rgb(hi_r, hi_g, hi_b)
        };
        let label_color = if row == label_row { hi_color } else { Color::DarkGray };
        spans.push(Span::styled(label, Style::default().fg(label_color)));
        spans.push(Span::styled("│", Style::default().fg(Color::DarkGray)));

        for col in 0..feat_w {
            let fill_frac = (col_max[col] as f64 + 1.0).ln() / log_max;
            let filled_rows = (fill_frac * levels as f64).round() as usize;
            let is_filled = if strand == '+' {
                row >= levels.saturating_sub(filled_rows)
            } else {
                row < filled_rows
            };
            if basic {
                let (ch, color) = if is_filled {
                    ('█', basic_hi)
                } else {
                    (' ', Color::Reset)
                };
                spans.push(Span::styled(ch.to_string(), Style::default().fg(color)));
            } else {
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
        }
        lines.push(Line::from(spans));
    }

    f.render_widget(Paragraph::new(lines), area);
}
