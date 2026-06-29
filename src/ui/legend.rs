use ratatui::{
    Frame,
    layout::Rect,
    style::{Color, Style},
    text::{Line, Span},
    widgets::Paragraph,
};

use crate::app::App;
use super::genomemap::skew_to_color;

pub(super) fn draw_legend(f: &mut Frame, _app: &App, area: Rect) {
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

    let col_p  = Color::Rgb(72,  210, 180);
    let col_m  = Color::Rgb(240, 148,  50);
    let col_gc = Color::Rgb(140, 140, 200);
    let col_g  = Color::Rgb(100, 210, 120);
    let col_c  = Color::Rgb(160,  80, 200);

    let mut lines: Vec<Line<'static>> = Vec::with_capacity(total_h);

    let blank_rows = total_h.saturating_sub(legend_h);
    for _ in 0..blank_rows { lines.push(blank()); }

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
    lines.push(Line::from(vec![
        Span::styled(" hi ", Style::default().fg(col_p).bg(bg)), gap(),
        Span::styled(" hi ", Style::default().fg(col_m).bg(bg)), gap(),
        Span::styled(" +G ", Style::default().fg(col_g).bg(bg)), tail(),
    ]));

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

    lines.push(Line::from(vec![
        Span::styled(" lo ", Style::default().fg(dim  ).bg(bg)), gap(),
        Span::styled(" lo ", Style::default().fg(dim  ).bg(bg)), gap(),
        Span::styled(" +C ", Style::default().fg(col_c).bg(bg)), tail(),
    ]));

    while lines.len() < total_h { lines.push(blank()); }

    f.render_widget(Paragraph::new(lines).style(Style::default().bg(bg)), area);
}
