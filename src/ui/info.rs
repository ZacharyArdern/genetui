use ratatui::{
    Frame,
    layout::Rect,
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::Paragraph,
};

use crate::app::App;

pub(super) fn draw_info_panel(f: &mut Frame, app: &App, area: Rect) {
    let bg   = if app.basic_mode { Color::Black } else { Color::Rgb(22, 22, 38) };
    let dim  = if app.basic_mode { Color::DarkGray } else { Color::Rgb(80, 80, 110) };
    let bold = if app.basic_mode { Color::White } else { Color::Rgb(205, 214, 244) };
    let hi   = if app.basic_mode { Color::White } else { Color::Rgb(160, 200, 240) };
    let shortcuts = Line::from(vec![
        Span::styled(" /:search  ", Style::default().fg(Color::White).bg(bg)),
        Span::styled("d:menu  ", Style::default().fg(Color::White).bg(bg)),
        Span::styled("w:browser app  ", Style::default().fg(Color::White).bg(bg)),
        Span::styled("q:quit", Style::default().fg(Color::White).bg(bg)),
    ]);

    let lines = if let Some(idx) = app.hovered {
        let feat = &app.active_features()[idx];
        let len_bp = feat.end.saturating_sub(feat.start) + 1;
        let kind = if feat.is_orf { "ORF" } else { "gene" };
        let locus = if feat.locus_tag.is_empty() { "—".to_string() } else { feat.locus_tag.clone() };

        let line1 = Line::from(vec![
            Span::styled(" name: ", Style::default().fg(dim).bg(bg)),
            Span::styled(feat.name.clone(), Style::default().fg(bold).bg(bg).add_modifier(Modifier::BOLD)),
            Span::styled("   locus: ", Style::default().fg(dim).bg(bg)),
            Span::styled(locus, Style::default().fg(hi).bg(bg)),
            Span::styled(format!("   {}", kind), Style::default().fg(dim).bg(bg)),
        ]);
        let line2 = Line::from(vec![
            Span::styled(" coords: ", Style::default().fg(dim).bg(bg)),
            Span::styled(format!("{}–{}", feat.start, feat.end), Style::default().fg(hi).bg(bg)),
            Span::styled("   strand: ", Style::default().fg(dim).bg(bg)),
            Span::styled(feat.strand.to_string(), Style::default().fg(hi).bg(bg)),
            Span::styled("   length: ", Style::default().fg(dim).bg(bg)),
            Span::styled(crate::core::format_bp(len_bp), Style::default().fg(hi).bg(bg)),
        ]);
        vec![line1, line2]
    } else if let Some(map_idx) = app.hovered_map {
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
                Span::styled(name, Style::default().fg(bold).bg(bg).add_modifier(Modifier::BOLD)),
            ]),
            Line::from(vec![
                Span::styled(
                    format!(" {}  ({})", if map_idx == 0 { "chromosome" } else { "plasmid" },
                        crate::core::format_bp(size)),
                    Style::default().fg(dim).bg(bg),
                ),
                Span::styled("  —  click to navigate", Style::default().fg(dim).bg(bg)),
            ]),
        ]
    } else if app.hovered_legend {
        vec![
            Line::from(Span::styled(
                " colour bars: gene density per 10 kb (+ strand, − strand) and GC skew",
                Style::default().fg(hi).bg(bg),
            )),
            shortcuts.clone(),
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
            shortcuts.clone(),
        ]
    };

    f.render_widget(
        Paragraph::new(lines).style(Style::default().bg(bg)),
        area,
    );
}

pub(super) fn draw_status(f: &mut Frame, app: &App, area: Rect) {
    let (msg, fg) = if app.blast_file_open {
        let target = if app.blast_target_idx == 0 { "6FT" } else { "GFF proteins" };
        (format!(" DIAMOND ({}) query: {}▌  (Tab: complete, Enter: run, Esc: cancel)",
            target, app.blast_file_path),
         Color::Rgb(230, 180, 60))
    } else if app.blast_running {
        let dots = match (app.anim_tick / 4) % 3 { 0 => ".", 1 => "..", _ => "..." };
        (format!(" DIAMOND running{}", dots), Color::Rgb(230, 180, 60))
    } else if app.search_mode {
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
