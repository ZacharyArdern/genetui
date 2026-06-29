use ratatui::{
    Frame,
    layout::Rect,
    style::{Color, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, Paragraph},
};

use crate::app::App;

pub(super) fn draw_search_popup(f: &mut Frame, app: &mut App) {
    let results = &app.search_results;
    if results.is_empty() { return; }

    let n = results.len();
    let max_rows: usize = 12;
    let visible = n.min(max_rows);
    let w: u16 = 46;
    let h: u16 = visible as u16 + 2;

    let total = f.area();
    let x = total.width.saturating_sub(w + 1);
    let y = total.height.saturating_sub(h + 4);
    let popup = Rect { x, y, width: w, height: h };
    app.search_popup_rect = popup;

    let feats = app.active_features();

    f.render_widget(Clear, popup);

    let bg  = Color::Rgb(20, 20, 38);
    let sel = Color::Rgb(60, 60, 110);
    let dim = Color::Rgb(90, 90, 130);
    let fg  = Color::Rgb(210, 215, 250);
    let hi  = Color::Rgb(137, 180, 250);

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

pub(super) fn draw_display_menu(f: &mut Frame, app: &mut App) {
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
