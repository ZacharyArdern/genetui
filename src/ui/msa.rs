use ratatui::{
    Frame,
    layout::Rect,
    style::{Color, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, Paragraph},
};

use crate::app::App;

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

pub(super) fn draw_msa_panel(f: &mut Frame, app: &mut App, area: Rect) {
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
