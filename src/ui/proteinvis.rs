use ratatui::{
    Frame,
    layout::Rect,
    style::{Color, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, Paragraph},
};

use crate::app::App;

pub(super) fn draw_protein_panel(f: &mut Frame, app: &App, area: Rect) {
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

    let n_res = panel.atoms.len();
    let avg_plddt = if n_res > 0 {
        panel.atoms.iter().map(|a| a.plddt as f64).sum::<f64>() / n_res as f64
    } else { 0.0 };
    let bg_rgb = Color::Rgb(18, 18, 30);

    let render_area = if inner.height > 1 {
        Rect { height: inner.height - 1, ..inner }
    } else {
        inner
    };

    if app.kitty_native {
        // Full Kitty terminal: render high-res image via unicode placeholder widget.
        f.render_widget(Paragraph::new("").style(Style::default().bg(bg_rgb)), inner);
        if let Some(ref img) = panel.img_cache {
            if let Some(widget) = crate::pv::kitty_png::KittyPngImage::new(img, render_area) {
                f.render_widget(widget, render_area);
            }
        }
    } else {
        // Braille rendering: works in any terminal including VSCode.
        f.render_widget(Paragraph::new("").style(Style::default().bg(Color::Black)), inner);
        let braille = crate::pv::braille::render_braille(
            &panel.atoms,
            &panel.camera,
            render_area,
        );
        f.render_widget(braille.style(Style::default().bg(Color::Black)), render_area);
    }

    if inner.height > 1 {
        let status = format!("  {} residues   avg pLDDT {:.0}   drag to rotate", n_res, avg_plddt);
        let st_area = Rect { y: inner.y + inner.height - 1, height: 1, ..inner };
        f.render_widget(
            Paragraph::new(status).style(Style::default().fg(Color::Rgb(70, 80, 120)).bg(bg_rgb)),
            st_area,
        );
    }
}
