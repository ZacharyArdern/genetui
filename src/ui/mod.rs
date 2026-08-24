use ratatui::{
    Frame,
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Style},
    widgets::{Block, Borders},
};

use crate::app::App;

mod genetracks;
mod genomemap;
mod coveragetracks;
mod proteinvis;
mod msa;
mod legend;
mod info;
mod search;

pub(super) const MINIMAP_W: u16 = 44;
pub(super) const MINIMAP_H: u16 = 22;

pub fn draw(f: &mut Frame, app: &mut App) {
    let size = f.area();

    const COV_H: u16 = 4;
    let cov_shown   = !app.coverages.is_empty() && app.display_opts.show_coverage && app.active_genome == 0;
    let n_cov       = if cov_shown { app.coverages.len() } else { 0 };
    let cov_extra   = if cov_shown { COV_H * 2 * n_cov as u16 } else { 0 };
    let fixed_rows  = genetracks::TRACK_PHYS as u16 + 1 + 2 + 1 + cov_extra;
    let minimap_h   = MINIMAP_H.min(size.height.saturating_sub(fixed_rows));
    let track_min   = genetracks::TRACK_PHYS as u16 + 1 + cov_extra.max(0);

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

    const LEGEND_W: u16 = 16;
    const PLASMID_W: u16 = 22;
    let maps_min_w = LEGEND_W + MINIMAP_W + PLASMID_W;
    let (protein_area_opt, maps_area) = if app.protein.is_some()
        && bottom_area.width > maps_min_w + 30
    {
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

    let legend_w  = if app.display_opts.show_legend  { LEGEND_W } else { 0 };
    let minimap_w = if app.display_opts.show_chr_map { MINIMAP_W.min(maps_area.width) } else { 0 };

    let legend_area  = Rect { x: maps_area.x, y: maps_area.y,
        width: legend_w, height: maps_area.height };

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

    // Paint the whole track area black first so any unused remainder rows don't ghost.
    f.render_widget(ratatui::widgets::Block::default().style(ratatui::style::Style::default().bg(Color::Black)), track_area);

    if cov_shown && app.display_opts.show_gene_tracks {
        // n_cov + tracks above genes, gene track, n_cov - tracks below
        let mut constraints: Vec<Constraint> = Vec::new();
        for _ in 0..n_cov { constraints.push(Constraint::Length(COV_H)); }
        constraints.push(Constraint::Length(genetracks::TRACK_PHYS as u16));
        for _ in 0..n_cov { constraints.push(Constraint::Length(COV_H)); }
        constraints.push(Constraint::Min(0));
        let track_split = Layout::default()
            .direction(Direction::Vertical)
            .constraints(constraints)
            .split(track_area);
        for i in 0..n_cov {
            coveragetracks::draw_coverage_track(f, app, track_split[i], '+', i);
        }
        genetracks::draw_gene_tracks(f, app, track_split[n_cov]);
        for i in 0..n_cov {
            coveragetracks::draw_coverage_track(f, app, track_split[n_cov + 1 + i], '-', i);
        }
    } else if cov_shown {
        let mut constraints: Vec<Constraint> = Vec::new();
        for _ in 0..n_cov { constraints.push(Constraint::Length(COV_H)); }
        for _ in 0..n_cov { constraints.push(Constraint::Length(COV_H)); }
        constraints.push(Constraint::Min(0));
        let track_split = Layout::default()
            .direction(Direction::Vertical)
            .constraints(constraints)
            .split(track_area);
        for i in 0..n_cov {
            coveragetracks::draw_coverage_track(f, app, track_split[i], '+', i);
        }
        for i in 0..n_cov {
            coveragetracks::draw_coverage_track(f, app, track_split[n_cov + i], '-', i);
        }
    } else if app.display_opts.show_gene_tracks {
        genetracks::draw_gene_tracks(f, app, track_area);
    }

    f.render_widget(Block::default().style(Style::default().bg(Color::Black)), bottom_area);

    if app.protein.is_some() && app.active_panel == crate::app::ActivePanel::Genome {
        let border_block = Block::default()
            .borders(Borders::ALL)
            .border_style(Style::default().fg(Color::White));
        f.render_widget(border_block, track_area);
    }

    if let Some(pa) = protein_area_opt { proteinvis::draw_protein_panel(f, app, pa); }
    if msa_open                           { msa::draw_msa_panel(f, app, msa_area); }
    if app.display_opts.show_chr_map      { genomemap::draw_minimap(f, app, minimap_area); }
    if app.display_opts.show_legend       { legend::draw_legend(f, app, legend_area); }
    if app.display_opts.show_plasmid_maps { genomemap::draw_plasmid_maps(f, app, plasmid_area); }
    info::draw_info_panel(f, app, info_area);
    info::draw_status(f, app, status_area);
    if app.search_popup_open   { search::draw_search_popup(f, app); }
    if app.display_menu_open   { search::draw_display_menu(f, app); }
    if app.search_menu_open    { search::draw_search_menu(f, app); }
    if app.blast_target_open   { search::draw_blast_target_menu(f, app); }
    if app.blast_file_open && !app.blast_completions.is_empty() {
        search::draw_blast_completions(f, app);
    }
    if app.upload_file_open {
        search::draw_upload_completions(f, app);
    }
}
