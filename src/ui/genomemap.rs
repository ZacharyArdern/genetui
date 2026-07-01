use std::f64::consts::PI;

use ratatui::{
    Frame,
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Style},
    text::{Line, Span},
    widgets::{Paragraph, canvas::{Canvas, Points}},
};

use crate::app::App;
use crate::core::nice_tick_spacing;

pub(super) fn skew_to_color(t: f64) -> Color {
    let t = t.clamp(-1.0, 1.0);
    if t >= 0.0 {
        let r = (80.0 - t * 50.0).round() as u8;
        let g = (80.0 + t * 150.0).round() as u8;
        let b = (80.0 - t * 40.0).round() as u8;
        Color::Rgb(r, g, b)
    } else {
        let s = -t;
        let r = (80.0 + s * 140.0).round() as u8;
        let g = (80.0 - s * 55.0).round() as u8;
        let b = (80.0 + s * 175.0).round() as u8;
        Color::Rgb(r, g, b)
    }
}

pub(super) fn draw_minimap(f: &mut Frame, app: &App, area: Rect) {
    let genome_size = app.genome_size as f64;
    let w = area.width  as f64;
    let h = area.height as f64;
    if w == 0.0 || h == 0.0 { return; }

    const CELL_RATIO: f64 = 2.15;
    let x_range = 1.1;
    let y_range = x_range * (CELL_RATIO * h) / w;

    const R_OUT:      f64 = 0.82;
    const R_IN:       f64 = R_OUT * 0.70;
    const R_MID:      f64 = (R_OUT + R_IN) * 0.5;
    const R_GC_INNER: f64 = R_OUT + 0.02;
    const R_GC_OUTER: f64 = R_OUT + 0.07;
    const GENE_OUT:   f64 = R_OUT + 0.14;

    let dot_size = x_range / (2.0 * w);
    let ang_step = (dot_size / R_OUT).max(0.004);
    let rad_step = dot_size.max(0.008);

    let pos_angle = |pos: f64| -> f64 { (pos / genome_size) * 2.0 * PI };

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

    let mut skew_map: std::collections::HashMap<(u8,u8,u8), Vec<(f64,f64)>> = Default::default();
    {
        let n_win = app.gc_skew.len();
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
            let t = if n_win > 0 {
                let win = ((frac * n_win as f64) as usize).min(n_win - 1);
                ((app.gc_skew[win] - gc_min) / gc_range * 2.0 - 1.0).clamp(-1.0, 1.0)
            } else {
                0.0
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

    let mut density_map: std::collections::HashMap<(u8,u8,u8), Vec<(f64,f64)>> = Default::default();
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

    let mut track_pts: Vec<(f64,f64)> = Vec::new();
    {
        let n_track = (2.0 * PI / ang_step).ceil() as usize;
        for i in 0..n_track {
            let a = i as f64 * ang_step;
            track_pts.push((GENE_OUT * a.sin(), GENE_OUT * a.cos()));
        }
    }

    let chromosome_active = app.active_genome == 0;

    let view_centre = (app.view_start as f64 + app.view_end as f64) / 2.0;
    let ma = pos_angle(view_centre.min(genome_size));
    let (msa, mca) = (ma.sin(), ma.cos());
    let mx = GENE_OUT * msa;
    let my = GENE_OUT * mca;
    let ddx = x_range / w;

    let char_w = 2.0 * x_range / w;
    let line_h = 2.0 * y_range / h;

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
                ctx.print(lbl.x, lbl.y, Line::styled(
                    lbl.text.clone(), Style::default().fg(Color::Rgb(170, 170, 170)),
                ));
            }
            ctx.draw(&Points { coords: &track_pts, color: Color::Rgb(60, 60, 80) });
            if chromosome_active {
                ctx.print(mx - char_w * 0.5, my + line_h * 0.35, Line::styled(
                    "◆", Style::default().fg(Color::Yellow),
                ));
                ctx.print(you_x, you_y1, Line::styled(
                    you_text, Style::default().fg(Color::Yellow),
                ));
                ctx.print(you_x, you_y2, Line::styled(
                    "here", Style::default().fg(Color::Yellow),
                ));
            }
        });

    f.render_widget(canvas, area);

    {
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

pub(super) fn draw_plasmid_maps(f: &mut Frame, app: &mut App, area: Rect) {
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

    let plasmid_rects_vec: Vec<Rect> = (0..n).map(|i| rects[i]).collect();
    app.plasmid_rects = plasmid_rects_vec.clone();

    for i in 0..n {
        let is_active = app.active_genome == i + 1;
        let rect = rects[i];
        draw_plasmid_minimap(f, app, i, rect, is_active);
    }
}

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
    skew_pts.sort_unstable_by_key(|(c,_)| match c { Color::Rgb(r,g,b) => (*r,*g,*b), _ => (0,0,0) });

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
        let n_track = (2.0 * PI / ang_step).ceil() as usize;
        for i in 0..n_track {
            let a = i as f64 * ang_step;
            p_track_pts.push((GENE_OUT_P * a.sin(), GENE_OUT_P * a.cos()));
        }
        let view_centre = (plasmid.view_start as f64 + plasmid.view_end as f64) / 2.0;
        let ma = (view_centre.min(plasmid.genome_size as f64) / plasmid.genome_size as f64) * 2.0 * PI;
        let (msa, mca) = (ma.sin(), ma.cos());
        p_mx = GENE_OUT_P * msa;
        p_my = GENE_OUT_P * mca;
        let you_w = "you are".len() as f64 * p_char_w;
        p_you_x = if p_mx + 3.0 * p_ddx + you_w <= x_range { p_mx + 3.0 * p_ddx } else { p_mx - you_w - 3.0 * p_ddx };
        let line_h_p = 2.0 * y_range / h;
        p_you_y1 = p_my + line_h_p * 0.5;
        p_you_y2 = p_my - line_h_p * 0.6;
    }

    let label = {
        let words: Vec<&str> = plasmid.name.split_whitespace().collect();
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

    const POS_TICK_IN:  f64 = GENE_OUT_P + 0.01;
    const POS_TICK_OUT: f64 = GENE_OUT_P + 0.08;
    const POS_LABEL_R:  f64 = GENE_OUT_P + 0.14;

    let gsize = plasmid.genome_size;
    let pos_tick_size = nice_tick_spacing(gsize, 60);

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
            if !pos_tick_pts.is_empty() {
                ctx.draw(&Points { coords: &pos_tick_pts, color: Color::Rgb(130, 130, 160) });
            }
            for lbl in &pos_lbls {
                let x = if lbl.right_align {
                    lbl.x - lbl.text.len() as f64 * p_char_w
                } else {
                    lbl.x
                };
                ctx.print(x, lbl.y, Line::styled(
                    lbl.text.clone(), Style::default().fg(Color::Rgb(130, 130, 160)),
                ));
            }
            if is_active {
                ctx.draw(&Points { coords: &p_track_pts, color: Color::Rgb(60, 60, 80) });
                let p_line_h = 2.0 * y_range / h;
                ctx.print(p_mx - p_char_w * 0.5, p_my + p_line_h * 0.35, Line::styled(
                    "◆", Style::default().fg(Color::Yellow),
                ));
                ctx.print(p_you_x, p_you_y1, Line::styled(
                    "you are", Style::default().fg(Color::Yellow),
                ));
                ctx.print(p_you_x, p_you_y2, Line::styled(
                    "here", Style::default().fg(Color::Yellow),
                ));
            }
            ctx.print(label_x, 0.0, Line::styled(
                label.clone(), Style::default().fg(label_color),
            ));
        });

    f.render_widget(canvas, area);
}
