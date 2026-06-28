use genetui::core;
use genetui::app;
use genetui::ui;
use genetui::server;

use std::collections::HashMap;
use std::io;
use std::time::Duration;

use anyhow::Result;
use clap::Parser;
use crossterm::{
    event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyCode, MouseEventKind},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{backend::CrosstermBackend, Terminal};
use tokio::process::Command as TokioCommand;

use app::{App, PlasmidData};
use core::{compute_alignment_coverage, compute_gc_skew, find_stop_codons, parse_fasta, parse_gff, translate};

#[derive(Parser, Debug)]
#[command(name = "genetui", about = "Terminal genome browser")]
struct Cli {
    /// Input files: GFF/GFF3, FASTA (.fa/.fasta/.fna), BAM/SAM/CRAM — in any order
    #[arg(required = true)]
    files: Vec<String>,

    /// Path to MiniFold-MLX repo directory (uses fold.py inside); falls back to `minifold` in PATH
    #[arg(long)]
    minifold_mlx: Option<String>,

    /// Directory to store/cache PDB files (default: current directory)
    #[arg(long)]
    pdbs: Option<String>,

    /// Alias for --pdbs (legacy)
    #[arg(long, hide = true)]
    out_dir: Option<String>,

    /// BAM/SAM/CRAM file for coverage (alternative to positional)
    #[arg(long)]
    bam: Option<String>,

    /// Path to DIAMOND database (.dmnd) for homolog search
    #[arg(long)]
    dmnd: Option<String>,

    /// Path to FAMSA binary (optional; falls back to PATH then ~/Science/Programs/FAMSA/famsa)
    #[arg(long)]
    famsa: Option<String>,
}


/// Returns the full path to a script if it is found in PATH.
fn which_script(name: &str) -> Option<String> {
    let out = std::process::Command::new("which").arg(name).output().ok()?;
    if out.status.success() {
        let path = String::from_utf8_lossy(&out.stdout).trim().to_string();
        if !path.is_empty() { Some(path) } else { None }
    } else {
        None
    }
}

fn find_famsa(explicit: Option<String>) -> Option<String> {
    if let Some(p) = explicit {
        return Some(p);
    }
    if std::process::Command::new("famsa").arg("--version").output().map(|o| o.status.success()).unwrap_or(false) {
        return Some("famsa".to_string());
    }
    let home = std::env::var("HOME").unwrap_or_default();
    let candidate = format!("{}/Science/Programs/FAMSA/famsa", home);
    if std::path::Path::new(&candidate).exists() {
        return Some(candidate);
    }
    None
}

fn classify_files(cli: &Cli) -> (Option<String>, Option<String>, Option<String>) {
    let mut fasta: Option<String> = None;
    let mut gff:   Option<String> = None;
    let mut bam:   Option<String> = cli.bam.clone();
    for path in &cli.files {
        let lower = path.to_ascii_lowercase();
        let ext = lower.trim_end_matches(".gz");
        if ext.ends_with(".gff") || ext.ends_with(".gff3") || ext.ends_with(".gtf") {
            gff = Some(path.clone());
        } else if ext.ends_with(".fa") || ext.ends_with(".fasta") || ext.ends_with(".fna")
               || ext.ends_with(".fas") || ext.ends_with(".faa") {
            fasta = Some(path.clone());
        } else if ext.ends_with(".bam") || ext.ends_with(".sam") || ext.ends_with(".cram") {
            bam = Some(path.clone());
        } else {
            eprintln!("Warning: unrecognised file extension, skipping: {path}");
        }
    }
    (fasta, gff, bam)
}

#[tokio::main]
async fn main() -> Result<()> {
    let cli = Cli::parse();
    let (fasta_path, gff_path, bam_path) = classify_files(&cli);
    let gff_path = gff_path.ok_or_else(|| anyhow::anyhow!("No GFF/GFF3 file provided"))?;

    // Build fold command: prefer --minifold_mlx path, fall back to `minifold` in PATH.
    let default_out_dir = || {
        cli.pdbs.clone()
            .or_else(|| cli.out_dir.clone())
            .unwrap_or_else(|| ".".to_string())
    };
    // Priority: --minifold_mlx <dir> (fold.py) → minifold in PATH → predict.py in PATH → disabled
    let out_dir_arg = default_out_dir();
    let (on_click_cmd, use_predict_py): (Option<String>, bool) =
        if let Some(ref dir) = cli.minifold_mlx {
            let fold_py = std::path::Path::new(dir).join("fold.py");
            (Some(format!("python {} {{protein_fasta}} --out_dir {}", fold_py.display(), out_dir_arg)), false)
        } else if which_script("minifold").is_some() {
            (Some(format!("minifold {{protein_fasta}} --out_dir {}", out_dir_arg)), false)
        } else if let Some(py) = which_script("predict.py") {
            (Some(format!("python {} {{protein_fasta}} --out_dir {}", py, out_dir_arg)), true)
        } else {
            (None, false)
        };

    let (genome_seq, genome_name, seqname_offsets, plasmid_entries) =
        if let Some(ref path) = fasta_path {
            let fasta = parse_fasta(path)?;
            (Some(fasta.main_seq), fasta.main_name, fasta.seqname_offsets, fasta.plasmids)
        } else {
            (None, String::new(), HashMap::new(), Vec::new())
        };

    let (features, gff_genome_size) = parse_gff(&gff_path, &seqname_offsets)?;

    let genome_size = genome_seq
        .as_ref()
        .map(|s| s.len() as u64)
        .unwrap_or(gff_genome_size);

    let stop_codons: HashMap<(char, u8), Vec<usize>> = if let Some(ref seq) = genome_seq {
        find_stop_codons(seq)
    } else {
        HashMap::new()
    };

    let gc_skew: Vec<f64> = if let Some(ref seq) = genome_seq {
        compute_gc_skew(seq, 500)
    } else {
        vec![]
    };

    let plasmids: Vec<PlasmidData> = plasmid_entries
        .iter()
        .map(|entry| {
            let local_features = features
                .iter()
                .filter(|f| f.seqname == entry.seqname)
                .map(|f| {
                    let mut lf = f.clone();
                    lf.start = f.start.saturating_sub(entry.offset_in_main);
                    lf.end   = f.end  .saturating_sub(entry.offset_in_main);
                    lf
                })
                .collect();
            let psize = entry.sequence.len() as u64;
            PlasmidData {
                name: entry.name.clone(),
                features: local_features,
                genome_size: psize,
                gc_skew: compute_gc_skew(&entry.sequence, 200),
                view_start: 1,
                view_end: psize.min(20_000).max(1),
                sequence: entry.sequence.clone(),
            }
        })
        .collect();

    let coverage = if let Some(ref bam_path) = bam_path {
        match compute_alignment_coverage(bam_path, fasta_path.as_deref(), &seqname_offsets, genome_size) {
            Ok(cov) => Some(cov),
            Err(e) => {
                eprintln!("Warning: could not read BAM file: {e}");
                None
            }
        }
    } else {
        None
    };

    let fold_out_dir: Option<String> = if cli.pdbs.is_some() || cli.out_dir.is_some() || on_click_cmd.is_some() {
        Some(default_out_dir())
    } else {
        None
    };

    let dmnd_db   = cli.dmnd.clone();
    let famsa_bin = find_famsa(cli.famsa.clone());

    let mut app = App::new(features, genome_size, genome_seq, genome_name, stop_codons, on_click_cmd, gc_skew, plasmids, coverage, fold_out_dir, dmnd_db, famsa_bin);
    app.fold_predict_py = use_predict_py && app.on_click_cmd.is_some();
    app.source_files = std::iter::once(Some(gff_path.clone()))
        .chain([fasta_path.clone(), bam_path.clone()])
        .filter_map(|p| p)
        .collect();

    const WEB_PORT: u16 = 7890;
    let web_server = server::WebServer::new();
    let web_server2 = web_server.clone();
    tokio::spawn(async move { web_server2.run(WEB_PORT).await; });
    push_browser_state(&app, &web_server);

    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen, EnableMouseCapture)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;

    let res = run_loop(&mut terminal, &mut app, &web_server).await;

    disable_raw_mode()?;
    execute!(
        terminal.backend_mut(),
        LeaveAlternateScreen,
        DisableMouseCapture
    )?;
    terminal.show_cursor()?;

    res
}

async fn run_loop(
    terminal: &mut ratatui::Terminal<CrosstermBackend<io::Stdout>>,
    app: &mut App,
    web: &std::sync::Arc<server::WebServer>,
) -> Result<()> {
    let mut fold_rx: Option<tokio::sync::oneshot::Receiver<anyhow::Result<String>>> = None;
    let mut msa_rx: Option<tokio::sync::oneshot::Receiver<anyhow::Result<Vec<(String, String)>>>> = None;

    loop {
        // Expire old status messages
        app.tick_status();

        // Poll for completed fold result
        let mut fold_completed = false;
        if let Some(ref mut rx) = fold_rx {
            match rx.try_recv() {
                Ok(result) => {
                    if let Some(ref mut panel) = app.protein {
                        match result {
                            Ok(pdb) => {
                                panel.atoms = genetui::protein::parse_pdb(&pdb);
                                panel.pdb_raw = pdb;
                                panel.folding = false;
                                panel.img_cache = None;
                                fold_completed = true;
                            }
                            Err(e) => {
                                panel.error = Some(format!("{}", e));
                                panel.folding = false;
                            }
                        }
                    }
                    fold_rx = None;
                }
                Err(tokio::sync::oneshot::error::TryRecvError::Empty) => {}
                Err(_) => { fold_rx = None; }
            }
        }
        if fold_completed {
            // Force full redraw to clear any text that minifold wrote to the terminal.
            terminal.clear()?;
            push_browser_state(app, web);
        }

        // Poll for completed MSA
        let mut msa_completed = false;
        if let Some(ref mut rx) = msa_rx {
            match rx.try_recv() {
                Ok(result) => {
                    if let Some(ref mut panel) = app.msa {
                        panel.loading = false;
                        match result {
                            Ok(seqs) => {
                                let n = seqs.len();
                                panel.sequences = seqs;
                                app.set_status(format!("MSA ready: {} sequences", n));
                            }
                            Err(e) => {
                                panel.error = Some(format!("{}", e));
                                app.set_status(format!("MSA error: {}", e));
                            }
                        }
                    }
                    msa_rx = None;
                    msa_completed = true;
                }
                Err(tokio::sync::oneshot::error::TryRecvError::Empty) => {}
                Err(_) => { msa_rx = None; }
            }
        }
        if msa_completed {
            terminal.clear()?;
        }

        // Apply viewport commands from browser
        if let Some((vs, ve)) = web.take_viewport_cmd() {
            if ve > vs && ve <= app.genome_size {
                app.view_start = vs;
                app.view_end   = ve;
            }
        }

        // Render protein image before draw so the frame is never blank
        if let Some(ref panel) = app.protein {
            if !panel.atoms.is_empty() && panel.img_cache.is_none() {
                let atoms_snap = panel.atoms.clone();
                let camera_snap = panel.camera.clone();
                let rect = app.protein_panel_rect;
                if rect.width > 4 && rect.height > 4 {
                    let iw = rect.width.saturating_sub(2);
                    let ih = rect.height.saturating_sub(3);
                    if iw > 0 && ih > 0 {
                        let img = genetui::render3d::render_protein_image(&atoms_snap, &camera_snap, iw, ih);
                        if let Some(ref mut panel_mut) = app.protein {
                            panel_mut.img_cache = Some(img);
                        }
                    }
                }
            }
        }

        terminal.draw(|f| ui::draw(f, app))?;

        if event::poll(Duration::from_millis(100))? {
            match event::read()? {
                Event::Key(key) => {
                    if app.display_menu_open {
                        match key.code {
                            KeyCode::Char('q') => { break; }
                            KeyCode::Esc | KeyCode::Char('d') => {
                                app.display_menu_open = false;
                            }
                            KeyCode::Up | KeyCode::Char('k') => {
                                let n = app.display_menu_item_count();
                                app.display_menu_idx = (app.display_menu_idx + n - 1) % n;
                            }
                            KeyCode::Down | KeyCode::Char('j') => {
                                let n = app.display_menu_item_count();
                                app.display_menu_idx = (app.display_menu_idx + 1) % n;
                            }
                            KeyCode::Char(' ') | KeyCode::Enter => {
                                app.toggle_display_opt();
                            }
                            _ => {}
                        }
                    } else if app.search_popup_open {
                        match key.code {
                            KeyCode::Esc => {
                                app.search_popup_open = false;
                            }
                            KeyCode::Up | KeyCode::Char('k') => {
                                if app.search_popup_idx > 0 {
                                    app.search_popup_idx -= 1;
                                    app.search_idx = app.search_popup_idx;
                                    app.jump_to_feature(app.search_results[app.search_popup_idx]);
                                }
                            }
                            KeyCode::Down | KeyCode::Char('j') => {
                                let n = app.search_results.len();
                                if app.search_popup_idx + 1 < n {
                                    app.search_popup_idx += 1;
                                    app.search_idx = app.search_popup_idx;
                                    app.jump_to_feature(app.search_results[app.search_popup_idx]);
                                }
                            }
                            KeyCode::Enter | KeyCode::Char(' ') => {
                                app.search_popup_select();
                            }
                            _ => {}
                        }
                    } else if app.search_mode {
                        match key.code {
                            KeyCode::Esc => {
                                app.search_mode = false;
                                app.search_query.clear();
                            }
                            KeyCode::Enter => {
                                app.search_mode = false;
                                app.do_search();
                            }
                            KeyCode::Backspace => {
                                app.search_query.pop();
                            }
                            KeyCode::Char(c) => {
                                app.search_query.push(c);
                            }
                            _ => {}
                        }
                    } else if app.active_panel == crate::app::ActivePanel::Msa {
                        if let Some(ref mut panel) = app.msa {
                            match key.code {
                                KeyCode::Esc => {
                                    app.msa = None;
                                    app.active_panel = crate::app::ActivePanel::Genome;
                                    msa_rx = None;
                                }
                                KeyCode::Up        | KeyCode::Char('k') => { if panel.viewport_row > 0 { panel.viewport_row -= 1; } }
                                KeyCode::Down      | KeyCode::Char('j') => { if panel.viewport_row + 1 < panel.sequences.len() { panel.viewport_row += 1; } }
                                KeyCode::Left      | KeyCode::Char('h') => { panel.viewport_col = panel.viewport_col.saturating_sub(1); }
                                KeyCode::Right     | KeyCode::Char('l') => { panel.viewport_col = panel.viewport_col.saturating_add(1); }
                                KeyCode::PageUp                         => { panel.viewport_col = panel.viewport_col.saturating_sub(20); }
                                KeyCode::PageDown                       => { panel.viewport_col = panel.viewport_col.saturating_add(20); }
                                KeyCode::Tab => {
                                    app.active_panel = if app.protein.is_some() {
                                        crate::app::ActivePanel::Protein
                                    } else {
                                        crate::app::ActivePanel::Genome
                                    };
                                }
                                _ => {}
                            }
                        }
                    } else {
                        match key.code {
                            KeyCode::Esc if app.protein.is_some() => {
                                genetui::kitty::kitty_delete();
                                app.protein = None;
                                app.active_panel = crate::app::ActivePanel::Genome;
                                let _ = terminal.clear();
                            }
                            KeyCode::Esc if app.msa.is_some() => {
                                app.msa = None;
                                msa_rx = None;
                                let _ = terminal.clear();
                            }
                            KeyCode::Char('m') => {
                                start_msa_search(app, &mut msa_rx).await;
                            }
                            KeyCode::Char('f') => {
                                start_fold(app, &mut fold_rx).await;
                            }
                            KeyCode::Tab => {
                                app.active_panel = match (app.protein.is_some(), app.msa.is_some(), &app.active_panel) {
                                    (true,  _,    crate::app::ActivePanel::Genome)  => crate::app::ActivePanel::Protein,
                                    (_,     true, crate::app::ActivePanel::Protein) => crate::app::ActivePanel::Msa,
                                    (false, true, crate::app::ActivePanel::Genome)  => crate::app::ActivePanel::Msa,
                                    _                                                => crate::app::ActivePanel::Genome,
                                };
                            }
                            KeyCode::Char('a') => {
                                if let Some(ref panel) = app.protein {
                                    if !panel.aa_seq.is_empty() {
                                        let home = std::env::var("HOME").unwrap_or_default();
                                        let path = format!("{}/Downloads/{}.fa", home, panel.gene_name);
                                        let content = format!(">{}\n{}\n", panel.gene_name, panel.aa_seq);
                                        match std::fs::write(&path, &content) {
                                            Ok(_)  => app.set_status(format!("saved: {}", path)),
                                            Err(e) => app.set_status(format!("save failed: {}", e)),
                                        }
                                    }
                                }
                            }
                            KeyCode::Char('n') if app.protein.is_some() => {
                                if let Some(ref panel) = app.protein {
                                    if !panel.nt_seq.is_empty() {
                                        let home = std::env::var("HOME").unwrap_or_default();
                                        let path = format!("{}/Downloads/{}.fna", home, panel.gene_name);
                                        let content = format!(">{}\n{}\n", panel.gene_name, panel.nt_seq);
                                        match std::fs::write(&path, &content) {
                                            Ok(_)  => app.set_status(format!("saved: {}", path)),
                                            Err(e) => app.set_status(format!("save failed: {}", e)),
                                        }
                                    }
                                }
                            }
                            KeyCode::Char('+') | KeyCode::Char('=')
                                if app.protein.as_ref().map(|p| !p.atoms.is_empty()).unwrap_or(false) =>
                            {
                                if let Some(ref mut p) = app.protein {
                                    p.camera.zoom *= 1.35;
                                    p.img_cache = None;
                                }
                            }
                            KeyCode::Char('-')
                                if app.protein.as_ref().map(|p| !p.atoms.is_empty()).unwrap_or(false) =>
                            {
                                if let Some(ref mut p) = app.protein {
                                    p.camera.zoom /= 1.35;
                                    p.img_cache = None;
                                }
                            }
                            KeyCode::Left | KeyCode::Char('h')
                                if app.active_panel == crate::app::ActivePanel::Protein
                                    && app.protein.as_ref().map(|p| !p.atoms.is_empty()).unwrap_or(false) =>
                            {
                                if let Some(ref mut p) = app.protein {
                                    p.camera.rot_y -= 0.15;
                                    p.img_cache = None;
                                }
                            }
                            KeyCode::Right | KeyCode::Char('l')
                                if app.active_panel == crate::app::ActivePanel::Protein
                                    && app.protein.as_ref().map(|p| !p.atoms.is_empty()).unwrap_or(false) =>
                            {
                                if let Some(ref mut p) = app.protein {
                                    p.camera.rot_y += 0.15;
                                    p.img_cache = None;
                                }
                            }
                            KeyCode::Up | KeyCode::Char('k')
                                if app.active_panel == crate::app::ActivePanel::Protein
                                    && app.protein.as_ref().map(|p| !p.atoms.is_empty()).unwrap_or(false) =>
                            {
                                if let Some(ref mut p) = app.protein {
                                    p.camera.rot_x += 0.15;
                                    p.img_cache = None;
                                }
                            }
                            KeyCode::Down | KeyCode::Char('j')
                                if app.active_panel == crate::app::ActivePanel::Protein
                                    && app.protein.as_ref().map(|p| !p.atoms.is_empty()).unwrap_or(false) =>
                            {
                                if let Some(ref mut p) = app.protein {
                                    p.camera.rot_x -= 0.15;
                                    p.img_cache = None;
                                }
                            }
                            KeyCode::Left | KeyCode::Char('h') => {
                                app.pan(-0.3);
                                push_browser_state(app, web);
                            }
                            KeyCode::Right | KeyCode::Char('l') => {
                                app.pan(0.3);
                                push_browser_state(app, web);
                            }
                            KeyCode::Char('+') | KeyCode::Char('=') => {
                                app.zoom(0.5);
                                push_browser_state(app, web);
                            }
                            KeyCode::Char('-') => {
                                app.zoom(2.0);
                                push_browser_state(app, web);
                            }
                            KeyCode::Char('/') => {
                                app.search_mode = true;
                                app.search_query.clear();
                                app.status_msg.clear();
                            }
                            KeyCode::Char('n') => {
                                app.search_next();
                                push_browser_state(app, web);
                            }
                            KeyCode::Char('w') => {
                                open_browser(WEB_PORT);
                            }
                            KeyCode::Char('d') => {
                                app.display_menu_open = true;
                            }
                            KeyCode::Char('q') => {
                                break;
                            }
                            _ => {}
                        }
                    }
                }
                Event::Mouse(mouse_event) => {
                    let col = mouse_event.column;
                    let row = mouse_event.row;
                    match mouse_event.kind {
                        MouseEventKind::Down(event::MouseButton::Left) => {
                            let r = row as usize; let c = col as usize;
                            if app.msa.is_some() && rect_contains(&app.msa_panel_rect, r, c) {
                                app.active_panel = crate::app::ActivePanel::Msa;
                            } else if app.protein.is_some() && rect_contains(&app.protein_panel_rect, r, c) {
                                app.active_panel = crate::app::ActivePanel::Protein;
                                if let Some(ref mut panel) = app.protein {
                                    panel.drag_last = Some((col, row));
                                }
                            } else {
                                app.active_panel = crate::app::ActivePanel::Genome;
                            }
                            handle_hover(app, row as usize, col as usize);
                            if handle_click(app, row as usize, col as usize).await {
                                fold_rx = None;
                                genetui::kitty::kitty_delete();
                                // Clear any text the previous fold subprocess wrote to the terminal.
                                let _ = terminal.clear();
                            }
                            push_browser_state(app, web);
                        }
                        MouseEventKind::Drag(event::MouseButton::Left) => {
                            if app.active_panel == crate::app::ActivePanel::Protein {
                                if let Some(ref mut panel) = app.protein {
                                    if let Some((lc, lr)) = panel.drag_last {
                                        panel.camera.rot_y += (col as i16 - lc as i16) as f64 * 0.05;
                                        panel.camera.rot_x += (row as i16 - lr as i16) as f64 * 0.05;
                                        panel.drag_last = Some((col, row));
                                        panel.img_cache = None;
                                    }
                                }
                            }
                            handle_hover(app, row as usize, col as usize);
                        }
                        MouseEventKind::Up(_) => {
                            if let Some(ref mut panel) = app.protein {
                                panel.drag_last = None;
                            }
                            handle_hover(app, row as usize, col as usize);
                        }
                        _ => {
                            handle_hover(app, row as usize, col as usize);
                        }
                    }
                }
                _ => {}
            }
        }
    }
    Ok(())
}

async fn start_fold(
    app: &mut App,
    fold_rx: &mut Option<tokio::sync::oneshot::Receiver<anyhow::Result<String>>>,
) {
    let panel = match app.protein.as_mut() {
        Some(p) if !p.folding && !p.aa_seq.is_empty() && fold_rx.is_none() => p,
        _ => return,
    };
    let aa        = panel.aa_seq.clone();
    let gene_name = panel.gene_name.clone();
    panel.folding = true;
    panel.error   = None;
    panel.atoms.clear();
    panel.img_cache = None;
    let (tx, rx) = tokio::sync::oneshot::channel();
    *fold_rx = Some(rx);

    let safe_id: String = gene_name.chars()
        .map(|c| if c.is_alphanumeric() || c == '_' { c } else { '_' })
        .collect();

    if let Some(ref out_dir) = app.fold_out_dir.clone() {
        let expected_pdb = format!("{}/{}.pdb", out_dir, safe_id);
        if let Ok(pdb) = std::fs::read_to_string(&expected_pdb) {
            if let Some(ref mut p) = app.protein { p.folding = false; }
            let _ = tx.send(Ok(pdb));
            return;
        }
    }

    if let (Some(ref cmd_tmpl), Some(ref out_dir)) =
        (app.on_click_cmd.clone(), app.fold_out_dir.clone())
    {
        let pid       = std::process::id();
        let tmp_fasta = format!("/tmp/{}_{}.fa", safe_id, pid);
        let fasta_base = format!("{}_{}", safe_id, pid);
        let expected_pdb = if app.fold_predict_py {
            // jwohlwend/minifold: output goes into a nested results subdirectory
            format!("{}/minifold_results_{}/{}.pdb", out_dir, fasta_base, safe_id)
        } else {
            format!("{}/{}.pdb", out_dir, safe_id)
        };
        let fasta_content = format!(">{}\n{}\n", safe_id, aa);
        let out_dir2 = out_dir.clone();
        if let Err(e) = std::fs::create_dir_all(out_dir).and_then(|_| std::fs::write(&tmp_fasta, &fasta_content)) {
            if let Some(ref mut p) = app.protein {
                p.error   = Some(format!("write fasta: {}", e));
                p.folding = false;
            }
            *fold_rx = None;
        } else {
            let cmd = cmd_tmpl.replace("{protein_fasta}", &tmp_fasta);
            let tmp_fasta2 = tmp_fasta.clone();
            tokio::spawn(async move {
                let res = TokioCommand::new("sh")
                    .arg("-c").arg(&cmd).output().await;
                let _ = std::fs::remove_file(&tmp_fasta2);
                let result = match res {
                    Ok(out) if out.status.success() => {
                        match std::fs::read_to_string(&expected_pdb) {
                            Ok(pdb) => Ok(pdb),
                            Err(_) => {
                                // Process exited 0 but no PDB — tool may have a length limit.
                                let stderr = String::from_utf8_lossy(&out.stderr);
                                let stdout = String::from_utf8_lossy(&out.stdout);
                                let extra = if !stderr.is_empty() {
                                    format!(": {}", &stderr[..stderr.len().min(200)])
                                } else if !stdout.is_empty() {
                                    format!(": {}", &stdout[..stdout.len().min(200)])
                                } else {
                                    " (no output — sequence may exceed the tool's length limit)".into()
                                };
                                Err(anyhow::anyhow!("No PDB produced{}", extra.trim_end()))
                            }
                        }
                    }
                    Ok(out) => {
                        let stderr = String::from_utf8_lossy(&out.stderr);
                        let stdout = String::from_utf8_lossy(&out.stdout);
                        let msg = if !stderr.is_empty() { &stderr[..stderr.len().min(250)] }
                                  else if !stdout.is_empty() { &stdout[..stdout.len().min(250)] }
                                  else { "exited with error (no output)" };
                        Err(anyhow::anyhow!("fold failed: {}", msg.trim()))
                    }
                    Err(e) => Err(anyhow::anyhow!("spawn: {}", e)),
                };
                let _ = tx.send(result);
            });
            let _ = out_dir2;
        }
    } else {
        if let Some(ref mut p) = app.protein {
            p.error   = Some("No local folding tool configured (use --minifold_mlx)".to_string());
            p.folding = false;
        }
        *fold_rx = None;
    }
}

async fn start_msa_search(
    app: &mut App,
    msa_rx: &mut Option<tokio::sync::oneshot::Receiver<anyhow::Result<Vec<(String, String)>>>>,
) {
    let (gene_name, aa_seq) = match app.protein.as_ref() {
        Some(p) if !p.aa_seq.is_empty() => (p.gene_name.clone(), p.aa_seq.clone()),
        _ => { app.set_status("Select a gene first (protein sequence needed for MSA)"); return; }
    };
    let dmnd_db = match app.dmnd_db.clone() {
        Some(db) => db,
        None => { app.set_status("No DIAMOND db specified — use --dmnd <path.dmnd>"); return; }
    };
    let famsa_bin = match app.famsa_bin.clone() {
        Some(b) => b,
        None => { app.set_status("FAMSA not found — use --famsa or add to PATH"); return; }
    };

    app.msa = Some(crate::app::MsaPanel {
        gene_name: gene_name.clone(),
        sequences: Vec::new(),
        loading: true,
        error: None,
        viewport_row: 0,
        viewport_col: 0,
    });
    app.active_panel = crate::app::ActivePanel::Msa;
    app.set_status(format!("MSA: searching homologs of {}…", gene_name));

    let pid = std::process::id();
    let query_fa = format!("/tmp/genetui_msa_{}_query.fa", pid);
    let hits_tsv = format!("/tmp/genetui_msa_{}_hits.tsv", pid);
    let seqs_fa  = format!("/tmp/genetui_msa_{}_seqs.fa",  pid);
    let aln_fa   = format!("/tmp/genetui_msa_{}_aln.fa",   pid);
    let safe_name = gene_name.replace([' ', '/', '|'], "_");

    if std::fs::write(&query_fa, format!(">{}\n{}\n", safe_name, aa_seq)).is_err() {
        if let Some(ref mut m) = app.msa { m.loading = false; m.error = Some("Failed to write temp FASTA".into()); }
        return;
    }

    let (tx, rx) = tokio::sync::oneshot::channel();
    *msa_rx = Some(rx);

    tokio::spawn(async move {
        let result: anyhow::Result<Vec<(String, String)>> = async {
            let out = TokioCommand::new("diamond")
                .args(["blastp", "--db", &dmnd_db, "--query", &query_fa,
                       "--out", &hits_tsv,
                       "--outfmt", "6", "qseqid", "sseqid", "evalue", "bitscore", "full_sseq",
                       "--evalue", "1e-10", "--max-target-seqs", "20",
                       "--threads", "4", "--quiet", "--masking", "0"])
                .output().await?;
            if !out.status.success() {
                let err = String::from_utf8_lossy(&out.stderr);
                anyhow::bail!("DIAMOND: {}", &err[..err.len().min(300)]);
            }

            let query_seq = aa_seq.clone();
            let mut seqs: Vec<(String, String)> = vec![(safe_name.clone(), query_seq)];
            for line in std::fs::read_to_string(&hits_tsv)?.lines() {
                let cols: Vec<&str> = line.split('\t').collect();
                if cols.len() >= 5 {
                    let id  = cols[1].to_string();
                    let seq = cols[4].replace('-', "");
                    if !seq.is_empty() && seq != "N/A" && !seqs.iter().any(|(i,_)| i == &id) {
                        seqs.push((id, seq));
                    }
                }
            }
            if seqs.len() < 2 { anyhow::bail!("No homologs found (e-value < 1e-10)"); }

            let fasta: String = seqs.iter().map(|(id, s)| format!(">{}\n{}\n", id, s)).collect();
            std::fs::write(&seqs_fa, fasta)?;

            let out = TokioCommand::new(&famsa_bin)
                .args([seqs_fa.as_str(), aln_fa.as_str(), "-t", "4"])
                .output().await?;
            if !out.status.success() {
                anyhow::bail!("FAMSA: {}", String::from_utf8_lossy(&out.stderr).chars().take(300).collect::<String>());
            }

            let aligned = std::fs::read_to_string(&aln_fa)?;
            let mut result: Vec<(String, String)> = Vec::new();
            let (mut cur_id, mut cur_seq) = (String::new(), String::new());
            for line in aligned.lines() {
                if let Some(rest) = line.strip_prefix('>') {
                    if !cur_id.is_empty() { result.push((cur_id.clone(), cur_seq.clone())); cur_seq.clear(); }
                    cur_id = rest.split_whitespace().next().unwrap_or("").to_string();
                } else { cur_seq.push_str(line.trim()); }
            }
            if !cur_id.is_empty() { result.push((cur_id, cur_seq)); }

            for f in [&query_fa, &hits_tsv, &seqs_fa, &aln_fa] { let _ = std::fs::remove_file(f); }
            Ok(result)
        }.await;
        let _ = tx.send(result);
    });
}

fn handle_hover(app: &mut App, row: usize, col: usize) {
    let hit = app
        .hit_map
        .iter()
        .find(|(r, cs, ce, _)| *r == row && col >= *cs && col <= *ce)
        .copied();
    app.hovered = hit.map(|(_, _, _, idx)| idx);

    if rect_contains(&app.minimap_rect, row, col) {
        app.hovered_map = Some(0);
    } else {
        let plasmid_hit = app.plasmid_rects.iter().enumerate()
            .find(|(_, r)| rect_contains(r, row, col))
            .map(|(i, _)| i + 1);
        app.hovered_map = plasmid_hit;
    }

    app.hovered_legend = rect_contains(&app.legend_rect, row, col);
}

fn rect_contains(rect: &ratatui::layout::Rect, row: usize, col: usize) -> bool {
    row >= rect.y as usize && row < (rect.y + rect.height) as usize
    && col >= rect.x as usize && col < (rect.x + rect.width) as usize
}

fn minimap_click_pos(area: ratatui::layout::Rect, row: usize, col: usize, genome_size: u64) -> Option<u64> {
    let w = area.width  as f64;
    let h = area.height as f64;
    if w == 0.0 || h == 0.0 { return None; }

    const X_RANGE:    f64 = 1.1;
    const CELL_RATIO: f64 = 2.15;
    const R_IN:       f64 = 0.82 * 0.70;
    const CLICK_OUTER: f64 = 1.06;

    let y_range = X_RANGE * (CELL_RATIO * h) / w;

    let cx = (col as f64 - area.x as f64 + 0.5) / w * (2.0 * X_RANGE) - X_RANGE;
    let cy = y_range - (row as f64 - area.y as f64 + 0.5) / h * (2.0 * y_range);

    let dist = (cx * cx + cy * cy).sqrt();

    if dist < R_IN {
        Some(1)
    } else if dist <= CLICK_OUTER {
        let angle = cx.atan2(cy);
        let frac  = (angle / (2.0 * std::f64::consts::PI)).rem_euclid(1.0);
        let pos   = ((frac * genome_size as f64) as u64).max(1).min(genome_size);
        Some(pos)
    } else {
        None
    }
}

/// Returns true if a new gene panel was opened (caller should reset fold_rx).
async fn handle_click(app: &mut App, row: usize, col: usize) -> bool {
    if app.search_popup_open && rect_contains(&app.search_popup_rect, row, col) {
        let item = row.saturating_sub(app.search_popup_rect.y as usize + 1);
        if item < app.search_results.len() {
            app.search_popup_idx = item;
            app.search_popup_select();
        }
        return false;
    }

    if app.display_menu_open && rect_contains(&app.display_menu_rect, row, col) {
        let menu_row = row.saturating_sub(app.display_menu_rect.y as usize + 1);
        if menu_row < app.display_menu_item_count() {
            app.display_menu_idx = menu_row;
            app.toggle_display_opt();
        }
        return false;
    }

    if rect_contains(&app.minimap_rect, row, col) {
        let genome_size = app.genome_size;
        let cur_span    = app.view_end.saturating_sub(app.view_start).max(1);
        app.switch_genome(0);
        if let Some(pos) = minimap_click_pos(app.minimap_rect, row, col, genome_size) {
            let half  = cur_span / 2;
            let start = pos.saturating_sub(half).max(1);
            let end   = (start + cur_span).min(genome_size);
            app.goto(start, end);
        }
        return false;
    }

    for (i, rect) in app.plasmid_rects.clone().iter().enumerate() {
        if rect_contains(rect, row, col) {
            let psize    = app.plasmids.get(i).map(|p| p.genome_size).unwrap_or(1);
            let cur_span = app.plasmids.get(i)
                .map(|p| p.view_end.saturating_sub(p.view_start).max(1))
                .unwrap_or(psize);
            app.switch_genome(i + 1);
            if let Some(pos) = minimap_click_pos(*rect, row, col, psize) {
                let half  = cur_span / 2;
                let start = pos.saturating_sub(half).max(1);
                let end   = (start + cur_span).min(psize);
                app.goto_active(start, end);
            }
            return false;
        }
    }

    let hit = app
        .hit_map
        .iter()
        .find(|(r, cs, ce, _)| *r == row && col >= *cs && col <= *ce)
        .copied();

    if let Some((_, _, _, feat_idx)) = hit {
        app.selected = Some(feat_idx);

        let feat = app.active_features()[feat_idx].clone();
        let feat_name   = feat.name.clone();
        let feat_locus  = feat.locus_tag.clone();
        let feat_start  = feat.start;
        let feat_end    = feat.end;
        let feat_strand = feat.strand;

        app.set_status(format!("selected: {} ({})", feat_name, feat_locus));

        let seq_source: Option<String> = if app.active_genome > 0 {
            app.plasmids.get(app.active_genome - 1).map(|p| p.sequence.clone())
        } else {
            app.genome_seq.clone()
        };

        if let Some(seq) = seq_source {
            let start_idx = (feat_start as usize).saturating_sub(1);
            let end_idx   = (feat_end   as usize).min(seq.len());
            let raw = &seq[start_idx..end_idx];
            let nt_seq = if feat_strand == '-' {
                crate::core::reverse_complement_str(raw).to_uppercase()
            } else {
                raw.to_uppercase()
            };
            let aa_seq = translate(&nt_seq);

            app.open_protein_panel(feat_name.clone(), nt_seq, aa_seq.clone());

            if let Some(ref cmd_template) = app.on_click_cmd.clone() {
                if aa_seq.is_empty() {
                    app.set_status(format!("{}: empty translation", feat_name));
                    return true;
                }
                let label = if feat_locus.is_empty() { &feat_name } else { &feat_locus };
                let tmp_path = format!("/tmp/{}_prot_{}.fa", label, std::process::id());
                let fasta_content = format!(
                    ">{name}|{lt}|{strand}|{start}-{end}\n{seq}\n",
                    name = feat_name, lt = feat_locus, strand = feat_strand,
                    start = feat_start, end = feat_end, seq = aa_seq
                );
                if std::fs::write(&tmp_path, &fasta_content).is_err() {
                    app.set_status(format!("{}: failed to write temp fasta", feat_name));
                    return true;
                }
                let cmd = cmd_template.replace("{protein_fasta}", &tmp_path);
                let short = if cmd.len() > 70 { format!("{}…", &cmd[..70]) } else { cmd.clone() };
                app.set_status(format!("▶ {} {}", feat_name, short));
                let tmp_path_clone = tmp_path.clone();
                let feat_name_clone = feat_name.clone();
                tokio::spawn(async move {
                    let result = TokioCommand::new("sh").arg("-c").arg(&cmd).output().await;
                    let _ = std::fs::remove_file(&tmp_path_clone);
                    match result {
                        Ok(out) if out.status.success() => {
                            let stdout = String::from_utf8_lossy(&out.stdout);
                            eprintln!("✓ {} {}", feat_name_clone, &stdout[..stdout.len().min(100)]);
                        }
                        Ok(out) => {
                            let stderr = String::from_utf8_lossy(&out.stderr);
                            eprintln!("✗ {} {}", feat_name_clone, &stderr[..stderr.len().min(100)]);
                        }
                        Err(e) => { eprintln!("✗ {} error: {}", feat_name_clone, e); }
                    }
                });
            }
            return true;
        } else {
            app.set_status(format!("{}: no FASTA loaded", feat_name));
        }
    }
    false
}

const WEB_PORT: u16 = 7890;

fn push_browser_state(app: &App, web: &server::WebServer) {
    const PLUS_COLORS:  &[(u8,u8,u8)] = &[
        (72,210,180),(80,170,230),(100,220,150),(55,225,200),(95,185,240),(120,215,170),
    ];
    const MINUS_COLORS: &[(u8,u8,u8)] = &[
        (240,148,50),(225,118,70),(245,175,45),(230,132,80),(238,160,55),(220,105,75),
    ];
    let features = app.active_features().iter()
        .map(|f| {
            let (r,g,b) = if f.strand == '+' {
                PLUS_COLORS[f.color_idx % PLUS_COLORS.len()]
            } else {
                MINUS_COLORS[f.color_idx % MINUS_COLORS.len()]
            };
            server::BrowserFeature {
                name:      if f.name.is_empty() { f.locus_tag.clone() } else { f.name.clone() },
                start:     f.start,
                end:       f.end,
                strand:    f.strand,
                color:     format!("#{:02x}{:02x}{:02x}", r, g, b),
                noncoding: f.noncoding,
            }
        })
        .collect();

    let active_seq: Option<&str> = if app.active_genome == 0 {
        app.genome_seq.as_deref()
    } else {
        app.plasmids.get(app.active_genome - 1).map(|p| p.sequence.as_str()).filter(|s| !s.is_empty())
    };
    let (sequence, seq_start) = active_seq.map(|s| {
        let center = (app.active_view_start() as usize + app.active_view_end() as usize) / 2;
        let half   = 500_000usize;
        let vs = center.saturating_sub(half);
        let ve = (center + half).min(s.len());
        if vs < ve { (s[vs..ve].to_string(), vs as u64) } else { (String::new(), 0u64) }
    }).unwrap_or_default();

    let (protein_pdb, protein_name) = app.protein.as_ref()
        .filter(|p| !p.pdb_raw.is_empty())
        .map(|p| (p.pdb_raw.clone(), p.gene_name.clone()))
        .unwrap_or_default();

    // Coverage: send viewport-context bins for chromosome only
    let (cov_plus, cov_minus, cov_bin_sz, cov_bin_start) =
        if app.active_genome == 0 {
            if let Some(ref cov) = app.coverage {
                let bs   = cov.bin_size;
                let vs   = app.active_view_start();
                let ve   = app.active_view_end();
                let span = ve.saturating_sub(vs);
                let ctx_s = vs.saturating_sub(span);
                let ctx_e = (ve + span).min(cov.genome_size);
                let fb = (ctx_s / bs) as usize;
                let lb = ((ctx_e / bs) as usize).min(cov.plus.len().saturating_sub(1));
                (cov.plus[fb..=lb].to_vec(), cov.minus[fb..=lb].to_vec(), bs, fb as u64 * bs)
            } else { (vec![], vec![], 1000, 0) }
        } else { (vec![], vec![], 1000, 0) };

    web.push(server::BrowserState {
        view_start:  app.active_view_start(),
        view_end:    app.active_view_end(),
        genome_size: app.active_genome_size(),
        genome_name: if app.active_genome == 0 {
            app.genome_name.clone()
        } else {
            app.plasmids.get(app.active_genome - 1)
                .map(|p| p.name.clone())
                .unwrap_or_else(|| app.genome_name.clone())
        },
        features,
        sequence,
        seq_start,
        protein_pdb,
        protein_name,
        coverage_plus:      cov_plus,
        coverage_minus:     cov_minus,
        coverage_bin_size:  cov_bin_sz,
        coverage_bin_start: cov_bin_start,
        input_files: app.source_files.clone(),
    });
}

fn open_browser(port: u16) {
    let url = format!("http://localhost:{}", port);
    #[cfg(target_os = "macos")]
    let _ = std::process::Command::new("open").arg(&url).spawn();
    #[cfg(target_os = "linux")]
    let _ = std::process::Command::new("xdg-open").arg(&url).spawn();
    #[cfg(target_os = "windows")]
    let _ = std::process::Command::new("cmd").args(["/c", "start", &url]).spawn();
}
