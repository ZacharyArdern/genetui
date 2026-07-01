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
use core::{compute_alignment_coverage, ensure_bam_index, fetch_reads_in_region, compute_gc_skew, find_stop_codons, parse_fasta, parse_gff, reverse_complement_str, translate};

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

    /// Basic colour mode for standard terminals (no truecolor required)
    #[arg(long)]
    basic: bool,
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
    let (fasta_path, gff_path, mut bam_path) = classify_files(&cli);
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

    // Ensure BAM is sorted+indexed; may return a different (sorted) path.
    if let Some(ref bp) = bam_path {
        match ensure_bam_index(bp) {
            Ok(actual) => { bam_path = Some(actual); }
            Err(e)     => { eprintln!("Warning: could not build BAM index: {e}"); }
        }
    }

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

    let dmnd_db = match cli.dmnd.clone() {
        Some(path) if path.ends_with(".fasta") || path.ends_with(".fasta.gz") || path.ends_with(".fa") || path.ends_with(".fa.gz") => {
            let db_path = format!("{}.dmnd", path);
            if !std::path::Path::new(&db_path).exists() {
                eprintln!("Building DIAMOND database from {} → {}", path, db_path);
                let status = std::process::Command::new("diamond")
                    .args(["makedb", "--in", &path, "--db", &db_path, "--quiet"])
                    .status();
                match status {
                    Ok(s) if s.success() => eprintln!("DIAMOND db built: {}", db_path),
                    Ok(s) => eprintln!("diamond makedb failed (exit {})", s),
                    Err(e) => eprintln!("Failed to run diamond: {}", e),
                }
            }
            Some(db_path)
        }
        other => other,
    };
    let famsa_bin = find_famsa(cli.famsa.clone());

    let mut app = App::new(features, genome_size, genome_seq, genome_name, stop_codons, on_click_cmd, gc_skew, plasmids, coverage, fold_out_dir, dmnd_db, famsa_bin);
    app.bam_path = bam_path.clone();
    app.seqname_offsets = seqname_offsets.clone();
    app.fold_predict_py = use_predict_py && app.on_click_cmd.is_some();
    app.basic_mode = cli.basic;
    if cli.basic {
        app.display_opts.show_legend = false;
    }
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
    let mut diamond_rx: Option<tokio::sync::oneshot::Receiver<anyhow::Result<Vec<core::Feature>>>> = None;

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
            push_browser_state(app, web);
        }

        // Poll for completed DIAMOND blast
        let mut diamond_completed = false;
        if let Some(ref mut rx) = diamond_rx {
            match rx.try_recv() {
                Ok(result) => {
                    app.blast_running = false;
                    match result {
                        Ok(feats) => {
                            let n = feats.len();
                            app.blast_features = feats;
                            app.set_status(format!("DIAMOND: {} hits", n));
                            diamond_completed = true;
                        }
                        Err(e) => {
                            app.blast_error = Some(format!("{}", e));
                            app.set_status(format!("DIAMOND error: {}", e));
                        }
                    }
                    diamond_rx = None;
                }
                Err(tokio::sync::oneshot::error::TryRecvError::Empty) => {}
                Err(_) => { diamond_rx = None; app.blast_running = false; }
            }
        }
        if diamond_completed {
            push_browser_state(app, web);
        }

        // Apply viewport commands from browser
        if let Some((vs, ve)) = web.take_viewport_cmd() {
            if ve > vs && ve <= app.genome_size {
                app.view_start = vs;
                app.view_end   = ve;
                push_browser_state(app, web);
            }
        }

        // Apply fold commands from browser (find gene, open protein panel, start fold)
        if let Some(gene_name) = web.take_fold_cmd() {
            let feat_opt = app.active_features().iter().find(|f| {
                f.name == gene_name || f.locus_tag == gene_name
            }).cloned();
            if let Some(feat) = feat_opt {
                let seq_source: Option<String> = if app.active_genome > 0 {
                    app.plasmids.get(app.active_genome - 1).map(|p| p.sequence.clone())
                } else {
                    app.genome_seq.clone()
                };
                if let Some(seq) = seq_source {
                    let start_idx = (feat.start as usize).saturating_sub(1);
                    let end_idx   = (feat.end as usize).min(seq.len());
                    let raw = &seq[start_idx..end_idx];
                    let nt_seq = if feat.strand == '-' {
                        crate::core::reverse_complement_str(raw).to_uppercase()
                    } else {
                        raw.to_uppercase()
                    };
                    let aa_seq = translate(&nt_seq);
                    app.open_protein_panel(feat.name.clone(), nt_seq, aa_seq);
                    start_fold(app, &mut fold_rx).await;
                } else {
                    app.set_status(format!("{}: no FASTA loaded", gene_name));
                }
            } else {
                app.set_status(format!("fold: gene '{}' not found", gene_name));
            }
            push_browser_state(app, web);
        }

        // Apply switch_genome commands from browser
        if let Some((genome_idx, vs, ve)) = web.take_switch_cmd() {
            app.switch_genome(genome_idx);
            let gsize = app.active_genome_size();
            if vs < ve && ve <= gsize {
                if app.active_genome == 0 {
                    app.view_start = vs;
                    app.view_end   = ve;
                } else if let Some(p) = app.plasmids.get_mut(app.active_genome.saturating_sub(1)) {
                    p.view_start = vs;
                    p.view_end   = ve;
                }
            }
            push_browser_state(app, web);
        }

        // Apply MSA commands from browser (find gene, open protein panel, start MSA)
        if let Some(gene_name) = web.take_msa_cmd() {
            let feat_opt = app.active_features().iter().find(|f| {
                f.name == gene_name || f.locus_tag == gene_name
            }).cloned();
            if let Some(feat) = feat_opt {
                let seq_source: Option<String> = if app.active_genome > 0 {
                    app.plasmids.get(app.active_genome - 1).map(|p| p.sequence.clone())
                } else {
                    app.genome_seq.clone()
                };
                if let Some(seq) = seq_source {
                    let start_idx = (feat.start as usize).saturating_sub(1);
                    let end_idx   = (feat.end as usize).min(seq.len());
                    let raw = &seq[start_idx..end_idx];
                    let nt_seq = if feat.strand == '-' {
                        crate::core::reverse_complement_str(raw).to_uppercase()
                    } else {
                        raw.to_uppercase()
                    };
                    let aa_seq = translate(&nt_seq);
                    app.open_protein_panel(feat.name.clone(), nt_seq, aa_seq);
                    start_msa_search(app, &mut msa_rx).await;
                } else {
                    app.set_status(format!("{}: no FASTA loaded", gene_name));
                }
            } else {
                app.set_status(format!("msa: gene '{}' not found", gene_name));
            }
            push_browser_state(app, web);
        }

        // Handle browser: set DIAMOND DB path
        if let Some(path) = web.take_set_dmnd_cmd() {
            let expanded = if path.starts_with("~/") {
                format!("{}{}", std::env::var("HOME").unwrap_or_default(), &path[1..])
            } else { path.clone() };
            app.dmnd_db = Some(expanded);
            app.needs_dmnd_db = false;
            app.path_completions.clear();
            app.set_status(format!("DIAMOND database set: {}", path));
            if let Some(gene_name) = app.pending_msa_gene.take() {
                let feat_opt = app.active_features().iter().find(|f| f.name == gene_name || f.locus_tag == gene_name).cloned();
                if let Some(feat) = feat_opt {
                    let seq_source = if app.active_genome > 0 {
                        app.plasmids.get(app.active_genome - 1).map(|p| p.sequence.clone())
                    } else { app.genome_seq.clone() };
                    if let Some(seq) = seq_source {
                        let start_idx = (feat.start as usize).saturating_sub(1);
                        let end_idx = (feat.end as usize).min(seq.len());
                        let raw = &seq[start_idx..end_idx];
                        let nt_seq = if feat.strand == '-' {
                            crate::core::reverse_complement_str(raw).to_uppercase()
                        } else { raw.to_uppercase() };
                        let aa_seq = translate(&nt_seq);
                        app.open_protein_panel(feat.name.clone(), nt_seq, aa_seq);
                        start_msa_search(app, &mut msa_rx).await;
                    }
                }
            }
            push_browser_state(app, web);
        }

        // Handle browser: path tab-completion request
        if let Some(prefix) = web.take_complete_path_cmd() {
            app.path_completions = blast_completions(&prefix);
            push_browser_state(app, web);
        }

        // Render protein image before draw so the frame is never blank.
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
                        app.kitty_image_dirty = true;
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
                    } else if app.search_menu_open {
                        match key.code {
                            KeyCode::Esc => { app.search_menu_open = false; }
                            KeyCode::Up | KeyCode::Char('k') => {
                                if app.search_menu_idx > 0 { app.search_menu_idx -= 1; }
                            }
                            KeyCode::Down | KeyCode::Char('j') => {
                                if app.search_menu_idx < 1 { app.search_menu_idx += 1; }
                            }
                            KeyCode::Enter | KeyCode::Char(' ') => {
                                app.search_menu_open = false;
                                if app.search_menu_idx == 0 {
                                    // Gene / Coordinate search
                                    app.search_mode = true;
                                    app.search_query.clear();
                                    app.status_msg.clear();
                                } else {
                                    // DIAMOND blast: open target menu
                                    app.blast_target_open = true;
                                    app.blast_target_idx = 0;
                                }
                            }
                            _ => {}
                        }
                    } else if app.blast_target_open {
                        match key.code {
                            KeyCode::Esc => { app.blast_target_open = false; }
                            KeyCode::Up | KeyCode::Char('k') => {
                                if app.blast_target_idx > 0 { app.blast_target_idx -= 1; }
                            }
                            KeyCode::Down | KeyCode::Char('j') => {
                                if app.blast_target_idx < 1 { app.blast_target_idx += 1; }
                            }
                            KeyCode::Enter | KeyCode::Char(' ') => {
                                app.blast_target_open = false;
                                app.blast_file_open = true;
                                app.blast_file_path.clear();
                                app.blast_completions.clear();
                                app.blast_completion_idx = 0;
                            }
                            _ => {}
                        }
                    } else if app.blast_file_open {
                        match key.code {
                            KeyCode::Esc => {
                                app.blast_file_open = false;
                                app.blast_file_path.clear();
                                app.blast_completions.clear();
                            }
                            KeyCode::Enter => {
                                app.blast_file_open = false;
                                app.blast_completions.clear();
                                start_diamond(app, &mut diamond_rx).await;
                            }
                            KeyCode::Backspace => {
                                app.blast_file_path.pop();
                                app.blast_completions.clear();
                            }
                            KeyCode::Tab => {
                                if app.blast_completions.is_empty() {
                                    app.blast_completions = blast_completions(&app.blast_file_path);
                                    app.blast_completion_idx = 0;
                                } else {
                                    app.blast_completion_idx = (app.blast_completion_idx + 1) % app.blast_completions.len();
                                }
                                if !app.blast_completions.is_empty() {
                                    app.blast_file_path = app.blast_completions[app.blast_completion_idx].clone();
                                }
                            }
                            KeyCode::Char(c) => {
                                app.blast_file_path.push(c);
                                app.blast_completions.clear();
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
                                    p.img_cache = None; app.kitty_image_dirty = true;
                                }
                            }
                            KeyCode::Char('-')
                                if app.protein.as_ref().map(|p| !p.atoms.is_empty()).unwrap_or(false) =>
                            {
                                if let Some(ref mut p) = app.protein {
                                    p.camera.zoom /= 1.35;
                                    p.img_cache = None; app.kitty_image_dirty = true;
                                }
                            }
                            KeyCode::Left | KeyCode::Char('h')
                                if app.active_panel == crate::app::ActivePanel::Protein
                                    && app.protein.as_ref().map(|p| !p.atoms.is_empty()).unwrap_or(false) =>
                            {
                                if let Some(ref mut p) = app.protein {
                                    p.camera.rot_y -= 0.15;
                                    p.img_cache = None; app.kitty_image_dirty = true;
                                }
                            }
                            KeyCode::Right | KeyCode::Char('l')
                                if app.active_panel == crate::app::ActivePanel::Protein
                                    && app.protein.as_ref().map(|p| !p.atoms.is_empty()).unwrap_or(false) =>
                            {
                                if let Some(ref mut p) = app.protein {
                                    p.camera.rot_y += 0.15;
                                    p.img_cache = None; app.kitty_image_dirty = true;
                                }
                            }
                            KeyCode::Up | KeyCode::Char('k')
                                if app.active_panel == crate::app::ActivePanel::Protein
                                    && app.protein.as_ref().map(|p| !p.atoms.is_empty()).unwrap_or(false) =>
                            {
                                if let Some(ref mut p) = app.protein {
                                    p.camera.rot_x += 0.15;
                                    p.img_cache = None; app.kitty_image_dirty = true;
                                }
                            }
                            KeyCode::Down | KeyCode::Char('j')
                                if app.active_panel == crate::app::ActivePanel::Protein
                                    && app.protein.as_ref().map(|p| !p.atoms.is_empty()).unwrap_or(false) =>
                            {
                                if let Some(ref mut p) = app.protein {
                                    p.camera.rot_x -= 0.15;
                                    p.img_cache = None; app.kitty_image_dirty = true;
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
                                app.search_menu_open = true;
                                app.search_menu_idx = 0;
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
                            // Track minimap drag start
                            if rect_contains(&app.minimap_rect, r, c) {
                                app.minimap_dragging = true;
                                app.plasmid_drag_idx = None;
                            } else {
                                let pi = app.plasmid_rects.iter().enumerate()
                                    .find(|(_, rect)| rect_contains(rect, r, c))
                                    .map(|(i, _)| i);
                                app.plasmid_drag_idx = pi;
                                app.minimap_dragging = pi.is_none();
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
                            } else if app.minimap_dragging {
                                let genome_size = app.genome_size;
                                app.switch_genome(0);
                                if let Some(pos) = minimap_click_pos(app.minimap_rect, row as usize, col as usize, genome_size) {
                                    let cur_span = app.view_end.saturating_sub(app.view_start).max(1);
                                    let half = cur_span / 2;
                                    let start = pos.saturating_sub(half);
                                    let end = (start + cur_span).min(genome_size);
                                    app.goto(start, end);
                                    push_browser_state(app, web);
                                }
                            } else if let Some(pi) = app.plasmid_drag_idx {
                                if let Some(psize) = app.plasmids.get(pi).map(|p| p.genome_size) {
                                    let rect = app.plasmid_rects.get(pi).copied().unwrap_or_default();
                                    app.switch_genome(pi + 1);
                                    if let Some(pos) = minimap_click_pos(rect, row as usize, col as usize, psize) {
                                        let cur_span = app.plasmids.get(pi)
                                            .map(|p| p.view_end.saturating_sub(p.view_start).max(1))
                                            .unwrap_or(psize);
                                        let half = cur_span / 2;
                                        let start = pos.saturating_sub(half);
                                        let end = (start + cur_span).min(psize);
                                        app.goto_active(start, end);
                                        push_browser_state(app, web);
                                    }
                                }
                            }
                            handle_hover(app, row as usize, col as usize);
                        }
                        MouseEventKind::Up(_) => {
                            if let Some(ref mut panel) = app.protein {
                                panel.drag_last = None;
                            }
                            app.minimap_dragging = false;
                            app.plasmid_drag_idx = None;
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
        None => {
            app.set_status("No DIAMOND database — specify path below");
            app.needs_dmnd_db = true;
            app.pending_msa_gene = app.protein.as_ref().map(|p| p.gene_name.clone());
            return;
        }
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

async fn start_diamond(
    app: &mut App,
    diamond_rx: &mut Option<tokio::sync::oneshot::Receiver<anyhow::Result<Vec<core::Feature>>>>,
) {
    if app.blast_running { return; }
    let query_fasta = app.blast_file_path.trim().to_string();
    if query_fasta.is_empty() {
        app.set_status("DIAMOND: no query file specified");
        return;
    }
    let expanded_path = if query_fasta.starts_with("~/") {
        let home = std::env::var("HOME").unwrap_or_default();
        format!("{}{}", home, &query_fasta[1..])
    } else {
        query_fasta.clone()
    };
    if !std::path::Path::new(&expanded_path).exists() {
        app.set_status(format!("DIAMOND: file not found: {}", query_fasta));
        return;
    }
    let genome_seq = match app.genome_seq.clone() {
        Some(s) => s,
        None => { app.set_status("DIAMOND: no genome sequence loaded"); return; }
    };
    let use_6ft = app.blast_target_idx == 0;
    let features_for_blastp = if !use_6ft { app.features.clone() } else { Vec::new() };
    let genome_size = app.genome_size;

    app.blast_running = true;
    app.blast_error = None;
    app.blast_features.clear();
    app.set_status("DIAMOND running\u{2026}");

    let (tx, rx) = tokio::sync::oneshot::channel();
    *diamond_rx = Some(rx);

    tokio::spawn(async move {
        let result = run_diamond(genome_seq, genome_size, expanded_path, use_6ft, features_for_blastp).await;
        let _ = tx.send(result);
    });
}

async fn run_diamond(
    genome_seq: String,
    genome_size: u64,
    query_fasta: String,
    use_6ft: bool,
    annotated_features: Vec<core::Feature>,
) -> anyhow::Result<Vec<core::Feature>> {
    use std::process::Stdio;
    use std::io::Write;

    let pid = std::process::id();
    let tmp = std::env::temp_dir();
    let db_path   = tmp.join(format!("genetui_db_{}", pid));
    let out_tsv   = tmp.join(format!("genetui_hits_{}.tsv", pid));

    let (query_file, feat_map): (std::path::PathBuf, HashMap<String, (u64, u64, char)>) = if use_6ft {
        let genome_fasta = tmp.join(format!("genetui_genome_{}.fna", pid));
        let mut f = std::fs::File::create(&genome_fasta)?;
        writeln!(f, ">genome")?;
        writeln!(f, "{}", genome_seq)?;
        (genome_fasta, HashMap::new())
    } else {
        let protein_fasta = tmp.join(format!("genetui_proteins_{}.faa", pid));
        let mut content = String::new();
        let mut feat_map: HashMap<String, (u64, u64, char)> = HashMap::new();
        for (i, feat) in annotated_features.iter().enumerate() {
            if feat.noncoding || feat.end <= feat.start { continue; }
            let s = (feat.start as usize).saturating_sub(1).min(genome_seq.len());
            let e = (feat.end as usize).min(genome_seq.len());
            if e <= s || e - s < 30 { continue; }
            let raw = &genome_seq[s..e];
            let to_translate = if feat.strand == '-' {
                reverse_complement_str(raw)
            } else {
                raw.to_string()
            };
            let aa = translate(&to_translate);
            if aa.len() < 10 { continue; }
            let id = format!("f{}", i);
            content.push_str(&format!(">{}\n{}\n", id, aa));
            feat_map.insert(id, (feat.start, feat.end, feat.strand));
        }
        if content.is_empty() {
            return Err(anyhow::anyhow!("No translatable CDS features found in GFF"));
        }
        std::fs::write(&protein_fasta, content)?;
        (protein_fasta, feat_map)
    };

    // Build DIAMOND DB from query FASTA
    let db_status = TokioCommand::new("diamond")
        .args(["makedb", "--in", query_fasta.as_str(), "--db", db_path.to_str().unwrap_or(""), "--quiet"])
        .stdout(Stdio::null()).stderr(Stdio::piped())
        .output().await
        .map_err(|e| anyhow::anyhow!("diamond not found on PATH: {}", e))?;
    if !db_status.status.success() {
        let e = String::from_utf8_lossy(&db_status.stderr);
        return Err(anyhow::anyhow!("diamond makedb failed: {}", &e[..e.len().min(300)]));
    }

    // Run blast
    let blast_cmd = if use_6ft { "blastx" } else { "blastp" };
    let db_dmnd = format!("{}.dmnd", db_path.display());
    let blast_out = TokioCommand::new("diamond")
        .args([
            blast_cmd,
            "--query", query_file.to_str().unwrap_or(""),
            "--db", &db_dmnd,
            "--outfmt", "6", "qseqid", "sseqid", "qstart", "qend", "pident", "evalue",
            "--out", out_tsv.to_str().unwrap_or(""),
            "--max-target-seqs", "1",
            "--quiet",
        ])
        .stdout(Stdio::null()).stderr(Stdio::piped())
        .output().await?;
    if !blast_out.status.success() {
        let e = String::from_utf8_lossy(&blast_out.stderr);
        return Err(anyhow::anyhow!("diamond {} failed: {}", blast_cmd, &e[..e.len().min(300)]));
    }

    // Parse results
    let tsv = std::fs::read_to_string(&out_tsv).unwrap_or_default();
    let mut results: Vec<core::Feature> = Vec::new();
    let mut seen: std::collections::HashSet<(u64, u64)> = std::collections::HashSet::new();

    if use_6ft {
        for line in tsv.lines() {
            let p: Vec<&str> = line.split('\t').collect();
            if p.len() < 6 { continue; }
            let qs: i64 = p[2].parse().unwrap_or(0);
            let qe: i64 = p[3].parse().unwrap_or(0);
            let (start, end, strand) = if qs <= qe {
                (qs as u64, qe as u64, '+')
            } else {
                (qe as u64, qs as u64, '-')
            };
            if start == 0 || end <= start || end > genome_size { continue; }
            if !seen.insert((start, end)) { continue; }
            results.push(core::Feature {
                start, end, strand,
                name: p[1].to_string(),
                locus_tag: String::new(),
                color_idx: results.len(),
                is_orf: false, noncoding: false,
                seqname: String::new(),
            });
        }
    } else {
        for line in tsv.lines() {
            let p: Vec<&str> = line.split('\t').collect();
            if p.len() < 6 { continue; }
            if let Some(&(start, end, strand)) = feat_map.get(p[0]) {
                if !seen.insert((start, end)) { continue; }
                results.push(core::Feature {
                    start, end, strand,
                    name: p[1].to_string(),
                    locus_tag: String::new(),
                    color_idx: results.len(),
                    is_orf: false, noncoding: false,
                    seqname: String::new(),
                });
            }
        }
    }

    // Cleanup temp files
    let _ = std::fs::remove_file(&query_file);
    let _ = std::fs::remove_file(&out_tsv);
    let _ = std::fs::remove_file(&db_dmnd);

    Ok(results)
}

fn blast_completions(prefix: &str) -> Vec<String> {
    let home = std::env::var("HOME").unwrap_or_default();
    let expanded = if prefix.starts_with("~/") {
        format!("{}{}", home, &prefix[1..])
    } else {
        prefix.to_string()
    };
    let (dir_part, file_prefix) = if expanded.ends_with('/') {
        (expanded.clone(), String::new())
    } else if let Some(pos) = expanded.rfind('/') {
        (expanded[..=pos].to_string(), expanded[pos + 1..].to_string())
    } else {
        (String::from("./"), expanded.clone())
    };
    // Reconstruct display prefix using original (with ~)
    let orig_dir = if prefix.starts_with("~/") {
        if let Some(pos) = prefix.rfind('/') { prefix[..=pos].to_string() }
        else { "~/".to_string() }
    } else {
        dir_part.clone()
    };

    let mut completions = Vec::new();
    if let Ok(entries) = std::fs::read_dir(&dir_part) {
        for entry in entries.flatten() {
            let fname = entry.file_name().to_string_lossy().to_string();
            if fname.starts_with('.') { continue; }
            if fname.to_ascii_lowercase().starts_with(&file_prefix.to_ascii_lowercase()) {
                let is_dir = entry.path().is_dir();
                let full = if is_dir {
                    format!("{}{}/", orig_dir, fname)
                } else {
                    format!("{}{}", orig_dir, fname)
                };
                completions.push(full);
            }
        }
    }
    completions.sort();
    completions
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

    let feat_to_browser = |feats: &[core::Feature]| -> Vec<server::BrowserFeature> {
        feats.iter().map(|f| {
            let (r,g,b) = if f.strand == '+' {
                PLUS_COLORS[f.color_idx % PLUS_COLORS.len()]
            } else {
                MINUS_COLORS[f.color_idx % MINUS_COLORS.len()]
            };
            server::BrowserFeature {
                name:      if f.name.is_empty() { f.locus_tag.clone() } else { f.name.clone() },
                locus_tag: f.locus_tag.clone(),
                start:     f.start,
                end:       f.end,
                strand:    f.strand,
                color:     format!("#{:02x}{:02x}{:02x}", r, g, b),
                noncoding: f.noncoding,
                is_orf:    f.is_orf,
            }
        }).collect()
    };

    let downsample_skew = |raw: &[f64]| -> Vec<f32> {
        if raw.is_empty() { return vec![]; }
        let target = 400.min(raw.len());
        (0..target).map(|i| raw[i * raw.len() / target] as f32).collect()
    };

    let main_features = feat_to_browser(&app.features);
    let main_gc_skew  = downsample_skew(&app.gc_skew);

    let plasmid_names:    Vec<String>              = app.plasmids.iter().map(|p| p.name.clone()).collect();
    let plasmid_sizes:    Vec<u64>                 = app.plasmids.iter().map(|p| p.genome_size).collect();
    let plasmid_features: Vec<Vec<server::BrowserFeature>> = app.plasmids.iter().map(|p| feat_to_browser(&p.features)).collect();
    let plasmid_gc_skew:  Vec<Vec<f32>>           = app.plasmids.iter().map(|p| downsample_skew(&p.gc_skew)).collect();

    let features = feat_to_browser(app.active_features());

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

    // Per-read data: fetch from indexed BAM when viewport is narrow enough
    const READ_FETCH_BP: u64 = 50_000;
    let (reads_p, reads_m) = if app.active_genome == 0 {
        if let Some(ref bam) = app.bam_path {
            let vs = app.active_view_start();
            let ve = app.active_view_end();
            if ve.saturating_sub(vs) <= READ_FETCH_BP {
                match fetch_reads_in_region(bam, &app.seqname_offsets, vs, ve, 8000) {
                    Ok(rds) => {
                        let rp = rds.iter().filter(|r| !r.2).map(|r| [r.0, r.1]).collect();
                        let rm = rds.iter().filter(|r|  r.2).map(|r| [r.0, r.1]).collect();
                        (rp, rm)
                    }
                    Err(_) => (vec![], vec![]),
                }
            } else { (vec![], vec![]) }
        } else { (vec![], vec![]) }
    } else { (vec![], vec![]) };

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

    // Stop codons: send viewport-context positions (±1 span) as 1-based u32, main genome only.
    // Plasmids fall back to JS sequence-window scanning (they're typically small).
    let (stop_codons_plus, stop_codons_minus) = if app.active_genome == 0 && !app.stop_codons.is_empty() {
        let vs = app.active_view_start() as usize;
        let ve = app.active_view_end() as usize;
        let span = ve.saturating_sub(vs);
        let ctx_s = vs.saturating_sub(span);
        let ctx_e = ve + span;
        let slice_for = |strand: char| -> Vec<Vec<u32>> {
            (0u8..3).map(|fr| {
                app.stop_codons.get(&(strand, fr))
                    .map(|v| {
                        let lo = v.partition_point(|&p| p < ctx_s);
                        let hi = v.partition_point(|&p| p <= ctx_e);
                        v[lo..hi].iter().map(|&p| p as u32).collect()
                    })
                    .unwrap_or_default()
            }).collect()
        };
        (slice_for('+'), slice_for('-'))
    } else {
        // Empty outer vec (length 0) signals JS to use client-side fallback scan.
        // vec![[],[],[]] has length 3 which the JS check misreads as valid server data.
        (vec![], vec![])
    };

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
        reads_plus:  reads_p,
        reads_minus: reads_m,
        gc_skew: vec![],  // redundant with main_gc_skew/plasmid_gc_skew; kept for compatibility
        input_files: app.source_files.clone(),
        msa_gene: app.msa.as_ref()
            .filter(|m| !m.loading && m.error.is_none() && !m.sequences.is_empty())
            .map(|m| m.gene_name.clone())
            .unwrap_or_default(),
        msa_sequences: app.msa.as_ref()
            .filter(|m| !m.loading && m.error.is_none())
            .map(|m| m.sequences.iter().map(|(id, seq)| [id.clone(), seq.clone()]).collect())
            .unwrap_or_default(),
        msa_error: app.msa.as_ref()
            .filter(|m| !m.loading)
            .and_then(|m| m.error.clone())
            .unwrap_or_default(),
        status_msg: if app.protein.as_ref().map(|p| p.folding).unwrap_or(false) {
            "Folding protein\u{2026}".into()
        } else if app.msa.as_ref().map(|m| m.loading).unwrap_or(false) {
            "Building MSA\u{2026}".into()
        } else {
            app.status_msg.clone()
        },
        needs_dmnd_db: app.needs_dmnd_db,
        path_completions: app.path_completions.clone(),
        active_genome: app.active_genome,
        main_name:     app.genome_name.clone(),
        main_size:     app.genome_size,
        main_features,
        main_gc_skew,
        plasmid_names,
        plasmid_sizes,
        plasmid_features,
        plasmid_gc_skew,
        stop_codons_plus,
        stop_codons_minus,
        blast_features: app.blast_features.iter().map(|f| {
            server::BrowserFeature {
                name:      if f.name.is_empty() { f.locus_tag.clone() } else { f.name.clone() },
                locus_tag: f.locus_tag.clone(),
                start:     f.start,
                end:       f.end,
                strand:    f.strand,
                color:     "#e040fb".to_string(),
                noncoding: false,
                is_orf:    false,
            }
        }).collect(),
        contig_names:  {
            let mut sv: Vec<(&String, &u64)> = app.seqname_offsets.iter().collect();
            sv.sort_by_key(|(_, &off)| off);
            sv.iter().map(|(n, _)| (*n).clone()).collect()
        },
        contig_starts: {
            let mut sv: Vec<(&String, &u64)> = app.seqname_offsets.iter().collect();
            sv.sort_by_key(|(_, &off)| off);
            sv.iter().map(|(_, &off)| off).collect()
        },
        contig_ends: {
            let mut sv: Vec<(&String, &u64)> = app.seqname_offsets.iter().collect();
            sv.sort_by_key(|(_, &off)| off);
            let mut ends = Vec::with_capacity(sv.len());
            for i in 0..sv.len() {
                if i + 1 < sv.len() {
                    // gap = 20k Ns before each subsequent contig
                    ends.push(sv[i + 1].1.saturating_sub(20_000));
                } else {
                    ends.push(app.genome_size);
                }
            }
            ends
        },
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
