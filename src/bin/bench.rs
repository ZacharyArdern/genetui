use std::time::Instant;
use genetui::core::{parse_gff, parse_fasta, find_stop_codons, Feature};
use std::collections::HashMap;

const WIDTH: usize = 200;
const LABEL_W: usize = 4;
const FEAT_W: usize = WIDTH - LABEL_W;
const N_ITERS: usize = 500;

fn bp_to_col(bp: u64, view_start: u64, span: f64) -> usize {
    let col = ((bp as f64 - view_start as f64) / span * (FEAT_W - 1) as f64).round() as i64;
    col.max(0).min(FEAT_W as i64 - 1) as usize
}

fn render_viewport(
    features: &[Feature],
    stop_codons: &HashMap<(char, u8), Vec<usize>>,
    view_start: u64,
    view_end: u64,
) {
    let span = (view_end - view_start) as f64;
    if span <= 0.0 { return; }

    // bucket visible features
    let mut buckets: HashMap<(char, u8), Vec<(usize, usize, usize)>> = HashMap::new();
    for (i, f) in features.iter().enumerate() {
        if f.end < view_start || f.start > view_end { continue; }
        let cs = bp_to_col(f.start.max(view_start), view_start, span);
        let ce = bp_to_col(f.end.min(view_end), view_start, span);
        let fr = if f.strand == '+' {
            ((f.start as usize).saturating_sub(1) % 3) as u8
        } else {
            ((f.end as usize).saturating_sub(1) % 3) as u8
        };
        buckets.entry((f.strand, fr)).or_default().push((i, cs, ce));
    }

    // stop codon columns
    let mut stop_cols: HashMap<(char, u8), Vec<usize>> = HashMap::new();
    for &s in &['+', '-'] {
        for fr in 0u8..3 {
            let cols: Vec<usize> = stop_codons.get(&(s, fr)).map(|v| {
                v.iter()
                    .filter(|&&p| p as u64 >= view_start && (p as u64) <= view_end)
                    .map(|&p| bp_to_col(p as u64, view_start, span))
                    .collect()
            }).unwrap_or_default();
            stop_cols.insert((s, fr), cols);
        }
    }

    // render 7 rows
    let rows: &[(Option<char>, u8)] = &[
        (Some('+'), 0), (Some('+'), 1), (Some('+'), 2),
        (None, 0), // ruler
        (Some('-'), 0), (Some('-'), 1), (Some('-'), 2),
    ];
    for (strand_opt, fr) in rows {
        let mut cells = vec![(' ', false); FEAT_W]; // (char, colored)
        if let Some(strand) = strand_opt {
            for &col in stop_cols.get(&(*strand, *fr)).unwrap_or(&vec![]) {
                if col < FEAT_W { cells[col] = ('│', false); }
            }
            for &(_fi, cs, ce) in buckets.get(&(*strand, *fr)).unwrap_or(&vec![]) {
                let fw = (ce + 1).saturating_sub(cs).max(1);
                let tip = if *strand == '+' { '❯' } else { '❮' };
                for j in 0..fw {
                    let col = cs + j;
                    if col < FEAT_W {
                        let ch = if j == if *strand == '+' { fw - 1 } else { 0 } { tip } else { ' ' };
                        cells[col] = (ch, true);
                    }
                }
            }
        } else {
            // ruler: just fill
            let _ruler: String = (0..FEAT_W).map(|_| '─').collect();
        }
        // simulate span building
        let _line: String = cells.iter().map(|(c, _)| c).collect();
    }
}

fn main() {
    let gff   = "/Users/zacharyardern/Science/Programs/genetui/GCF_000005845.2_ASM584v2_genomic.gff.gz";
    let fasta = "/Users/zacharyardern/Science/Programs/genetui/GCF_000005845.2_ASM584v2_genomic.fna";

    println!("Loading data...");
    let parsed_fasta = parse_fasta(fasta).expect("fasta");
    let (features, genome_size) = parse_gff(gff, &parsed_fasta.seqname_offsets).expect("gff");
    let stop_codons = find_stop_codons(&parsed_fasta.main_seq);
    println!("  {} features, genome {} bp", features.len(), genome_size);

    // warmup
    for _ in 0..10 {
        render_viewport(&features, &stop_codons, 1, 50_000);
    }

    // benchmark
    let step = genome_size / N_ITERS as u64;
    let mut times: Vec<f64> = Vec::with_capacity(N_ITERS);
    for i in 0..N_ITERS {
        let vs = 1 + i as u64 * step;
        let ve = vs + 50_000;
        let t0 = Instant::now();
        render_viewport(&features, &stop_codons, vs, ve);
        times.push(t0.elapsed().as_secs_f64() * 1000.0);
    }

    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mean = times.iter().sum::<f64>() / times.len() as f64;
    let p50  = times[times.len() / 2];
    let p99  = times[(times.len() as f64 * 0.99) as usize];
    let fps  = 1000.0 / mean;

    println!("\nRust render benchmark ({N_ITERS} viewports, width={WIDTH}):");
    println!("  mean:  {mean:.3} ms");
    println!("  p50:   {p50:.3} ms");
    println!("  p99:   {p99:.3} ms");
    println!("  → max theoretical scroll rate: {fps:.0} fps");
}
