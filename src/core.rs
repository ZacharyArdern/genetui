use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::Result;
use flate2::read::GzDecoder;

/// A genomic feature (gene or ORF) parsed from a GFF3 file.
#[derive(Debug, Clone)]
pub struct Feature {
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub name: String,
    pub locus_tag: String,
    pub color_idx: usize,
    pub is_orf: bool,
    pub noncoding: bool,
    pub seqname: String,
}

/// One contig identified as a plasmid (contains "plasmid" in header).
pub struct PlasmidEntry {
    pub name: String,
    pub seqname: String,
    pub sequence: String,
    pub offset_in_main: u64,
}

/// Result of parsing a (possibly multi-contig) FASTA file.
pub struct ParsedFasta {
    /// All contigs concatenated with 10 kb N gaps.
    pub main_seq: String,
    /// Display name from the first non-plasmid header (or first header).
    pub main_name: String,
    /// seqname (first word of header) → byte offset in main_seq.
    pub seqname_offsets: HashMap<String, u64>,
    /// Plasmid contigs, stored separately for their own circle maps.
    pub plasmids: Vec<PlasmidEntry>,
    pub total_len: u64,
}

/// Open a file as a boxed BufRead, handling optional .gz compression.
fn open_reader(path: &str) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if path.ends_with(".gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Parse a GFF3 file (plain or .gz). Returns (features, genome_size).
/// Only "gene" and "ORF" type rows are returned.
/// When `seqname_offsets` is non-empty, feature coordinates are shifted by the
/// per-contig offset so they are global within the concatenated sequence; features
/// on seqnames not present in the map are skipped.
/// When `seqname_offsets` is empty, no offset is applied and genome_size is taken
/// from GFF `region` tags (single-contig / no-FASTA mode).
pub fn parse_gff(path: &str, seqname_offsets: &HashMap<String, u64>) -> Result<(Vec<Feature>, u64)> {
    let reader = open_reader(path)?;
    let mut features = Vec::new();
    let mut genome_size: u64 = 0;
    let mut color_counter: usize = 0;
    let use_offsets = !seqname_offsets.is_empty();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.splitn(9, '\t').collect();
        if parts.len() < 9 {
            continue;
        }
        let seqname = parts[0];
        let ftype = parts[2];
        let local_start: u64 = parts[3].parse().unwrap_or(0);
        let local_end:   u64 = parts[4].parse().unwrap_or(0);
        let strand = parts[6].chars().next().unwrap_or('+');
        let attrs = parts[8];

        if ftype == "region" {
            if !use_offsets && local_end > genome_size {
                genome_size = local_end;
            }
            continue;
        }
        if ftype != "gene" && ftype != "ORF" {
            continue;
        }

        let (start, end) = if use_offsets {
            match seqname_offsets.get(seqname) {
                Some(&off) => (local_start + off, local_end + off),
                None => continue,
            }
        } else {
            (local_start, local_end)
        };

        if end > genome_size {
            genome_size = end;
        }

        let is_orf = ftype == "ORF";
        let mut name = String::new();
        let mut locus_tag = String::new();
        let mut gene_biotype = String::new();

        // Detect GTF vs GFF3 by attribute format: GTF uses key "value"; GFF3 uses key=value
        let is_gtf = attrs.contains('"');
        if is_gtf {
            for attr in attrs.split(';') {
                let attr = attr.trim();
                let mut parts = attr.splitn(2, ' ');
                let key = parts.next().unwrap_or("").trim();
                let val = parts.next().unwrap_or("").trim().trim_matches('"');
                match key {
                    "gene_name" | "Name"     => name       = val.to_string(),
                    "gene_id"  | "locus_tag" => locus_tag  = val.to_string(),
                    "gene_biotype" | "gene_type" => gene_biotype = val.to_string(),
                    _ => {}
                }
            }
        } else {
            for attr in attrs.split(';') {
                if let Some(val) = attr.strip_prefix("Name=") {
                    name = val.to_string();
                } else if let Some(val) = attr.strip_prefix("locus_tag=") {
                    locus_tag = val.to_string();
                } else if let Some(val) = attr.strip_prefix("gene_biotype=") {
                    gene_biotype = val.to_string();
                }
            }
        }
        if name.is_empty() {
            name = locus_tag.clone();
        }

        const NC_BIOTYPES: &[&str] = &[
            "rRNA", "tRNA", "tmRNA", "ncRNA", "antisense_RNA",
            "lnc_RNA", "misc_RNA", "ribozyme", "sRNA",
        ];
        let noncoding = !is_orf && NC_BIOTYPES.iter().any(|&b| gene_biotype == b);

        features.push(Feature {
            start,
            end,
            strand,
            name,
            locus_tag,
            color_idx: color_counter,
            is_orf,
            noncoding,
            seqname: seqname.to_string(),
        });
        color_counter += 1;
    }

    Ok((features, genome_size))
}

/// Parse a FASTA file (plain or .gz), handling multiple contigs.
///
/// All contigs are concatenated into `main_seq` with 10 kb N gaps between them.
/// Contigs whose header contains "plasmid" (case-insensitive) are additionally
/// stored in `plasmids` so the UI can draw a separate circle map for each.
pub fn parse_fasta(path: &str) -> Result<ParsedFasta> {
    let reader = open_reader(path)?;

    // Collect all (header, sequence) entries
    let mut entries: Vec<(String, String)> = Vec::new();
    let mut cur_header = String::new();
    let mut cur_seq    = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !cur_header.is_empty() {
                entries.push((cur_header.clone(), std::mem::take(&mut cur_seq)));
            }
            cur_header = line[1..].trim().to_string();
        } else {
            cur_seq.push_str(line.trim());
        }
    }
    if !cur_header.is_empty() {
        entries.push((cur_header, cur_seq));
    }

    if entries.is_empty() {
        return Ok(ParsedFasta {
            main_seq: String::new(),
            main_name: String::new(),
            seqname_offsets: HashMap::new(),
            plasmids: Vec::new(),
            total_len: 0,
        });
    }

    // Single entry: fast path (no gaps, no plasmid split)
    if entries.len() == 1 {
        let (header, seq) = entries.into_iter().next().unwrap();
        let seq = seq.to_uppercase();
        let seqname = header.split_whitespace().next().unwrap_or("").to_string();
        let total_len = seq.len() as u64;
        let mut seqname_offsets = HashMap::new();
        seqname_offsets.insert(seqname, 0u64);
        return Ok(ParsedFasta {
            main_seq: seq,
            main_name: header,
            seqname_offsets,
            plasmids: Vec::new(),
            total_len,
        });
    }

    // Multiple entries: concatenate with 20 kb N gaps
    const GAP: usize = 20_000;
    let gap_str = "N".repeat(GAP);

    let mut main_seq = String::new();
    let mut seqname_offsets: HashMap<String, u64> = HashMap::new();
    let mut plasmids: Vec<PlasmidEntry> = Vec::new();
    let mut main_name = String::new();

    for (i, (header, seq)) in entries.iter().enumerate() {
        if i > 0 {
            main_seq.push_str(&gap_str);
        }
        let offset = main_seq.len() as u64;
        let seqname = header.split_whitespace().next().unwrap_or("").to_string();
        seqname_offsets.insert(seqname.clone(), offset);

        let seq_upper = seq.to_uppercase();
        let is_plasmid = header.to_ascii_lowercase().contains("plasmid");

        if is_plasmid {
            plasmids.push(PlasmidEntry {
                name: header.clone(),
                seqname: seqname.clone(),
                sequence: seq_upper.clone(),
                offset_in_main: offset,
            });
        } else if main_name.is_empty() {
            main_name = header.clone();
        }

        main_seq.push_str(&seq_upper);
    }

    if main_name.is_empty() {
        main_name = entries[0].0.clone();
    }

    let total_len = main_seq.len() as u64;
    Ok(ParsedFasta { main_seq, main_name, seqname_offsets, plasmids, total_len })
}

/// Reverse-complement a DNA string (public for use in main).
pub fn reverse_complement_str(seq: &str) -> String {
    reverse_complement(seq)
}

/// Reverse-complement a DNA string.
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'a' => 't',
            't' => 'a',
            'g' => 'c',
            'c' => 'g',
            other => other,
        })
        .collect()
}

/// Returns a map of (strand, frame) → sorted 1-based positions of stop codons.
/// Stop codons: TAA, TAG, TGA.
/// For + strand: frame = pos0 % 3, position = pos0 + 1 (1-based).
/// For - strand: frame derived from rightmost coordinate, position = leftmost 1-based.
pub fn find_stop_codons(seq: &str) -> HashMap<(char, u8), Vec<usize>> {
    let mut result: HashMap<(char, u8), Vec<usize>> = HashMap::new();
    for &s in &['+', '-'] {
        for fr in 0u8..3 {
            result.insert((s, fr), Vec::new());
        }
    }

    let upper: Vec<u8> = seq.to_uppercase().bytes().collect();
    let n = upper.len();

    // + strand
    let mut i = 0;
    while i + 2 < n {
        let codon = &upper[i..i + 3];
        if matches!(codon, b"TAA" | b"TAG" | b"TGA") {
            let frame = (i % 3) as u8;
            result.get_mut(&('+', frame)).unwrap().push(i + 1);
        }
        i += 1;
    }

    // - strand: search reverse complement
    let rc = reverse_complement(seq).to_uppercase();
    let rc_bytes: Vec<u8> = rc.bytes().collect();
    let mut i = 0;
    while i + 2 < n {
        let codon = &rc_bytes[i..i + 3];
        if matches!(codon, b"TAA" | b"TAG" | b"TGA") {
            // frame on minus strand: based on position from right end
            let fr = ((n - i - 1) % 3) as u8;
            // leftmost 1-based coordinate on + strand
            let pos = n - i - 2;
            result.get_mut(&('-', fr)).unwrap().push(pos);
        }
        i += 1;
    }

    // Sort all
    for v in result.values_mut() {
        v.sort_unstable();
    }

    result
}

/// Translate a DNA sequence (5'→3') to protein using the standard genetic code.
/// Stops at the first stop codon.
pub fn translate(dna: &str) -> String {
    let dna = dna.to_uppercase();
    let bytes = dna.as_bytes();
    let mut aa = String::new();
    let mut i = 0;
    while i + 2 < bytes.len() {
        let codon = &bytes[i..i + 3];
        let residue = codon_to_aa(codon);
        if residue == '*' {
            break;
        }
        aa.push(residue);
        i += 3;
    }
    aa
}

fn codon_to_aa(codon: &[u8]) -> char {
    match codon {
        b"TTT" | b"TTC" => 'F',
        b"TTA" | b"TTG" => 'L',
        b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        b"ATT" | b"ATC" | b"ATA" => 'I',
        b"ATG" => 'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" => 'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        b"TAT" | b"TAC" => 'Y',
        b"TAA" | b"TAG" => '*',
        b"CAT" | b"CAC" => 'H',
        b"CAA" | b"CAG" => 'Q',
        b"AAT" | b"AAC" => 'N',
        b"AAA" | b"AAG" => 'K',
        b"GAT" | b"GAC" => 'D',
        b"GAA" | b"GAG" => 'E',
        b"TGT" | b"TGC" => 'C',
        b"TGA" => '*',
        b"TGG" => 'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" => 'R',
        b"AGT" | b"AGC" => 'S',
        b"AGA" | b"AGG" => 'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
        _ => 'X',
    }
}

/// Compute GC skew over n_windows equal windows covering the genome.
/// Returns (G - C) / (G + C) per window. Returns 0.0 if G+C == 0.
pub fn compute_gc_skew(seq: &str, n_windows: usize) -> Vec<f64> {
    let seq = seq.as_bytes();
    let len = seq.len();
    if len == 0 || n_windows == 0 {
        return vec![0.0; n_windows];
    }
    let win = (len / n_windows).max(1);
    (0..n_windows)
        .map(|i| {
            let start = i * win;
            let end = ((i + 1) * win).min(len);
            let slice = &seq[start..end];
            let g = slice.iter().filter(|&&b| b == b'G' || b == b'g').count() as f64;
            let c = slice.iter().filter(|&&b| b == b'C' || b == b'c').count() as f64;
            if g + c == 0.0 { 0.0 } else { (g - c) / (g + c) }
        })
        .collect()
}

/// Format bp position as human-readable string.
pub fn format_bp(bp: u64) -> String {
    if bp >= 1_000_000 {
        format!("{:.2}M", bp as f64 / 1_000_000.0)
    } else if bp >= 1_000 {
        format!("{:.1}k", bp as f64 / 1_000.0)
    } else {
        bp.to_string()
    }
}

/// Like format_bp but uses enough decimal places so adjacent ticks `tick` apart are distinct.
pub fn format_bp_tick(bp: u64, tick: u64) -> String {
    if bp >= 1_000_000 {
        let prec: usize = match tick {
            t if t >= 1_000_000 => 0,
            t if t >= 100_000   => 1,
            t if t >= 10_000    => 2,
            t if t >= 1_000     => 3,
            _                   => 4,
        };
        format!("{:.prec$}M", bp as f64 / 1_000_000.0, prec = prec)
    } else if bp >= 1_000 {
        let prec: usize = match tick {
            t if t >= 1_000 => 0,
            t if t >= 100   => 1,
            _               => 2,
        };
        format!("{:.prec$}k", bp as f64 / 1_000.0, prec = prec)
    } else {
        bp.to_string()
    }
}

/// Per-strand read coverage in equal-width bins across the genome.
pub struct StrandCoverage {
    pub plus:  Vec<u32>,
    pub minus: Vec<u32>,
    pub bin_size: u64,
    pub genome_size: u64,
}


/// Build a BAM index (.bai) if one does not already exist alongside the BAM.
/// Ensure an indexed BAM exists and return the path to use.
/// If the BAM is unsorted, sorts it first (creating a .sorted.bam beside it).
/// Returns the path that should actually be opened (may differ from input).
pub fn ensure_bam_index(bam_path: &str) -> Result<String> {
    let bai = format!("{}.bai", bam_path);
    let csi = format!("{}.csi", bam_path);
    if std::path::Path::new(&bai).exists() || std::path::Path::new(&csi).exists() {
        return Ok(bam_path.to_string());
    }

    eprintln!("Building BAM index for {} ...", bam_path);

    // Try indexing directly (works if BAM is already sorted)
    let indexed = std::process::Command::new("samtools")
        .args(["index", bam_path])
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    if indexed && (std::path::Path::new(&bai).exists() || std::path::Path::new(&csi).exists()) {
        eprintln!("BAM index built.");
        return Ok(bam_path.to_string());
    }

    // BAM is likely unsorted — sort first, then index
    let sorted_path = if bam_path.ends_with(".bam") {
        format!("{}.sorted.bam", &bam_path[..bam_path.len()-4])
    } else {
        format!("{}.sorted.bam", bam_path)
    };

    if !std::path::Path::new(&sorted_path).exists() {
        eprintln!("BAM appears unsorted, sorting to {} ...", sorted_path);
        let ok = std::process::Command::new("samtools")
            .args(["sort", "-o", &sorted_path, bam_path])
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if !ok {
            anyhow::bail!("samtools sort failed for {}", bam_path);
        }
        eprintln!("BAM sorted.");
    }

    let sorted_bai = format!("{}.bai", sorted_path);
    let sorted_csi = format!("{}.csi", sorted_path);
    if !std::path::Path::new(&sorted_bai).exists() && !std::path::Path::new(&sorted_csi).exists() {
        let ok = std::process::Command::new("samtools")
            .args(["index", &sorted_path])
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if !ok {
            anyhow::bail!("samtools index failed for {}", sorted_path);
        }
        eprintln!("Sorted BAM indexed.");
    }

    Ok(sorted_path)
}

/// Fetch individual reads overlapping [vs, ve) from an indexed BAM.
/// Returns (global_start, global_end, is_reverse) tuples sorted by start, capped at max_reads.
pub fn fetch_reads_in_region(
    bam_path: &str,
    seqname_offsets: &HashMap<String, u64>,
    vs: u64,
    ve: u64,
    max_reads: usize,
) -> Result<Vec<(u32, u32, bool)>> {
    use rust_htslib::bam::{IndexedReader, Read};
    // Avoid htslib printing "could not retrieve index file" to stderr
    let bai = format!("{}.bai", bam_path);
    let csi = format!("{}.csi", bam_path);
    if !std::path::Path::new(&bai).exists() && !std::path::Path::new(&csi).exists() {
        return Ok(vec![]);
    }
    let mut reader = IndexedReader::from_path(bam_path)
        .map_err(|e| anyhow::anyhow!("IndexedReader: {}", e))?;
    let header = rust_htslib::bam::Header::from_template(reader.header());
    let hview  = rust_htslib::bam::HeaderView::from_header(&header);
    let n_refs  = hview.target_count() as usize;

    let ref_offsets: Vec<Option<u64>> = (0..n_refs).map(|i| {
        let name = std::str::from_utf8(hview.tid2name(i as u32)).ok()?;
        if seqname_offsets.is_empty() { Some(0u64) }
        else { seqname_offsets.get(name).copied() }
    }).collect();

    let mut reads: Vec<(u32, u32, bool)> = Vec::new();

    for (tid, off_opt) in ref_offsets.iter().enumerate() {
        let offset = match off_opt { Some(o) => *o, None => continue };
        let local_s = if vs > offset { (vs - offset) as i64 } else { 0i64 };
        let local_e = if ve > offset { (ve - offset) as i64 } else { continue };
        if local_e <= local_s { continue; }
        let _ = reader.fetch((tid as u32, local_s, local_e));
        for result in reader.records() {
            if reads.len() >= max_reads { break; }
            let rec = match result { Ok(r) => r, Err(_) => continue };
            if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() { continue; }
            let pos = rec.pos(); if pos < 0 { continue; }
            let gs = offset + pos as u64;
            let ge = gs + rec.seq_len() as u64;
            if ge < vs || gs >= ve { continue; }
            reads.push((gs as u32, ge as u32, rec.is_reverse()));
        }
        if reads.len() >= max_reads { break; }
    }
    reads.sort_by_key(|r| r.0);
    Ok(reads)
}

/// Compute BAM/SAM/CRAM coverage using rust-htslib (htslib C bindings).
pub fn compute_alignment_coverage(
    path: &str,
    ref_path: Option<&str>,
    seqname_offsets: &HashMap<String, u64>,
    genome_size: u64,
) -> Result<StrandCoverage> {
    use rust_htslib::bam::{Read, Reader};

    const BIN: u64 = 1000;
    let n = ((genome_size + BIN - 1) / BIN) as usize;
    let mut plus  = vec![0u32; n];
    let mut minus = vec![0u32; n];

    let n_threads = std::thread::available_parallelism()
        .map(|p| p.get()).unwrap_or(1).saturating_sub(1);

    // Prevent htslib from going to the network (EBI) to fetch CRAM references.
    // We only need alignment positions, not sequence data, so no reference is needed.
    // Setting REF_PATH to "." restricts htslib to local files only.
    unsafe { std::env::set_var("REF_PATH", "."); }

    let mut reader = Reader::from_path(path)?;
    if let Some(r) = ref_path {
        // Only pass the reference if it has a usable .fai index alongside it.
        // Plain gzip files (.gz without .fai) cannot be used by htslib as CRAM references.
        let fai = format!("{}.fai", r);
        if std::path::Path::new(&fai).exists() {
            let _ = reader.set_reference(r);
        } else if r.ends_with(".gz") {
            eprintln!(
                "Note: skipping gzip-compressed FASTA as CRAM reference (no .fai index). \
                 To enable: gunzip '{}' && samtools faidx '{}'",
                r, r.trim_end_matches(".gz")
            );
        } else {
            let _ = reader.set_reference(r);
        }
    }
    if n_threads > 0 {
        reader.set_threads(n_threads)?;
    }

    let header = rust_htslib::bam::Header::from_template(reader.header());
    let hview  = rust_htslib::bam::HeaderView::from_header(&header);

    // Pre-build tid → global offset Vec.
    let n_refs = hview.target_count() as usize;
    let ref_offsets: Vec<Option<u64>> = (0..n_refs)
        .map(|i| {
            let name = std::str::from_utf8(hview.tid2name(i as u32)).ok()?;
            if seqname_offsets.is_empty() { Some(0u64) }
            else { seqname_offsets.get(name).copied() }
        })
        .collect();

    if !seqname_offsets.is_empty() && ref_offsets.iter().all(|o| o.is_none()) {
        let bam_refs: Vec<_> = (0..n_refs)
            .filter_map(|i| std::str::from_utf8(hview.tid2name(i as u32)).ok().map(|s| s.to_string()))
            .collect();
        eprintln!(
            "Warning: BAM reference names ({}) don't match FASTA sequence names — coverage will be empty.",
            bam_refs.join(", ")
        );
    }

    let mut decode_errors = 0u32;
    for result in reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => { decode_errors += 1; continue; }
        };
        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            continue;
        }
        let tid = record.tid();
        if tid < 0 { continue; }
        let offset = match ref_offsets.get(tid as usize).and_then(|o| *o) {
            Some(o) => o,
            None => continue,
        };
        let pos = record.pos();
        if pos < 0 { continue; }
        let bin = ((offset + pos as u64) / BIN) as usize;
        if bin >= n { continue; }
        if record.is_reverse() {
            minus[bin] = minus[bin].saturating_add(1);
        } else {
            plus[bin] = plus[bin].saturating_add(1);
        }
    }

    let total_reads = plus.iter().map(|&v| v as u64).sum::<u64>()
        + minus.iter().map(|&v| v as u64).sum::<u64>();
    if decode_errors > 0 {
        eprintln!(
            "Warning: {decode_errors} CRAM record(s) failed to decode (reference mismatch?). \
             {} reads counted.",
            total_reads
        );
    }
    if total_reads == 0 && decode_errors == 0 {
        eprintln!("Warning: coverage is all-zero — BAM/CRAM may be empty or reference names may not match.");
    }

    Ok(StrandCoverage { plus, minus, bin_size: BIN, genome_size })
}

/// Calculate a nice tick spacing for a ruler.
pub fn nice_tick_spacing(span: u64, width: usize) -> u64 {
    let raw = span as f64 / (width as f64 / 12.0).max(1.0);
    let mag = 10u64.pow(raw.max(1.0).log10().floor() as u32);
    for &m in &[1u64, 2, 5, 10] {
        if mag * m >= raw as u64 {
            return (mag * m).max(1);
        }
    }
    (mag * 10).max(1)
}
