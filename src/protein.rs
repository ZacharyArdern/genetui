use std::collections::HashMap;
use ratatui::style::Color;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Secondary { Helix, Sheet, Coil }

#[derive(Clone, Debug)]
pub struct ProteinAtom {
    pub residue: u32,
    /// Normalised coords: centred at origin, max radius ≈ 0.85.
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub plddt: f32,
    pub secondary: Secondary,
}

/// Parse Cα atoms from a PDB string (ESMFold/minifold output).
/// Secondary structure is taken from HELIX/SHEET records when present,
/// otherwise assigned geometrically from Cα coordinates (P-SEA-style).
pub fn parse_pdb(pdb: &str) -> Vec<ProteinAtom> {
    let pdb_ss = parse_secondary(pdb);

    let mut raw: Vec<(u32, f64, f64, f64, f32)> = Vec::new();
    for line in pdb.lines() {
        if !line.starts_with("ATOM") || line.len() < 66 { continue; }
        if line.get(12..16).map(|s| s.trim()) != Some("CA") { continue; }
        let res:   u32 = parse_col(line, 22, 26);
        let x:     f64 = parse_col_f(line, 30, 38);
        let y:     f64 = parse_col_f(line, 38, 46);
        let z:     f64 = parse_col_f(line, 46, 54);
        let plddt: f32 = parse_col_f(line, 60, 66) as f32;
        raw.push((res, x, y, z, plddt));
    }
    if raw.is_empty() { return Vec::new(); }

    // Assign SS: prefer PDB records, fall back to geometry when absent.
    let geo_ss = if pdb_ss.is_empty() { assign_secondary_ca(&raw) } else { vec![] };
    let ss_for = |i: usize, res: u32| -> Secondary {
        if !pdb_ss.is_empty() {
            pdb_ss.get(&res).copied().unwrap_or(Secondary::Coil)
        } else {
            geo_ss.get(i).copied().unwrap_or(Secondary::Coil)
        }
    };

    let n = raw.len() as f64;
    let cx = raw.iter().map(|a| a.1).sum::<f64>() / n;
    let cy = raw.iter().map(|a| a.2).sum::<f64>() / n;
    let cz = raw.iter().map(|a| a.3).sum::<f64>() / n;
    let max_r = raw.iter()
        .map(|a| ((a.1-cx).powi(2) + (a.2-cy).powi(2) + (a.3-cz).powi(2)).sqrt())
        .fold(1.0f64, f64::max);
    let scale = 0.85 / max_r;

    raw.into_iter().enumerate().map(|(i, (res, x, y, z, plddt))| ProteinAtom {
        residue: res,
        x: (x - cx) * scale,
        y: (y - cy) * scale,
        z: (z - cz) * scale,
        plddt,
        secondary: ss_for(i, res),
    }).collect()
}

/// Cα-only secondary structure assignment (P-SEA-style).
///
/// Helix:  d(Cα_i, Cα_{i+4}) ∈ [4.2, 5.8] Å — characteristic of α-helix (~5.1 Å ideal).
///         Requires ≥ 5 consecutive residues in the detected span.
/// Sheet:  d(Cα_i, Cα_{i+2}) > 6.0 Å — extended conformation (helix gives ~5.5 Å).
///         Requires ≥ 3 consecutive residues; helices take priority.
fn assign_secondary_ca(raw: &[(u32, f64, f64, f64, f32)]) -> Vec<Secondary> {
    let n = raw.len();
    let mut ss = vec![Secondary::Coil; n];

    let d = |i: usize, j: usize| -> f64 {
        let (_, xi, yi, zi, _) = raw[i];
        let (_, xj, yj, zj, _) = raw[j];
        ((xi-xj).powi(2) + (yi-yj).powi(2) + (zi-zj).powi(2)).sqrt()
    };

    // --- Helix: use d(i, i+3) as primary criterion.
    // Minifold produces helices where d(i,i+4) ≈ 6.0–6.3 Å (above the ideal
    // 5.1 Å), but d(i,i+3) sits cleanly at 5.0–5.6 Å vs >7.5 Å outside.
    // Threshold [4.5, 6.2] covers α-helices and 3₁₀-helices robustly.
    let mut hcand = vec![false; n];
    for i in 0..n.saturating_sub(3) {
        let v = d(i, i + 3);
        if v >= 4.5 && v <= 6.2 {
            for k in i..=(i + 3) { hcand[k] = true; }
        }
    }
    // Require minimum run of 4 (a 3-residue i+3 span touches 4 atoms)
    let mut i = 0;
    while i < n {
        if hcand[i] {
            let start = i;
            while i < n && hcand[i] { i += 1; }
            if i - start >= 4 {
                for k in start..i { ss[k] = Secondary::Helix; }
            }
        } else { i += 1; }
    }

    // --- Sheet: extended d(i, i+2) > 6.0 Å, not already helix ---
    let mut scand = vec![false; n];
    for i in 0..n.saturating_sub(2) {
        if ss[i] == Secondary::Coil && ss[i+2] == Secondary::Coil && d(i, i+2) > 6.0 {
            scand[i] = true; scand[i+1] = true; scand[i+2] = true;
        }
    }
    let mut i = 0;
    while i < n {
        if scand[i] {
            let start = i;
            while i < n && scand[i] { i += 1; }
            if i - start >= 3 {
                for k in start..i {
                    if ss[k] == Secondary::Coil { ss[k] = Secondary::Sheet; }
                }
            }
        } else { i += 1; }
    }

    ss
}

fn parse_col<T: std::str::FromStr>(line: &str, s: usize, e: usize) -> T
where T: Default { line.get(s..e).and_then(|x| x.trim().parse().ok()).unwrap_or_default() }

fn parse_col_f(line: &str, s: usize, e: usize) -> f64 {
    line.get(s..e).and_then(|x| x.trim().parse().ok()).unwrap_or(0.0)
}

fn parse_secondary(pdb: &str) -> HashMap<u32, Secondary> {
    let mut map = HashMap::new();
    for line in pdb.lines() {
        if line.starts_with("HELIX") && line.len() >= 37 {
            let s: u32 = parse_col(line, 21, 25);
            let e: u32 = parse_col(line, 33, 37);
            for r in s..=e { map.insert(r, Secondary::Helix); }
        } else if line.starts_with("SHEET") && line.len() >= 37 {
            let s: u32 = parse_col(line, 22, 26);
            let e: u32 = parse_col(line, 33, 37);
            for r in s..=e { map.insert(r, Secondary::Sheet); }
        }
    }
    map
}

/// Rotate by rot_y (around Y), then rot_x (around X). Returns (screen_x, screen_y, depth).
pub fn project(x: f64, y: f64, z: f64, rot_x: f64, rot_y: f64) -> (f64, f64, f64) {
    let (sy, cy) = (rot_y.sin(), rot_y.cos());
    let x1 =  x * cy + z * sy;
    let z1 = -x * sy + z * cy;
    let (sx, cx) = (rot_x.sin(), rot_x.cos());
    let y2 = y * cx - z1 * sx;
    let z2 = y * sx + z1 * cx;
    (x1, y2, z2)
}

/// AlphaFold-style pLDDT colour scheme.
pub fn plddt_color(plddt: f32) -> Color {
    if plddt >= 90.0      { Color::Rgb(  0,  53, 214) }
    else if plddt >= 70.0 { Color::Rgb(101, 203, 243) }
    else if plddt >= 50.0 { Color::Rgb(255, 219,  19) }
    else                  { Color::Rgb(255, 125,  69) }
}

/// Pre-project all atoms into colour-batched screen points for the Canvas.
/// `y_scale` corrects for terminal cell aspect ratio (cell_h / cell_w ≈ 2.1).
pub fn build_point_batches(
    atoms: &[ProteinAtom],
    rot_x: f64, rot_y: f64,
    y_scale: f64,
) -> Vec<(Color, Vec<(f64, f64)>)> {
    let mut map: HashMap<(u8, u8, u8), Vec<(f64, f64)>> = HashMap::new();

    for pair in atoms.windows(2) {
        let a = &pair[0];
        let b = &pair[1];
        // Skip chain breaks (ESMFold may produce multi-chain outputs)
        if b.residue.saturating_sub(a.residue) > 1 { continue; }

        let steps = 10usize;
        for i in 0..=steps {
            let t  = i as f64 / steps as f64;
            let x  = a.x + t * (b.x - a.x);
            let y  = a.y + t * (b.y - a.y);
            let z  = a.z + t * (b.z - a.z);
            let pl = a.plddt + t as f32 * (b.plddt - a.plddt);
            let ss = if t < 0.5 { a.secondary } else { b.secondary };

            let (sx, sy, _) = project(x, y, z, rot_x, rot_y);
            let sy = sy * y_scale;

            let c = plddt_color(pl);
            let (r, g, bl) = match c { Color::Rgb(r,g,b) => (r,g,b), _ => (128,128,128) };
            let pts = map.entry((r, g, bl)).or_default();

            pts.push((sx, sy));

            // Extra lateral dots: helices → wide blob; sheets → flat bar
            match ss {
                Secondary::Helix => {
                    let d = 0.03;
                    pts.push((sx + d, sy));
                    pts.push((sx - d, sy));
                    pts.push((sx,     sy + d * 0.5));
                    pts.push((sx,     sy - d * 0.5));
                }
                Secondary::Sheet => {
                    let d = 0.018;
                    pts.push((sx + d, sy));
                    pts.push((sx - d, sy));
                }
                Secondary::Coil => {}
            }
        }
    }

    map.into_iter().map(|((r,g,b), pts)| (Color::Rgb(r,g,b), pts)).collect()
}

