//! Ribbon/cartoon geometry generation for HD protein rendering.
//!
//! Generates triangle meshes from protein backbone data by:
//! 1. Extracting C-alpha positions per chain
//! 2. Fitting a Catmull-Rom spline through backbone positions
//! 3. Computing Frenet-Serret local coordinate frames along the spline
//! 4. Extruding secondary-structure-dependent cross-sections
//! 5. Connecting consecutive cross-sections into triangle strips
//!
//! The output mesh is in world space; the caller projects through the camera.

use crate::pv::{MoleculeType, Protein, SecondaryStructure};
use crate::pv::color::{ColorScheme, color_to_rgb};

// ---------------------------------------------------------------------------
// Constants & LOD configuration
// ---------------------------------------------------------------------------

/// Default number of spline subdivisions between each pair of C-alpha atoms.
const DEFAULT_SPLINE_SUBDIVISIONS: usize = 14;

/// Default number of vertices around the coil/turn tube cross-section.
const DEFAULT_COIL_SEGMENTS: usize = 12;

/// Reduced spline subdivisions for large structures (>5000 residues).
const LARGE_SPLINE_SUBDIVISIONS: usize = 4;

/// Reduced coil segments for large structures (>5000 residues).
const LARGE_COIL_SEGMENTS: usize = 6;

/// Level-of-detail configuration for ribbon mesh generation.
/// Large structures use reduced subdivision counts to cut triangle count
/// with zero visible difference at terminal resolution.
#[derive(Debug, Clone, Copy)]
struct LodConfig {
    spline_subdivisions: usize,
    coil_segments: usize,
}

impl LodConfig {
    fn normal() -> Self {
        Self {
            spline_subdivisions: DEFAULT_SPLINE_SUBDIVISIONS,
            coil_segments: DEFAULT_COIL_SEGMENTS,
        }
    }

    fn large() -> Self {
        Self {
            spline_subdivisions: LARGE_SPLINE_SUBDIVISIONS,
            coil_segments: LARGE_COIL_SEGMENTS,
        }
    }

    /// Pick the appropriate LOD based on residue count.
    fn for_residue_count(residue_count: usize) -> Self {
        if residue_count > crate::pv::LARGE_STRUCTURE_THRESHOLD {
            Self::large()
        } else {
            Self::normal()
        }
    }
}

/// Cross-section dimensions (in Angstroms).
const HELIX_HALF_WIDTH: f64 = 1.17;
const HELIX_HALF_HEIGHT: f64 = 0.40;

const SHEET_HALF_WIDTH: f64 = 1.50;
const SHEET_HALF_HEIGHT: f64 = 0.20;

const SHEET_ARROW_HALF_WIDTH: f64 = 2.20;

const COIL_RADIUS: f64 = 0.40;

// ---------------------------------------------------------------------------
// Output type
// ---------------------------------------------------------------------------

/// A single triangle in the ribbon mesh, ready for rasterization.
#[derive(Debug, Clone)]
pub struct RibbonTriangle {
    /// Three vertices in 3D world space, each `[x, y, z]`.
    pub verts: [[f64; 3]; 3],
    /// Base RGB color of this face.
    pub color: [u8; 3],
    /// Outward-facing unit normal of the triangle.
    pub normal: [f64; 3],
}

// ---------------------------------------------------------------------------
// 3-component vector helpers (no external crate)
// ---------------------------------------------------------------------------

type V3 = [f64; 3];

#[inline]
fn v3_add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn v3_sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn v3_scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn v3_dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn v3_cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn v3_len(a: V3) -> f64 {
    v3_dot(a, a).sqrt()
}

#[inline]
fn v3_normalize(a: V3) -> V3 {
    let l = v3_len(a);
    if l < 1e-12 {
        [0.0, 1.0, 0.0] // fallback up vector
    } else {
        v3_scale(a, 1.0 / l)
    }
}

// ---------------------------------------------------------------------------
// Sheet frame guide helpers
// ---------------------------------------------------------------------------

use crate::pv::Residue;

/// Compute the carbonyl C→O direction for a protein residue.
///
/// Looks for backbone atoms named "C" (carbonyl carbon) and "O" (carbonyl
/// oxygen), computes the normalized direction vector from C to O.  As a
/// fallback, tries CA→O if C is missing.  Returns `None` for nucleic acid
/// residues or when the required atoms are absent.
fn residue_frame_hint(residue: &Residue) -> Option<V3> {
    let find_atom = |name: &str| -> Option<V3> {
        residue
            .atoms
            .iter()
            .find(|a| a.name.trim() == name)
            .map(|a| [a.x, a.y, a.z])
    };

    let o_pos = find_atom("O")?;

    // Prefer backbone C atom; fall back to CA.
    let origin = find_atom("C").or_else(|| find_atom("CA"))?;

    let dir = v3_sub(o_pos, origin);
    let len = v3_len(dir);
    if len < 1e-12 {
        return None;
    }
    Some(v3_scale(dir, 1.0 / len))
}

/// Project vector `v` perpendicular to `onto` and normalize the result.
/// Returns `None` if the projection is degenerate (v is parallel to onto).
fn project_perp_normalize(v: V3, onto: V3) -> Option<V3> {
    let d = v3_dot(v, onto);
    let perp = v3_sub(v, v3_scale(onto, d));
    let len = v3_len(perp);
    if len < 1e-12 {
        None
    } else {
        Some(v3_scale(perp, 1.0 / len))
    }
}

/// Adjust binormal/normal frames for sheet regions using carbonyl direction
/// hints (frame guides).
fn apply_sheet_frame_guides(spline_points: &mut [SplinePoint]) {
    let n = spline_points.len();
    if n == 0 {
        return;
    }

    // Identify contiguous sheet runs.
    let mut i = 0;
    while i < n {
        if spline_points[i].ss != SecondaryStructure::Sheet {
            i += 1;
            continue;
        }

        // Found the start of a sheet run.
        let run_start = i;
        while i < n && spline_points[i].ss == SecondaryStructure::Sheet {
            i += 1;
        }
        let run_end = i; // exclusive

        // Process this sheet run: project hints, align signs, blend.
        let mut prev_hint_perp: Option<V3> = None;

        for sp in spline_points[run_start..run_end].iter_mut() {
            let hint = match sp.frame_hint {
                Some(h) => h,
                None => continue,
            };

            let tangent = sp.tangent;

            // Project hint perpendicular to tangent.
            let mut hint_perp = match project_perp_normalize(hint, tangent) {
                Some(h) => h,
                None => continue,
            };

            // Align sign with previous hint to avoid flipping.
            if let Some(prev) = prev_hint_perp {
                if v3_dot(hint_perp, prev) < 0.0 {
                    hint_perp = v3_scale(hint_perp, -1.0);
                }
            }
            prev_hint_perp = Some(hint_perp);

            // Blend with existing binormal: 65% hint, 35% existing.
            let existing_b = sp.binormal;
            let blended = v3_add(v3_scale(hint_perp, 0.65), v3_scale(existing_b, 0.35));
            let new_binormal = v3_normalize(blended);

            // Recompute normal from tangent x binormal.
            let new_normal = v3_normalize(v3_cross(tangent, new_binormal));

            sp.binormal = new_binormal;
            sp.normal = new_normal;
        }
    }
}

/// Flatten sheet ribbon frames by averaging binormals across each sheet run.
///
/// Without full backbone atoms the carbonyl hints are absent, so Frenet frames
/// can twist freely and make flat β-sheets look bumpy.  This function collects
/// the Frenet binormals for each contiguous sheet run, sign-corrects them so
/// they all point the same way, averages them into a consensus direction, then
/// re-projects that consensus onto the local tangent plane so the ribbon lies
/// flat along the whole strand.
fn smooth_sheet_frames(spline_points: &mut [SplinePoint]) {
    let n = spline_points.len();
    let mut i = 0;
    while i < n {
        if spline_points[i].ss != SecondaryStructure::Sheet {
            i += 1;
            continue;
        }
        let start = i;
        while i < n && spline_points[i].ss == SecondaryStructure::Sheet { i += 1; }
        let end = i;

        // Walk forward, flipping each binormal to match its predecessor.
        let mut corrected: Vec<V3> = Vec::with_capacity(end - start);
        let mut prev = spline_points[start].binormal;
        corrected.push(prev);
        for sp in &spline_points[start + 1..end] {
            let b = if v3_dot(sp.binormal, prev) < 0.0 {
                v3_scale(sp.binormal, -1.0)
            } else {
                sp.binormal
            };
            corrected.push(b);
            prev = b;
        }

        // Average into a consensus binormal direction.
        let mut sum = [0.0_f64; 3];
        for b in &corrected { sum = v3_add(sum, *b); }
        let avg_b = v3_normalize(sum);

        // Re-project consensus onto each tangent plane and apply.
        for sp in spline_points[start..end].iter_mut() {
            let t = sp.tangent;
            let b_perp = v3_sub(avg_b, v3_scale(t, v3_dot(avg_b, t)));
            let bl = v3_len(b_perp);
            if bl < 1e-12 { continue; }
            let new_b = v3_scale(b_perp, 1.0 / bl);
            sp.binormal = new_b;
            sp.normal   = v3_normalize(v3_cross(t, new_b));
        }
    }
}

// ---------------------------------------------------------------------------
// Catmull-Rom spline
// ---------------------------------------------------------------------------

/// Evaluate VMD's modified Catmull-Rom spline (slope 1.25 vs standard 2.0).
/// Slope 1.25 gives fuller, rounder helix loops — matches molar_vis / VMD NewCartoon.
fn catmull_rom(p0: V3, p1: V3, p2: V3, p3: V3, t: f64) -> V3 {
    // is = 1/slope = 0.8; standard CR uses is = 0.5 (slope=2).
    const IS: f64 = 1.0 / 1.25;
    let mut out = [0.0; 3];
    for k in 0..3 {
        let q0 = -IS*p0[k] + (2.0-IS)*p1[k] + (IS-2.0)*p2[k] + IS*p3[k];
        let q1 = 2.0*IS*p0[k] + (IS-3.0)*p1[k] + (3.0-2.0*IS)*p2[k] - IS*p3[k];
        let q2 = -IS*p0[k] + IS*p2[k];
        out[k] = ((q0*t + q1)*t + q2)*t + p1[k];
    }
    out
}

// ---------------------------------------------------------------------------
// Spline point with metadata
// ---------------------------------------------------------------------------

/// A single point along the backbone spline, carrying the local coordinate
/// frame and secondary-structure annotation.
struct SplinePoint {
    pos: V3,
    tangent: V3,
    normal: V3,
    binormal: V3,
    /// Carbonyl direction hint propagated from the nearest control point.
    /// Used by `apply_sheet_frame_guides` to orient sheet ribbons.
    frame_hint: Option<V3>,
    ss: SecondaryStructure,
    color: [u8; 3],
    /// True when this point lies within the arrowhead region at the end of a
    /// sheet run.  The arrowhead linearly widens over the last two original
    /// residue spans (2 * spline_subdivisions points).
    arrow_t: Option<f64>, // 0.0 = start of arrow, 1.0 = tip
}

// ---------------------------------------------------------------------------
// Cross-section generation
// ---------------------------------------------------------------------------

/// Build the cross-section ring for a given spline point.  Returns the
/// world-space positions of the cross-section vertices.
fn cross_section(sp: &SplinePoint, lod: &LodConfig) -> Vec<V3> {
    match sp.ss {
        SecondaryStructure::Helix => ribbon_cross_section(sp, HELIX_HALF_WIDTH, HELIX_HALF_HEIGHT),
        SecondaryStructure::Sheet => {
            let hw = if let Some(t) = sp.arrow_t {
                // arrow_t=1 at arrow base (full barb width), 0 at strand tip (point).
                SHEET_ARROW_HALF_WIDTH * t
            } else {
                SHEET_HALF_WIDTH
            };
            ribbon_cross_section(sp, hw.max(0.01), SHEET_HALF_HEIGHT)
        }
        SecondaryStructure::Turn | SecondaryStructure::Coil => coil_cross_section(sp, lod),
    }
}

/// Flat ribbon cross-section with 4 vertices (top-left, top-right,
/// bottom-right, bottom-left) forming a thin rectangular profile.
fn ribbon_cross_section(sp: &SplinePoint, half_w: f64, half_h: f64) -> Vec<V3> {
    let n = sp.normal;
    let b = sp.binormal;

    // Four corners of the rectangular cross-section.
    let tl = v3_add(sp.pos, v3_add(v3_scale(b, -half_w), v3_scale(n, half_h)));
    let tr = v3_add(sp.pos, v3_add(v3_scale(b, half_w), v3_scale(n, half_h)));
    let br = v3_add(sp.pos, v3_add(v3_scale(b, half_w), v3_scale(n, -half_h)));
    let bl = v3_add(sp.pos, v3_add(v3_scale(b, -half_w), v3_scale(n, -half_h)));

    vec![tl, tr, br, bl]
}

/// Circular tube cross-section for coil/turn regions.
fn coil_cross_section(sp: &SplinePoint, lod: &LodConfig) -> Vec<V3> {
    let n = sp.normal;
    let b = sp.binormal;
    let segs = lod.coil_segments;
    let mut pts = Vec::with_capacity(segs);
    for i in 0..segs {
        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (segs as f64);
        let (sin_a, cos_a) = angle.sin_cos();
        let offset = v3_add(
            v3_scale(n, cos_a * COIL_RADIUS),
            v3_scale(b, sin_a * COIL_RADIUS),
        );
        pts.push(v3_add(sp.pos, offset));
    }
    pts
}

// ---------------------------------------------------------------------------
// Triangle normal
// ---------------------------------------------------------------------------

fn triangle_normal(v0: V3, v1: V3, v2: V3) -> V3 {
    let e1 = v3_sub(v1, v0);
    let e2 = v3_sub(v2, v0);
    v3_normalize(v3_cross(e1, e2))
}

// ---------------------------------------------------------------------------
// Emit triangle strip between two cross-sections
// ---------------------------------------------------------------------------

/// Connect two consecutive cross-section rings with a triangle strip.
/// Both rings must have the same number of vertices.
fn emit_strip(ring_a: &[V3], ring_b: &[V3], color: [u8; 3], out: &mut Vec<RibbonTriangle>) {
    let n = ring_a.len();
    debug_assert_eq!(n, ring_b.len());
    if n == 0 {
        return;
    }

    for i in 0..n {
        let j = (i + 1) % n;

        let a0 = ring_a[i];
        let a1 = ring_a[j];
        let b0 = ring_b[i];
        let b1 = ring_b[j];

        // Quad (a0, a1, b1, b0) -> two triangles.
        // Triangle 1: a0, a1, b0
        let n1 = triangle_normal(a0, a1, b0);
        out.push(RibbonTriangle {
            verts: [a0, a1, b0],
            color,
            normal: n1,
        });

        // Triangle 2: a1, b1, b0
        let n2 = triangle_normal(a1, b1, b0);
        out.push(RibbonTriangle {
            verts: [a1, b1, b0],
            color,
            normal: n2,
        });
    }
}

/// Emit a cap (disc) to close the end of a tube or ribbon.
fn emit_cap(
    ring: &[V3],
    center: V3,
    color: [u8; 3],
    facing_forward: bool,
    out: &mut Vec<RibbonTriangle>,
) {
    let n = ring.len();
    if n < 3 {
        return;
    }
    for i in 0..n {
        let j = (i + 1) % n;
        let (v0, v1) = if facing_forward {
            (ring[i], ring[j])
        } else {
            (ring[j], ring[i])
        };
        let norm = triangle_normal(center, v0, v1);
        out.push(RibbonTriangle {
            verts: [center, v0, v1],
            color,
            normal: norm,
        });
    }
}

// ---------------------------------------------------------------------------
// Transition cross-sections between different secondary structure types
// ---------------------------------------------------------------------------

fn resample_ring(ring: &[V3], target_count: usize) -> Vec<V3> {
    let n = ring.len();
    if n == target_count {
        return ring.to_vec();
    }
    if n == 0 || target_count == 0 {
        return Vec::new();
    }

    let mut out = Vec::with_capacity(target_count);
    for i in 0..target_count {
        let frac = (i as f64) / (target_count as f64) * (n as f64);
        let idx = frac as usize;
        let t = frac - idx as f64;
        let a = ring[idx % n];
        let b = ring[(idx + 1) % n];
        out.push(v3_add(v3_scale(a, 1.0 - t), v3_scale(b, t)));
    }
    out
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Generate the complete ribbon/cartoon triangle mesh for a protein.
///
/// The returned triangles are in world space.  The caller should project each
/// vertex through the camera and then rasterize.
pub fn generate_ribbon_mesh(protein: &Protein, color_scheme: &ColorScheme) -> Vec<RibbonTriangle> {
    let mut triangles: Vec<RibbonTriangle> = Vec::new();
    let lod = LodConfig::for_residue_count(protein.residue_count());

    for chain in &protein.chains {
        match chain.molecule_type {
            MoleculeType::Protein => {
                generate_chain_ribbon(chain, color_scheme, &lod, &mut triangles);
            }
        }
    }

    triangles
}

// ---------------------------------------------------------------------------
// Per-chain generation
// ---------------------------------------------------------------------------

/// C-alpha record extracted from a residue.
struct CaRecord {
    pos: V3,
    ss: SecondaryStructure,
    color: [u8; 3],
    /// Carbonyl C→O direction hint for sheet frame orientation.
    frame_hint: Option<V3>,
}

/// Demote isolated single-residue helix/sheet runs to coil so they don't
/// generate spurious stubs or arrows (matches molar_vis / VMD behaviour).
fn demote_singletons(ss: &mut [SecondaryStructure]) {
    let n = ss.len();
    let mut i = 0;
    while i < n {
        let c = ss[i];
        let mut j = i + 1;
        while j < n && ss[j] == c { j += 1; }
        if c != SecondaryStructure::Coil && c != SecondaryStructure::Turn && j - i < 2 {
            for k in i..j { ss[k] = SecondaryStructure::Coil; }
        }
        i = j;
    }
}

fn generate_chain_ribbon(
    chain: &crate::pv::Chain,
    color_scheme: &ColorScheme,
    lod: &LodConfig,
    out: &mut Vec<RibbonTriangle>,
) {
    // 1. Collect C-alpha positions, SS types, colors, and frame hints.
    let cas: Vec<CaRecord> = chain
        .residues
        .iter()
        .filter_map(|res| {
            let ca = res.atoms.iter().find(|a| a.is_backbone)?;
            let color = color_to_rgb(color_scheme.residue_color(res, chain));
            Some(CaRecord {
                pos: [ca.x, ca.y, ca.z],
                ss: res.secondary_structure,
                color,
                frame_hint: residue_frame_hint(res),
            })
        })
        .collect();

    let n = cas.len();
    if n < 2 {
        return;
    }

    // 1b. Demote isolated single-residue SS runs to coil.
    let mut ss_classes: Vec<SecondaryStructure> = cas.iter().map(|r| r.ss).collect();
    demote_singletons(&mut ss_classes);

    // 1c. Pre-smooth β-strand Cα positions with weighted neighbour average
    // `(2·CAᵢ + CAᵢ₋₁ + CAᵢ₊₁) / 4` (VMD / molar_vis). This cancels the
    // alternating zig-zag pleat of a β-strand before splining, which is the
    // primary cause of bumpy-looking sheet ribbons.
    let smoothed: Vec<V3> = (0..n).map(|i| {
        if ss_classes[i] == SecondaryStructure::Sheet {
            let prev = cas[i.saturating_sub(1)].pos;
            let next = cas[(i + 1).min(n - 1)].pos;
            v3_add(v3_scale(cas[i].pos, 0.5), v3_scale(v3_add(prev, next), 0.25))
        } else {
            cas[i].pos
        }
    }).collect();

    // 2. Identify arrowhead regions (strand end → taper to a point).
    let mut arrow_regions: Vec<(f64, f64)> = Vec::new(); // (base, tip) in residue coords
    {
        let mut i = 0;
        while i < n {
            if ss_classes[i] == SecondaryStructure::Sheet {
                let run_start = i;
                while i < n && ss_classes[i] == SecondaryStructure::Sheet { i += 1; }
                let run_end = i; // exclusive
                // Arrow: last 1.6 residues of strand taper from arrow_base width → 0.
                let tip  = (run_end - 1) as f64;
                let base = (tip - 1.6_f64).max(run_start as f64);
                arrow_regions.push((base, tip));
            } else {
                i += 1;
            }
        }
    }

    // 3. Generate spline points with Catmull-Rom interpolation.
    let mut spline_points: Vec<SplinePoint> = Vec::new();

    for seg in 0..n - 1 {
        // Indices for the four control points, clamping at endpoints.
        let i0 = if seg == 0 { 0 } else { seg - 1 };
        let i1 = seg;
        let i2 = seg + 1;
        let i3 = if seg + 2 >= n { n - 1 } else { seg + 2 };

        let p0 = smoothed[i0];
        let p1 = smoothed[i1];
        let p2 = smoothed[i2];
        let p3 = smoothed[i3];

        let subdivs = if seg == n - 2 {
            lod.spline_subdivisions + 1 // include the last point
        } else {
            lod.spline_subdivisions
        };

        for sub in 0..subdivs {
            let t = sub as f64 / lod.spline_subdivisions as f64;
            let pos = catmull_rom(p0, p1, p2, p3, t);

            // Interpolated secondary structure: use the nearer residue (demoted).
            let ss = if t < 0.5 { ss_classes[i1] } else { ss_classes[i2] };
            // Interpolated color: use the nearer residue.
            let color = if t < 0.5 { cas[i1].color } else { cas[i2].color };

            // Arrow_t: fraction along arrowhead (0=base/wide, 1=tip/point).
            let global_pos = seg as f64 + t;
            let arrow_t = arrow_regions.iter().find_map(|&(base, tip)| {
                if global_pos >= base && global_pos <= tip + 1.0 {
                    // Inside arrow: taper from arrow_base → 0 as r → tip.
                    Some(((tip - global_pos) / (tip - base).max(1e-3)).clamp(0.0, 1.0))
                } else {
                    None
                }
            });

            // Nearest-neighbor frame hint: snap to the closer control point.
            let frame_hint = if t < 0.5 {
                cas[i1].frame_hint
            } else {
                cas[i2].frame_hint
            };

            spline_points.push(SplinePoint {
                pos,
                tangent: [0.0, 0.0, 0.0], // computed below
                normal: [0.0, 1.0, 0.0],
                binormal: [0.0, 0.0, 1.0],
                frame_hint,
                ss,
                color,
                arrow_t,
            });
        }
    }

    if spline_points.len() < 2 {
        return;
    }

    // 4. Compute tangents via finite differences.
    let sp_len = spline_points.len();
    for i in 0..sp_len {
        let prev = if i == 0 { 0 } else { i - 1 };
        let next = if i == sp_len - 1 { sp_len - 1 } else { i + 1 };
        let t = v3_normalize(v3_sub(spline_points[next].pos, spline_points[prev].pos));
        spline_points[i].tangent = t;
    }

    // 5. Compute Frenet-Serret frames with a propagated reference normal.
    {
        // Choose initial normal perpendicular to first tangent.
        let t0 = spline_points[0].tangent;
        let arbitrary = if t0[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let mut prev_normal = v3_normalize(v3_cross(t0, arbitrary));

        for sp in spline_points.iter_mut() {
            let t = sp.tangent;
            // Project previous normal onto plane perpendicular to current tangent.
            let proj = v3_scale(t, v3_dot(prev_normal, t));
            let mut n = v3_sub(prev_normal, proj);
            let nl = v3_len(n);
            if nl < 1e-12 {
                // Degenerate: pick a new arbitrary normal.
                let arb = if t[0].abs() < 0.9 {
                    [1.0, 0.0, 0.0]
                } else {
                    [0.0, 1.0, 0.0]
                };
                n = v3_normalize(v3_cross(t, arb));
            } else {
                n = v3_scale(n, 1.0 / nl);
            }
            let b = v3_normalize(v3_cross(t, n));

            sp.normal = n;
            sp.binormal = b;
            prev_normal = n;
        }
    }

    // 5b. Apply sheet frame guides (carbonyl-direction-based orientation).
    apply_sheet_frame_guides(&mut spline_points);

    // 5c. Smooth sheet frames: average binormal across each sheet run so the
    // ribbon lies flat even when carbonyl hints are absent (Cα-only models).
    smooth_sheet_frames(&mut spline_points);

    // 6. Build cross-sections and emit triangle strips.
    let mut prev_ring = cross_section(&spline_points[0], lod);

    for sp in spline_points.iter().skip(1) {
        let mut curr_ring = cross_section(sp, lod);
        let color = sp.color;

        // Handle cross-section vertex count mismatch at SS transitions.
        if prev_ring.len() != curr_ring.len() {
            let target = prev_ring.len().max(curr_ring.len());
            if prev_ring.len() != target {
                prev_ring = resample_ring(&prev_ring, target);
            }
            if curr_ring.len() != target {
                curr_ring = resample_ring(&curr_ring, target);
            }
        }

        emit_strip(&prev_ring, &curr_ring, color, out);
        prev_ring = curr_ring;
    }

    // 7. Cap the ends of the ribbon.
    let first_ring = cross_section(&spline_points[0], lod);
    emit_cap(
        &first_ring,
        spline_points[0].pos,
        spline_points[0].color,
        false,
        out,
    );

    let last_ring = cross_section(spline_points.last().unwrap(), lod);
    emit_cap(
        &last_ring,
        spline_points.last().unwrap().pos,
        spline_points.last().unwrap().color,
        true,
        out,
    );
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pv::{Atom, Residue};

    /// Helper to build a minimal atom.
    fn make_atom(name: &str, x: f64, y: f64, z: f64) -> Atom {
        Atom {
            name: name.to_string(),
            element: "C".to_string(),
            x,
            y,
            z,
            b_factor: 0.0,
            is_backbone: name == "CA",
            is_hetero: false,
        }
    }

    #[test]
    fn test_residue_frame_hint_with_co() {
        // Residue with C at origin and O along the +X axis.
        let res = Residue {
            name: "ALA".to_string(),
            seq_num: 1,
            atoms: vec![
                make_atom("CA", 0.0, 0.0, 0.0),
                make_atom("C", 1.0, 0.0, 0.0),
                make_atom("O", 4.0, 0.0, 0.0),
            ],
            secondary_structure: SecondaryStructure::Sheet,
        };
        let hint = residue_frame_hint(&res).expect("should produce a hint");
        // Direction should be [1, 0, 0] (normalized C→O).
        assert!((hint[0] - 1.0).abs() < 1e-9);
        assert!(hint[1].abs() < 1e-9);
        assert!(hint[2].abs() < 1e-9);
    }

    #[test]
    fn test_residue_frame_hint_missing_o() {
        // Residue with C but no O.
        let res = Residue {
            name: "ALA".to_string(),
            seq_num: 1,
            atoms: vec![
                make_atom("CA", 0.0, 0.0, 0.0),
                make_atom("C", 1.0, 0.0, 0.0),
                make_atom("N", 2.0, 0.0, 0.0),
            ],
            secondary_structure: SecondaryStructure::Sheet,
        };
        assert!(residue_frame_hint(&res).is_none());
    }

    #[test]
    fn test_residue_frame_hint_fallback_ca_o() {
        // Residue with CA and O but no C — should use CA→O fallback.
        let res = Residue {
            name: "ALA".to_string(),
            seq_num: 1,
            atoms: vec![
                make_atom("CA", 0.0, 0.0, 0.0),
                make_atom("O", 0.0, 3.0, 0.0),
            ],
            secondary_structure: SecondaryStructure::Sheet,
        };
        let hint = residue_frame_hint(&res).expect("should fallback to CA->O");
        assert!(hint[0].abs() < 1e-9);
        assert!((hint[1] - 1.0).abs() < 1e-9);
        assert!(hint[2].abs() < 1e-9);
    }

    #[test]
    fn test_apply_sheet_guides_flipped_hints() {
        let mut points: Vec<SplinePoint> = Vec::new();
        for i in 0..6 {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            points.push(SplinePoint {
                pos: [i as f64, 0.0, 0.0],
                tangent: [1.0, 0.0, 0.0],
                normal: [0.0, 1.0, 0.0],
                binormal: [0.0, 0.0, 1.0],
                frame_hint: Some([0.0, sign * 0.5, 0.5]),
                ss: SecondaryStructure::Sheet,
                color: [255, 255, 255],
                arrow_t: None,
            });
        }

        apply_sheet_frame_guides(&mut points);

        for i in 1..points.len() {
            let d = v3_dot(points[i].binormal, points[i - 1].binormal);
            assert!(
                d > 0.0,
                "Binormals at {} and {} should agree in sign, got dot={}",
                i - 1,
                i,
                d
            );
        }
    }

    #[test]
    fn test_coil_points_unaffected() {
        let original_binormal: V3 = [0.0, 0.0, 1.0];
        let original_normal: V3 = [0.0, 1.0, 0.0];

        let mut points: Vec<SplinePoint> = Vec::new();
        for i in 0..4 {
            points.push(SplinePoint {
                pos: [i as f64, 0.0, 0.0],
                tangent: [1.0, 0.0, 0.0],
                normal: original_normal,
                binormal: original_binormal,
                frame_hint: Some([0.0, 0.7, 0.7]),
                ss: SecondaryStructure::Coil,
                color: [128, 128, 128],
                arrow_t: None,
            });
        }

        apply_sheet_frame_guides(&mut points);

        for (i, sp) in points.iter().enumerate() {
            assert_eq!(
                sp.binormal, original_binormal,
                "Coil point {} binormal should be unchanged",
                i
            );
            assert_eq!(
                sp.normal, original_normal,
                "Coil point {} normal should be unchanged",
                i
            );
        }
    }

    #[test]
    fn test_helix_points_unaffected() {
        let original_binormal: V3 = [0.0, 0.0, 1.0];
        let original_normal: V3 = [0.0, 1.0, 0.0];

        let mut points: Vec<SplinePoint> = Vec::new();
        for i in 0..4 {
            points.push(SplinePoint {
                pos: [i as f64, 0.0, 0.0],
                tangent: [1.0, 0.0, 0.0],
                normal: original_normal,
                binormal: original_binormal,
                frame_hint: Some([0.0, 0.7, 0.7]),
                ss: SecondaryStructure::Helix,
                color: [128, 128, 128],
                arrow_t: None,
            });
        }

        apply_sheet_frame_guides(&mut points);

        for (i, sp) in points.iter().enumerate() {
            assert_eq!(
                sp.binormal, original_binormal,
                "Helix point {} binormal should be unchanged",
                i
            );
            assert_eq!(
                sp.normal, original_normal,
                "Helix point {} normal should be unchanged",
                i
            );
        }
    }
}
