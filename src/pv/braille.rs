use ratatui::layout::Rect;
use ratatui::widgets::Paragraph;

use crate::protein::{plddt_color, ProteinAtom};
use crate::pv::camera::Camera;
use crate::pv::framebuffer::{framebuffer_to_braille_widget, Framebuffer};
use ratatui::style::Color;

/// Render protein backbone as a braille `Paragraph` that fits `area`.
///
/// Each terminal cell maps to a 2×4 pixel block in the framebuffer, giving
/// 4× the spatial resolution of half-block rendering.  Backbone Cα–Cα bonds
/// are drawn as thick lines coloured by pLDDT confidence.
pub fn render_braille<'a>(
    atoms: &[ProteinAtom],
    camera: &Camera,
    area: Rect,
) -> Paragraph<'static> {
    if atoms.is_empty() || area.width < 2 || area.height < 2 {
        return Paragraph::new("");
    }

    let w = area.width as usize;
    let h = area.height as usize;

    // Braille resolution: 2 px wide × 4 px tall per terminal cell.
    let fb_w = w * 2;
    let fb_h = h * 4;
    let mut fb = Framebuffer::new(fb_w, fb_h);

    let cache = camera.projection_cache();

    // Scale projected coords (range ~[-1, 1] at zoom=1) to pixel space.
    // Protein is normalised so max radius ≈ 0.85; zoom=1 fills the panel.
    let scale_x = (fb_w as f64) * 0.45;
    let scale_y = (fb_h as f64) * 0.45;
    let cx = fb_w as f64 / 2.0;
    let cy = fb_h as f64 / 2.0;

    for pair in atoms.windows(2) {
        let a = &pair[0];
        let b = &pair[1];
        // Skip chain breaks
        if b.residue.saturating_sub(a.residue) > 1 {
            continue;
        }

        let pa = cache.project(a.x, a.y, a.z);
        let pb = cache.project(b.x, b.y, b.z);

        let ax = cx + pa.x * scale_x;
        let ay = cy - pa.y * scale_y; // screen Y is inverted
        let bx = cx + pb.x * scale_x;
        let by_ = cy - pb.y * scale_y;

        let color = color_rgb(plddt_color(a.plddt));

        // Thick line (4 px) for backbone visibility at braille resolution.
        fb.draw_thick_line_3d(
            [ax, ay, pa.z as f64],
            [bx, by_, pb.z as f64],
            color,
            4.0,
        );
    }

    framebuffer_to_braille_widget(&fb)
}

#[inline]
fn color_rgb(c: Color) -> [u8; 3] {
    match c {
        Color::Rgb(r, g, b) => [r, g, b],
        _ => [128, 128, 128],
    }
}
