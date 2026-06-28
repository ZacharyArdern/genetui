use image::DynamicImage;
use crate::pv;
use crate::pv::camera::Camera;
use crate::pv::color::ColorScheme;
use crate::pv::ribbon::generate_ribbon_mesh;
use crate::pv::hd::render_hd_framebuffer;

/// Render a protein structure to a DynamicImage.
///
/// Atoms are centered and normalized to radius=1.0 in `from_protein_atoms`.
/// `camera.zoom` is treated as a user multiplier on top of the auto-fit zoom,
/// so zoom=1.0 (default) always fills ~90% of the panel.
pub fn render_protein_image(
    atoms: &[crate::protein::ProteinAtom],
    camera: &Camera,
    cols: u16,
    rows: u16,
) -> DynamicImage {
    let pw = (cols as u32 * 9).clamp(400, 1200) as f64;
    let ph = (rows as u32 * 17).clamp(200, 600) as f64;

    let protein = pv::from_protein_atoms(atoms);
    let color_scheme = ColorScheme;
    let mesh = generate_ribbon_mesh(&protein, &color_scheme);

    // Auto-fit zoom: match ProteinView's formula so the protein fills 90% of
    // the smaller pixel dimension at camera.zoom == 1.0 (user default).
    let radius = protein.bounding_radius().max(1.0);
    let auto_zoom = 0.9 * pw.min(ph) / (2.0 * radius);
    let mut fitted_camera = camera.clone();
    fitted_camera.zoom = camera.zoom * auto_zoom;

    let fb = render_hd_framebuffer(&protein, &fitted_camera, &color_scheme, pw, ph, &mesh);
    DynamicImage::ImageRgba8(fb.to_rgba_image())
}
