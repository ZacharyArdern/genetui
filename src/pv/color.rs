use crate::pv::{Chain, Residue, SecondaryStructure};
use ratatui::style::Color;

#[derive(Debug, Clone, Copy)]
pub struct ColorScheme;

impl ColorScheme {
    pub fn residue_color(&self, residue: &Residue, _chain: &Chain) -> Color {
        match residue.secondary_structure {
            SecondaryStructure::Helix => Color::Rgb(100, 180, 255),
            SecondaryStructure::Sheet => Color::Rgb(255, 200, 0),
            SecondaryStructure::Turn | SecondaryStructure::Coil => Color::Rgb(0, 204, 0),
        }
    }
}

pub fn color_to_rgb(color: Color) -> [u8; 3] {
    match color {
        Color::Rgb(r, g, b) => [r, g, b],
        _ => [180, 180, 180],
    }
}
