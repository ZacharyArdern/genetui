// The modules in this directory (camera, framebuffer, ribbon, braille, etc.) are
// derived from ProteinView by Tristan Farmer: https://github.com/001TMF/ProteinView
pub mod braille;
pub mod camera;
pub mod color;
pub mod framebuffer;
pub mod ribbon;
pub mod hd;
pub mod kitty_png;

pub(crate) const LARGE_STRUCTURE_THRESHOLD: usize = 5000;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SecondaryStructure { Helix, Sheet, Turn, Coil }

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MoleculeType { Protein }

#[derive(Debug, Clone)]
pub struct Atom {
    pub name: String,
    pub element: String,
    pub x: f64, pub y: f64, pub z: f64,
    pub b_factor: f64,
    pub is_backbone: bool,
    pub is_hetero: bool,
}

#[derive(Debug, Clone)]
pub struct Residue {
    pub name: String,
    pub seq_num: i32,
    pub atoms: Vec<Atom>,
    pub secondary_structure: SecondaryStructure,
}

#[derive(Debug, Clone)]
pub struct Chain {
    pub id: String,
    pub residues: Vec<Residue>,
    pub molecule_type: MoleculeType,
}

#[derive(Debug, Clone)]
pub struct Protein {
    pub chains: Vec<Chain>,
}

impl Protein {
    pub fn residue_count(&self) -> usize {
        self.chains.iter().map(|c| c.residues.len()).sum()
    }

    pub fn bounding_radius(&self) -> f64 {
        self.chains.iter()
            .flat_map(|c| &c.residues)
            .flat_map(|r| &r.atoms)
            .map(|a| (a.x * a.x + a.y * a.y + a.z * a.z).sqrt())
            .fold(0.0_f64, f64::max)
    }
}

/// Convert our ProteinAtom slice into a Protein model for rendering.
/// Centers the protein at the origin and normalizes to bounding_radius = 1.0
/// so that auto-zoom in render_protein_image always fits the view.
pub fn from_protein_atoms(atoms: &[crate::protein::ProteinAtom]) -> Protein {
    if atoms.is_empty() {
        return Protein { chains: vec![] };
    }

    // parse_pdb already centered + normalized coords to radius ≈ 0.85.
    // Ribbon constants (HELIX_HALF_WIDTH=1.30, COIL_RADIUS=0.40) are in Ångstroms
    // and designed for raw PDB scale (~15-30 Å radius).
    // Rescale to TARGET_RADIUS so cross-sections are ~5-10% of protein radius.
    const TARGET_RADIUS: f64 = 15.0;

    let n = atoms.len() as f64;
    let cx = atoms.iter().map(|a| a.x).sum::<f64>() / n;
    let cy = atoms.iter().map(|a| a.y).sum::<f64>() / n;
    let cz = atoms.iter().map(|a| a.z).sum::<f64>() / n;

    let current_r = atoms.iter()
        .map(|a| ((a.x-cx).powi(2) + (a.y-cy).powi(2) + (a.z-cz).powi(2)).sqrt())
        .fold(0.0_f64, f64::max)
        .max(1e-6);

    let scale = TARGET_RADIUS / current_r;

    // Group by residue index
    let mut residues: Vec<Residue> = Vec::new();
    let mut current_res_idx = u32::MAX;

    for atom in atoms {
        if atom.residue != current_res_idx {
            current_res_idx = atom.residue;
            let ss = match atom.secondary {
                crate::protein::Secondary::Helix => SecondaryStructure::Helix,
                crate::protein::Secondary::Sheet => SecondaryStructure::Sheet,
                crate::protein::Secondary::Coil  => SecondaryStructure::Coil,
            };
            residues.push(Residue {
                name: "ALA".to_string(),
                seq_num: residues.len() as i32 + 1,
                atoms: Vec::new(),
                secondary_structure: ss,
            });
        }
        if let Some(res) = residues.last_mut() {
            let is_first = res.atoms.is_empty();
            res.atoms.push(Atom {
                name: if is_first { "CA".to_string() } else { "CB".to_string() },
                element: "C".to_string(),
                x: (atom.x - cx) * scale,
                y: (atom.y - cy) * scale,
                z: (atom.z - cz) * scale,
                b_factor: atom.plddt as f64,
                is_backbone: is_first,
                is_hetero: false,
            });
        }
    }

    Protein {
        chains: vec![Chain {
            id: "A".to_string(),
            residues,
            molecule_type: MoleculeType::Protein,
        }],
    }
}
