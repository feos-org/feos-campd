use serde::{Deserialize, Serialize};

mod comt_camd;
mod polynomial;
mod supermolecule;
pub use comt_camd::{CoMTCAMD, GroupCount};
pub use supermolecule::{SegmentAndBondCount, SuperMolecule};

/// A generic molecular representation to be used in an [OptimizationProblem](super::OptimizationProblem).
pub trait MolecularRepresentation {
    type ChemicalRecord;
    fn build(&self, y: Vec<f64>) -> Self::ChemicalRecord;

    fn variables(&self) -> usize;
}

/// Molecular representation for a fixed molecule.
#[derive(Clone, Serialize, Deserialize)]
pub struct FixedMolecule;

impl MolecularRepresentation for FixedMolecule {
    type ChemicalRecord = ();

    fn variables(&self) -> usize {
        0
    }

    fn build(&self, _: Vec<f64>) {}
}
