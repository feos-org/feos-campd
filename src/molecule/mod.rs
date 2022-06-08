use serde::{Deserialize, Serialize};

mod polynomial;
mod supermolecule;
pub use supermolecule::{SegmentAndBondCount, SuperMolecule};

pub trait MolecularRepresentation {
    type ChemicalRecord;
    fn build(&self, y: Vec<f64>) -> Self::ChemicalRecord;

    fn variables(&self) -> usize;
}

#[derive(Clone, Serialize, Deserialize)]
pub struct FixedMolecule;

impl MolecularRepresentation for FixedMolecule {
    type ChemicalRecord = ();

    fn variables(&self) -> usize {
        0
    }

    fn build(&self, _: Vec<f64>) {}
}
