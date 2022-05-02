use feos_core::parameter::ChemicalRecord;
use serde::{Deserialize, Serialize};

mod polynomial;
mod supermolecule;
pub use supermolecule::SuperMolecule;

pub trait MolecularRepresentation {
    fn build(&self, y: Vec<f64>) -> ChemicalRecord;

    fn variables(&self) -> usize;
}

#[derive(Clone, Serialize, Deserialize)]
pub struct FixedMolecule(pub ChemicalRecord);

impl MolecularRepresentation for FixedMolecule {
    fn variables(&self) -> usize {
        0
    }

    fn build(&self, _: Vec<f64>) -> ChemicalRecord {
        self.0.clone()
    }
}
