use feos_core::parameter::{ChemicalRecord, Identifier};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

mod polynomial;
mod supermolecule;
pub use supermolecule::SuperMolecule;

pub trait MolecularRepresentation {
    fn build(&self, y: Vec<f64>) -> ChemicalRecord;

    fn variables(&self) -> usize;
}

#[derive(Clone, Serialize, Deserialize)]
pub struct FixedMolecule;

impl MolecularRepresentation for FixedMolecule {
    fn variables(&self) -> usize {
        0
    }

    fn build(&self, _: Vec<f64>) -> ChemicalRecord {
        let identifier = Identifier::new("", Some("FixedMolecule"), None, None, None, None);
        ChemicalRecord::new_count(identifier, HashMap::new(), None)
    }
}
