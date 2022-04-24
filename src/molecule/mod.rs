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
pub struct FixedMolecule(pub String);

impl MolecularRepresentation for FixedMolecule {
    fn variables(&self) -> usize {
        0
    }

    fn build(&self, _: Vec<f64>) -> ChemicalRecord {
        let segments = HashMap::from([(self.0.clone(), 1.0)]);
        let identifier = Identifier::new("", Some(&self.0), None, None, None, None);
        ChemicalRecord::new_count(identifier, segments, None)
    }
}
