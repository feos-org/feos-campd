use feos_core::parameter::{ChemicalRecord, Identifier};
use serde::{Deserialize, Serialize};

mod supermolecule;
pub use supermolecule::SuperMolecule;

#[derive(Clone, Serialize, Deserialize)]
pub enum MolecularRepresentation {
    SuperMolecule(SuperMolecule),
    FixedMolecule(FixedMolecule),
}

#[derive(Clone, Serialize, Deserialize)]
pub struct FixedMolecule {
    components: Vec<String>,
}

impl FixedMolecule {
    pub fn new(components: Vec<String>) -> Self {
        Self { components }
    }

    pub fn variables(&self) -> usize {
        self.components.len()
    }

    pub fn build(&self, y: Vec<f64>) -> ChemicalRecord {
        let segments = self.components.iter().cloned().zip(y.into_iter()).collect();
        let identifier = Identifier::new("", Some("HomoGc"), None, None, None, None);
        ChemicalRecord::new_count(identifier, segments, None)
    }
}
