use feos_core::parameter::{ChemicalRecord, Identifier};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

mod supermolecule;
pub use supermolecule::{Molecule, SuperMolecule};

#[derive(Clone, Serialize, Deserialize)]
pub enum MolecularRepresentation {
    SuperMolecule(SuperMolecule),
    HomoGc(HomoGc),
    FixedMolecule(FixedMolecule),
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct FunctionalGroup {
    group: String,
    smiles: String,
}

impl FunctionalGroup {
    pub fn new(group: String, smiles: String) -> Self {
        Self { smiles, group }
    }
}

impl From<(String, String)> for FunctionalGroup {
    fn from(tuple: (String, String)) -> Self {
        let (group, smiles) = tuple;
        Self { group, smiles }
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct HomoGc {
    pub size: usize,
    functional_groups: Vec<FunctionalGroup>,
}

impl HomoGc {
    pub fn new(size: usize, functional_groups: Vec<FunctionalGroup>) -> Self {
        Self {
            size,
            functional_groups,
        }
    }

    pub fn variables(&self) -> usize {
        4 + self.functional_groups.len()
    }

    pub fn build(&self, y: Vec<f64>) -> ChemicalRecord {
        let mut segments = HashMap::new();
        let mut y_iter = y.into_iter();

        [">C<", ">CH", "CH2", "CH3"].iter().for_each(|&c| {
            segments.insert(c.to_string(), y_iter.next().unwrap());
        });
        self.functional_groups.iter().for_each(|f| {
            *segments.entry(f.group.clone()).or_insert(0.0) += y_iter.next().unwrap();
        });

        let identifier = Identifier::new("", Some("HomoGc"), None, None, None, None);
        ChemicalRecord::new_count(identifier, segments, None)
    }
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
