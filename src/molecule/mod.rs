use serde::{Deserialize, Serialize};

mod comt_camd;
mod polynomial;
mod supermolecule;
mod supermolecule_disjunct;
pub use comt_camd::{CoMTCAMD, GroupCount};
pub use supermolecule::{SegmentAndBondCount, SuperMolecule};
pub use supermolecule_disjunct::SuperMoleculeDisjunct;

/// A generic molecular representation to be used in an [OptimizationProblem](super::OptimizationProblem).
pub trait MolecularRepresentation<const N: usize> {
    type ChemicalRecord;
    fn build(&self, y: &[f64]) -> [Self::ChemicalRecord; N];

    fn variables(&self) -> usize;

    fn smiles(&self, y: &[usize]) -> [String; N];
}

/// Molecular representation for a fixed molecule.
#[derive(Clone, Serialize, Deserialize)]
pub struct FixedMolecule<C> {
    smiles: String,
    chemical_record: C,
}

impl<C> FixedMolecule<C> {
    pub fn new(smiles: String, chemical_record: C) -> Self {
        Self {
            smiles,
            chemical_record,
        }
    }
}

impl<C: Clone> MolecularRepresentation<1> for FixedMolecule<C> {
    type ChemicalRecord = C;

    fn variables(&self) -> usize {
        0
    }

    fn build(&self, _: &[f64]) -> [C; 1] {
        [self.chemical_record.clone()]
    }

    fn smiles(&self, _: &[usize]) -> [String; 1] {
        [self.smiles.clone()]
    }
}

impl<
        C,
        M1: MolecularRepresentation<1, ChemicalRecord = C>,
        M2: MolecularRepresentation<1, ChemicalRecord = C>,
    > MolecularRepresentation<2> for (M1, M2)
{
    type ChemicalRecord = C;

    fn build(&self, y: &[f64]) -> [Self::ChemicalRecord; 2] {
        let (y1, y2) = y.split_at(self.0.variables());
        let [cr1] = self.0.build(y1);
        let [cr2] = self.1.build(y2);
        [cr1, cr2]
    }

    fn variables(&self) -> usize {
        self.0.variables() + self.1.variables()
    }

    fn smiles(&self, y: &[usize]) -> [String; 2] {
        let (y1, y2) = y.split_at(self.0.variables());
        let [smiles1] = self.0.smiles(y1);
        let [smiles2] = self.1.smiles(y2);
        [smiles1, smiles2]
    }
}
