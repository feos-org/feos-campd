use crate::ChemicalRecord;
use good_lp::{Constraint, Variable};
use num_dual::DualNum;

mod comt_camd;
mod disjunction;
mod polynomial;
mod supermolecule;
pub use comt_camd::CoMTCAMD;
pub use disjunction::Disjunction;
pub use supermolecule::SuperMolecule;

pub trait MolecularRepresentation<const N_Y: usize> {
    fn structure_variables(&self) -> [(i32, i32); N_Y];

    fn constraints(&self, y: [Variable; N_Y]) -> Vec<Constraint>;

    fn smiles(&self, y: [f64; N_Y]) -> String;

    fn build_molecule<D: DualNum<f64> + Copy>(
        &self,
        feature_variables: [D; N_Y],
    ) -> ChemicalRecord<D>;
}
