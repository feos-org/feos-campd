use crate::{solver::GeneralConstraint, ChemicalRecord};
use feos::core::{FeosResult, Total};
use nalgebra::Const;
use num_dual::DualNum;

#[cfg(test)]
mod orc;
#[cfg(test)]
pub use orc::OrganicRankineCycle;

/// Information about bounds and initial values of continuous variables.
#[derive(Clone, Copy)]
pub struct ContinuousVariable {
    pub lobnd: f64,
    pub upbnd: f64,
    pub init: f64,
}

impl ContinuousVariable {
    pub fn new(lobnd: f64, upbnd: f64, init: f64) -> Self {
        Self { lobnd, upbnd, init }
    }
}

/// Generic process model to be used in an [IntegratedDesign](../IntegratedDesign).
pub trait ProcessModel<const N_X: usize, const N: usize> {
    fn variables(&self) -> [ContinuousVariable; N_X];

    fn constraints(&self) -> Vec<GeneralConstraint>;

    fn evaluate<E: Total<Const<N>, D>, D: DualNum<f64> + Copy>(
        &self,
        eos: &E,
        chemical_records: [&ChemicalRecord<D>; N],
        x: [D; N_X],
    ) -> FeosResult<(D, Vec<D>)>;
}
