use crate::{solver::GeneralConstraint, ChemicalRecord};
use feos_ad::{HelmholtzEnergyWrapper, ParametersAD};
use feos_core::EosResult;
use num_dual::DualNum;

#[cfg(test)]
mod orc;
#[cfg(test)]
pub use orc::OrganicRankineCycle;

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

pub trait ProcessModel<E: ParametersAD, const N_X: usize, const N: usize> {
    fn variables(&self) -> [ContinuousVariable; N_X];

    fn constraints(&self) -> Vec<GeneralConstraint>;

    fn evaluate<D: DualNum<f64> + Copy>(
        &self,
        eos: &HelmholtzEnergyWrapper<E, D, N>,
        chemical_records: [&ChemicalRecord<D>; N],
        x: [D; N_X],
    ) -> EosResult<(D, Vec<D>)>;
}
