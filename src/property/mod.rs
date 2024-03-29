use crate::variables::ExplicitVariable;
use crate::Constraint;
use feos::core::{IdealGas, Residual};
use indexmap::IndexMap;
#[cfg(feature = "knitro_rs")]
use knitro_rs::{Knitro, KnitroError};
use std::array;
use std::sync::Arc;

mod gc_pcsaft;
mod ideal_gas;
mod pcsaft;
pub use gc_pcsaft::GcPcSaftPropertyModel;
use ideal_gas::IdealGasRecord;
pub use pcsaft::PcSaftPropertyModel;

/// A generic property model to be used in an [OptimizationProblem](super::OptimizationProblem).
pub trait PropertyModel<const N: usize> {
    type EquationOfState: Residual + IdealGas;

    fn parameter_variables(
        &self,
        index_feature_vars: &[IndexMap<String, i32>; N],
    ) -> Vec<ExplicitVariable>;

    fn evaluate_parameter_variables(&self, feature_vars: &[IndexMap<String, f64>; N]) -> Vec<f64> {
        let mut vars = Vec::new();
        let mut index = 0;
        let index_feature_vars = array::from_fn(|i| {
            let (index_feature_vars, v): (IndexMap<_, _>, Vec<_>) = feature_vars[i]
                .iter()
                .enumerate()
                .map(|(i, (k, &v))| ((k.clone(), (index + i) as i32), v))
                .unzip();
            index += v.len();
            vars.extend(v);
            index_feature_vars
        });
        self.parameter_variables(&index_feature_vars)
            .into_iter()
            .map(|v| v.evaluate(&vars))
            .collect()
    }

    fn symmetry_constraints(
        &self,
        index_feature_vars: &[IndexMap<String, i32>; N],
        index_parameter_vars: &[i32],
    ) -> Vec<Constraint>;

    fn build_eos(&self, parameters: &[f64]) -> Arc<Self::EquationOfState>;

    #[cfg(feature = "knitro_rs")]
    fn setup_knitro(
        &self,
        kc: &Knitro,
        index_feature_vars: &[IndexMap<String, i32>; N],
    ) -> Result<Vec<i32>, KnitroError> {
        // define parameter variables
        let index_parameter_vars = self
            .parameter_variables(index_feature_vars)
            .into_iter()
            .map(|f| f.setup_knitro(kc))
            .collect::<Result<Vec<_>, _>>()?;

        // Add symmetry constraint
        self.symmetry_constraints(index_feature_vars, &index_parameter_vars)
            .into_iter()
            .try_for_each(|c| c.setup_knitro(kc))?;

        Ok(index_parameter_vars)
    }
}
