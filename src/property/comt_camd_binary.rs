use super::comt_camd::build_records;
use super::CoMTCAMDPropertyModel;
use crate::variables::{Constraint, ExplicitVariable};
use crate::PropertyModel;
use feos::core::parameter::{BinaryRecord, Parameter};
use feos::core::EquationOfState;
use feos::ideal_gas::{IdealGasModel, Joback};
use feos::pcsaft::{PcSaft, PcSaftBinaryRecord, PcSaftParameters};
use indexmap::IndexMap;
use std::iter;
use std::sync::Arc;

#[derive(Clone)]
pub struct CoMTCAMDBinaryPropertyModel {
    molecules: [CoMTCAMDPropertyModel; 2],
    k_ij: Option<IndexMap<[String; 2], f64>>,
    symmetry_constraints: bool,
    viscosity: bool,
}

impl CoMTCAMDBinaryPropertyModel {
    pub fn new(
        molecules: [CoMTCAMDPropertyModel; 2],
        k_ij: Option<Vec<BinaryRecord<String, f64>>>,
        symmetry_constraints: bool,
        viscosity: bool,
    ) -> Self {
        let k_ij = k_ij.map(|k_ij| {
            k_ij.iter()
                .flat_map(|r| {
                    [
                        ([r.id1.clone(), r.id2.clone()], r.model_record),
                        ([r.id2.clone(), r.id1.clone()], r.model_record),
                    ]
                })
                .collect()
        });
        Self {
            molecules,
            k_ij,
            symmetry_constraints,
            viscosity,
        }
    }
}

impl PropertyModel<2> for CoMTCAMDBinaryPropertyModel {
    type EquationOfState = EquationOfState<IdealGasModel, PcSaft>;

    fn parameter_variables(
        &self,
        [index_feature_vars1, index_feature_vars2]: &[IndexMap<String, i32>; 2],
    ) -> Vec<ExplicitVariable> {
        let mut variables = self.molecules[0].parameter_variables(&[index_feature_vars1.clone()]);
        variables.extend(self.molecules[1].parameter_variables(&[index_feature_vars2.clone()]));

        // k_ij
        if let Some(k_ij) = &self.k_ij {
            let (vars, qcoefs): (Vec<_>, Vec<_>) = index_feature_vars1
                .iter()
                .flat_map(|(name1, &var1)| {
                    index_feature_vars2.iter().map(move |(name2, &var2)| {
                        (
                            (var1, var2),
                            k_ij.get(&[name1.clone(), name2.clone()])
                                .copied()
                                .unwrap_or_default(),
                        )
                    })
                })
                .unzip();
            let (qvars1, qvars2) = vars.into_iter().unzip();
            variables.push(
                ExplicitVariable::new("k_ij".into()).quadratic_struct(qvars1, qvars2, qcoefs),
            );
        }

        variables
    }

    fn symmetry_constraints(
        &self,
        [index_feature_vars1, index_feature_vars2]: &[IndexMap<String, i32>; 2],
        index_parameter_vars: &[i32],
    ) -> Vec<Constraint> {
        let mut constraints = Vec::new();
        if self.symmetry_constraints {
            // Exclude A/A systems
            let (vars, coefs): (Vec<_>, _) = index_feature_vars1
                .values()
                .chain(index_feature_vars2.values())
                .map(|&i| (i, i))
                .zip(iter::repeat(1.0))
                .chain(
                    index_feature_vars1
                        .values()
                        .copied()
                        .zip(index_feature_vars2.values().copied())
                        .zip(iter::repeat(-2.0)),
                )
                .unzip();
            let (vars1, vars2) = vars.into_iter().unzip();
            constraints.push(
                Constraint::new()
                    .quadratic_struct(vars1, vars2, coefs)
                    .lobnd(1.0),
            );

            // Exclude A/B - B/A symmetry
            let n_params = index_parameter_vars.len() / 2;
            let vars = vec![index_parameter_vars[3], index_parameter_vars[n_params + 3]];
            let coefs = vec![1.0, -1.0];
            constraints.push(Constraint::new().linear_struct(vars, coefs).upbnd(0.0))
        }
        println!("{constraints:?}");
        constraints
    }

    fn build_eos(&self, parameters: &[f64]) -> Arc<Self::EquationOfState> {
        let (p1, p2) = parameters.split_at(self.molecules[0].parameters.len());

        let k_ij = self
            .k_ij
            .is_some()
            .then(|| parameters[parameters.len() - 1]);

        let (pcsaft, joback): (Vec<_>, Vec<_>) = [p1, p2]
            .into_iter()
            .map(|p| build_records(p, self.viscosity))
            .unzip();
        let kij = PcSaftBinaryRecord::new(k_ij, None, None);
        let pcsaft = Arc::new(PcSaftParameters::new_binary(pcsaft, Some(kij)).unwrap());
        let joback = Arc::new(Joback::from_records(joback, None).unwrap());
        Arc::new(EquationOfState::new(
            Arc::new(IdealGasModel::Joback(joback)),
            Arc::new(PcSaft::new(pcsaft.clone())),
        ))
    }
}
