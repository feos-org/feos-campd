use super::comt_camd::build_records;
use super::{CoMTCAMD, CoMTCAMDPropertyModel, MolecularRepresentation};
use crate::variables::{Constraint, ExplicitVariable, StructureVariables};
use crate::PropertyModel;
use feos::core::parameter::{BinaryRecord, Parameter};
use feos::core::EquationOfState;
use feos::ideal_gas::{IdealGasModel, Joback};
use feos::pcsaft::{PcSaft, PcSaftBinaryRecord, PcSaftParameters};
use std::collections::HashMap;
use std::iter;
use std::sync::Arc;

#[derive(Clone)]
pub struct CoMTCAMDBinary {
    molecules: [CoMTCAMD; 2],
    symmetry_constraints: bool,
}

impl CoMTCAMDBinary {
    pub fn new(molecules: [CoMTCAMD; 2], symmetry_constraints: bool) -> Self {
        Self {
            molecules,
            symmetry_constraints,
        }
    }

    pub fn get_initial_values(
        &self,
        structure: [&str; 2],
        groups: [&HashMap<&str, usize>; 2],
    ) -> Vec<f64> {
        let y0 = self.molecules[0].get_initial_values(structure[0], groups[0]);
        let y1 = self.molecules[1].get_initial_values(structure[1], groups[1]);
        [y0, y1].concat()
    }
}

impl MolecularRepresentation for CoMTCAMDBinary {
    fn structure_variables(&self) -> StructureVariables {
        let mut variables = self.molecules[0].structure_variables();
        variables.append(&mut self.molecules[1].structure_variables());
        variables
    }

    fn feature_variables(&self, index_structure_vars: &[i32]) -> HashMap<String, ExplicitVariable> {
        let (index_structure_vars1, index_structure_vars2) =
            index_structure_vars.split_at(self.molecules[0].structure_variables().len());
        self.molecules[0]
            .feature_variables(index_structure_vars1)
            .into_iter()
            .zip(iter::repeat(0))
            .chain(
                self.molecules[1]
                    .feature_variables(index_structure_vars2)
                    .into_iter()
                    .zip(iter::repeat(1)),
            )
            .map(|((name, var), i)| (format!("{name}{i}"), var))
            .collect()
    }

    fn constraints(&self, index_structure_vars: &[i32]) -> Vec<Constraint> {
        let (index_structure_vars1, index_structure_vars2) =
            index_structure_vars.split_at(self.molecules[0].structure_variables().len());
        let mut constraints = self.molecules[0].constraints(index_structure_vars1);
        constraints.append(&mut self.molecules[1].constraints(index_structure_vars2));

        if self.symmetry_constraints {
            // Exclude A/A systems
            let (vars, coefs): (Vec<_>, _) = index_structure_vars1
                .iter()
                .chain(index_structure_vars2)
                .map(|&i| (i, i))
                .zip(iter::repeat(1.0))
                .chain(
                    index_structure_vars1
                        .iter()
                        .copied()
                        .zip(index_structure_vars2.iter().copied())
                        .zip(iter::repeat(-2.0)),
                )
                .unzip();
            let (vars1, vars2) = vars.into_iter().unzip();
            constraints.push(
                Constraint::new()
                    .quadratic_struct(vars1, vars2, coefs)
                    .lobnd(1.0),
            );
        }

        constraints
    }

    fn smiles(&self, y: &[usize]) -> Vec<String> {
        let (y1, y2) = y.split_at(self.molecules[0].structure_variables().len());
        let smiles1 = self.molecules[0].smiles(y1);
        let smiles2 = self.molecules[1].smiles(y2);
        [smiles1, smiles2].concat()
    }
}

#[derive(Clone)]
pub struct CoMTCAMDBinaryPropertyModel {
    molecules: [CoMTCAMDPropertyModel; 2],
    k_ij: Option<HashMap<[String; 2], f64>>,
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

impl PropertyModel for CoMTCAMDBinaryPropertyModel {
    type EquationOfState = EquationOfState<IdealGasModel, PcSaft>;

    fn parameter_variables(
        &self,
        index_feature_vars: &HashMap<String, i32>,
    ) -> Vec<ExplicitVariable> {
        let mut index_feature_vars1 = HashMap::new();
        let mut index_feature_vars2 = HashMap::new();
        index_feature_vars.iter().for_each(|(k, &v)| {
            let (k, i) = k.split_at(k.len() - 1);
            if i == "0" {
                index_feature_vars1.insert(k.to_string(), v);
            } else {
                index_feature_vars2.insert(k.to_string(), v);
            }
        });

        let mut variables = self.molecules[0].parameter_variables(&index_feature_vars1);
        variables.append(&mut self.molecules[1].parameter_variables(&index_feature_vars2));

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

    fn symmetry_constraint(&self, index_parameter_vars: &[i32]) -> Option<Constraint> {
        self.symmetry_constraints.then(|| {
            let n_params = index_parameter_vars.len() / 2;
            let vars = vec![index_parameter_vars[3], index_parameter_vars[n_params + 3]];
            let coefs = vec![1.0, -1.0];
            Constraint::new().linear_struct(vars, coefs).upbnd(0.0)
        })
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
