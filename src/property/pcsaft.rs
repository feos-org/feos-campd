use super::ideal_gas::JOBACK_PARAMETER_OFFSETS;
use super::IdealGasRecord;
use crate::variables::{Constraint, ExplicitVariable};
use crate::PropertyModel;
use feos::core::parameter::{
    BinaryRecord, Identifier, IdentifierOption, Parameter, ParameterError, PureRecord,
    SegmentRecord,
};
use feos::core::EquationOfState;
use feos::ideal_gas::{Joback, JobackRecord};
use feos::pcsaft::{PcSaft, PcSaftBinaryRecord, PcSaftParameters, PcSaftRecord};
use indexmap::IndexMap;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::iter;
use std::path::Path;
use std::sync::Arc;

struct GroupRecord {
    identifier: String,
    molarweight: f64,
    pcsaft_record: PcSaftRecord,
    ideal_gas_record: IdealGasRecord,
}

impl From<(SegmentRecord<PcSaftRecord>, IdealGasRecord)> for GroupRecord {
    fn from((pcsaft, ideal_gas_record): (SegmentRecord<PcSaftRecord>, IdealGasRecord)) -> Self {
        Self {
            identifier: pcsaft.identifier,
            molarweight: pcsaft.molarweight,
            pcsaft_record: pcsaft.model_record,
            ideal_gas_record,
        }
    }
}

/// The homosegmented group contribution method for the PC-SAFT equation of state.
#[derive(Clone)]
pub struct PcSaftPropertyModel {
    parameter_names: Vec<String>,
    parameters: Vec<IndexMap<String, f64>>,
    parameter_offsets: Vec<f64>,
    viscosity: bool,
    k_ij: Option<IndexMap<[String; 2], f64>>,
    symmetry_constraints: bool,
}

impl PcSaftPropertyModel {
    fn new(
        groups: Vec<GroupRecord>,
        k_ij: Option<Vec<BinaryRecord<String, f64>>>,
        symmetry_constraints: bool,
    ) -> Self {
        // Extract PC-SAFT parameters
        let mut parameters = Vec::new();
        parameters.push(groups.iter().map(|s| s.molarweight).collect());
        let pcsaft: Vec<_> = groups.iter().map(|s| &s.pcsaft_record).collect();
        parameters.push(pcsaft.iter().map(|_| 1.0).collect());
        parameters.push(pcsaft.iter().map(|p| p.m).collect());
        parameters.push(pcsaft.iter().map(|p| p.m * p.sigma.powi(3)).collect());
        parameters.push(pcsaft.iter().map(|p| p.m * p.epsilon_k).collect());
        parameters.push(pcsaft.iter().map(|p| p.mu.unwrap_or_default()).collect());
        let mut parameter_offsets = vec![0.0; 6];
        let mut parameter_names = vec![
            "molarweight".into(),
            "s".into(),
            "m".into(),
            "m_sigma3".into(),
            "m_epsilon".into(),
            "mu".into(),
        ];

        // Extract viscosity parameters
        let viscosity_records: Option<Vec<_>> = pcsaft.iter().map(|p| p.viscosity).collect();
        let viscosity = viscosity_records.is_some();
        if let Some(viscosity) = viscosity_records {
            let [mut a_vec, mut b_vec, mut c_vec, mut d_vec] =
                [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
            pcsaft
                .iter()
                .zip(viscosity.iter())
                .for_each(|(&p, &[a, b, c, d])| {
                    a_vec.push(p.m * p.sigma.powi(3) * a);
                    b_vec.push(p.m * p.sigma.powi(3) * b);
                    c_vec.push(c);
                    d_vec.push(d);
                });
            parameters.extend_from_slice(&[a_vec, b_vec, c_vec, d_vec]);
            parameter_offsets.extend_from_slice(&[0.0; 4]);
            parameter_names.extend_from_slice(&[
                "viscosity_a".into(),
                "viscosity_b".into(),
                "viscosity_c".into(),
                "viscosity_d".into(),
            ])
        }

        // Extract ideal_gas parameters
        let ig: Vec<_> = groups
            .iter()
            .map(|s| s.ideal_gas_record.joback_params())
            .collect();
        parameters.push(ig.iter().map(|[a, _, _, _, _]| *a).collect());
        parameters.push(ig.iter().map(|[_, b, _, _, _]| *b).collect());
        parameters.push(ig.iter().map(|[_, _, c, _, _]| *c).collect());
        parameters.push(ig.iter().map(|[_, _, _, d, _]| *d).collect());
        parameters.push(ig.iter().map(|[_, _, _, _, e]| *e).collect());
        parameter_offsets.extend_from_slice(&JOBACK_PARAMETER_OFFSETS);
        parameter_names.extend_from_slice(&[
            "joback_a".into(),
            "joback_b".into(),
            "joback_c".into(),
            "joback_d".into(),
            "joback_e".into(),
        ]);

        let groups: Vec<_> = groups.iter().map(|r| r.identifier.clone()).collect();
        let parameters = parameters
            .into_iter()
            .map(|p| groups.iter().cloned().zip(p).collect())
            .collect();

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
            parameter_names,
            parameters,
            parameter_offsets,
            viscosity,
            k_ij,
            symmetry_constraints,
        }
    }

    pub fn from_records_groups(
        pcsaft: Vec<SegmentRecord<PcSaftRecord>>,
        ideal_gas: Vec<SegmentRecord<IdealGasRecord>>,
        binary: Option<Vec<BinaryRecord<String, f64>>>,
        symmetry_constraints: bool,
    ) -> Self {
        let mut ideal_gas: HashMap<_, _> = ideal_gas
            .into_iter()
            .map(|record| (record.identifier, record.model_record))
            .collect();
        let groups: Vec<_> = pcsaft
            .into_iter()
            .filter_map(|record| ideal_gas.remove(&record.identifier).map(|ig| (record, ig)))
            .map(GroupRecord::from)
            .collect();
        Self::new(groups, binary, symmetry_constraints)
    }

    pub fn from_records_molecules(
        pcsaft: Vec<PureRecord<PcSaftRecord>>,
        ideal_gas: Vec<PureRecord<IdealGasRecord>>,
        binary: Option<Vec<BinaryRecord<Identifier, PcSaftBinaryRecord>>>,
        identifier_option: IdentifierOption,
        symmetry_constraints: bool,
    ) -> Self {
        let pcsaft = pcsaft
            .into_iter()
            .filter_map(|record| {
                record
                    .identifier
                    .as_string(identifier_option)
                    .map(|identifier| {
                        SegmentRecord::new(identifier, record.molarweight, record.model_record)
                    })
            })
            .collect();
        let ideal_gas = ideal_gas
            .into_iter()
            .filter_map(|record| {
                record
                    .identifier
                    .as_string(identifier_option)
                    .map(|identifier| {
                        SegmentRecord::new(identifier, record.molarweight, record.model_record)
                    })
            })
            .collect();
        let binary = binary.map(|binary| {
            binary
                .into_iter()
                .filter_map(|record| {
                    record.id1.as_string(identifier_option).and_then(|id1| {
                        record
                            .id2
                            .as_string(identifier_option)
                            .map(|id2| BinaryRecord::new(id1, id2, record.model_record.k_ij))
                    })
                })
                .collect()
        });
        Self::from_records_groups(pcsaft, ideal_gas, binary, symmetry_constraints)
    }

    pub fn from_json_groups<FP: AsRef<Path>>(
        file_pcsaft: FP,
        file_ideal_gas: FP,
        file_binary: Option<FP>,
        symmetry_constraints: bool,
    ) -> Result<Self, ParameterError> {
        let pcsaft = SegmentRecord::from_json(file_pcsaft)?;
        let ideal_gas = SegmentRecord::from_json(file_ideal_gas)?;
        let binary = file_binary
            .map(|f| BinaryRecord::from_json(f))
            .transpose()?;
        Ok(Self::from_records_groups(
            pcsaft,
            ideal_gas,
            binary,
            symmetry_constraints,
        ))
    }

    pub fn from_json_molecules<FP: AsRef<Path>>(
        file_pcsaft: FP,
        file_ideal_gas: FP,
        file_binary: Option<FP>,
        identifier_option: IdentifierOption,
        symmetry_constraints: bool,
    ) -> Result<Self, ParameterError> {
        let pcsaft = serde_json::from_reader(BufReader::new(File::open(file_pcsaft)?))?;
        let ideal_gas = serde_json::from_reader(BufReader::new(File::open(file_ideal_gas)?))?;
        let binary = file_binary
            .map(|f| {
                Ok::<_, ParameterError>(serde_json::from_reader(BufReader::new(File::open(f)?))?)
            })
            .transpose()?;
        Ok(Self::from_records_molecules(
            pcsaft,
            ideal_gas,
            binary,
            identifier_option,
            symmetry_constraints,
        ))
    }
}

impl PropertyModel<1> for PcSaftPropertyModel {
    type EquationOfState = EquationOfState<Joback, PcSaft>;

    fn parameter_variables(
        &self,
        [index_feature_vars]: &[IndexMap<String, i32>; 1],
    ) -> Vec<ExplicitVariable> {
        self.parameters
            .iter()
            .zip(&self.parameter_names)
            .zip(&self.parameter_offsets)
            .map(|((pars, names), &offset)| {
                let (vars, coefs): (Vec<_>, Vec<_>) = index_feature_vars
                    .iter()
                    .filter(|&(name, _)| (!name.contains('-')))
                    .map(|(name, vars)| (vars, pars[name]))
                    .unzip();
                ExplicitVariable::new(names.clone())
                    .linear_struct(vars, coefs)
                    .cons(offset)
            })
            .collect()
    }

    fn symmetry_constraints(&self, _: &[IndexMap<String, i32>; 1], _: &[i32]) -> Vec<Constraint> {
        vec![]
    }

    fn build_eos(&self, parameters: &[f64]) -> Arc<Self::EquationOfState> {
        let (pcsaft, joback) = build_records(parameters, self.viscosity);
        let pcsaft = Arc::new(PcSaftParameters::new_pure(pcsaft).unwrap());
        let joback = Arc::new(Joback::new_pure(joback).unwrap());

        Arc::new(EquationOfState::new(
            joback,
            Arc::new(PcSaft::new(pcsaft.clone())),
        ))
    }
}

pub(super) fn build_records(
    parameters: &[f64],
    viscosity: bool,
) -> (PureRecord<PcSaftRecord>, PureRecord<JobackRecord>) {
    let pcsaft = parameters[..6].try_into().unwrap();
    let mut k = 6;
    let viscosity = viscosity.then(|| {
        let viscosity = parameters[k..k + 4].try_into().unwrap();
        k += 4;
        viscosity
    });
    let pcsaft = build_pcsaft_record(pcsaft, viscosity);

    let [a, b, c, d, e] = parameters[k..k + 5].try_into().unwrap();
    let joback = PureRecord::new(Default::default(), 0.0, JobackRecord::new(a, b, c, d, e));

    (pcsaft, joback)
}

fn build_pcsaft_record(
    parameters: [f64; 6],
    viscosity: Option<[f64; 4]>,
) -> PureRecord<PcSaftRecord> {
    let [mw, _, m, m_sigma_3, m_epsilon_k, mu] = parameters;
    let viscosity = viscosity.map(|viscosity| {
        let [mut a, mut b, c, d] = viscosity;
        a -= 0.5 * m.ln();
        b /= m_sigma_3.powf(0.45);
        [a, b, c, d]
    });
    PureRecord::new(
        Default::default(),
        mw,
        PcSaftRecord::new(
            m,
            (m_sigma_3 / m).cbrt(),
            m_epsilon_k / m,
            Some(mu),
            None,
            None,
            None,
            None,
            None,
            None,
            viscosity,
            None,
            None,
        ),
    )
}

impl PropertyModel<2> for PcSaftPropertyModel {
    type EquationOfState = EquationOfState<Joback, PcSaft>;

    fn parameter_variables(
        &self,
        [index_feature_vars1, index_feature_vars2]: &[IndexMap<String, i32>; 2],
    ) -> Vec<ExplicitVariable> {
        let mut variables = self.parameter_variables(&[index_feature_vars1.clone()]);
        variables.extend(self.parameter_variables(&[index_feature_vars2.clone()]));

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
            let (qvars1, qvars2): (Vec<_>, Vec<_>) = vars.into_iter().unzip();
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

        constraints
    }

    fn build_eos(&self, parameters: &[f64]) -> Arc<Self::EquationOfState> {
        let (p1, p2) = parameters.split_at(self.parameters.len());

        let k_ij = self
            .k_ij
            .is_some()
            .then(|| parameters[parameters.len() - 1] / (p1[1] * p2[1]));

        let (pcsaft, joback): (Vec<_>, Vec<_>) = [p1, p2]
            .into_iter()
            .map(|p| build_records(p, self.viscosity))
            .unzip();
        let kij = PcSaftBinaryRecord::new(k_ij, None, None);
        let pcsaft = Arc::new(PcSaftParameters::new_binary(pcsaft, Some(kij)).unwrap());
        let joback = Arc::new(Joback::from_records(joback, None).unwrap());
        Arc::new(EquationOfState::new(
            joback,
            Arc::new(PcSaft::new(pcsaft.clone())),
        ))
    }
}
