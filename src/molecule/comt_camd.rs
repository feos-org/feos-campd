use super::MolecularRepresentation;
use crate::variables::{Constraint, ExplicitVariable, StructureVariables, Variable};
use crate::PropertyModel;
use feos::core::parameter::{Parameter, ParameterError, PureRecord};
use feos::core::EquationOfState;
use feos::ideal_gas::{DipprRecord, IdealGasModel, Joback, JobackRecord};
use feos::pcsaft::{PcSaft, PcSaftParameters, PcSaftRecord};
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::BufReader;
use std::iter;
use std::path::Path;
use std::sync::Arc;

#[derive(Clone, Serialize, Deserialize)]
struct Structure {
    name: String,
    groups: Option<Vec<String>>,
    count: usize,
    rings: usize,
}

#[derive(Clone, Serialize, Deserialize)]
struct Group {
    name: String,
    open_bonds: usize,
    n_max: usize,
    molarweight: f64,
    pcsaft_parameters: PcSaftRecord,
    ideal_gas_parameters: IdealGasRecord,
}

#[derive(Clone, Serialize, Deserialize)]
enum IdealGasRecord {
    Joback(JobackRecord),
    Dippr(DipprRecord),
}

const JOBACK_PARAMETER_OFFSETS: [f64; 5] = [-37.93, 0.21, -3.91e-4, 2.06e-7, 0.0];

#[derive(Clone, Serialize, Deserialize)]
pub struct CoMTCAMD {
    structures: Vec<Structure>,
    groups: Vec<Group>,
}

impl CoMTCAMD {
    pub fn from_json<FP: AsRef<Path>>(file: FP) -> Result<Self, ParameterError> {
        Ok(serde_json::from_reader(BufReader::new(File::open(file)?))?)
    }

    pub fn get_initial_values(&self, structure: &str, groups: &IndexMap<&str, usize>) -> Vec<f64> {
        self.structures
            .iter()
            .map(|s| if s.name == structure { 1.0 } else { 0.0 })
            .chain(self.groups.iter().flat_map(|g| {
                let count = *groups.get(&g.name as &str).unwrap_or(&0);
                std::iter::repeat(1.0)
                    .take(count)
                    .chain(std::iter::repeat(0.0).take(g.n_max - count))
            }))
            .collect()
    }
}

impl MolecularRepresentation for CoMTCAMD {
    fn structure_variables(&self) -> StructureVariables {
        self.structures
            .iter()
            .map(|s| (s.name.clone(), Variable::binary()))
            .chain(self.groups.iter().flat_map(|g| {
                (0..g.n_max).map(|i| {
                    let name = if g.n_max == 1 {
                        g.name.to_string()
                    } else {
                        format!("{}_{i}", g.name)
                    };
                    (name, Variable::binary())
                })
            }))
            .collect()
    }

    fn feature_variables(
        &self,
        index_structure_vars: &[i32],
    ) -> IndexMap<String, ExplicitVariable> {
        let mut index = index_structure_vars[self.structures.len()];
        self.groups
            .iter()
            .map(|g| {
                let mut vars = Vec::new();
                let mut coefs = Vec::new();
                for _ in 0..g.n_max {
                    vars.push(index);
                    coefs.push(1.0);
                    index += 1;
                }
                (
                    g.name.clone(),
                    ExplicitVariable::new(g.name.clone()).linear_struct(vars, coefs),
                )
            })
            .collect()
    }

    fn constraints(&self, index_structure_vars: &[i32]) -> Vec<Constraint> {
        let mut constraints = Vec::new();

        let (y_vars, n_vars) = index_structure_vars.split_at(self.structures.len());

        let mut n: IndexMap<_, _> = IndexMap::new();
        self.groups
            .iter()
            .flat_map(|g| std::iter::repeat(&g.name as &str).take(g.n_max))
            .zip(n_vars.to_vec())
            .for_each(|(g, n_var)| n.entry(g).or_insert(vec![]).push(n_var));

        // Exactly one structure can be active
        constraints.push(
            Constraint::new()
                .linear_struct(y_vars.to_vec(), vec![1.0; y_vars.len()])
                .eqbnd(1.0),
        );

        // Connect molecular structures and segments
        for (structure, &y) in self.structures.iter().zip(y_vars) {
            if let Some(groups) = &structure.groups {
                if !groups.is_empty() {
                    let (vars, coefs): (Vec<_>, Vec<_>) = groups
                        .iter()
                        .flat_map(|g| n[g as &str].iter().copied())
                        .zip(std::iter::repeat(1.0))
                        .chain(std::iter::once((y, -(structure.count as f64))))
                        .unzip();
                    constraints.push(Constraint::new().linear_struct(vars, coefs).eqbnd(0.0));
                }
            } else {
                let coefs = vec![1.0; n_vars.len()];
                constraints.push(
                    Constraint::new()
                        .linear_struct(n_vars.to_vec(), coefs)
                        .eqbnd(1.0),
                );
            };
        }

        // octett rule
        let n_open_bonds = self
            .groups
            .iter()
            .flat_map(|g| iter::repeat(g.open_bonds).take(g.n_max));
        let (vars, coefs): (Vec<_>, Vec<_>) = n_vars
            .iter()
            .zip(n_open_bonds)
            .map(|(&v, o)| (v, 2.0 - o as f64))
            .chain(
                self.structures
                    .iter()
                    .zip(y_vars)
                    .map(|(s, &y)| (y, 2.0 * (s.rings as f64 - 1.0))),
            )
            .unzip();
        constraints.push(Constraint::new().linear_struct(vars, coefs).eqbnd(0.0));

        // // single molecule
        // for &n in n_vars {
        //     let (vars, coefs): (Vec<_>, Vec<_>) = n_vars
        //         .iter()
        //         .zip(&self.n_open_bonds)
        //         .map(|(&v, &o)| (v, if v == n { 2.0 - o as f64 } else { 1.0 }))
        //         .chain([y["alkane"], y["alkene"], y["P"], y["M"]].map(|i| (i, -2.0)))
        //         .unzip();
        //     constraints.push(Constraint::new(vars, coefs).lobnd(0.0));
        // }

        // order groups
        for vars in n.values() {
            for vars in vars.windows(2) {
                constraints.push(
                    Constraint::new()
                        .linear_struct(vars.to_vec(), vec![1.0, -1.0])
                        .lobnd(0.0),
                )
            }
        }

        constraints
    }

    fn smiles(&self, y: &[usize]) -> Vec<String> {
        let (_, n) = y.split_at(self.structures.len());
        let groups: Vec<_> = self
            .groups
            .iter()
            .flat_map(|g| std::iter::repeat(&g.name as &str).take(g.n_max))
            .zip(n)
            .filter_map(|(g, &n)| (n == 1).then_some(g))
            .collect();
        vec![groups.join(",")]
    }
}

#[derive(Clone)]
pub struct CoMTCAMDPropertyModel {
    parameter_names: Vec<String>,
    pub parameters: Vec<IndexMap<String, f64>>,
    parameter_offsets: Vec<f64>,
    viscosity: bool,
}

impl CoMTCAMDPropertyModel {
    pub fn from_json<FP: AsRef<Path>>(file: FP) -> Result<Self, ParameterError> {
        let record: CoMTCAMD = serde_json::from_reader(BufReader::new(File::open(file)?))?;

        // Extract PC-SAFT parameters
        let mut parameters = Vec::new();
        parameters.push(record.groups.iter().map(|s| s.molarweight).collect());
        let pcsaft: Vec<_> = record.groups.iter().map(|s| &s.pcsaft_parameters).collect();
        parameters.push(pcsaft.iter().map(|p| p.m).collect());
        parameters.push(pcsaft.iter().map(|p| p.m * p.sigma.powi(3)).collect());
        parameters.push(pcsaft.iter().map(|p| p.m * p.epsilon_k).collect());
        parameters.push(pcsaft.iter().map(|p| p.mu.unwrap_or_default()).collect());
        let mut parameter_offsets = vec![0.0; 5];
        let mut parameter_names = vec![
            "molarweight".into(),
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
        let ig: Vec<_> = record
            .groups
            .iter()
            .map(|s| match &s.ideal_gas_parameters {
                IdealGasRecord::Joback(JobackRecord { a, b, c, d, e }) => (*a, *b, *c, *d, *e),
                IdealGasRecord::Dippr(dippr) => {
                    let DipprRecord::DIPPR100(v) = &dippr else {
                        panic!("All components have to use Eq. 100!")
                    };
                    let &[a, b, c, d, e, ..] = &v[..] else {
                        panic!("Eq. 100 needs to be provided with 5 coefficients!")
                    };
                    let [a0, b0, c0, d0, e0] = JOBACK_PARAMETER_OFFSETS;
                    (
                        a * 0.001 - a0,
                        b * 0.001 - b0,
                        c * 0.001 - c0,
                        d * 0.001 - d0,
                        e * 0.001 - e0,
                    )
                }
            })
            .collect();
        parameters.push(ig.iter().map(|(a, _, _, _, _)| *a).collect());
        parameters.push(ig.iter().map(|(_, b, _, _, _)| *b).collect());
        parameters.push(ig.iter().map(|(_, _, c, _, _)| *c).collect());
        parameters.push(ig.iter().map(|(_, _, _, d, _)| *d).collect());
        parameters.push(ig.iter().map(|(_, _, _, _, e)| *e).collect());
        parameter_offsets.extend_from_slice(&JOBACK_PARAMETER_OFFSETS);
        parameter_names.extend_from_slice(&[
            "joback_a".into(),
            "joback_b".into(),
            "joback_c".into(),
            "joback_d".into(),
            "joback_e".into(),
        ]);

        let groups: Vec<_> = record.groups.into_iter().map(|r| r.name).collect();
        let parameters = parameters
            .into_iter()
            .map(|p| groups.iter().cloned().zip(p).collect())
            .collect();

        Ok(Self {
            parameter_names,
            parameters,
            parameter_offsets,
            viscosity,
        })
    }
}

impl PropertyModel for CoMTCAMDPropertyModel {
    type EquationOfState = EquationOfState<IdealGasModel, PcSaft>;

    fn parameter_variables(
        &self,
        index_feature_vars: &IndexMap<String, i32>,
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

    fn symmetry_constraint(&self, _: &[i32]) -> Option<Constraint> {
        None
    }

    fn build_eos(&self, parameters: &[f64]) -> Arc<Self::EquationOfState> {
        let (pcsaft, joback) = build_records(parameters, self.viscosity);
        let pcsaft = Arc::new(PcSaftParameters::new_pure(pcsaft).unwrap());
        let joback = Arc::new(Joback::new_pure(joback).unwrap());

        Arc::new(EquationOfState::new(
            Arc::new(IdealGasModel::Joback(joback)),
            Arc::new(PcSaft::new(pcsaft.clone())),
        ))
    }
}

pub(super) fn build_records(
    parameters: &[f64],
    viscosity: bool,
) -> (PureRecord<PcSaftRecord>, PureRecord<JobackRecord>) {
    let pcsaft = parameters[..5].try_into().unwrap();
    let mut k = 5;
    let viscosity = viscosity.then(|| {
        k += 5;
        parameters[k..k + 4].try_into().unwrap()
    });
    let pcsaft = build_pcsaft_record(pcsaft, viscosity);

    let [a, b, c, d, e] = parameters[k..k + 5].try_into().unwrap();
    let joback = PureRecord::new(Default::default(), 0.0, JobackRecord::new(a, b, c, d, e));

    (pcsaft, joback)
}

fn build_pcsaft_record(
    parameters: [f64; 5],
    viscosity: Option<[f64; 4]>,
) -> PureRecord<PcSaftRecord> {
    let [mw, m, m_sigma_3, m_epsilon_k, mu] = parameters;
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
