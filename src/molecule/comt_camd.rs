use super::MolecularRepresentation;
use crate::variables::{ContinuousVariables, DiscreteVariables, LinearConstraint, Variable};
use crate::PropertyModel;
use feos::core::parameter::{Parameter, ParameterError, PureRecord};
use feos::core::EquationOfState;
use feos::ideal_gas::{DipprRecord, Joback, JobackRecord};
use feos::pcsaft::{PcSaft, PcSaftParameters, PcSaftRecord};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
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

const PARAMETER_OFFSETS: [f64; 10] = [
    0.0, 0.0, 0.0, 0.0, 0.0, -37.93, 0.21, -3.91e-4, 2.06e-7, 0.0,
];

#[derive(Clone, Serialize, Deserialize)]
pub struct CoMTCAMD {
    structures: Vec<Structure>,
    groups: Vec<Group>,
    #[serde(skip)]
    n_open_bonds: Vec<usize>,
    #[serde(skip)]
    parameters: [Vec<f64>; 10],
}

impl CoMTCAMD {
    pub fn from_json<FP: AsRef<Path>>(file: FP) -> Result<Self, ParameterError> {
        let mut result: Self = serde_json::from_reader(BufReader::new(File::open(file)?))?;

        // Replicate groups
        let segments: Vec<_> = result
            .groups
            .iter()
            .flat_map(|g| std::iter::repeat(g).take(g.n_max))
            .cloned()
            .collect();

        // Extract PC-SAFT parameters
        result.parameters[0].extend(segments.iter().map(|s| s.molarweight));
        let pcsaft: Vec<_> = segments.iter().map(|s| &s.pcsaft_parameters).collect();
        result.parameters[1].extend(pcsaft.iter().map(|p| p.m));
        result.parameters[2].extend(pcsaft.iter().map(|p| p.m * p.sigma.powi(3)));
        result.parameters[3].extend(pcsaft.iter().map(|p| p.m * p.epsilon_k));
        result.parameters[4].extend(pcsaft.iter().map(|p| p.mu.unwrap_or_default()));

        // Extract valence
        result
            .n_open_bonds
            .extend(segments.iter().map(|s| s.open_bonds));

        // Extract ideal_gas parameters
        let ig: Vec<_> = segments
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
                    let [_, _, _, _, _, a0, b0, c0, d0, e0] = PARAMETER_OFFSETS;
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
        result.parameters[5].extend(ig.iter().map(|(a, _, _, _, _)| *a));
        result.parameters[6].extend(ig.iter().map(|(_, b, _, _, _)| *b));
        result.parameters[7].extend(ig.iter().map(|(_, _, c, _, _)| *c));
        result.parameters[8].extend(ig.iter().map(|(_, _, _, d, _)| *d));
        result.parameters[9].extend(ig.iter().map(|(_, _, _, _, e)| *e));

        Ok(result)
    }
}

impl CoMTCAMD {
    pub fn get_initial_values(&self, structure: &str, groups: &HashMap<&str, usize>) -> Vec<f64> {
        self.groups
            .iter()
            .flat_map(|g| {
                let count = *groups.get(&g.name as &str).unwrap_or(&0);
                std::iter::repeat(1.0)
                    .take(count)
                    .chain(std::iter::repeat(0.0).take(g.n_max - count))
            })
            .chain(
                self.structures
                    .iter()
                    .map(|s| if s.name == structure { 1.0 } else { 0.0 }),
            )
            .collect()
    }
}

impl MolecularRepresentation<1> for CoMTCAMD {
    type ChemicalRecord = Vec<f64>;

    fn build(&self, _: &[f64], p: &[f64]) -> [Self::ChemicalRecord; 1] {
        [p.to_vec()]
    }

    fn structure_variables(&self) -> DiscreteVariables {
        vec![Variable::binary(); self.parameters[0].len() + self.structures.len()].into()
    }

    fn parameter_variables(&self) -> ContinuousVariables {
        vec![Variable::free(); 10].into()
    }

    fn determine_parameters(&self, y: &[f64]) -> Vec<f64> {
        self.parameters
            .iter()
            .zip(PARAMETER_OFFSETS)
            .map(|(pars, offset)| offset + pars.iter().zip(y).map(|(p, y)| p * y).sum::<f64>())
            .collect()
    }

    fn constraints(
        &self,
        index_structure_vars: &[i32],
        index_parameter_vars: Option<&[i32]>,
    ) -> Vec<LinearConstraint> {
        let mut constraints = Vec::new();

        let (n_vars, y_vars) = index_structure_vars.split_at(self.parameters[0].len());

        let mut n: HashMap<_, _> = HashMap::new();
        self.groups
            .iter()
            .flat_map(|g| std::iter::repeat(&g.name as &str).take(g.n_max))
            .zip(n_vars.to_vec())
            .for_each(|(g, n_var)| n.entry(g).or_insert(vec![]).push(n_var));

        // Exactly one structure can be active
        constraints
            .push(LinearConstraint::new(y_vars.to_vec(), vec![1.0; y_vars.len()]).eqbnd(1.0));

        // Connect molecular structures and segments
        for (structure, &y) in self.structures.iter().zip(y_vars) {
            if let Some(groups) = &structure.groups {
                let (vars, coefs): (Vec<_>, Vec<_>) = groups
                    .iter()
                    .flat_map(|g| n[g as &str].iter().copied())
                    .zip(std::iter::repeat(1.0))
                    .chain(std::iter::once((y, -(structure.count as f64))))
                    .unzip();
                constraints.push(LinearConstraint::new(vars, coefs).eqbnd(0.0));
            } else {
                let coefs = vec![1.0; n_vars.len()];
                constraints.push(LinearConstraint::new(n_vars.to_vec(), coefs).eqbnd(1.0));
            };
        }

        // octett rule
        let (vars, coefs): (Vec<_>, Vec<_>) = n_vars
            .iter()
            .zip(&self.n_open_bonds)
            .map(|(&v, &o)| (v, 2.0 - o as f64))
            .chain(
                self.structures
                    .iter()
                    .zip(y_vars)
                    .map(|(s, &y)| (y, 2.0 * (s.rings as f64 - 1.0))),
            )
            .unzip();
        constraints.push(LinearConstraint::new(vars, coefs).eqbnd(0.0));

        // // single molecule
        // for &n in n_vars {
        //     let (vars, coefs): (Vec<_>, Vec<_>) = n_vars
        //         .iter()
        //         .zip(&self.n_open_bonds)
        //         .map(|(&v, &o)| (v, if v == n { 2.0 - o as f64 } else { 1.0 }))
        //         .chain([y["alkane"], y["alkene"], y["P"], y["M"]].map(|i| (i, -2.0)))
        //         .unzip();
        //     constraints.push(LinearConstraint::new(vars, coefs).lobnd(0.0));
        // }

        // order groups
        for vars in n.values() {
            for vars in vars.windows(2) {
                constraints.push(LinearConstraint::new(vars.to_vec(), vec![1.0, -1.0]).lobnd(0.0))
            }
        }

        // parameter constraints
        if let Some(p) = index_parameter_vars {
            constraints.extend(p.iter().zip(&self.parameters).zip(&PARAMETER_OFFSETS).map(
                |((&p, pars), &offset)| {
                    let (vars, coefs): (Vec<_>, Vec<_>) = n_vars
                        .iter()
                        .zip(pars)
                        .chain(std::iter::once((&p, &-1.0)))
                        .unzip();
                    LinearConstraint::new(vars, coefs).eqbnd(-offset)
                },
            ));
        }

        constraints
    }

    fn smiles(&self, y: &[usize]) -> [String; 1] {
        let (n, _) = y.split_at(self.parameters[0].len());
        let groups: Vec<_> = self
            .groups
            .iter()
            .flat_map(|g| std::iter::repeat(&g.name as &str).take(g.n_max))
            .zip(n)
            .filter_map(|(g, &n)| (n == 1).then_some(g))
            .collect();
        [groups.join(",")]
    }
}

pub struct CoMTCAMDPropertyModel;

impl PropertyModel<Vec<f64>> for CoMTCAMDPropertyModel {
    type EquationOfState = EquationOfState<Joback, PcSaft>;

    fn build_eos(
        &self,
        chemical_records: &[Vec<f64>],
    ) -> Result<Arc<Self::EquationOfState>, ParameterError> {
        let (pcsaft, joback): (Vec<_>, Vec<_>) = chemical_records
            .iter()
            .map(|p| {
                let &[mw, m, m_sigma_3, m_epsilon_k, mu] = &p[..5] else {
                    unreachable!()
                };
                let pcsaft = PureRecord::new(
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
                        None,
                        None,
                        None,
                    ),
                );
                let &[a, b, c, d, e] = &p[5..] else {
                    unreachable!()
                };
                let joback =
                    PureRecord::new(Default::default(), 0.0, JobackRecord::new(a, b, c, d, e));
                (pcsaft, joback)
            })
            .unzip();
        let pcsaft = Arc::new(PcSaftParameters::from_records(pcsaft, None).unwrap());
        let joback = Arc::new(Joback::from_records(joback, None).unwrap());
        Ok(Arc::new(EquationOfState::new(
            joback,
            Arc::new(PcSaft::new(pcsaft.clone())),
        )))
    }
}
