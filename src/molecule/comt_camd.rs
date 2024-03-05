use super::MolecularRepresentation;
use crate::variables::{Constraint, ParameterVariables, StructureVariables, Variable};
use crate::PropertyModel;
use feos::core::parameter::{BinaryRecord, Parameter, ParameterError, PureRecord};
use feos::core::EquationOfState;
use feos::ideal_gas::{DipprRecord, IdealGasModel, Joback, JobackRecord};
use feos::pcsaft::{PcSaft, PcSaftBinaryRecord, PcSaftParameters, PcSaftRecord};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
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
    #[serde(skip)]
    n_open_bonds: Vec<usize>,
    #[serde(skip)]
    parameters: Vec<Vec<f64>>,
    #[serde(skip)]
    parameter_offsets: Vec<f64>,
    #[serde(skip)]
    viscosity: bool,
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
        let p = &mut result.parameters;
        p.push(segments.iter().map(|s| s.molarweight).collect());
        let pcsaft: Vec<_> = segments.iter().map(|s| &s.pcsaft_parameters).collect();
        p.push(pcsaft.iter().map(|p| p.m).collect());
        p.push(pcsaft.iter().map(|p| p.m * p.sigma.powi(3)).collect());
        p.push(pcsaft.iter().map(|p| p.m * p.epsilon_k).collect());
        p.push(pcsaft.iter().map(|p| p.mu.unwrap_or_default()).collect());
        result.parameter_offsets.extend_from_slice(&[0.0; 5]);

        // Extract viscosity parameters
        let viscosity: Option<Vec<_>> = pcsaft.iter().map(|p| p.viscosity).collect();
        if let Some(viscosity) = viscosity {
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
            p.extend_from_slice(&[a_vec, b_vec, c_vec, d_vec]);
            result.parameter_offsets.extend_from_slice(&[0.0; 4]);
            result.viscosity = true;
        }

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
        p.push(ig.iter().map(|(a, _, _, _, _)| *a).collect());
        p.push(ig.iter().map(|(_, b, _, _, _)| *b).collect());
        p.push(ig.iter().map(|(_, _, c, _, _)| *c).collect());
        p.push(ig.iter().map(|(_, _, _, d, _)| *d).collect());
        p.push(ig.iter().map(|(_, _, _, _, e)| *e).collect());
        result
            .parameter_offsets
            .extend_from_slice(&JOBACK_PARAMETER_OFFSETS);

        Ok(result)
    }

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

impl MolecularRepresentation for CoMTCAMD {
    type ChemicalRecord = Vec<f64>;

    fn build(&self, _: &[f64], p: &[f64]) -> Self::ChemicalRecord {
        p.to_vec()
    }

    fn structure_variables(&self) -> StructureVariables {
        vec![Variable::binary(); self.parameters[0].len() + self.structures.len()].into()
    }

    fn parameter_variables(&self) -> ParameterVariables {
        vec![Variable::free(); self.parameters.len()].into()
    }

    fn determine_parameters(&self, y: &[f64]) -> Vec<f64> {
        self.parameters
            .iter()
            .zip(&self.parameter_offsets)
            .map(|(pars, offset)| offset + pars.iter().zip(y).map(|(p, y)| p * y).sum::<f64>())
            .collect()
    }

    fn constraints(
        &self,
        index_structure_vars: &[i32],
        index_parameter_vars: Option<&[i32]>,
    ) -> Vec<Constraint> {
        let mut constraints = Vec::new();

        let (n_vars, y_vars) = index_structure_vars.split_at(self.parameters[0].len());

        let mut n: HashMap<_, _> = HashMap::new();
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

        // parameter constraints
        if let Some(p) = index_parameter_vars {
            constraints.extend(
                p.iter()
                    .zip(&self.parameters)
                    .zip(&self.parameter_offsets)
                    .map(|((&p, pars), &offset)| {
                        let (vars, coefs): (Vec<_>, Vec<_>) = n_vars
                            .iter()
                            .zip(pars)
                            .chain(std::iter::once((&p, &-1.0)))
                            .unzip();
                        Constraint::new().linear_struct(vars, coefs).eqbnd(-offset)
                    }),
            );
        }

        constraints
    }

    fn smiles(&self, y: &[usize]) -> Vec<String> {
        let (n, _) = y.split_at(self.parameters[0].len());
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
pub struct CoMTCAMDBinary {
    molecules: [CoMTCAMD; 2],
    k_ij: Option<Vec<Vec<f64>>>,
    symmetry_constraints: bool,
    viscosity: bool,
}

impl CoMTCAMDBinary {
    pub fn new(
        molecules: [CoMTCAMD; 2],
        symmetry_constraints: bool,
        k_ij: Option<Vec<BinaryRecord<String, f64>>>,
    ) -> Self {
        let viscosity = if molecules[0].viscosity == molecules[1].viscosity {
            molecules[0].viscosity
        } else {
            panic!("Inconsistent viscosity data!");
        };
        let k_ij = k_ij.map(|k_ij| {
            let record_map: HashMap<_, _> = k_ij
                .iter()
                .flat_map(|r| {
                    [
                        ((&r.id1, &r.id2), r.model_record),
                        ((&r.id2, &r.id1), r.model_record),
                    ]
                })
                .collect();
            let segments = molecules.clone().map(|m| {
                m.groups
                    .iter()
                    .flat_map(|g| std::iter::repeat(g).take(g.n_max))
                    .cloned()
                    .collect::<Vec<_>>()
            });
            segments[0]
                .iter()
                .map(|s1| {
                    segments[1]
                        .iter()
                        .map(|s2| {
                            record_map
                                .get(&(&s1.name, &s2.name))
                                .copied()
                                .unwrap_or(0.0)
                        })
                        .collect()
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

impl MolecularRepresentation for CoMTCAMDBinary {
    type ChemicalRecord = ([Vec<f64>; 2], f64);

    fn build(&self, y: &[f64], p: &[f64]) -> Self::ChemicalRecord {
        let (y1, y2) = y.split_at(self.molecules[0].structure_variables().len());
        let (p1, p2) = p.split_at(self.molecules[0].parameter_variables().len());
        let (p2, k_ij) = p2.split_at(self.molecules[1].parameter_variables().len());
        let k_ij = match *k_ij {
            [k_ij] => k_ij,
            [] => 0.0,
            _ => unreachable!(),
        };
        let cr1 = self.molecules[0].build(y1, p1);
        let cr2 = self.molecules[1].build(y2, p2);
        ([cr1, cr2], k_ij)
    }

    fn structure_variables(&self) -> StructureVariables {
        let mut variables = self.molecules[0].structure_variables();
        variables.append(&mut self.molecules[1].structure_variables());
        variables
    }

    fn parameter_variables(&self) -> ParameterVariables {
        let mut variables = self.molecules[0].parameter_variables();
        variables.append(&mut self.molecules[1].parameter_variables());
        // k_ij
        if self.k_ij.is_some() {
            variables.push(Variable::free());
        }
        variables
    }

    fn determine_parameters(&self, y: &[f64]) -> Vec<f64> {
        let (y1, y2) = y.split_at(self.molecules[0].structure_variables().len());
        let p1 = self.molecules[0].determine_parameters(y1);
        let p2 = self.molecules[0].determine_parameters(y2);
        let k_ij = if let Some(k_ij) = &self.k_ij {
            k_ij.iter()
                .zip(y1)
                .flat_map(|(k_i, &n1)| k_i.iter().zip(y2).map(move |(&k, &n2)| k * n1 * n2))
                .sum::<f64>()
        } else {
            0.0
        };
        [p1, p2, vec![k_ij]].concat()
    }

    fn constraints(
        &self,
        index_structure_vars: &[i32],
        index_parameter_vars: Option<&[i32]>,
    ) -> Vec<Constraint> {
        let (index_structure_vars1, index_structure_vars2) =
            index_structure_vars.split_at(self.molecules[0].structure_variables().len());
        let (index_parameter_vars1, index_parameter_vars2, index_k_ij) =
            if let Some(index_parameter_vars) = index_parameter_vars {
                let (p1, p2) =
                    index_parameter_vars.split_at(self.molecules[0].parameter_variables().len());
                let (p2, k_ij) = p2.split_at(self.molecules[1].parameter_variables().len());
                (Some(p1), Some(p2), k_ij)
            } else {
                (None, None, [].as_slice())
            };
        let mut constraints =
            self.molecules[0].constraints(index_structure_vars1, index_parameter_vars1);
        constraints.append(
            &mut self.molecules[1].constraints(index_structure_vars2, index_parameter_vars2),
        );

        // sum rule for k_ij
        if let Some(k_ij) = &self.k_ij {
            let (vars, coefs): (Vec<_>, Vec<_>) = k_ij
                .iter()
                .zip(index_structure_vars1)
                .flat_map(|(k_i, &n1)| {
                    k_i.iter()
                        .zip(index_structure_vars2)
                        .map(move |(&k, &n2)| ((n1, n2), k))
                })
                .unzip();
            let (vars1, vars2): (Vec<_>, Vec<_>) = vars.into_iter().unzip();
            constraints.push(
                Constraint::new()
                    .quadratic_struct(vars1, vars2, coefs)
                    .linear_struct(index_k_ij.to_vec(), vec![-1.0])
                    .eqbnd(0.0),
            );
        }

        if self.symmetry_constraints {
            // Exclude A/B - B/A symmetries
            let coefs1 = &self.molecules[0].parameters[3];
            let coefs2 = &self.molecules[1].parameters[3];
            let (vars, coefs): (Vec<_>, Vec<_>) = index_structure_vars1
                .iter()
                .zip(coefs1.iter().copied())
                .chain(index_structure_vars2.iter().zip(coefs2.iter().map(|c| -c)))
                .unzip();
            constraints.push(Constraint::new().linear_struct(vars, coefs).upbnd(0.0));

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

impl PropertyModel<Vec<f64>> for CoMTCAMD {
    type EquationOfState = EquationOfState<IdealGasModel, PcSaft>;

    fn build_eos(
        &self,
        parameters: Vec<f64>,
    ) -> Result<Arc<Self::EquationOfState>, ParameterError> {
        let (pcsaft, joback) = build_records(parameters, self.viscosity);
        let pcsaft = Arc::new(PcSaftParameters::new_pure(pcsaft).unwrap());
        let joback = Arc::new(Joback::new_pure(joback).unwrap());

        Ok(Arc::new(EquationOfState::new(
            Arc::new(IdealGasModel::Joback(joback)),
            Arc::new(PcSaft::new(pcsaft.clone())),
        )))
    }
}

impl PropertyModel<([Vec<f64>; 2], f64)> for CoMTCAMDBinary {
    type EquationOfState = EquationOfState<IdealGasModel, PcSaft>;

    fn build_eos(
        &self,
        (parameters, kij): ([Vec<f64>; 2], f64),
    ) -> Result<Arc<Self::EquationOfState>, ParameterError> {
        let (pcsaft, joback): (Vec<_>, Vec<_>) = parameters
            .into_iter()
            .map(|p| build_records(p, self.viscosity))
            .unzip();
        let kij = PcSaftBinaryRecord::new(Some(kij), None, None);
        let pcsaft = Arc::new(PcSaftParameters::new_binary(pcsaft, Some(kij)).unwrap());
        let joback = Arc::new(Joback::from_records(joback, None).unwrap());
        Ok(Arc::new(EquationOfState::new(
            Arc::new(IdealGasModel::Joback(joback)),
            Arc::new(PcSaft::new(pcsaft.clone())),
        )))
    }
}

fn build_records(
    parameters: Vec<f64>,
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
