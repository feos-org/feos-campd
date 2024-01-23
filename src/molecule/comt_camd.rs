use super::MolecularRepresentation;
use crate::variables::{ContinuousVariables, DiscreteVariables, LinearConstraint, Variable};
use crate::PropertyModel;
use feos::core::parameter::{
    IdentifierOption, Parameter, ParameterError, PureRecord, SegmentRecord,
};
use feos::core::EquationOfState;
use feos::ideal_gas::{Dippr, DipprRecord, JobackRecord};
use feos::pcsaft::{PcSaft, PcSaftParameters, PcSaftRecord};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::Arc;

const STRUCTURES: [&str; 7] = ["alkane", "alkene", "arom", "hex", "pent", "P", "M"];

const GROUPS: [(&str, usize, usize); 20] = [
    ("CH3", 1, 10),
    ("CH2", 2, 10),
    (">CH", 3, 10),
    (">C<", 4, 10),
    ("=CH2", 1, 1),
    ("=CH", 2, 2),
    ("=C<", 3, 2),
    ("C≡CH", 1, 1),
    ("CH2_hex", 2, 6),
    ("CH_hex", 3, 6),
    ("CH2_pent", 2, 5),
    ("CH_pent", 3, 5),
    ("CH_arom", 2, 6),
    ("C_arom", 3, 6),
    ("CH=O", 1, 1),
    (">C=O", 2, 1),
    ("OCH3", 1, 1),
    ("OCH2", 2, 1),
    ("HCOO", 1, 1),
    ("COO", 2, 1),
    // ("OH", 1, 1),
    // ("NH2", 1, 1),
];

const PARAMETER_OFFSETS: [f64; 10] = [
    0.0, 0.0, 0.0, 0.0, 0.0, -37.93, 0.21, -3.91e-4, 2.06e-7, 0.0,
];

#[derive(Clone, Serialize, Deserialize)]
pub struct CoMTCAMD {
    n_open_bonds: Vec<usize>,
    n_max: Vec<usize>,
    mw: Vec<f64>,
    m: Vec<f64>,
    m_sigma_3: Vec<f64>,
    m_epsilon_k: Vec<f64>,
    mu: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
    e: Vec<f64>,
}

impl CoMTCAMD {
    pub fn new(fixed_molecules: &[&str]) -> Result<Self, ParameterError> {
        // Read group parameters
        let segment_records = SegmentRecord::from_json("tests/sauer2014_homo.json")?;
        let mut segment_records: HashMap<_, _> = segment_records
            .into_iter()
            .map(|r| (r.identifier, (r.molarweight, r.model_record)))
            .collect();
        let segments = GROUPS
            .iter()
            .map(|&(g, v, n)| (segment_records.remove(g).unwrap(), v, n));

        // Read molecule parameters
        let molecules = fixed_molecules.iter().map(|s| {
            let record = &PcSaftParameters::from_json(
                vec![s],
                "tests/esper2023.json",
                None,
                IdentifierOption::Name,
            )
            .unwrap()
            .pure_records[0];
            ((record.molarweight, record.model_record.clone()), 0, 1)
        });

        // Extract PC-SAFT parameters
        let pcsaft: Vec<_> = segments.chain(molecules).collect();
        let mw: Vec<_> = pcsaft.iter().map(|((mw, _), _, _)| *mw).collect();
        let m: Vec<_> = pcsaft.iter().map(|((_, r), _, _)| r.m).collect();
        let m_sigma_3: Vec<_> = pcsaft
            .iter()
            .map(|((_, r), _, _)| r.m * r.sigma.powi(3))
            .collect();
        let m_epsilon_k: Vec<_> = pcsaft
            .iter()
            .map(|((_, r), _, _)| r.m * r.epsilon_k)
            .collect();
        let mu: Vec<_> = pcsaft
            .iter()
            .map(|((_, r), _, _)| r.mu.unwrap_or_default())
            .collect();

        // Extract valence and maximum count
        let (n_open_bonds, n_max): (Vec<_>, Vec<_>) =
            pcsaft.iter().map(|(_, v, n)| (*v, *n)).unzip();

        // Read group ideal gas parameters
        let segment_records = SegmentRecord::from_json("tests/joback1987.json")?;
        let mut segment_records: HashMap<_, _> = segment_records
            .into_iter()
            .map(|r| (r.identifier, r.model_record))
            .collect();
        let segments = GROUPS.iter().map(|&(g, _, _)| {
            let JobackRecord { a, b, c, d, e } = segment_records.remove(g).unwrap();
            (a, b, c, d, e)
        });

        // Read molecule ideal gas parameters
        let molecules = fixed_molecules.iter().map(|s| {
            let dippr = Dippr::from_json(
                vec![s],
                "tests/poling2000.json",
                None,
                IdentifierOption::Name,
            )
            .unwrap();
            let (dippr, _) = dippr.records();
            let DipprRecord::DIPPR100(v) = &dippr[0].model_record else {
                panic!("All components have to use Eq. 100!")
            };
            let &[a, b, c, d, e, ..] = &v[..] else {
                panic!("Eq. 100 needs to be provided with 5 coefficients!")
            };
            let [_, _, _, _, _, a0, b0, c0, d0, e0] = PARAMETER_OFFSETS;
            (a - a0, b - b0, c - c0, d - d0, e - e0)
        });

        // Extract ideal_gas parameters
        let ideal_gas: Vec<_> = segments.chain(molecules).collect();
        let a: Vec<_> = ideal_gas.iter().map(|(a, _, _, _, _)| *a).collect();
        let b: Vec<_> = ideal_gas.iter().map(|(_, b, _, _, _)| *b).collect();
        let c: Vec<_> = ideal_gas.iter().map(|(_, _, c, _, _)| *c).collect();
        let d: Vec<_> = ideal_gas.iter().map(|(_, _, _, d, _)| *d).collect();
        let e: Vec<_> = ideal_gas.iter().map(|(_, _, _, _, e)| *e).collect();

        Ok(Self {
            n_open_bonds,
            n_max,
            mw,
            m,
            m_sigma_3,
            m_epsilon_k,
            mu,
            a,
            b,
            c,
            d,
            e,
        })
    }
}

impl MolecularRepresentation<1> for CoMTCAMD {
    type ChemicalRecord = Vec<f64>;

    fn build(&self, _: &[f64], p: &[f64]) -> [Self::ChemicalRecord; 1] {
        [p.to_vec()]
    }

    fn structure_variables(&self) -> DiscreteVariables {
        self.n_max
            .iter()
            .map(|&n| Variable::integer(0, n))
            .chain([Variable::binary(); STRUCTURES.len()])
            .collect()
    }

    fn parameter_variables(&self) -> ContinuousVariables {
        vec![Variable::continuous(None, None); 10].into()
    }

    fn determine_parameters(&self, y: &[f64]) -> Vec<f64> {
        [
            &self.mw,
            &self.m,
            &self.m_sigma_3,
            &self.m_epsilon_k,
            &self.mu,
            &self.a,
            &self.b,
            &self.c,
            &self.d,
            &self.e,
        ]
        .iter()
        .map(|pars| pars.iter().zip(y).map(|(p, y)| p * y).sum())
        .collect()
    }

    fn constraints(
        &self,
        index_structure_vars: &[i32],
        index_parameter_vars: Option<&[i32]>,
    ) -> Vec<LinearConstraint> {
        let mut constraints = Vec::new();

        let (n_vars, y_vars) = index_structure_vars.split_at(self.m.len());

        let y: HashMap<_, _> = STRUCTURES.iter().copied().zip(y_vars.to_vec()).collect();
        let n: HashMap<_, _> = GROUPS
            .iter()
            .map(|&(g, _, _)| g)
            .zip(n_vars.to_vec())
            .collect();

        // Exactly one structure can be active
        constraints.push(LinearConstraint::new(y_vars.to_vec(), vec![1.0; y.len()]).eqbnd(1.0));

        // Connect molecular structures and segments
        // Alkenes
        let vars = vec![n["=CH2"], n["=CH"], n["=C<"], y["alkene"]];
        constraints.push(LinearConstraint::new(vars, vec![1.0, 1.0, 1.0, -2.0]).eqbnd(0.0));

        // Aromatics
        let vars = vec![n["CH_arom"], n["C_arom"], y["arom"]];
        constraints.push(LinearConstraint::new(vars, vec![1.0, 1.0, -6.0]).eqbnd(0.0));

        // Cyclohexanes
        let vars = vec![n["CH2_hex"], n["CH_hex"], y["hex"]];
        constraints.push(LinearConstraint::new(vars, vec![1.0, 1.0, -6.0]).eqbnd(0.0));

        // Cyclopentanes
        let vars = vec![n["CH2_pent"], n["CH_pent"], y["pent"]];
        constraints.push(LinearConstraint::new(vars, vec![1.0, 1.0, -5.0]).eqbnd(0.0));

        // Polar groups
        let vars = ["C≡CH", "CH=O", ">C=O", "OCH3", "OCH2", "HCOO", "COO"];
        let (vars, coefs): (Vec<_>, Vec<_>) = vars
            .into_iter()
            .map(|g| (n[g], 1.0))
            .chain(std::iter::once((y["P"], -1.0)))
            .unzip();
        constraints.push(LinearConstraint::new(vars, coefs).eqbnd(0.0));

        // Molecules
        let (vars, coefs): (Vec<_>, Vec<_>) = n_vars[GROUPS.len()..]
            .iter()
            .copied()
            .zip(std::iter::repeat(1.0))
            .chain(std::iter::once((y["M"], -1.0)))
            .unzip();
        constraints.push(LinearConstraint::new(vars, coefs).eqbnd(0.0));

        // octett rule
        let (vars, coefs): (Vec<_>, Vec<_>) = GROUPS
            .iter()
            .map(|&(g, o, _)| (n[g], 2.0 - o as f64))
            .chain([y["alkane"], y["alkene"], y["P"], y["M"]].map(|i| (i, -2.0)))
            .unzip();
        constraints.push(LinearConstraint::new(vars, coefs).eqbnd(0.0));

        // single molecule
        for (g, _, _) in GROUPS {
            let (vars, coefs): (Vec<_>, Vec<_>) = GROUPS
                .iter()
                .map(|&(j, o, _)| (n[j], if j == g { 2.0 - o as f64 } else { 1.0 }))
                .chain([y["alkane"], y["alkene"], y["P"], y["M"]].map(|i| (i, -2.0)))
                .unzip();
            constraints.push(LinearConstraint::new(vars, coefs).lobnd(0.0));
        }

        // PC-SAFT parameters
        if let Some(p) = index_parameter_vars {
            constraints.extend(
                p.iter()
                    .zip([
                        &self.mw,
                        &self.m,
                        &self.m_sigma_3,
                        &self.m_epsilon_k,
                        &self.mu,
                        &self.a,
                        &self.b,
                        &self.c,
                        &self.d,
                        &self.e,
                    ])
                    .zip(&PARAMETER_OFFSETS)
                    .map(|((&p, pars), &offset)| {
                        let (vars, coefs): (Vec<_>, Vec<_>) = n_vars
                            .iter()
                            .zip(pars)
                            .chain(std::iter::once((&p, &-1.0)))
                            .unzip();
                        LinearConstraint::new(vars, coefs).eqbnd(-offset)
                    }),
            );
        }

        constraints
    }

    fn smiles(&self, y: &[usize]) -> [String; 1] {
        let groups: Vec<_> = GROUPS
            .iter()
            .zip(y)
            .flat_map(|(&(g, _, _), &y)| vec![g; y])
            .collect();
        [groups.join(",")]
    }
}

pub struct CoMTCAMDPropertyModel;

impl PropertyModel<Vec<f64>> for CoMTCAMDPropertyModel {
    type EquationOfState = EquationOfState<Dippr, PcSaft>;

    fn build_eos(
        &self,
        chemical_records: &[Vec<f64>],
    ) -> Result<Arc<Self::EquationOfState>, ParameterError> {
        let (pcsaft, dippr): (Vec<_>, Vec<_>) = chemical_records
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
                let dippr = PureRecord::new(
                    Default::default(),
                    0.0,
                    DipprRecord::eq100(&[
                        a * 1000.0,
                        b * 1000.0,
                        c * 1000.0,
                        d * 1000.0,
                        e * 1000.0,
                    ]),
                );
                (pcsaft, dippr)
            })
            .unzip();
        let pcsaft = Arc::new(PcSaftParameters::from_records(pcsaft, None).unwrap());
        let dippr = Arc::new(Dippr::from_records(dippr, None).unwrap());
        Ok(Arc::new(EquationOfState::new(
            dippr,
            Arc::new(PcSaft::new(pcsaft.clone())),
        )))
    }
}
