use crate::molecule::IdealGasRecord;
use crate::variables::{Constraint, ExplicitVariable};
use crate::{CoMTCAMD, PropertyModel};
use feos::core::parameter::{Parameter, ParameterError, PureRecord};
use feos::core::EquationOfState;
use feos::ideal_gas::{DipprRecord, IdealGasModel, Joback, JobackRecord};
use feos::pcsaft::{PcSaft, PcSaftParameters, PcSaftRecord};
use indexmap::IndexMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;

const JOBACK_PARAMETER_OFFSETS: [f64; 5] = [-37.93, 0.21, -3.91e-4, 2.06e-7, 0.0];

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

impl PropertyModel<1> for CoMTCAMDPropertyModel {
    type EquationOfState = EquationOfState<IdealGasModel, PcSaft>;

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
