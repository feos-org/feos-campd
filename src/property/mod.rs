use crate::ChemicalRecord;
use feos_ad::eos::{GcPcSaft, GcPcSaftParameters, Joback, PcSaftBinary, PcSaftPure};
use feos_ad::{EquationOfStateAD, HelmholtzEnergyWrapper, ParametersAD, TotalHelmholtzEnergy};
use feos_core::parameter::{BinaryRecord, ParameterError, SegmentRecord};
use num_dual::DualNum;
use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;

pub trait PropertyModel<const N: usize> {
    type EquationOfState: TotalHelmholtzEnergy<N>;

    fn build_eos<D: DualNum<f64> + Copy>(
        &self,
        chemical_records: [&ChemicalRecord<D>; N],
    ) -> HelmholtzEnergyWrapper<Self::EquationOfState, D, N>;
}

#[derive(Clone, Copy)]
pub struct GcPcSaftPropertyModel;

impl<const N: usize> PropertyModel<N> for GcPcSaftPropertyModel {
    type EquationOfState = EquationOfStateAD<Joback, GcPcSaft<N>, N>;

    fn build_eos<D: DualNum<f64> + Copy>(
        &self,
        chemical_records: [&ChemicalRecord<D>; N],
    ) -> HelmholtzEnergyWrapper<EquationOfStateAD<Joback, GcPcSaft<N>, N>, D, N> {
        let groups = chemical_records.map(|cr| &cr.groups);
        let bonds = chemical_records.map(|cr| &cr.bonds);
        let gc_pcsaft_params = GcPcSaftParameters::from_groups(groups, bonds);
        let gc_pcsaft = GcPcSaft(gc_pcsaft_params.re());

        let joback_params = groups.map(Joback::from_group_counts);
        let joback = joback_params.map(|j| Joback(j.map(|j| j.re())));

        EquationOfStateAD::new(joback, gc_pcsaft)
            .wrap()
            .derivatives((joback_params, gc_pcsaft_params))
    }
}

#[derive(Clone)]
pub struct PcSaftPropertyModel<const ASSOC: bool> {
    parameters: HashMap<String, PcSaftRecord>,
    binary_parameters: Option<HashMap<[String; 2], f64>>,
}

#[derive(Deserialize, Clone)]
struct PcSaftRecord {
    m: f64,
    sigma: f64,
    epsilon_k: f64,
    #[serde(default)]
    mu: f64,
    #[serde(default)]
    kappa_ab: f64,
    #[serde(default)]
    epsilon_k_ab: f64,
    #[serde(default)]
    na: f64,
    #[serde(default)]
    nb: f64,
}

impl PcSaftPropertyModel<true> {
    pub fn full<P: AsRef<Path>>(file: P, file_binary: Option<P>) -> Result<Self, ParameterError> {
        Self::new(file, file_binary)
    }
}

impl PcSaftPropertyModel<false> {
    pub fn non_associating<P: AsRef<Path>>(
        file: P,
        file_binary: Option<P>,
    ) -> Result<Self, ParameterError> {
        Self::new(file, file_binary)
    }
}

impl<const ASSOC: bool> PcSaftPropertyModel<ASSOC> {
    fn new<P: AsRef<Path>>(file: P, file_binary: Option<P>) -> Result<Self, ParameterError> {
        let records: Vec<SegmentRecord<PcSaftRecord>> = SegmentRecord::from_json(file)?;
        let parameters = records
            .into_iter()
            .map(|s| (s.identifier, s.model_record))
            .filter(|(_, r)| r.na == 0.0 || ASSOC)
            .collect();

        let binary_records: Option<Vec<BinaryRecord<String, f64>>> = file_binary
            .map(|f| BinaryRecord::from_json(f))
            .transpose()?;

        let binary_parameters = binary_records.map(|br| {
            let mut binary_parameters = HashMap::new();
            br.into_iter().for_each(|br| {
                binary_parameters.insert([br.id1.clone(), br.id2.clone()], br.model_record);
                binary_parameters.insert([br.id2, br.id1], br.model_record);
            });
            binary_parameters
        });

        Ok(Self {
            parameters,
            binary_parameters,
        })
    }
}

impl PcSaftPropertyModel<true> {
    pub fn from_group_counts<D: DualNum<f64> + Copy>(
        &self,
        groups: &HashMap<&'static str, D>,
    ) -> [D; 8] {
        let mut m = D::zero();
        let mut m_sigma3 = D::zero();
        let mut m_epsilon_k = D::zero();
        let mut mu = D::zero();
        let mut kappa_ab = D::zero();
        let mut epsilon_k_ab = D::zero();
        let mut na = D::zero();
        let mut nb = D::zero();

        groups.iter().for_each(|(&g, &n)| {
            let p = &self.parameters[g];
            m += n * p.m;
            m_sigma3 += n * p.m * p.sigma.powi(3);
            m_epsilon_k += n * p.m * p.epsilon_k;
            mu += n * p.mu;
            kappa_ab += n * p.kappa_ab;
            epsilon_k_ab += n * p.epsilon_k_ab;
            na += n * p.na;
            nb += n * p.nb;
        });
        let sigma = (m_sigma3 / m).cbrt();
        let epsilon_k = m_epsilon_k / m;

        [m, sigma, epsilon_k, mu, kappa_ab, epsilon_k_ab, na, nb]
    }
}

impl PcSaftPropertyModel<false> {
    pub fn from_group_counts<D: DualNum<f64> + Copy>(
        &self,
        groups: &HashMap<&'static str, D>,
    ) -> [D; 4] {
        let mut m = D::zero();
        let mut m_sigma3 = D::zero();
        let mut m_epsilon_k = D::zero();
        let mut mu = D::zero();

        groups.iter().for_each(|(&g, &n)| {
            let p = &self.parameters[g];
            m += n * p.m;
            m_sigma3 += n * p.m * p.sigma.powi(3);
            m_epsilon_k += n * p.m * p.epsilon_k;
            mu += n * p.mu;
        });
        let sigma = (m_sigma3 / m).cbrt();
        let epsilon_k = m_epsilon_k / m;

        [m, sigma, epsilon_k, mu]
    }
}

impl<const ASSOC: bool> PcSaftPropertyModel<ASSOC> {
    pub fn kij_from_group_counts<D: DualNum<f64> + Copy>(
        &self,
        groups: [&HashMap<&'static str, D>; 2],
    ) -> D {
        let Some(binary_parameters) = &self.binary_parameters else {
            return D::zero();
        };
        let mut kij = D::zero();
        let mut n = D::zero();
        for (&g1, &c1) in groups[0].iter() {
            for (&g2, &c2) in groups[1].iter() {
                let kab = *binary_parameters
                    .get(&[g1.to_string(), g2.to_string()])
                    .unwrap_or(&0.0);
                n += c1 * c2;
                kij += c1 * c2 * kab;
            }
        }

        kij / n
    }
}

impl PropertyModel<1> for PcSaftPropertyModel<true> {
    type EquationOfState = EquationOfStateAD<Joback, PcSaftPure<8>, 1>;

    fn build_eos<D: DualNum<f64> + Copy>(
        &self,
        [cr]: [&ChemicalRecord<D>; 1],
    ) -> HelmholtzEnergyWrapper<Self::EquationOfState, D, 1> {
        let joback_params = Joback::from_group_counts(&cr.groups);
        let joback = Joback(joback_params.map(|x| x.re()));

        let pcsaft_params = self.from_group_counts(&cr.groups);
        let pcsaft = PcSaftPure(pcsaft_params.map(|x| x.re()));

        EquationOfStateAD::new([joback], pcsaft)
            .wrap()
            .derivatives(([joback_params], pcsaft_params))
    }
}

impl PropertyModel<1> for PcSaftPropertyModel<false> {
    type EquationOfState = EquationOfStateAD<Joback, PcSaftPure<4>, 1>;

    fn build_eos<D: DualNum<f64> + Copy>(
        &self,
        [cr]: [&ChemicalRecord<D>; 1],
    ) -> HelmholtzEnergyWrapper<Self::EquationOfState, D, 1> {
        let joback_params = Joback::from_group_counts(&cr.groups);
        let joback = Joback(joback_params.map(|x| x.re()));

        let pcsaft_params = self.from_group_counts(&cr.groups);
        let pcsaft = PcSaftPure(pcsaft_params.map(|x| x.re()));

        EquationOfStateAD::new([joback], pcsaft)
            .wrap()
            .derivatives(([joback_params], pcsaft_params))
    }
}

impl PropertyModel<2> for PcSaftPropertyModel<true> {
    type EquationOfState = EquationOfStateAD<Joback, PcSaftBinary<8>, 2>;

    fn build_eos<D: DualNum<f64> + Copy>(
        &self,
        chemical_records: [&ChemicalRecord<D>; 2],
    ) -> HelmholtzEnergyWrapper<Self::EquationOfState, D, 2> {
        let joback_params = chemical_records.map(|cr| Joback::from_group_counts(&cr.groups));
        let joback = joback_params.map(|j| Joback(j.map(|x| x.re())));

        let kij = self.kij_from_group_counts(chemical_records.map(|cr| &cr.groups));
        let pcsaft_params = chemical_records.map(|cr| self.from_group_counts(&cr.groups));
        let pcsaft = PcSaftBinary::new(pcsaft_params.map(|p| p.map(|p| p.re())), kij.re());

        EquationOfStateAD::new(joback, pcsaft)
            .wrap()
            .derivatives((joback_params, (pcsaft_params, kij)))
    }
}

impl PropertyModel<2> for PcSaftPropertyModel<false> {
    type EquationOfState = EquationOfStateAD<Joback, PcSaftBinary<4>, 2>;

    fn build_eos<D: DualNum<f64> + Copy>(
        &self,
        chemical_records: [&ChemicalRecord<D>; 2],
    ) -> HelmholtzEnergyWrapper<Self::EquationOfState, D, 2> {
        let joback_params = chemical_records.map(|cr| Joback::from_group_counts(&cr.groups));
        let joback = joback_params.map(|j| Joback(j.map(|x| x.re())));

        let kij = self.kij_from_group_counts(chemical_records.map(|cr| &cr.groups));
        let pcsaft_params = chemical_records.map(|cr| self.from_group_counts(&cr.groups));
        let pcsaft = PcSaftBinary::new(pcsaft_params.map(|p| p.map(|p| p.re())), kij.re());

        EquationOfStateAD::new(joback, pcsaft)
            .wrap()
            .derivatives((joback_params, (pcsaft_params, kij)))
    }
}
