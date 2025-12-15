use crate::ChemicalRecord;
use feos::core::parameter::{BinaryRecord, SegmentRecord};
use feos::core::{EquationOfState, FeosResult, Total};
use feos::gc_pcsaft::{GcPcSaftAD, GcPcSaftADParameters};
use feos::ideal_gas::Joback;
use feos::pcsaft::{PcSaftAssociationRecord, PcSaftBinary, PcSaftPure};
use nalgebra::Const;
use num_dual::DualNum;
use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;

/// A generic property model to be used in an [IntegratedDesign](../IntegratedDesign)
pub trait PropertyModel<const N: usize> {
    type EquationOfState<D: DualNum<f64> + Copy>: Total<Const<N>, D>;

    fn build_eos<D: DualNum<f64> + Copy>(
        &self,
        chemical_records: [&ChemicalRecord<D>; N],
    ) -> Self::EquationOfState<D>;
}

/// The heterosegmented gc-PC-SAFT equation of state.
#[derive(Clone, Copy)]
pub struct GcPcSaftPropertyModel;

impl<const N: usize> PropertyModel<N> for GcPcSaftPropertyModel {
    type EquationOfState<D: DualNum<f64> + Copy> =
        EquationOfState<[Joback<D>; N], GcPcSaftAD<D, N>>;

    fn build_eos<D: DualNum<f64> + Copy>(
        &self,
        chemical_records: [&ChemicalRecord<D>; N],
    ) -> Self::EquationOfState<D> {
        let groups = chemical_records.map(|cr| &cr.groups);
        let bonds = chemical_records.map(|cr| &cr.bonds);
        let gc_pcsaft_params = GcPcSaftADParameters::from_groups(groups, bonds);
        let gc_pcsaft = GcPcSaftAD(gc_pcsaft_params);

        let joback = groups.map(Joback::from_group_counts);

        EquationOfState::new(joback, gc_pcsaft)
    }
}

/// The homosegmented group contribution method for the PC-SAFT equation of state.
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
    pub fn full<P: AsRef<Path>>(file: P, file_binary: Option<P>) -> FeosResult<Self> {
        Self::new(file, file_binary)
    }
}

impl PcSaftPropertyModel<false> {
    pub fn non_associating<P: AsRef<Path>>(file: P, file_binary: Option<P>) -> FeosResult<Self> {
        Self::new(file, file_binary)
    }
}

impl<const ASSOC: bool> PcSaftPropertyModel<ASSOC> {
    fn new<P: AsRef<Path>>(file: P, file_binary: Option<P>) -> FeosResult<Self> {
        let records: Vec<SegmentRecord<PcSaftRecord, PcSaftAssociationRecord>> =
            SegmentRecord::from_json(file)?;
        let parameters = records
            .into_iter()
            .map(|s| (s.identifier, s.model_record))
            .filter(|(_, r)| r.na == 0.0 || ASSOC)
            .collect();

        let binary_records: Option<Vec<BinaryRecord<String, f64, PcSaftAssociationRecord>>> =
            file_binary
                .map(|f| BinaryRecord::from_json(f))
                .transpose()?;

        let binary_parameters = binary_records.map(|br| {
            let mut binary_parameters = HashMap::new();
            br.into_iter().for_each(|br| {
                if let Some(m) = br.model_record {
                    binary_parameters.insert([br.id1.clone(), br.id2.clone()], m);
                    binary_parameters.insert([br.id2, br.id1], m);
                }
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
    type EquationOfState<D: DualNum<f64> + Copy> =
        EquationOfState<[Joback<D>; 1], PcSaftPure<D, 8>>;

    fn build_eos<D: DualNum<f64> + Copy>(
        &self,
        [cr]: [&ChemicalRecord<D>; 1],
    ) -> Self::EquationOfState<D> {
        let joback = Joback::from_group_counts(&cr.groups);

        let pcsaft_params = self.from_group_counts(&cr.groups);
        let pcsaft = PcSaftPure(pcsaft_params);

        EquationOfState::new([joback], pcsaft)
    }
}

impl PropertyModel<1> for PcSaftPropertyModel<false> {
    type EquationOfState<D: DualNum<f64> + Copy> =
        EquationOfState<[Joback<D>; 1], PcSaftPure<D, 4>>;

    fn build_eos<D: DualNum<f64> + Copy>(
        &self,
        [cr]: [&ChemicalRecord<D>; 1],
    ) -> Self::EquationOfState<D> {
        let joback = Joback::from_group_counts(&cr.groups);

        let pcsaft_params = self.from_group_counts(&cr.groups);
        let pcsaft = PcSaftPure(pcsaft_params);

        EquationOfState::new([joback], pcsaft)
    }
}

impl PropertyModel<2> for PcSaftPropertyModel<true> {
    type EquationOfState<D: DualNum<f64> + Copy> =
        EquationOfState<[Joback<D>; 2], PcSaftBinary<D, 8>>;

    fn build_eos<D: DualNum<f64> + Copy>(
        &self,
        chemical_records: [&ChemicalRecord<D>; 2],
    ) -> Self::EquationOfState<D> {
        let joback = chemical_records.map(|cr| Joback::from_group_counts(&cr.groups));

        let kij = self.kij_from_group_counts(chemical_records.map(|cr| &cr.groups));
        let pcsaft_params = chemical_records.map(|cr| self.from_group_counts(&cr.groups));
        let pcsaft = PcSaftBinary::new(pcsaft_params, kij);

        EquationOfState::new(joback, pcsaft)
    }
}

impl PropertyModel<2> for PcSaftPropertyModel<false> {
    type EquationOfState<D: DualNum<f64> + Copy> =
        EquationOfState<[Joback<D>; 2], PcSaftBinary<D, 4>>;

    fn build_eos<D: DualNum<f64> + Copy>(
        &self,
        chemical_records: [&ChemicalRecord<D>; 2],
    ) -> Self::EquationOfState<D> {
        let joback = chemical_records.map(|cr| Joback::from_group_counts(&cr.groups));

        let kij = self.kij_from_group_counts(chemical_records.map(|cr| &cr.groups));
        let pcsaft_params = chemical_records.map(|cr| self.from_group_counts(&cr.groups));
        let pcsaft = PcSaftBinary::new(pcsaft_params, kij);

        EquationOfState::new(joback, pcsaft)
    }
}
