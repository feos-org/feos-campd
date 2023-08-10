use super::MolecularRepresentation;
use feos_core::parameter::{Identifier, SegmentCount};
use std::{borrow::Cow, collections::HashMap};

const GROUPS: [&str; 22] = [
    "CH3", "CH2", ">CH", ">C<", "=CH2", "=CH", "=C<", "C≡CH", "CH2_hex", "CH_hex", "CH2_pent",
    "CH_pent", "CH_arom", "C_arom", "CH=O", ">C=O", "OCH3", "OCH2", "HCOO", "COO", "OH", "NH2",
];

#[derive(Clone)]
pub struct GroupCount(HashMap<String, f64>);

impl SegmentCount for GroupCount {
    type Count = f64;

    fn identifier(&self) -> Cow<Identifier> {
        Cow::Owned(Identifier::default())
    }

    fn segment_count(&self) -> Cow<HashMap<String, Self::Count>> {
        Cow::Borrowed(&self.0)
    }
}

pub struct CoMTCAMD;

impl MolecularRepresentation for CoMTCAMD {
    type ChemicalRecord = GroupCount;

    fn build(&self, y: Vec<f64>) -> Self::ChemicalRecord {
        GroupCount(
            GROUPS
                .iter()
                .map(|g| g.to_string())
                .zip(y.iter().copied())
                .collect(),
        )
    }

    fn variables(&self) -> usize {
        GROUPS.len() + 6
    }
}

#[cfg(feature = "knitro_rs")]
mod knitro {
    use super::*;
    use crate::knitro::MolecularRepresentationKnitro;
    use knitro_rs::{Knitro, KnitroError, KN_VARTYPE_BINARY, KN_VARTYPE_INTEGER};

    const STRUCTURES: [&str; 6] = ["alkane", "alkene", "arom", "hex", "pent", "P"];
    const OPEN_BONDS: [usize; 22] = [
        1, 2, 3, 4, 1, 2, 3, 1, 2, 3, 2, 3, 2, 3, 1, 2, 1, 2, 1, 2, 1, 1,
    ];
    const NMAX: [usize; 22] = [
        10, 10, 10, 10, 1, 2, 2, 1, 6, 6, 5, 5, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1,
    ];

    impl MolecularRepresentationKnitro for CoMTCAMD {
        fn setup_knitro(&self, kc: &Knitro) -> Result<(), KnitroError> {
            // group counts
            let n = kc.add_vars(GROUPS.len(), Some(KN_VARTYPE_INTEGER))?;
            for (&i, &nmax) in n.iter().zip(NMAX.iter()) {
                kc.set_var_lobnd(i, 0.0)?;
                kc.set_var_upbnd(i, nmax as f64)?;
            }
            let n: HashMap<_, _> = GROUPS.iter().copied().zip(n.into_iter()).collect();

            // structure variables
            let y = kc.add_vars(STRUCTURES.len(), Some(KN_VARTYPE_BINARY))?;

            // Exactly one structure can be active
            add_eq_constraint(kc, &y, &vec![1.0; y.len()], 1.0)?;
            let y: HashMap<_, _> = STRUCTURES.iter().copied().zip(y.into_iter()).collect();

            // Connect molecular structures and segments
            // Alkenes
            let vars = [n["=CH2"], n["=CH"], n["=C<"], y["alkene"]];
            add_eq_constraint(kc, &vars, &[1.0, 1.0, 1.0, -2.0], 0.0)?;

            // Aromatics
            let vars = [n["CH_arom"], n["C_arom"], y["arom"]];
            add_eq_constraint(kc, &vars, &[1.0, 1.0, -6.0], 0.0)?;

            // Cyclohexanes
            let vars = [n["CH2_hex"], n["CH_hex"], y["hex"]];
            add_eq_constraint(kc, &vars, &[1.0, 1.0, -6.0], 0.0)?;

            // Cyclopentanes
            let vars = [n["CH2_pent"], n["CH_pent"], y["pent"]];
            add_eq_constraint(kc, &vars, &[1.0, 1.0, -5.0], 0.0)?;

            // Polar groups
            let vars = [
                "C≡CH", "CH=O", ">C=O", "OCH3", "OCH2", "HCOO", "COO", "OH", "NH2",
            ];
            let vars: Vec<_> = vars
                .into_iter()
                .map(|g| n[g])
                .chain(std::iter::once(y["P"]))
                .collect();
            let coefs = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0];
            add_eq_constraint(kc, &vars, &coefs, 0.0)?;

            // octett rule
            let (vars, coefs): (Vec<_>, Vec<_>) = GROUPS
                .iter()
                .zip(OPEN_BONDS)
                .map(|(&g, o)| (n[g], 2.0 - o as f64))
                .chain([y["alkane"], y["alkene"], y["P"]].map(|i| (i, -2.0)))
                .unzip();
            add_eq_constraint(kc, &vars, &coefs, 0.0)?;

            // single molecule
            for g in GROUPS {
                let c = kc.add_con()?;
                let (vars, coefs): (Vec<_>, Vec<_>) = GROUPS
                    .iter()
                    .zip(OPEN_BONDS)
                    .map(|(&j, o)| (n[j], if j == g { 2.0 - o as f64 } else { 1.0 }))
                    .chain([y["alkane"], y["alkene"], y["P"]].map(|i| (i, -2.0)))
                    .unzip();
                kc.add_con_linear_struct_one(c, &vars, &coefs)?;
                kc.set_con_lobnd(c, 0.0)?;
            }

            Ok(())
        }
    }
    fn add_eq_constraint(
        kc: &Knitro,
        vars: &[i32],
        coefs: &[f64],
        eq: f64,
    ) -> Result<(), KnitroError> {
        let c = kc.add_con()?;
        kc.add_con_linear_struct_one(c, vars, coefs)?;
        kc.set_con_eqbnd(c, eq)?;
        Ok(())
    }
}
