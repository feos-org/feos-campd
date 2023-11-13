use super::{LinearConstraint, MolecularRepresentation, Variable};
use feos_core::parameter::{Identifier, SegmentCount};
use serde::{Deserialize, Serialize};
use std::{borrow::Cow, collections::HashMap};

const STRUCTURES: [&str; 6] = ["alkane", "alkene", "arom", "hex", "pent", "P"];

const GROUPS: [&str; 22] = [
    "CH3", "CH2", ">CH", ">C<", "=CH2", "=CH", "=C<", "C≡CH", "CH2_hex", "CH_hex", "CH2_pent",
    "CH_pent", "CH_arom", "C_arom", "CH=O", ">C=O", "OCH3", "OCH2", "HCOO", "COO", "OH", "NH2",
];

const OPEN_BONDS: [usize; 22] = [
    1, 2, 3, 4, 1, 2, 3, 1, 2, 3, 2, 3, 2, 3, 1, 2, 1, 2, 1, 2, 1, 1,
];
const NMAX: [usize; 22] = [
    10, 10, 10, 10, 1, 2, 2, 1, 6, 6, 5, 5, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1,
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

#[derive(Clone, Copy, Serialize, Deserialize)]
pub struct CoMTCAMD;

impl MolecularRepresentation<1> for CoMTCAMD {
    type ChemicalRecord = GroupCount;

    fn build(&self, y: &[f64]) -> [Self::ChemicalRecord; 1] {
        [GroupCount(
            GROUPS
                .iter()
                .map(|g| g.to_string())
                .zip(y.iter().copied())
                .collect(),
        )]
    }

    fn variables(&self) -> Vec<Variable> {
        NMAX.iter()
            .map(|&n| Variable::integer(None, 0, n))
            .chain([Variable::binary(None); 6])
            .collect()
    }

    fn constraints(&self, index_vars: &[i32]) -> Vec<LinearConstraint> {
        let mut constraints = Vec::new();

        let y_vars = index_vars[GROUPS.len()..].to_vec();
        let y: HashMap<_, _> = STRUCTURES.iter().copied().zip(y_vars.clone()).collect();
        let n: HashMap<_, _> = GROUPS.iter().copied().zip(index_vars.to_vec()).collect();

        // Exactly one structure can be active
        constraints.push(LinearConstraint::new(y_vars, vec![1.0; y.len()]).eqbnd(1.0));

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
        let vars = [
            "C≡CH", "CH=O", ">C=O", "OCH3", "OCH2", "HCOO", "COO", "OH", "NH2",
        ];
        let vars: Vec<_> = vars
            .into_iter()
            .map(|g| n[g])
            .chain(std::iter::once(y["P"]))
            .collect();
        let coefs = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0];
        constraints.push(LinearConstraint::new(vars, coefs).eqbnd(0.0));

        // octett rule
        let (vars, coefs): (Vec<_>, Vec<_>) = GROUPS
            .iter()
            .zip(OPEN_BONDS)
            .map(|(&g, o)| (n[g], 2.0 - o as f64))
            .chain([y["alkane"], y["alkene"], y["P"]].map(|i| (i, -2.0)))
            .unzip();
        constraints.push(LinearConstraint::new(vars, coefs).eqbnd(0.0));

        // single molecule
        for g in GROUPS {
            let (vars, coefs): (Vec<_>, Vec<_>) = GROUPS
                .iter()
                .zip(OPEN_BONDS)
                .map(|(&j, o)| (n[j], if j == g { 2.0 - o as f64 } else { 1.0 }))
                .chain([y["alkane"], y["alkene"], y["P"]].map(|i| (i, -2.0)))
                .unzip();
            constraints.push(LinearConstraint::new(vars, coefs).lobnd(0.0));
        }

        constraints
    }

    fn smiles(&self, _: &[usize]) -> [String; 1] {
        [String::new()]
    }
}
