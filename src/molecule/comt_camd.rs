use super::MolecularRepresentation;
use crate::ChemicalRecord;
use good_lp::{constraint, Constraint, Expression, Variable};
use num_dual::DualNum;
use quantity::MolarWeight;
use std::collections::HashMap;
use std::sync::LazyLock;

pub(super) struct Group {
    pub name: &'static str,
    open_bonds: usize,
    n_max: i32,
    heavy_atoms: i32,
    pub molar_weight: f64,
}

impl Group {
    const fn new(name: &'static str, open_bonds: usize, n_max: i32, hcno: [i32; 4]) -> Self {
        let [h, c, n, o] = hcno;
        let heavy_atoms = c + n + o;
        // values from PubChem https://pubchem.ncbi.nlm.nih.gov/ptable/atomic-mass/
        let molar_weight =
            h as f64 * 1.008 + c as f64 * 12.011 + n as f64 * 14.007 + o as f64 * 15.999;
        Self {
            name,
            open_bonds,
            n_max,
            heavy_atoms,
            molar_weight,
        }
    }
}

const N_GROUPS: usize = 22;
pub(super) const GROUPS: [Group; N_GROUPS] = [
    Group::new("CH3", 1, 10, [3, 1, 0, 0]),
    Group::new("CH2", 2, 10, [2, 1, 0, 0]),
    Group::new(">CH", 3, 10, [1, 1, 0, 0]),
    Group::new(">C<", 4, 10, [0, 1, 0, 0]),
    Group::new("=CH2", 1, 1, [2, 1, 0, 0]),
    Group::new("=CH", 2, 2, [1, 1, 0, 0]),
    Group::new("=C<", 3, 2, [0, 1, 0, 0]),
    Group::new("C≡CH", 1, 1, [1, 2, 0, 0]),
    Group::new("CH2_hex", 2, 6, [2, 1, 0, 0]),
    Group::new("CH_hex", 3, 6, [1, 1, 0, 0]),
    Group::new("CH2_pent", 2, 5, [2, 1, 0, 0]),
    Group::new("CH_pent", 3, 5, [1, 1, 0, 0]),
    Group::new("CH_arom", 2, 6, [1, 1, 0, 0]),
    Group::new("C_arom", 3, 6, [0, 1, 0, 0]),
    Group::new("CH=O", 1, 1, [1, 1, 0, 1]),
    Group::new(">C=O", 2, 1, [0, 1, 0, 1]),
    Group::new("OCH3", 1, 1, [3, 1, 0, 1]),
    Group::new("OCH2", 2, 1, [2, 1, 0, 1]),
    Group::new("HCOO", 1, 1, [1, 1, 0, 2]),
    Group::new("COO", 2, 1, [0, 1, 0, 2]),
    Group::new("OH", 1, 1, [1, 0, 0, 1]),
    Group::new("NH2", 1, 1, [2, 0, 1, 0]),
];

struct Structure {
    groups: Vec<&'static str>,
    count: usize,
    rings: usize,
}

impl Structure {
    fn new(groups: Vec<&'static str>, count: usize, rings: usize) -> Self {
        Self {
            groups,
            count,
            rings,
        }
    }
}

const POLAR_GROUPS: [&str; 9] = [
    "C≡CH", "CH=O", ">C=O", "OCH3", "OCH2", "HCOO", "COO", "OH", "NH2",
];
const N_STRUCTURES: usize = 6;
static STRUCTURES: LazyLock<[Structure; N_STRUCTURES]> = LazyLock::new(|| {
    [
        Structure::new(vec![], 0, 0),                      // alkanes
        Structure::new(vec!["=CH2", "=CH", "=C<"], 2, 0),  // alkene
        Structure::new(vec!["CH_arom", "C_arom"], 6, 1),   // aromatics
        Structure::new(vec!["CH2_hex", "CH_hex"], 6, 1),   // cyclohexanes
        Structure::new(vec!["CH2_pent", "CH_pent"], 5, 1), // cyclopentanes
        Structure::new(POLAR_GROUPS.to_vec(), 1, 0),       // polars
    ]
});

const N_Y: usize = N_GROUPS + N_STRUCTURES;

/// Molecular representation based on group counts, as used in 1-stage CoMT-CAMD.
pub struct CoMTCAMD;

impl MolecularRepresentation<N_Y> for CoMTCAMD {
    fn structure_variables(&self) -> [(i32, i32); N_Y] {
        let mut y_max = [(0, 0); N_Y];
        y_max[..N_STRUCTURES].fill((0, 1));
        y_max[N_STRUCTURES..].copy_from_slice(&GROUPS.map(|g| (0, g.n_max)));
        y_max
    }

    fn constraints(&self, y: [Variable; N_Y]) -> Vec<Constraint> {
        let (s, n) = y.split_at(STRUCTURES.len());
        let n_dict: HashMap<_, _> = GROUPS
            .into_iter()
            .zip(n)
            .map(|(g, v)| (g.name, v))
            .collect();

        let mut constraints = Vec::new();

        // exactly one structure can be active
        constraints.push(constraint!(s.iter().sum::<Expression>() == 1.0));

        // at least three heavy atoms
        let expr = n.iter().zip(GROUPS).map(|(&n, g)| n * g.heavy_atoms as f64);
        constraints.push(constraint!(expr.sum::<Expression>() >= 3.0));

        // Connect molecular structure and segments
        for (structure, &y) in STRUCTURES.iter().zip(s) {
            if !structure.groups.is_empty() {
                let expr: Expression = structure.groups.iter().map(|g| n_dict[g as &str]).sum();
                constraints.push(constraint!(expr == (structure.count as f64) * y));
            }
        }

        // octett rule
        let expr: Expression = n
            .iter()
            .zip(GROUPS)
            .map(|(&v, g)| v * (2.0 - g.open_bonds as f64))
            .chain(
                STRUCTURES
                    .iter()
                    .zip(s)
                    .map(|(s, &y)| 2.0 * (s.rings as f64 - 1.0) * y),
            )
            .sum();
        constraints.push(constraint!(expr == 0.0));

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

        constraints
    }

    fn smiles(&self, y: [f64; N_Y]) -> String {
        let (_, n) = y.split_at(STRUCTURES.len());
        let groups: Vec<_> = n
            .iter()
            .zip(GROUPS)
            .filter(|&(&c, _)| c >= 1.0)
            .map(|(&c, g)| format!("{}: {c}", g.name))
            .collect();
        format!("[{}]", groups.join(", "))
    }

    fn build_molecule<D: DualNum<f64> + Copy>(
        &self,
        feature_variables: [D; N_Y],
    ) -> ChemicalRecord<D> {
        let mut group_counts = [D::zero(); N_GROUPS];
        group_counts.copy_from_slice(&feature_variables[N_STRUCTURES..]);
        let mut molar_weight = D::from(0.0);
        let groups: HashMap<_, _> = GROUPS
            .iter()
            .zip(group_counts)
            .map(|(g, c)| {
                molar_weight += c * g.molar_weight;
                (g.name, c)
            })
            .collect();
        ChemicalRecord::new(
            groups,
            HashMap::new(),
            MolarWeight::new(molar_weight * 1e-3),
        )
    }
}
